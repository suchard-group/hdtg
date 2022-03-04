#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;


std::pair<VectorXd, VectorXd> Hamiltonian(const VectorXd position, 
                                          const VectorXd momentum, 
                                          const double time) {
  return std::make_pair(momentum * sin(time) + position * cos(time), 
                        momentum * cos(time) - position * sin(time));
}


VectorXd ReflectMomentum(const VectorXd position,
                         const VectorXd momentum,
                         const Map<MatrixXd> constraint_direc,
                         const Map<VectorXd> constraint_row_normsq,
                         const int bounce_idx) {
  return momentum - 2 * momentum.dot(constraint_direc.row(bounce_idx-1)) /  
    constraint_row_normsq(bounce_idx-1) * constraint_direc.row(bounce_idx-1).transpose();
}


std::pair<double, int> BounceTime(const VectorXd position,
                                  const VectorXd momentum,
                                  const Map<MatrixXd> constraint_direc,
                                  const Map<VectorXd> constraint_bound) {
  
  ArrayXd fa = (constraint_direc * momentum).array();
  ArrayXd fb = (constraint_direc * position).array();
  ArrayXd U = (fa.square() + fb.square()).sqrt();
  ArrayXd phi = -fa.binaryExpr(fb, [] (double a, double b) { return std::atan2(a,b);} );
  double min_time = std::numeric_limits<double>::infinity();
  int constraint_idx = -1;
  for (int i=0; i<constraint_bound.size(); ++i){
    if (U[i]>abs(constraint_bound[i])) {
      double bounce_time = -phi[i] + std::acos(-constraint_bound[i]/ U[i]);
      if (bounce_time < min_time) {
        min_time = bounce_time;
        constraint_idx = i+1;
      }
    }
  }
  return std::make_pair(min_time, constraint_idx);
}


// [[Rcpp::export]]
VectorXd GenerateWhitenedSample(const Map<VectorXd> initial_position,
                                const Map<VectorXd> initial_momentum,
                                const Map<MatrixXd> constraint_direc,
                                const Map<VectorXd> constraint_row_normsq,
                                const Map<VectorXd> constraint_bound,
                                float total_time){
  VectorXd position = initial_position;
  VectorXd momentum = initial_momentum;
  double travelled_time = 0;
  while (true) {
    std::pair<double, int> bounce = BounceTime(position,
                                               momentum,
                                               constraint_direc,
                                               constraint_bound);
    if (bounce.first < total_time - travelled_time) {
      double bounce_time = bounce.first;
      std::pair<VectorXd, VectorXd> hamiltonian = Hamiltonian(position, momentum, bounce_time);
      position = hamiltonian.first;
      momentum = ReflectMomentum(
        position,
        hamiltonian.second,
        constraint_direc,
        constraint_row_normsq,
        bounce.second
      );
      travelled_time += bounce_time;
    } else {
      double bounce_time = total_time - travelled_time;
      std::pair<VectorXd, VectorXd> hamiltonian = Hamiltonian(position, momentum, bounce_time);
      return hamiltonian.first;
    }
  }
}


// [[Rcpp::export]]
VectorXd WhitenPosition(const Map<VectorXd> position,
                        const Map<MatrixXd> constraint_direc,
                        const Map<VectorXd> constraint_bound,
                        const Map<MatrixXd> cholesky,
                        const Map<VectorXd> mean,
                        bool precision) {
  
  if (precision) {
    return cholesky * (position - mean);
  } else {
    return cholesky.transpose().triangularView<Eigen::Lower>().solve(position-mean);
  }
}


// [[Rcpp::export]]
VectorXd UnwhitenPosition(const Map<VectorXd> position,
                          const Map<MatrixXd> cholesky,
                          const Map<VectorXd> mean,
                          bool precision) {
  
  if (precision) {
    return cholesky.triangularView<Eigen::Upper>().solve(position) + mean;
  } else {
    return cholesky.transpose() * position + mean;
  }
}

                        