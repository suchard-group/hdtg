#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;


// [[Rcpp::export]]
Rcpp::List WhitenConstraints(const Map<MatrixXd> constraint_direc,
                             const Map<VectorXd> constraint_bound,
                             const Map<MatrixXd> cholesky,
                             const Map<VectorXd> mean,
                             bool precision) {
  if (precision) {
    ArrayXXd direc =  cholesky.transpose().triangularView<Eigen::Lower>().solve(
      constraint_direc.transpose()).transpose().array();
    return Rcpp::List::create(
      Rcpp::_["direc"] = direc, 
      Rcpp::_["direc_rownorm_sq"]= direc.square().rowwise().sum(),
      Rcpp::_["bound"] = constraint_bound + constraint_direc * mean);
  } else {
    ArrayXXd direc =  constraint_direc * cholesky.transpose();
    return Rcpp::List::create(
      Rcpp::_["direc"] = direc, 
      Rcpp::_["direc_rownorm_sq"]=direc.square().rowwise().sum(),
      Rcpp::_["bound"] = constraint_bound + constraint_direc * mean);
  }
}


//' Compute hamiltonian after specified time.
//'
//' @param position starting position
//' @param momentum starting momentum
//' @param time amount of time the system is run for
//' @return pair of new (position, momentum)
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
  // Eigen doesn't have an atan2 function
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


VectorXd UnwhitenPosition(const VectorXd position,
                          const Map<MatrixXd> cholesky,
                          const Map<VectorXd> mean,
                          bool precision) {
  
  if (precision) {
    return cholesky.triangularView<Eigen::Upper>().solve(position) + mean;
  } else {
    return cholesky.transpose() * position + mean;
  }
}


VectorXd GenerateWhitenedSample(const VectorXd initial_position,
                                const Map<VectorXd> initial_momentum,
                                const Map<MatrixXd> constraint_direc,
                                const Map<VectorXd> constraint_row_normsq,
                                const Map<VectorXd> constraint_bound,
                                float total_time){
  double bounce_time;
  double travelled_time = 0;
  VectorXd position = initial_position;
  VectorXd momentum = initial_momentum;
  while (true) {
    std::pair<double, int> bounce = BounceTime(position,
                                               momentum,
                                               constraint_direc,
                                               constraint_bound);
    if (bounce.first < total_time - travelled_time) {
      bounce_time = bounce.first;
      std::pair<VectorXd, VectorXd> hamiltonian = Hamiltonian(position, momentum, bounce_time);
      position = hamiltonian.first;
      momentum = ReflectMomentum(position, 
                                 hamiltonian.second, 
                                 constraint_direc, 
                                 constraint_row_normsq, 
                                 bounce.second);
      travelled_time += bounce_time;
    } else {
      bounce_time = total_time - travelled_time;
      std::pair<VectorXd, VectorXd> hamiltonian = Hamiltonian(position, momentum, bounce_time);
      return hamiltonian.first;
    }
  }
}


// [[Rcpp::export]]
VectorXd GenerateSample(const Map<VectorXd> initial_position,
                        const Map<VectorXd> initial_momentum,
                        const Map<MatrixXd> constraint_direc,
                        const Map<VectorXd> constraint_row_normsq,
                        const Map<VectorXd> constraint_bound,
                        const Map<MatrixXd> cholesky,
                        const Map<VectorXd> mean,
                        float total_time,
                        bool precision){
  VectorXd sample = WhitenPosition(initial_position,
                                   constraint_direc,
                                   constraint_bound,
                                   cholesky,
                                   mean,
                                   precision);
  sample = GenerateWhitenedSample(sample,
                                  initial_momentum,
                                  constraint_direc,
                                  constraint_row_normsq,
                                  constraint_bound,
                                  total_time);
  return UnwhitenPosition(sample, cholesky, mean, precision);
}

