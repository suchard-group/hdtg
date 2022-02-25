#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::ArrayXd;
using std::atan2;


// Unused at the moment
// // [[Rcpp::export]]
// VectorXd RcppMatVec(Map<MatrixXd> A, Map<VectorXd> b) {
// return A * b;
// }


void UpdateHamiltonian(VectorXd &position, 
                       VectorXd &momentum, 
                       const double time) {
  VectorXd new_position = momentum * sin(time) + position * cos(time);
  momentum = momentum * cos(time) - position * sin(time);
  position = new_position;
}


void ReflectMomentum(const VectorXd &position,
                     VectorXd &momentum,
                     const Map<MatrixXd> constraint_direc,
                     const Map<VectorXd> constraint_row_normsq,
                     const int bounce_idx) {
  momentum = momentum - 2 * momentum.dot(constraint_direc.row(bounce_idx-1)) /  
    constraint_row_normsq(bounce_idx-1) * constraint_direc.row(bounce_idx-1).transpose();
}


std::pair<double, int> RcppBounceTime(const VectorXd &position,
                                      const VectorXd &momentum,
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
  //return Rcpp::List::create(Rcpp::_["bounce_time"]=min_time, 
  //                          Rcpp::_["constraint_idx"]=constraint_idx);
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
  //Map<VectorXd> momentum = Rcpp::rnorm(constraint_direc.cols());  // Need to figure out how to do type conversion
  double travelled_time = 0;
  double bounce_time;
  while (true) {
    std::pair<double, int> bounce = RcppBounceTime(position,
                                                   momentum,
                                                   constraint_direc,
                                                   constraint_bound);
    if (bounce.first < total_time - travelled_time) {
      bounce_time = bounce.first;
      UpdateHamiltonian(position, momentum, bounce_time);
      ReflectMomentum(
        position,
        momentum,
        constraint_direc,
        constraint_row_normsq,
        bounce.second
      );
      travelled_time += bounce_time;
    } else {
      bounce_time = total_time - travelled_time;
      UpdateHamiltonian(position, momentum, bounce_time);
      return position;
    }
  }
}

