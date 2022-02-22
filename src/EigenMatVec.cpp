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

// [[Rcpp::export]]
Rcpp::List RcppHamiltonian(const Map<VectorXd> position, const Map<VectorXd> momentum, const double time) {
  return Rcpp::List::create(Rcpp::_["position"]=momentum * sin(time) + position * cos(time), 
                            Rcpp::_["momentum"]=momentum * cos(time) - position * sin(time));
}

// [[Rcpp::export]]
VectorXd RcppBounceMomentum(const Map<VectorXd> position,
                            const Map<VectorXd> momentum,
                            const Map<MatrixXd> constraint_direc,
                            const Map<MatrixXd> constraint_row_normsq,
                            const int bounce_idx) {
  return momentum - 2 * momentum.dot(constraint_direc.row(bounce_idx-1)) /  constraint_row_normsq(bounce_idx-1) * constraint_direc.row(bounce_idx-1).transpose();
}

// [[Rcpp::export]]
Rcpp::List RcppBounceTime(const Map<VectorXd> position,
                          const Map<VectorXd> momentum,
                          const Map<MatrixXd> constraint_direc,
                          const Map<VectorXd> constraint_bound) {
  
  ArrayXd fa = (constraint_direc * momentum).array();
  ArrayXd fb = (constraint_direc * position).array();
  ArrayXd U = (fa.square() + fb.square()).sqrt();
  ArrayXd phi = -fa.binaryExpr(fb, [] (double a, double b) { return std::atan2(a,b);} );
  double bounce_time = std::numeric_limits<double>::infinity();
  int constraint_idx = -1;
  for (int i=0; i<constraint_bound.size(); ++i){
    if (U[i]>abs(constraint_bound[i])) {
      double time = -phi[i] + std::acos(-constraint_bound[i]/ U[i]);
      if (time < bounce_time) {
        bounce_time = time;
        constraint_idx = i+1;
      }
    }
  }
  return Rcpp::List::create(Rcpp::_["bounce_time"]=bounce_time, 
                            Rcpp::_["constraint_idx"]=constraint_idx);
}
