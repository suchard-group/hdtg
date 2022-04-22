#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

//' Compute Cholesky decomposition of a matrix.
//'
//' @param A matrix to decompose
//' @return upper triangular matrix R such that A = R'R.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cholesky(const Eigen::Map<Eigen::MatrixXd> A) {
  return A.llt().matrixU();
}

//' Solve XA=B for two matrices A and B.
//'
//' Eigen does not have builtin methods to solve XA=B, only AX=B, so we
//' solve A'X'= B'
//'
//' @param A matrix
//' @param B matrix
//' @return X matrix
Eigen::MatrixXd solveFromRight(const Eigen::Map<Eigen::MatrixXd> A,
                               const Eigen::Map<Eigen::MatrixXd> B) {
  return A.transpose()
      .triangularView<Eigen::Lower>()
      .solve(B.transpose())
      .transpose();
}

//' Whiten constraints for use in generateUnwhitenedSample
//'
//' Transforms constraints of the form Fx+g >= 0 for a target normal 
//' distribution into the corresponding constraints for a standard normal.
//'
//' @param constraint_direc F matrix (k-by-d matrix where k is the number of 
//' linear constraints)
//' @param constraint_bound g vector (k dimensional)
//' @param cholesky_factor upper triangular matrix R from cholesky decomposition
//'  of precision or covariance matrix into R^TR
//' @param unconstrained_mean mean of unconstrained Gaussian
//' @param prec_parametrized boolean for whether parametrization is by precision
//'  (true) or covariance matrix (false)
//' @return List of new constraint directions, the squared row norms of those 
//' constraints (for computational efficiency later), and new bounds
//' @export
// [[Rcpp::export]]
Rcpp::List applyWhitenTransform(
    const Eigen::Map<Eigen::MatrixXd> constraint_direc,
    const Eigen::Map<Eigen::VectorXd> constraint_bound,
    const Eigen::Map<Eigen::MatrixXd> cholesky_factor,
    const Eigen::Map<Eigen::VectorXd> unconstrained_mean,
    bool prec_parametrized) {
  Eigen::ArrayXXd direc;
  if (prec_parametrized) {
    direc = solveFromRight(cholesky_factor, constraint_direc).array();
  } else {
    direc = constraint_direc * cholesky_factor.transpose();
  }
  return Rcpp::List::create(
      Rcpp::_["direc"] = direc,
      Rcpp::_["direc_row_norm_sq"] = direc.square().rowwise().sum(),
      Rcpp::_["bound"] =
          constraint_bound + constraint_direc * unconstrained_mean);
}

//' Compute Hamiltonian dynamics after specified time.
//'
//' @param position starting position
//' @param momentum starting momentum
//' @param time amount of time the system is run for
//' @return pair of new (position, momentum)
std::pair<Eigen::VectorXd, Eigen::VectorXd> advanceWhitenedDynamics(
    const Eigen::VectorXd position, const Eigen::VectorXd momentum,
    const double time) {
  return std::make_pair(momentum * sin(time) + position * cos(time),
                        momentum * cos(time) - position * sin(time));
}

//' Reflect momentum off of a constraint boundary.
//'
//' Given a constraint boundary, calculate the momentum as if that boundary
//' was a wall and there is an elastic collision, and the angle of incidence
//' equals the angle of reflection.
//'
//' @param momentum starting momentum
//' @param constraint_direc F matrix (k-by-d matrix where k is the number of
//' linear constraints)
//' @param constraint_row_normsq vector of squared row norms of constraint_direc
//' @param bounce_idx integer index of which constraint is being bounced off of
//' @param time amount of time the system is run for
//' @return momentum after bouncing
Eigen::VectorXd reflectMomentum(
    const Eigen::VectorXd momentum,
    const Eigen::Map<Eigen::MatrixXd> constraint_direc,
    const Eigen::Map<Eigen::VectorXd> constraint_row_norm_sq,
    const int bounce_idx) {
  return momentum - 2 * momentum.dot(constraint_direc.row(bounce_idx - 1)) /
                        constraint_row_norm_sq(bounce_idx - 1) *
                        constraint_direc.row(bounce_idx - 1).transpose();
}

//' Compute when the next bounce occurs and which constraint it occurs on.
//'
//' @param position starting position
//' @param momentum starting momentum
//' @param constraint_direc F matrix (k-by-d matrix where k is the number of
//' linear constraints)
//' @param constraint_bound g vector (k dimensional)
//' @return pair of new (time until bounce, constraint index corresponding to
//' bounce)
std::pair<double, int> computeNextBounce(
    const Eigen::VectorXd position, const Eigen::VectorXd momentum,
    const Eigen::Map<Eigen::MatrixXd> constraint_direc,
    const Eigen::Map<Eigen::VectorXd> constraint_bound) {
  Eigen::ArrayXd fa = (constraint_direc * momentum).array();
  Eigen::ArrayXd fb = (constraint_direc * position).array();
  Eigen::ArrayXd U = (fa.square() + fb.square()).sqrt();
  // Eigen doesn't have an atan2 function
  Eigen::ArrayXd phi =
      -fa.binaryExpr(fb, [](double a, double b) { return std::atan2(a, b); });
  double min_time = std::numeric_limits<double>::infinity();
  int constraint_idx = -1;
  for (int i = 0; i < constraint_bound.size(); ++i) {
    if (U[i] > abs(constraint_bound[i])) {
      double bounce_time = -phi[i] + std::acos(-constraint_bound[i] / U[i]);
      if (bounce_time < min_time) {
        min_time = bounce_time;
        constraint_idx = i + 1;
      }
    }
  }
  return std::make_pair(min_time, constraint_idx);
}

//' Whiten a given position into the standard normal frame.
//'
//' @param position starting position
//' @param momentum starting momentum
//' @param constraint_direc F matrix (k-by-d matrix where k is the number of
//' linear constraints)
//' @param constraint_bound g vector (k dimensional)
//' @param cholesky_factor upper triangular matrix R from cholesky decomposition
//' of precision or covariance matrix into R^TR
//' @param unconstrained_mean mean of unconstrained Gaussian
//' @param prec_parametrized boolean for whether parametrization is by
//' precision (true)
//' or covariance matrix (false)
//' @return vector of position in standard normal frame
// [[Rcpp::export]]
Eigen::VectorXd whitenPosition(
    const Eigen::Map<Eigen::VectorXd> position,
    const Eigen::Map<Eigen::MatrixXd> constraint_direc,
    const Eigen::Map<Eigen::VectorXd> constraint_bound,
    const Eigen::Map<Eigen::MatrixXd> cholesky_factor,
    const Eigen::Map<Eigen::VectorXd> unconstrained_mean,
    bool prec_parametrized) {
  if (prec_parametrized) {
    return cholesky_factor * (position - unconstrained_mean);
  } else {
    return cholesky_factor.transpose().triangularView<Eigen::Lower>().solve(
        position - unconstrained_mean);
  }
}

//' Convert a position from standard normal frame back to original frame.
//'
//' @param position starting position
//' @param cholesky_factor upper triangular matrix R from cholesky decomposition
//' of precision or covariance matrix into R^TR
//' @param unconstrained_mean mean of unconstrained Gaussian
//' @param prec_parametrized boolean for whether parametrization is by
//' precision (true)
//' or covariance matrix (false)
//' @return vector of position in original frame
// [[Rcpp::export]]
Eigen::VectorXd unwhitenPosition(
    const Eigen::VectorXd position,
    const Eigen::Map<Eigen::MatrixXd> cholesky_factor,
    const Eigen::Map<Eigen::VectorXd> unconstrained_mean,
    bool prec_parametrized) {
  if (prec_parametrized) {
    return cholesky_factor.triangularView<Eigen::Upper>().solve(position) +
           unconstrained_mean;
  } else {
    return cholesky_factor.transpose() * position + unconstrained_mean;
  }
}

//' Simulate bouncing particle in whitened frame.
//'
//' @param initial_position starting position
//' @param initial_momentum starting momentum
//' @param constraint_direc F matrix (k-by-d matrix where k is the number of
//' linear constraints)
//' @param constraint_row_normsq vector of squared row norms of constraint_direc
//' @param constraint_bound g vector (k dimensional)
//' @param total_time total time the particle will bounce for
//' @param param diagnostic_mode boolean for whether to return the bounce
//' distances for each sample
//' @return vector of position in standard normal frame
// [[Rcpp::export]]
Rcpp::List simulateWhitenedDynamics(
    const Eigen::Map<Eigen::VectorXd> initial_position,
    const Eigen::Map<Eigen::VectorXd> initial_momentum,
    const Eigen::Map<Eigen::MatrixXd> constraint_direc,
    const Eigen::Map<Eigen::VectorXd> constraint_row_norm_sq,
    const Eigen::Map<Eigen::VectorXd> constraint_bound, double total_time,
    bool diagnostic_mode) {
  int bounce_constraint;
  double bounce_time;
  double bounced_distance;
  Eigen::VectorXd new_position;
  Eigen::VectorXd position = initial_position;
  Eigen::VectorXd momentum = initial_momentum;
  Eigen::VectorXd bounce_distances;
  if (diagnostic_mode) {
    bounce_distances = Eigen::VectorXd(constraint_direc.cols());
  }
  int num_bounces = 0;
  double travelled_time = 0;
  while (true) {
    std::tie(bounce_time, bounce_constraint) = computeNextBounce(
        position, momentum, constraint_direc, constraint_bound);
    if (bounce_time < total_time - travelled_time) {
      if (diagnostic_mode) {
        std::tie(new_position, momentum) =
            advanceWhitenedDynamics(position, momentum, bounce_time);
        bounced_distance = (new_position - position).norm();
        if (num_bounces > bounce_distances.size()) {
          bounce_distances.conservativeResize(2 * num_bounces);
        }
        bounce_distances(num_bounces) = bounced_distance;
        num_bounces++;
        position = new_position;
      } else {
        std::tie(position, momentum) =
            advanceWhitenedDynamics(position, momentum, bounce_time);
      }
      momentum = reflectMomentum(momentum, constraint_direc,
                                 constraint_row_norm_sq, bounce_constraint);
      travelled_time += bounce_time;
    } else {
      bounce_time = total_time - travelled_time;
      std::tie(position, momentum) =
          advanceWhitenedDynamics(position, momentum, bounce_time);
      return Rcpp::List::create(
          Rcpp::Named("sample") = position,
          Rcpp::Named("bounce_distances") = bounce_distances.head(num_bounces));
    }
  }
}
