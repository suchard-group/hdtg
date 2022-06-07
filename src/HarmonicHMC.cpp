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

// ' Solve XA=B for two matrices A and B.
// '

// '
// ' @param A matrix
// ' @param B matrix
// ' @return X matrix
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
//' @param constraintDirec F matrix (k-by-d matrix where k is the number of
//' linear constraints)
//' @param constraintBound g vector (k dimensional)
//' @param choleskyFactor upper triangular matrix R from cholesky decomposition
//'  of precision or covariance matrix into R^TR
//' @param unconstrainedMean mean of unconstrained Gaussian
//' @param precParametrized boolean for whether parametrization is by precision
//'  (true) or covariance matrix (false)
//' @return List of new constraint directions, the squared row norms of those
//' constraints (for computational efficiency later), and new bounds
//' @export
// [[Rcpp::export]]
Rcpp::List applyWhitenTransform(
    const Eigen::Map<Eigen::MatrixXd> constraintDirec,
    const Eigen::Map<Eigen::VectorXd> constraintBound,
    const Eigen::Map<Eigen::MatrixXd> choleskyFactor,
    const Eigen::Map<Eigen::VectorXd> unconstrainedMean,
    bool precParametrized) {
  Eigen::ArrayXXd direc;
  if (precParametrized) {
    direc = solveFromRight(choleskyFactor, constraintDirec).array();
  } else {
    direc = constraintDirec * choleskyFactor.transpose();
  }
  return Rcpp::List::create(
      Rcpp::_["direc"] = direc,
      Rcpp::_["direcRowNormSq"] = direc.square().rowwise().sum(),
      Rcpp::_["bound"] = constraintBound + constraintDirec * unconstrainedMean);
}

// ' Compute Hamiltonian dynamics after specified time.
// '
// ' @param position starting position
// ' @param momentum starting momentum
// ' @param time amount of time the system is run for
// ' @return pair of new (position, momentum)
std::pair<Eigen::VectorXd, Eigen::VectorXd> advanceWhitenedDynamics(
    const Eigen::VectorXd position, const Eigen::VectorXd momentum,
    const double time) {
  return std::make_pair(momentum * sin(time) + position * cos(time),
                        momentum * cos(time) - position * sin(time));
}

// ' Reflect momentum off of a constraint boundary.
// '
// ' Given a constraint boundary, calculate the momentum as if that boundary
// ' was a wall and there is an elastic collision, and the angle of incidence
// ' equals the angle of reflection.
// '
// ' @param momentum starting momentum
// ' @param constraintDirec F matrix (k-by-d matrix where k is the number of
// ' linear constraints)
// ' @param constraintRowNormSq vector of squared row norms of constraintDirec
// ' @param bounceIdx integer index of which constraint is being bounced off of
// ' @param time amount of time the system is run for
// ' @return momentum after bouncing
Eigen::VectorXd reflectMomentum(
    const Eigen::VectorXd momentum,
    const Eigen::Map<Eigen::MatrixXd> constraintDirec,
    const Eigen::Map<Eigen::VectorXd> constraintRowNormSq,
    const int bounceIdx) {
  return momentum - 2 * momentum.dot(constraintDirec.row(bounceIdx - 1)) /
                        constraintRowNormSq(bounceIdx - 1) *
                        constraintDirec.row(bounceIdx - 1).transpose();
}

//  ' Compute when the next bounce occurs and which constraint it occurs on.
// '
// ' @param position starting position
// ' @param momentum starting momentum
// ' @param constraintDirec F matrix (k-by-d matrix where k is the number of
// ' linear constraints)
// ' @param constraintBound g vector (k dimensional)
// ' @return pair of new (time until bounce, constraint index corresponding to
// ' bounce)
std::pair<double, int> computeNextBounce(
    const Eigen::VectorXd position, const Eigen::VectorXd momentum,
    const Eigen::Map<Eigen::MatrixXd> constraintDirec,
    const Eigen::Map<Eigen::VectorXd> constraintBound) {
  Eigen::ArrayXd fa = (constraintDirec * momentum).array();
  Eigen::ArrayXd fb = (constraintDirec * position).array();
  Eigen::ArrayXd U = (fa.square() + fb.square()).sqrt();
  // Eigen doesn't have an atan2 function
  Eigen::ArrayXd phi =
      -fa.binaryExpr(fb, [](double a, double b) { return std::atan2(a, b); });
  double minTime = std::numeric_limits<double>::infinity();
  int constraintIdx = -1;
  for (int i = 0; i < constraintBound.size(); ++i) {
    if (U[i] > abs(constraintBound[i])) {
      double bounceTime = -phi[i] + std::acos(-constraintBound[i] / U[i]);
      if (bounceTime < minTime) {
        minTime = bounceTime;
        constraintIdx = i + 1;
      }
    }
  }
  return std::make_pair(minTime, constraintIdx);
}

//' Whiten a given position into the standard normal frame.
//'
//' @param position starting position
//' @param constraintDirec F matrix (k-by-d matrix where k is the number of
//' linear constraints)
//' @param constraintBound g vector (k dimensional)
//' @param choleskyFactor upper triangular matrix R from cholesky decomposition
//' of precision or covariance matrix into R^TR
//' @param unconstrainedMean mean of unconstrained Gaussian
//' @param precParametrized boolean for whether parametrization is by
//' precision (true)
//' or covariance matrix (false)
//' @return vector of position in standard normal frame
// [[Rcpp::export]]
Eigen::VectorXd whitenPosition(
    const Eigen::Map<Eigen::VectorXd> position,
    const Eigen::Map<Eigen::MatrixXd> constraintDirec,
    const Eigen::Map<Eigen::VectorXd> constraintBound,
    const Eigen::Map<Eigen::MatrixXd> choleskyFactor,
    const Eigen::Map<Eigen::VectorXd> unconstrainedMean,
    bool precParametrized) {
  if (precParametrized) {
    return choleskyFactor * (position - unconstrainedMean);
  } else {
    return choleskyFactor.transpose().triangularView<Eigen::Lower>().solve(
        position - unconstrainedMean);
  }
}

//' Convert a position from standard normal frame back to original frame.
//'
//' @param position starting position
//' @param choleskyFactor upper triangular matrix R from cholesky decomposition
//' of precision or covariance matrix into R^TR
//' @param unconstrainedMean mean of unconstrained Gaussian
//' @param precParametrized boolean for whether parametrization is by
//' precision (true)
//' or covariance matrix (false)
//' @return vector of position in original frame
// [[Rcpp::export]]
Eigen::VectorXd unwhitenPosition(
    const Eigen::VectorXd position,
    const Eigen::Map<Eigen::MatrixXd> choleskyFactor,
    const Eigen::Map<Eigen::VectorXd> unconstrainedMean,
    bool precParametrized) {
  if (precParametrized) {
    return choleskyFactor.triangularView<Eigen::Upper>().solve(position) +
           unconstrainedMean;
  } else {
    return choleskyFactor.transpose() * position + unconstrainedMean;
  }
}

//' Simulate bouncing particle in whitened frame.
//'
//' @param initialPosition starting position
//' @param initialMomentum starting momentum
//' @param constraintDirec F matrix (k-by-d matrix where k is the number of
//' linear constraints)
//' @param constraintRowNormSq vector of squared row norms of constraintDirec
//' @param constraintBound g vector (k dimensional)
//' @param integrationTime total time the particle will travel for
//' @param diagnosticMode boolean for whether to return the bounce
//' distances for each sample
//' @return vector of position in standard normal frame
// [[Rcpp::export]]
Rcpp::List simulateWhitenedDynamics(
    const Eigen::Map<Eigen::VectorXd> initialPosition,
    const Eigen::Map<Eigen::VectorXd> initialMomentum,
    const Eigen::Map<Eigen::MatrixXd> constraintDirec,
    const Eigen::Map<Eigen::VectorXd> constraintRowNormSq,
    const Eigen::Map<Eigen::VectorXd> constraintBound, double integrationTime,
    bool diagnosticMode) {
  int bounceConstraint;
  double bounceTime;
  double bouncedDistance;
  Eigen::VectorXd newPosition;
  Eigen::VectorXd position = initialPosition;
  Eigen::VectorXd momentum = initialMomentum;
  Eigen::VectorXd bounceDistances;
  if (diagnosticMode) {
    bounceDistances = Eigen::VectorXd(constraintDirec.cols());
  }
  int numBounces = 0;
  double travelledTime = 0;
  while (true) {
    std::tie(bounceTime, bounceConstraint) =
        computeNextBounce(position, momentum, constraintDirec, constraintBound);
    if (bounceTime < integrationTime - travelledTime) {
      if (diagnosticMode) {
        std::tie(newPosition, momentum) =
            advanceWhitenedDynamics(position, momentum, bounceTime);
        bouncedDistance = (newPosition - position).norm();
        if (numBounces >= bounceDistances.size()) {
          bounceDistances.conservativeResize(2 * numBounces);
        }
        bounceDistances(numBounces) = bouncedDistance;
        numBounces++;
        position = newPosition;
      } else {
        std::tie(position, momentum) =
            advanceWhitenedDynamics(position, momentum, bounceTime);
      }
      momentum = reflectMomentum(momentum, constraintDirec, constraintRowNormSq,
                                 bounceConstraint);
      travelledTime += bounceTime;
    } else {
      bounceTime = integrationTime - travelledTime;
      std::tie(position, momentum) =
          advanceWhitenedDynamics(position, momentum, bounceTime);
      return Rcpp::List::create(
          Rcpp::Named("position") = position,
          Rcpp::Named("bounceDistances") = bounceDistances.head(numBounces));
    }
  }
}
