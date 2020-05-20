
#include <unordered_map>

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel,RcppXsimd)]]
#include <RcppParallel.h>

#include "ZigZag.hpp"

//' @export
// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}

using ZigZagSharedPtr = std::shared_ptr<zz::AbstractZigZag>;

class ZigZagWrapper {
private:
  ZigZagSharedPtr zigZag;
  
public:
  ZigZagWrapper(ZigZagSharedPtr zigZag) : zigZag(zigZag) { }
  
  ZigZagSharedPtr& get() {
    return zigZag;
  }
};

using XPtrZigZagWrapper = Rcpp::XPtr<ZigZagWrapper>;

ZigZagSharedPtr& parsePtr(SEXP sexp) {
  XPtrZigZagWrapper ptr(sexp);
  if (!ptr) {
    Rcpp::stop("External pointer is uninitialized");
  }
  return ptr->get();
}

//' Create ZigZag engine object
//'
//' Helper function creates zigZag engine object with given latent dimension, location count and various
//' implementation details. Called by \code{MassivezigZag::engineInitial()}.
//'
//' @param embeddingDimension Dimension of latent locations.
//' @param locationCount Number of locations and size of distance matrix.
//' @param tbb Number of CPU cores to be used.
//' @param simd For CPU implementation: no SIMD (\code{0}), SSE (\code{1}) or AVX (\code{2}).
//' @param truncation Likelihood includes truncation term? Defaults to \code{TRUE}.
//' @param gpu Which GPU to use? If only 1 available, use \code{gpu=1}. Defaults to \code{0}, no GPU.
//' @param single Set \code{single=1} if your GPU does not accommodate doubles.
//' @return zigZag engine object.
//'
//' @export
// [[Rcpp::export(createEngine)]]
Rcpp::List createEngine(int dimension, 
                        std::vector<double>& mask, 
                        std::vector<double>& observed, 
                        long flags, long info, long seed) {
  
  auto zigZag = new ZigZagWrapper(
    zz::dispatch(dimension, mask.data(), observed.data(), flags, info, seed));
  
  XPtrZigZagWrapper engine(zigZag);
  
  Rcpp::List list = Rcpp::List::create(
    Rcpp::Named("engine") = engine,
    Rcpp::Named("dimension") = dimension,
    Rcpp::Named("mask") = mask,
    Rcpp::Named("observed") = observed,
 //   Rcpp::Named("dataInitialzied") = false,
//    Rcpp::Named("locationsInitialized") = false,
 //   Rcpp::Named("threads") = threads,
 //   Rcpp::Named("deviceNumber") = deviceNumber,
    Rcpp::Named("flags") = flags,
    Rcpp::Named("info") = info
  );
  
  return list;
}

// [[Rcpp::export(.doSomething)]]
void doSomething(SEXP sexp,
                     std::vector<double>& data) {
  auto ptr = parsePtr(sexp);
  //ptr->doSomething(data.data(), data.size());
}

// [[Rcpp::export(.getNextEvent)]]
Rcpp::List getNextEvent(SEXP sexp, 
                        NumericVector& position,
                        NumericVector& velocity,
                        NumericVector& action,
                        NumericVector& gradient,
                        NumericVector& momentum) {
  
  auto ptr = parsePtr(sexp);
  auto firstBounce =  ptr->getNextBounce(
    zz::DblSpan(position.begin(), position.end()),
    zz::DblSpan(velocity.begin(), velocity.end()),
    zz::DblSpan(action.begin(), action.end()),
    zz::DblSpan(gradient.begin(), gradient.end()),
    zz::DblSpan(momentum.begin(), momentum.end()));
  
  Rcpp::List list = Rcpp::List::create(
    Rcpp::Named("type") = firstBounce.type,
    Rcpp::Named("index") = firstBounce.index,
    Rcpp::Named("time") = firstBounce.time);
  
  return list;  
}


class RCallback : public zz::PrecisionColumnCallback {
public:
  RCallback(Function callback) : PrecisionColumnCallback(), callback(callback), column(-1) { }
  
  ~RCallback() { releaseColumn(); }
  
  const double* getColumn(int index) override {
    
    auto it = map.find(index);
    
    if (it == map.end()) {
      auto newElement = map.emplace(index, callback(index));
      return newElement.first->second.begin();
    } else {
      return it->second.begin();
    }
  }
  
  void releaseColumn() override {
    column = -1;
  }
  
private:
  Function callback;
  int column;
  std::unordered_map<int, NumericVector> map;
  NumericVector result;
};

// [[Rcpp::export(.operate)]]
Rcpp::List operate(SEXP sexp,
                   Function rCallback,
                   NumericVector& position,
                   NumericVector& velocity,
                   NumericVector& action,
                   NumericVector& gradient,
                   NumericVector& momentum,
                   double time) {

  RCallback callback(rCallback);

  auto ptr = parsePtr(sexp);
  auto returnValue =  ptr->operate(
    zz::DblSpan(position.begin(), position.end()),
    zz::DblSpan(velocity.begin(), velocity.end()),
    zz::DblSpan(action.begin(), action.end()),
    zz::DblSpan(gradient.begin(), gradient.end()),
    zz::DblSpan(momentum.begin(), momentum.end()),
    time, callback);

  Rcpp::List list = Rcpp::List::create(
    Rcpp::Named("returnValue") = returnValue,
    Rcpp::Named("position") = position);

  return list;
}