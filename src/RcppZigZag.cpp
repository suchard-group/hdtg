#include <unordered_map>

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppParallel,RcppXsimd)]]
#include <RcppParallel.h>

#include "ZigZag.h"
#include "NoUTurn.h"

using ZigZagSharedPtr = std::shared_ptr<zz::AbstractZigZag>;

class ZigZagWrapper {
private:
    ZigZagSharedPtr zigZag;

public:
    ZigZagWrapper(ZigZagSharedPtr zigZag) : zigZag(zigZag) {}

    ZigZagSharedPtr &get() {
        return zigZag;
    }
};

using XPtrZigZagWrapper = Rcpp::XPtr<ZigZagWrapper>;

ZigZagSharedPtr &parsePtr(SEXP sexp) {
    XPtrZigZagWrapper ptr(sexp);
    if (!ptr) {
        Rcpp::stop("External pointer is uninitialized");
    }
    return ptr->get();
}
// copied XPtr for sse. todo: make a template
using SseSharedPtr = std::shared_ptr<zz::ZigZag<zz::DoubleSseTypeInfo>>;

class SseWrapper {
private:
    SseSharedPtr Sse;

public:
    SseWrapper(SseSharedPtr Sse) : Sse(Sse) {}

    SseSharedPtr &get() {
        return Sse;
    }
};

using XPtrSseWrapper = Rcpp::XPtr<SseWrapper>;

SseSharedPtr &parsePtrSse(SEXP sexp) {
    XPtrSseWrapper ptr(sexp);
    if (!ptr) {
        Rcpp::stop("External pointer is uninitialized");
    }
    return ptr->get();
}

// copied XPtr for Nuts. todo: make a template

using NutsSharedPtr = std::shared_ptr<nuts::NoUTurn>;

class NutsWrapper {
private:
    NutsSharedPtr Nuts;

public:
    NutsWrapper(NutsSharedPtr Nuts) : Nuts(Nuts) {}

    NutsSharedPtr &get() {
        return Nuts;
    }
};

using XPtrNutsWrapper = Rcpp::XPtr<NutsWrapper>;

NutsSharedPtr &parsePtrNuts(SEXP sexp) {
    XPtrNutsWrapper ptr(sexp);
    if (!ptr) {
        Rcpp::stop("External pointer is uninitialized");
    }
    return ptr->get();
}

//' Create a Zigzag-HMC engine object
//' 
//' Create the C++ object to set up SIMD vectorization for speeding up calculations
//' for Zigzag-HMC ("Zigzag-HMC engine"). 
//'
//' @param dimension the dimension of MTN.
//' @param lowerBounds a vector specifying the lower bounds.
//' @param upperBounds a vector specifying the upper bounds.
//' @param seed random seed.
//' @param mean the mean vector.
//' @param precision the precision matrix.
//' @param flags which SIMD instruction set to use. 128 = SSE, 256 = AVX.
//' @return a list whose only element is the Zigzag-HMC engine object.
//' @export
// [[Rcpp::export(createEngine)]]
Rcpp::List createEngine(int dimension,
                        std::vector<double> &lowerBounds,
                        std::vector<double> &upperBounds,
                        long seed, 
                        NumericVector &mean,
                        NumericVector &precision,
                        long flags = 128) {
    std::vector<double> mask(dimension, 1);
    auto zigZag = new ZigZagWrapper(
      zz::dispatch(dimension, mask.data(), lowerBounds.data(), upperBounds.data(),
                         flags, 1, seed));

    XPtrZigZagWrapper engine(zigZag);

    auto ptr = parsePtrSse(engine);
    ptr->setMean(zz::DblSpan(mean.begin(), mean.end()));
    ptr->setPrecision(zz::DblSpan(precision.begin(), precision.end()));

    Rcpp::List list = Rcpp::List::create(Rcpp::Named("ptr") = engine);

    return list;
}

//' Create a Zigzag-NUTS engine object
//' 
//' Create the C++ object to set up SIMD vectorization for speeding up calculations
//' for Zigzag-NUTS ("Zigzag-NUTS engine"). 
//'
//' @param dimension the dimension of MTN.
//' @param lowerBounds a vector specifying the lower bounds.
//' @param upperBounds a vector specifying the upper bounds.
//' @param seed random seed.
//' @param stepSize the base step size for Zigzag-NUTS.
//' @param mean the mean vector.
//' @param precision the precision matrix.
//' @param flags which SIMD instruction set to use. 128 = SSE, 256 = AVX.
//' @return a list whose only element is the Zigzag-NUTS engine object.
//' @export
// [[Rcpp::export(createNutsEngine)]]
Rcpp::List createNutsEngine(int dimension,
                            std::vector<double> &lowerBounds,
                            std::vector<double> &upperBounds,
                            long seed,
                            double stepSize,
                            NumericVector &mean,
                            NumericVector &precision,
                            long flags = 128) {
    std::vector<double> mask(dimension, 1); 
    auto zigZag = new ZigZagWrapper(
            zz::dispatch(dimension, mask.data(), lowerBounds.data(), upperBounds.data(), flags, 1, seed));
    XPtrZigZagWrapper engineZZ(zigZag);

    // ptr to a zigzag obj
    auto ptr = parsePtrSse(engineZZ);
    ptr->setMean(zz::DblSpan(mean.begin(), mean.end()));
    ptr->setPrecision(zz::DblSpan(precision.begin(), precision.end()));

    // create a NUTS obj:
    auto nuts = new NutsWrapper(nuts::dispatchNuts(100, 10, seed, stepSize, ptr));
    XPtrNutsWrapper engineNuts(nuts);

    Rcpp::List list = Rcpp::List::create(Rcpp::Named("ptr") = engineNuts);

    return list;
}


//' Set the mean for the target MTN
//'
//' Set the mean vector for a given Zigzag-HMC engine object.
//'
//' @param sexp pointer to a Zigzag-HMC engine object.
//' @param mean the mean vector.
//' @export
// [[Rcpp::export(setMean)]]
void setMean(SEXP sexp, NumericVector &mean) {
    auto ptr = parsePtr(sexp);
    ptr->setMean(zz::DblSpan(mean.begin(), mean.end()));
}

//' Set the precision matrix for the target MTN
//' 
//' Set the precision matrix for a given Zigzag-HMC engine object.
//'
//' @param sexp pointer to a Zigzag-HMC engine object.
//' @param precision the precision matrix.
//' @export
// [[Rcpp::export(setPrecision)]]
void setPrecision(SEXP sexp, NumericVector &precision) {
    auto ptr = parsePtr(sexp);
    ptr->setPrecision(zz::DblSpan(precision.begin(), precision.end()));
}

// [[Rcpp::export(getNextEvent)]]
Rcpp::List getNextEvent(SEXP sexp,
                        NumericVector &position,
                        NumericVector &velocity,
                        NumericVector &action,
                        NumericVector &logpdfGradient,
                        NumericVector &momentum) {

    auto ptr = parsePtr(sexp);
    // Rcout << "action";
    // Rf_PrintValue(action);
    // Rcout << "\n";
    //
    // Rcout << "logpdfGradient";
    // Rf_PrintValue(logpdfGradient);
    // Rcout << "\n";
    //
    // Rcout << "momentum";
    // Rf_PrintValue(momentum);
    // Rcout << "\n";
    auto firstBounce = ptr->getNextBounce(
            zz::DblSpan(position.begin(), position.end()),
            zz::DblSpan(velocity.begin(), velocity.end()),
            zz::DblSpan(action.begin(), action.end()),
            zz::DblSpan(logpdfGradient.begin(), logpdfGradient.end()),
            zz::DblSpan(momentum.begin(), momentum.end()));

    Rcpp::List list = Rcpp::List::create(
            Rcpp::Named("type") = firstBounce.type,
            Rcpp::Named("index") = firstBounce.index,
            Rcpp::Named("time") = firstBounce.time);

    return list;
}


class RCallback : public zz::PrecisionColumnCallback {
public:
    RCallback(Function callback) : PrecisionColumnCallback(), callback(callback), column(-1) {}

    ~RCallback() { releaseColumn(); }

    const double *getColumn(int index) override {

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

// [[Rcpp::export(.oneIteration)]]
Rcpp::List oneIteration(SEXP sexp,
                        NumericVector &position,
                        NumericVector &momentum,
                        double time) {
    auto ptr = parsePtr(sexp);

    auto returnValue = ptr->operate(
            zz::DblSpan(position.begin(), position.end()),
            zz::DblSpan(momentum.begin(), momentum.end()),
            time
    );
    Rcpp::List list = Rcpp::List::create(
            Rcpp::Named("returnValue") = returnValue,
            Rcpp::Named("position") = position);

    return list;
}

// [[Rcpp::export(.oneIrreversibleIteration)]]
Rcpp::List oneIrreversibleIteration(SEXP sexp,
                                    NumericVector &position,
                                    NumericVector &velocity,
                                    double time) {
    auto ptr = parsePtr(sexp);

    auto returnValue = ptr->operateIrreversible(
            zz::DblSpan(position.begin(), position.end()),
            zz::DblSpan(velocity.begin(), velocity.end()),
            time
    );
    Rcpp::List list = Rcpp::List::create(
            Rcpp::Named("returnValue") = returnValue,
            Rcpp::Named("position") = position,
            Rcpp::Named("velocity") = velocity);

    return list;
}

// [[Rcpp::export(.oneNutsIteration)]]
Rcpp::List oneNutsIteration(SEXP sexp,
                            NumericVector &position,
                            NumericVector &momentum) {
    auto ptrNuts = parsePtrNuts(sexp);

    auto returnValue = ptrNuts->generateNextState(zz::DblSpan(position.begin(), position.end()),
                                            zz::DblSpan(momentum.begin(), momentum.end()));
    Rcpp::List list = Rcpp::List::create(Rcpp::Named("position") = returnValue);
    return list;
}
