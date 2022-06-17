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


// ' Create ZigZag engine object
// '
// ' @param dimension the dimension of MTN (d)
// ' @param a d-dimensional vector specifying the lower bounds.
// ' @param a d-dimensional vector specifying the upper bounds.
// ' @param flags 128
// ' @param info 1
// ' @param seed random seed
// ' @return a zigzag engine object.
// '
// [[Rcpp::export(createEngine)]]
Rcpp::List createEngine(int dimension,
                        std::vector<double> &lowerBounds,
                        std::vector<double> &upperBounds,
                        long flags, long info, long seed) {
    std::vector<double> mask(dimension, 1);
    auto zigZag = new ZigZagWrapper(
            zz::dispatch(dimension, mask.data(), lowerBounds.data(), upperBounds.data(),
                         flags, info, seed));

    XPtrZigZagWrapper engine(zigZag);

    Rcpp::List list = Rcpp::List::create(
            Rcpp::Named("engine") = engine,
            Rcpp::Named("dimension") = dimension,
            Rcpp::Named("mask") = mask,
            Rcpp::Named("flags") = flags,
            Rcpp::Named("info") = info
    );

    return list;
}

// ' Create ZigZag nuts engine object
// '
// ' Helper function creates zigZag nuts engine object with given latent dimension, location count and various
// ' implementation details. 
// '
// ' @param dimension the dimension of MTN (d).
// ' @param lowerBounds the d-dimensional lower bound.
// ' @param upperBounds the d-dimensional upper bound.
// ' @param flags 128.
// ' @param info 1.
// ' @param seed random seed.
// ' @param randomFlg set to TRUE.
// ' @param stepSize step size.
// ' @param mean mean vector.
// ' @param precision precision matrix.
// ' @return a zigzag-nuts engine object.
//'
// [[Rcpp::export(createNutsEngine)]]
Rcpp::List createNutsEngine(int dimension,
                            std::vector<double> &lowerBounds,
                            std::vector<double> &upperBounds,
                            long flags, long info, long seed,
                            bool randomFlg,
                            double stepSize,
                            NumericVector &mean,
                            NumericVector &precision) {
    std::vector<double> mask(dimension, 1); 
    auto zigZag = new ZigZagWrapper(
            zz::dispatch(dimension, mask.data(), lowerBounds.data(), upperBounds.data(), flags, info, seed));
    XPtrZigZagWrapper engineZZ(zigZag);

    // ptr to a zigzag obj
    auto ptr = parsePtrSse(engineZZ);
    ptr->setMean(zz::DblSpan(mean.begin(), mean.end()));
    ptr->setPrecision(zz::DblSpan(precision.begin(), precision.end()));

    // create a NUTS obj:
    auto nuts = new NutsWrapper(nuts::dispatchNuts(100, 10, seed, randomFlg, stepSize, ptr));
    XPtrNutsWrapper engineNuts(nuts);

    Rcpp::List list = Rcpp::List::create(Rcpp::Named("engine") = engineNuts);

    return list;
}

// ' Set mean for MTN
// '
// ' @param sexp pointer to zigzag object
// ' @param mean a numeric vector containing the MTN mean
// [[Rcpp::export(setMean)]]
void setMean(SEXP sexp, NumericVector &mean) {
    auto ptr = parsePtr(sexp);
    try {
        ptr->setMean(zz::DblSpan(mean.begin(), mean.end()));
    }

    catch (Rcpp::internal::InterruptedException &e) {
        Rcout << "Caught an interrupt!" << std::endl;
    }
}

// ' Set the precision matrix for the target MTN
// '
// ' @param sexp pointer to zigzag object
// ' @param precision the MTN precision matrix
// [[Rcpp::export(setPrecision)]]
void setPrecision(SEXP sexp, NumericVector &precision) {
    auto ptr = parsePtr(sexp);
    try {
        ptr->setPrecision(zz::DblSpan(precision.begin(), precision.end()));
    }

    catch (Rcpp::internal::InterruptedException &e) {
        Rcout << "Caught an interrupt!" << std::endl;
    }
}

// [[Rcpp::export(.doSomething)]]
void doSomething(SEXP sexp,
                 std::vector<double> &data) {
    auto ptr = parsePtr(sexp);
    //ptr->doSomething(data.data(), data.size());
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

// [[Rcpp::export(.oneNutsIteration)]]
Rcpp::List oneNutsIteration(SEXP sexp,
                            NumericVector &position,
                            NumericVector &momentum) {
    auto ptrNuts = parsePtrNuts(sexp);

    auto returnValue = ptrNuts->takeOneStep(zz::DblSpan(position.begin(), position.end()),
                                            zz::DblSpan(momentum.begin(), momentum.end()));
    Rcpp::List list = Rcpp::List::create(Rcpp::Named("position") = returnValue);
    return list;
}
