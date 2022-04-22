//
// Created by Marc Suchard on 2019-12-09.
//

#ifndef ZIG_ZAG_ABSTRACT_ZIG_ZAG_HPP
#define ZIG_ZAG_ABSTRACT_ZIG_ZAG_HPP

//#define TCB_SPAN_NAMESPACE_NAME std
//#define TCB_SPAN_NO_CONTRACT_CHECKING
#include "span.hpp"
#include "PrecisionColumn.hpp"

//namespace std {
//    template <typename T>
//    using span = tcb::span<T>;
//}

namespace zz {

    using DblSpan = tcb::span<double>;

    enum Flags {
        DOUBLE = 1u << 1u,
        FLOAT = 1u << 2u,
        TBB = 1u << 3u,
        OPENCL = 1u << 4u,
        SSE = 1u << 7u,
        AVX = 1u << 8u,
        AVX512 = 1u << 9u
    };

//    struct CpuAccumulate { };

#ifdef USE_TBB
    struct TbbAccumulate{ };
#endif

    class MinTravelInfo;

    class AbstractZigZag {
    public:

        AbstractZigZag() = default;

        virtual ~AbstractZigZag() = default;

        virtual double operate(DblSpan initialPosition,
                               DblSpan initialVelocity,
                               DblSpan initialAction,
                               DblSpan initialGradient,
                               DblSpan initialMomentum,
                               double time) = 0;

        virtual double operate(DblSpan initialPosition,
                               DblSpan initialMomentum,
                               double time) = 0;

        virtual void setMean(DblSpan meanVec) = 0;

        virtual void setPrecision(DblSpan precisionMat) = 0;

        virtual MinTravelInfo getNextBounce(DblSpan position,
                                            DblSpan velocity,
                                            DblSpan action,
                                            DblSpan gradient,
                                            DblSpan momentum) = 0;

//        virtual MinTravelInfo getNextBounceIrreversible(DblSpan position,
//                                                        DblSpan velocity,
//                                                        DblSpan action,
//                                                        DblSpan gradient) = 0;

        virtual void innerBounce(DblSpan position,
                                 DblSpan velocity,
                                 DblSpan action,
                                 DblSpan gradient,
                                 DblSpan momentum,
                                 double time, int index, int type) = 0;

        virtual void updateDynamics(DblSpan position,
                                    DblSpan velocity,
                                    DblSpan action,
                                    DblSpan gradient,
                                    DblSpan momentum,
                                    DblSpan column,
                                    double time, int index) = 0;
    };

    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args&&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }
}

#endif //ZIG_ZAG_ABSTRACT_ZIG_ZAG_HPP
