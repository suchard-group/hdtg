//
// Created by Marc Suchard on 2019-12-09.
//

#ifndef ZIG_ZAG_ABSTRACTZIGZAG_HPP
#define ZIG_ZAG_ABSTRACTZIGZAG_HPP

#define TCB_SPAN_NAMESPACE_NAME std
#define TCB_SPAN_NO_CONTRACT_CHECKING
#include "span.hpp"

namespace zz {

    enum Flags {
        FLOAT = 1 << 2,
        TBB = 1 << 3,
        OPENCL = 1 << 4,
        SSE = 1 << 7,
        AVX = 1 << 8,
        AVX512 = 1 << 9
    };

    struct CpuAccumulate { };

#ifdef USE_TBB
    struct TbbAccumulate{ };
#endif

    class MinTravelInfo;

    class AbstractZigZag {
    public:

        AbstractZigZag() = default;

        virtual ~AbstractZigZag() = default;

        virtual double operate(std::span<double> initialPosition,
                               std::span<double> initialVelocity,
                               std::span<double> initialAction,
                               std::span<double> initialGradient,
                               std::span<double> initialMomentum,
                               double time) = 0;

        virtual MinTravelInfo getNextBounce(std::span<double> position,
                                            std::span<double> velocity,
                                            std::span<double> action,
                                            std::span<double> gradient,
                                            std::span<double> momentum) = 0;
    };

    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args&&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }
}

#endif //ZIG_ZAG_ABSTRACTZIGZAG_HPP
