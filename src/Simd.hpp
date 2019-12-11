//
// Created by Marc Suchard on 2019-12-08.
//

#ifndef ZIG_ZAG_SIMD_HPP
#define ZIG_ZAG_SIMD_HPP


#include "xsimd/xsimd.hpp"
#include "ZigZag.hpp"

namespace zz {

    enum BounceType : int {
        NONE = 0,
        BOUNDARY,
        GRADIENT
    };

    class BounceInfo {
    public:

        BounceInfo() : type(BounceType::NONE),
                       index(-1),
                       time(std::numeric_limits<double>::infinity()) { }

        BounceInfo(int type, int index, double time) :
                type(type), index(index), time(time) { }

        ~BounceInfo() = default;

        bool isTimeRemaining() {
            return time > 0.0;
        }

        int getTypeInt() {
            if (type == BounceType::NONE) {
                return 0;
            } else if (type == BounceType::BOUNDARY) {
                return 1;
            } else if (type == BounceType::GRADIENT) {
                return 2;
            } else {
                exit(-1);
            }
        }

        int type;
        int index;
        double time;
    };

    class BounceState : public BounceInfo {
    public:

        BounceState() : BounceInfo(BounceType::NONE, -1, 0.0) { }

        explicit BounceState(double remainingTime) :
                BounceInfo(BounceType::NONE, -1, remainingTime) { }

        BounceState(int type, int index, double time) : BounceInfo(type, index, time) { }

        friend std::ostream& operator<<(std::ostream& os, const BounceState& rhs) {
            os << "remainingTime : " << rhs.time <<
            " lastBounceType: " << rhs.type << " in dim: " << rhs.index;
            return os;
        }
    };

    class MinTravelInfo : public BounceInfo {
    public:
        MinTravelInfo() : BounceInfo() { }

        MinTravelInfo(int minimumIndex, double minimumTime) :
                BounceInfo(BounceType::NONE, minimumIndex, minimumTime) { }

        MinTravelInfo(int type, int minimumIndex, double minimumTime) :
                BounceInfo(type, minimumIndex, minimumTime) { }

        friend std::ostream& operator<<(std::ostream& os, const MinTravelInfo& rhs) {
            os << "time = " << rhs.time << " @ " << rhs.index;
            return os;
        }

    };

    using D2 = xsimd::batch<double, 2>;
    using D2Bool = xsimd::batch_bool<double, 2>;
    using D2Index = xsimd::batch<int64_t, 2>;

//#define MERGE
#ifdef MERGE

    union pack_type {
        double d;
        int i[2];
    };

    constexpr double makePack() {
        pack_type pack;
        pack.i[0] = -1;
        pack.i[1] = BounceType::NONE;
        return pack.d;
    }

    class DoubleSseMinTravelInfo {
    public:
        DoubleSseMinTravelInfo() : time(std::numeric_limits<double>::infinity()) { }

        D2 time;
        D2 pack;
    };

#else

    class DoubleSseMinTravelInfo {
    public:
        DoubleSseMinTravelInfo() : type(BounceType::NONE),
                                   index(-1),
                                   time(std::numeric_limits<double>::infinity()) { }

        ~DoubleSseMinTravelInfo() = default;

        D2Index type;
        D2Index index;
        D2 time;
    };
#endif

    struct DoubleNoSimdTypeInfo {
        using BaseType = double;
        using SimdType = double;
        using BoolType = bool;
        using InfoType = MinTravelInfo;
        using IndexType = int;
        static const int SimdSize = 1;
    };

    struct DoubleSseTypeInfo {
        using BaseType = double;
        using SimdType = xsimd::batch<double, 2>;
        using BoolType = xsimd::batch_bool<double, 2>;
        using InfoType = DoubleSseMinTravelInfo;
        using IndexType = xsimd::batch<int64_t, 2>;
        static const int SimdSize = 2;
    };


//    template <typename T, size_t N>
//    using SimdBatch = xsimd::batch<T, N>;
//
//    template <typename T, size_t N>
//    using SimdBatchBool = xsimd::batch_bool<T, N>;

    template<typename SimdType, typename RealType>
    class SimdHelper {
    public:
        static inline SimdType get(const RealType *iterator);

        static inline void put(SimdType x, RealType *iterator);

        static inline SimdType infinity() noexcept; // TODO Make constexpr?

//        static const SimdType infinity = std::numeric_limits<double>::infinity();
    };

    template<>
    inline D2 SimdHelper<D2, D2::value_type>::get(const double *iterator) {
        return D2(iterator, xsimd::unaligned_mode());
//        return D2(iterator, xsimd::aligned_mode()); // TODO Check nearer end of development
    }

    template<>
    inline void SimdHelper<D2, D2::value_type>::put(D2 x, double *iterator) {
        x.store_aligned(iterator);
    }

    template<>
    inline D2 SimdHelper<D2, D2::value_type>::infinity() noexcept {
        return D2(std::numeric_limits<D2::value_type>::infinity(),
                  std::numeric_limits<D2::value_type>::infinity());
    }

    template<>
    inline double SimdHelper<double, double>::get(const double *iterator) {
        return *iterator;
    }

    template<>
    inline void SimdHelper<double, double>::put(double x, double *iterator) {
        *iterator = x;
    }

    template<>
    inline double SimdHelper<double, double>::infinity() noexcept {
        return std::numeric_limits<D2::value_type>::infinity();
    }

    template <typename T, typename B>
    inline T select(const B test, const T lhs, const T rhs);

    template<>
    inline double select(const bool test, const double lhs, const double rhs) {
        return test ? lhs : rhs;
    }

    template<>
    inline int select(const bool test, const int lhs, const int rhs) {
        return test ? lhs : rhs;
    }

    template<>
    inline D2 select(const D2Bool test, const D2 lhs, const D2 rhs) {
        return xsimd::select(test, lhs, rhs);
    }

    template<typename B>
    inline D2Index select(const B test, const D2Index lhs, const D2Index rhs) {
        return xsimd::select(test, lhs, rhs);
    }

    template <typename T>
    inline T infinity() noexcept;

    template<>
    inline double infinity() noexcept {
        return std::numeric_limits<double>::infinity();
    }

    template<>
    inline D2 infinity() noexcept {
        return D2(std::numeric_limits<double>::infinity());
    }

    template <typename T>
    inline T makeSimdIndex(int index) noexcept;

    template <>
    inline int makeSimdIndex(int index) noexcept {
        return index;
    }

    template <>
    inline D2Index makeSimdIndex(int index) noexcept {
        return D2Index(index, index + 1);
    }
}

#endif //ZIG_ZAG_SIMD_HPP
