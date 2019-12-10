//
// Created by Marc Suchard on 2019-12-03.
//

#ifndef ZIG_ZAG_ZIGZAG_HPP
#define ZIG_ZAG_ZIGZAG_HPP

#include <vector>
#include <cmath>

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/global_control.h"

#define TIMING

#ifdef TIMING
#include <map>
#include <iomanip>
#include "Timing.h"
#endif // TIMING

#include "dr_evomodel_operators_NativeZigZag.h"
#include "MemoryManagement.hpp"
#include "Simd.hpp"
#include "AbstractZigZag.hpp"

namespace zz {

    template <typename TypeInfo>
    class ZigZag : public AbstractZigZag {
    public:
        using RealType = typename TypeInfo::BaseType;
        using SimdType = typename TypeInfo::SimdType;
        using InfoType = typename TypeInfo::InfoType;
        using IndexType = typename TypeInfo::IndexType;
        static const int SimdSize = TypeInfo::SimdSize;

        using MaskType = double;

        ZigZag(size_t dimension,
               double *rawMask,
               double *rawObserved,
               jmethodID id,
               long flags,
               int nThreads) : AbstractZigZag(),
                               dimension(dimension),
                               mask(constructMask(rawMask, dimension)),
                               observed(constructMask(rawObserved, dimension)),
                               providerMethodId(id),
                               mmPosition(dimension),
                               mmVelocity(dimension),
                               mmAction(dimension),
                               mmGradient(dimension),
                               mmMomentum(dimension),
                               flags(flags),
                               nThreads(nThreads) {
            std::cerr << "c'tor ZigZag" << std::endl;

            if (flags & zz::Flags::TBB) {
                if (nThreads <= 0) {
                    nThreads = tbb::task_scheduler_init::default_num_threads();
                }

                std::cout << "Using " << nThreads << " threads" << std::endl;

                control = std::make_shared<tbb::global_control>(tbb::global_control::max_allowed_parallelism, nThreads);
            }
        }

        virtual ~ZigZag() {

#ifdef TIMING
            std::cerr << std::endl;
            for (auto& d : duration) {
                std::cerr << d.first << " " << std::scientific <<
                static_cast<double>(d.second) * 0.001 << std::endl;
            }
#endif
            std::cerr << "d'tor ZigZag" << std::endl;
        };

        double operate(std::span<double> position,
                       std::span<double> velocity,
                       std::span<double> action,
                       std::span<double> gradient,
                       std::span<double> momentum,
                       double time) {

            BounceState bounceState(time);

            while (bounceState.isTimeRemaining()) {

                const auto firstBounce = getNextBounce(
                        position, velocity,
                        action, gradient, momentum);
                
                bounceState = doBounce(bounceState, firstBounce,
                                       position, velocity, action, gradient, momentum);
            }

            return 0.0;
        }

        template <typename T>
        struct Dynamics {

            template <typename V, typename W>
            Dynamics(V& position,
                     V& velocity,
                     V& action,
                     V& gradient,
                     V& momentum,
                     const W& observed) : position(position.data()),
                                          velocity(velocity.data()),
                                          action(action.data()),
                                          gradient(gradient.data()),
                                          momentum(momentum.data()),
                                          observed(observed.data()) { }
            ~Dynamics() = default;

            T* position;
            T* velocity;
            T* action;
            T* gradient;
            T* momentum;
            const T* observed;
        };

        MinTravelInfo getNextBounce(std::span<double> position,
                                    std::span<double> velocity,
                                    std::span<double> action,
                                    std::span<double> gradient,
                                    std::span<double> momentum) {

#if 0
            std::vector<double> b;

            auto buffer = [&b](std::span<double>& in, mm::MemoryManager<double>& out) {
                mm::bufferedCopy(std::begin(in), std::end(in), std::begin(out), b);
            };

            buffer(position, mmPosition);
            buffer(velocity, mmVelocity);
            buffer(action, mmAction);
            buffer(gradient, mmGradient);
            buffer(momentum, mmMomentum);

            return getNextBounceImpl(
                    Dynamics<double>(mmPosition, mmVelocity, mmAction, mmGradient, mmMomentum, observed));
#else
            return getNextBounceImpl(
                    Dynamics<double>(position, velocity, action, gradient, momentum, observed));
#endif
        }

        template <typename R>
        MinTravelInfo getNextBounceImpl(const Dynamics<R>& dynamics) {

#ifdef TIMING
            auto start = zz::chrono::steady_clock::now();
#endif

            auto task = [&](const size_t begin, const size_t end) -> MinTravelInfo {

                const auto length = end - begin;
                const auto vectorCount = length - length % SimdSize;

                MinTravelInfo travel = vectorized_transform<SimdType, SimdSize>(begin, begin + vectorCount, dynamics,
                                                                                InfoType());

                if (begin + vectorCount < length) { // Edge-case
                    travel = vectorized_transform<RealType, 1>(begin + vectorCount, length, dynamics, travel);
                }

                return travel;
            };

            MinTravelInfo travel = (nThreads <= 1) ?
                                   task(size_t(0), dimension) :
                                   parallel_task(size_t(0), dimension, MinTravelInfo(),
                                                 task,
                                                 [](MinTravelInfo lhs, MinTravelInfo rhs) {
                                                     return (lhs.time < rhs.time) ? lhs : rhs;
                                                 });

#ifdef TIMING
            auto end = zz::chrono::steady_clock::now();
            duration["getNextBounce"] += zz::chrono::duration_cast<chrono::TimingUnits>(end - start).count();
#endif

            return travel;
        }

        template <typename T, typename F, typename G>
        inline T parallel_task(const size_t begin, const size_t end, T sum, F transform, G reduce) {

            return tbb::parallel_reduce(
                    tbb::blocked_range<size_t>(begin, end, (end - begin) / nThreads),
                    sum,
                    [transform, reduce](const tbb::blocked_range<size_t>& r, T sum) { // TODO Test &transform, &reduce
                        return reduce(sum, transform(r.begin(), r.end()));
                    },
                    reduce
            );
        }

    protected:

    private:

        template <typename S, int SimdSize, typename R, typename I, typename Int>
        MinTravelInfo vectorized_transform(Int i, const Int end,
                                           const Dynamics<R>& dynamics, I result) {

            const auto *position = dynamics.position;
            const auto *velocity = dynamics.velocity;
            const auto *action = dynamics.action;
            const auto *gradient = dynamics.gradient;
            const auto *momentum = dynamics.momentum;
            const auto *observed = dynamics.observed;

            for ( ; i < end; i += SimdSize) {

                const auto boundaryTime = findBoundaryTime(
                        SimdHelper<S, R>::get(position + i),
                        SimdHelper<S, R>::get(velocity + i),
                        SimdHelper<S, R>::get(observed + i)
                );

                reduce_min(result, boundaryTime, i, BounceType::BOUNDARY); // TODO Try: result = reduce_min(result, ...)

                const auto gradientTime = minimumPositiveRoot(
                        -SimdHelper<S, R>::get(action + i) / 2,
                         SimdHelper<S, R>::get(gradient + i),
                         SimdHelper<S, R>::get(momentum + i)
                );

                reduce_min(result, gradientTime, i, BounceType::GRADIENT);
            }

            return horizontal_min(result);
        };


        template <typename Vector>
        BounceState doBounce(BounceState bounceState, MinTravelInfo bounceInformation,
                Vector position,
                Vector velocity,
                Vector action,
                Vector gradient,
                Vector momentum) {

            return BounceState(0.0);
        }

        static inline void reduce_min(MinTravelInfo& result,
                                      const double time, const int index, const int type) {
            if (time < result.time) {
                result.time = time;
                result.index = index;
                result.type = type;
            }
        }

        static inline MinTravelInfo horizontal_min(MinTravelInfo result) {
            return result;
        }

        static inline void reduce_min(DoubleSseMinTravelInfo& result,
                const D2 time, const int index, const int type) {
            const auto lessThan = time < result.time;
            if (xsimd::any(lessThan)) {
                result.time = select(lessThan, time, result.time);
                const auto mask = static_cast<__m128i>(lessThan);
                result.index = select(mask, makeSimdIndex<D2Index>(index), result.index); // TODO Merge into single register?
                result.type = select(mask, D2Index(type), result.type);
            }
        }

        static inline MinTravelInfo horizontal_min(DoubleSseMinTravelInfo vector) {
            return (vector.time[0] < vector.time[1]) ?
                MinTravelInfo(vector.type[0], vector.index[0], vector.time[0]) :
                MinTravelInfo(vector.type[1], vector.index[1], vector.time[1]);
        }

        template <typename T>
        static inline T findBoundaryTime(const T position,
                                         const T velocity,
                                         const T observed) {

            return select(headingTowardsBoundary(position, velocity, observed),
                          abs(position / velocity),
                          infinity<T>());
        }

        template <typename T>
        static inline auto headingTowardsBoundary(const T position,
                                                  const T velocity,
                                                  const T observed)
        -> decltype(T(1.0) > T(0.0)) {
            return observed * position * velocity < T(0.0);
        }

        template <typename T>
        static inline T minimumPositiveRoot(const T a, const T b, const T c) {

            const auto discriminant = b * b - 4 * a * c;
            const auto sqrtDiscriminant = select(c == T(0.0), b, sqrt(abs(discriminant)));

            auto root1 = (-b - sqrtDiscriminant) / (2 * a);
            auto root2 = (-b + sqrtDiscriminant) / (2 * a);

            root1 = select(root1 > T(0.0), root1, infinity<T>());
            root2 = select(root2 > T(0.0), root2, infinity<T>());

            const auto root = select(root1 < root2, root1, root2);
            return select(discriminant < T(0.0), infinity<T>(), root);
        }

        template <typename T>
        static inline T sign(const T x) {
            const auto zero = T(0.0);
            return select(x == zero,
                          x,
                          select(x < zero,
                                 T(-1.0),
                                 T(+1.0))
            );
        }

        static mm::MemoryManager<MaskType> constructMask(double *raw, size_t length) {

            mm::MemoryManager<MaskType> mask;
            mask.reserve(length);

            std::transform(raw, raw + length, std::back_inserter(mask),
                           [](double x) {
                               return (x == 1.0) ? MaskType(1.0) : MaskType(0.0);
                           });

            return mask;
        }

        static double doMask(MaskType mask, double x) {
            return mask * x;
        }

        size_t dimension;
        mm::MemoryManager<MaskType> mask;
        mm::MemoryManager<MaskType> observed;
        jmethodID providerMethodId;

        mm::MemoryManager<double> mmPosition;
        mm::MemoryManager<double> mmVelocity;
        mm::MemoryManager<double> mmAction;
        mm::MemoryManager<double> mmGradient;
        mm::MemoryManager<double> mmMomentum;

        long flags;
        int nThreads;

        std::shared_ptr<tbb::global_control> control;

#ifdef TIMING
	std::map<std::string,long long> duration;
#endif
    };

    template <typename Integer, typename Transform, typename Reduce, typename Output>
    inline Output transform_reduce(Integer begin, Integer end,
                                   Output result, Transform transform, Reduce reduce) {
        for (; begin != end; ++begin) {
            result = reduce(result, transform(begin));
        }
        return result;
    }

//    template <typename Integer, typename Function>
//    inline void for_each(Integer begin, const Integer end, Function function, TbbAccumulate) {
//
//        tbb::parallel_for(
//                tbb::blocked_range<size_t>(begin, end
//                        //, 200
//                ),
//                [function](const tbb::blocked_range<size_t>& r) -> void {
//                    const auto end = r.end();
//                    for (auto i = r.begin(); i != end; ++i) {
//                        function(i);
//                    }
//                }
//        );
//    };
}

#endif //ZIG_ZAG_ZIGZAG_HPP
