//
// Created by Marc Suchard on 2019-12-03.
//

#ifndef ZIG_ZAG_ZIGZAG_HPP
#define ZIG_ZAG_ZIGZAG_HPP

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedMacroInspection" // Turn off warning for TBB_PREVIEW_GLOBAL_CONTROL

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

#include "threefry.h"
//#include "dr_evomodel_operators_NativeZigZag.h"
#include "MemoryManagement.hpp"
#include "Simd.hpp"
#include "AbstractZigZag.hpp"

namespace zz {

    template<typename TypeInfo>
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
               double *rawParameterSign,
               long flags,
               int nThreads,
               long seed,
               DblSpan precision = DblSpan()) : AbstractZigZag(),
                            dimension(dimension),
                            mask(constructMask(rawMask, dimension, true)),
                            observed(constructMask(rawObserved, dimension, true)),
                            parameterSign(constructMask(rawParameterSign, dimension, false)),
                            mmPosition(dimension),
                            mmVelocity(dimension),
                            mmAction(dimension),
                            mmGradient(dimension),
                            mmMomentum(dimension),
                            flags(flags),
                            nThreads(nThreads),
                            seed(seed),
                            fixedPrecision(precision){
            std::cerr << "ZigZag constructed" << std::endl;
            std::cout << '\n';
            if (flags & zz::Flags::TBB) {
                if (nThreads <= 0) {
                    nThreads = tbb::task_scheduler_init::default_num_threads();
                }

                std::cout << "Using " << nThreads << " threads" << std::endl;

                control = std::make_shared<tbb::global_control>(tbb::global_control::max_allowed_parallelism, nThreads);
            }

            rng.resize(static_cast<std::size_t>(nThreads));
            for (int i = 0; i < nThreads; ++i) {
                rng[i].seed(static_cast<std::uint64_t>(seed + i));
            }
        }

        virtual ~ZigZag() {
#ifdef TIMING
            std::cerr << std::endl;
            for (auto &d: duration) {
                std::cerr << d.first << " " << std::scientific <<
                          static_cast<double>(d.second) * 0.001 << std::endl;
            }
#endif
        };

        template<typename T>
        struct Dynamics {

            template<typename V, typename W>
            Dynamics(V &position,
                     V &velocity,
                     V &action,
                     V &gradient,
                     V &momentum,
                     const W &observed,
                     const W &parameterSign) : position(position.data()),
                                               velocity(velocity.data()),
                                               action(action.data()),
                                               gradient(gradient.data()),
                                               momentum(momentum.data()),
                                               observed(observed.data()),
                                               parameterSign(parameterSign.data()),
                                               column(nullptr) {}

            template<typename V, typename W>
            Dynamics(V &position,
                     V &velocity,
                     V &action,
                     V &gradient,
                     std::nullptr_t,
                     const W &observed,
                     const W &parameterSign) : position(position.data()),
                                               velocity(velocity.data()),
                                               action(action.data()),
                                               gradient(gradient.data()),
                                               momentum(nullptr),
                                               observed(observed.data()),
                                               parameterSign(parameterSign.data()),
                                               column(nullptr) {}

            template<typename V, typename W>
            Dynamics(V &position,
                     V &velocity,
                     V &action,
                     V &gradient,
                     V &momentum,
                     const W &observed,
                     const W &parameterSign,
                     V &column) :position(position.data()),
                                 velocity(velocity.data()),
                                 action(action.data()),
                                 gradient(gradient.data()),
                                 momentum(momentum.data()),
                                 observed(observed.data()),
                                 parameterSign(parameterSign.data()),
                                 column(column.data()) {}

            ~Dynamics() = default;

            T *position;
            T *velocity;
            T *action;
            T *gradient;
            T *momentum;
            const T *observed;
            const T *parameterSign;
            T *column;
        };

        double operate(DblSpan position,
                       DblSpan velocity,
                       DblSpan action,
                       DblSpan gradient,
                       DblSpan momentum,
                       double time,
                       int dimension) {
            Dynamics<double> dynamics(position, velocity, action, gradient, momentum, observed, parameterSign);
            return operateImpl(dynamics, time, dimension);
        }

        void innerBounce(DblSpan position,
                         DblSpan velocity,
                         DblSpan action,
                         DblSpan gradient,
                         DblSpan momentum,
                         double time, int index, int type) {
#ifdef TIMING
            auto start = zz::chrono::steady_clock::now();
#endif

            Dynamics<double> dynamics(position, velocity, action, gradient, momentum, observed, parameterSign);
            innerBounceImpl(dynamics, time, index, type);

#ifdef TIMING
            auto end = zz::chrono::steady_clock::now();
            duration["innerBounce"] += zz::chrono::duration_cast<chrono::TimingUnits>(end - start).count();
#endif
        }

        void updateDynamics(DblSpan position,
                            DblSpan velocity,
                            DblSpan action,
                            DblSpan gradient,
                            DblSpan momentum,
                            DblSpan column,
                            double time, int index) {

#ifdef TIMING
            auto start = zz::chrono::steady_clock::now();
#endif

            Dynamics<double> dynamics(position, velocity, action, gradient, momentum, observed, parameterSign, column);
            updateDynamicsImpl<SimdType, SimdSize>(dynamics, time, index);


#ifdef TIMING
            auto end = zz::chrono::steady_clock::now();
            duration["updateDynamics"] += zz::chrono::duration_cast<chrono::TimingUnits>(end - start).count();
#endif
        }

        template<typename T>
        double operateImpl(Dynamics<T> &dynamics, double time, int dimension) {

#ifdef TIMING
            auto start = zz::chrono::steady_clock::now();
#endif

            BounceState bounceState(BounceType::NONE, -1, time);

            while (bounceState.isTimeRemaining()) {

                const auto firstBounce = getNextBounce(dynamics);

                bounceState = doBounce(bounceState, firstBounce, dynamics, dimension);
            }

#ifdef TIMING
            auto end = zz::chrono::steady_clock::now();
            duration["operateImpl"] += zz::chrono::duration_cast<chrono::TimingUnits>(end - start).count();
#endif

            return 0.0;
        }

        MinTravelInfo getNextBounce(DblSpan position,
                                    DblSpan velocity,
                                    DblSpan action,
                                    DblSpan gradient,
                                    DblSpan momentum) {

#if 0
            std::vector<double> b;

            auto buffer = [&b](DblSpan& in, mm::MemoryManager<double>& out) {
                mm::bufferedCopy(std::begin(in), std::end(in), std::begin(out), b);
            };

            buffer(position, mmPosition);
            buffer(velocity, mmVelocity);
            buffer(action, mmAction);
            buffer(gradient, mmGradient);
            buffer(momentum, mmMomentum);

            return getNextBounce(
                    Dynamics<double>(mmPosition, mmVelocity, mmAction, mmGradient, mmMomentum, observed, parameterSign, dimension));
#else
            return getNextBounce(
                    Dynamics<double>(position, velocity, action, gradient, momentum, observed, parameterSign));
#endif
        }

        MinTravelInfo getNextBounceIrreversible(DblSpan position,
                                                DblSpan velocity,
                                                DblSpan action,
                                                DblSpan gradient) {

            return getNextBounceIrreversible(
                    Dynamics<double>(position, velocity, action, gradient, nullptr, observed, parameterSign));
        }

        template<typename R>
        MinTravelInfo getNextBounceIrreversible(const Dynamics<R> &dynamics) {

#ifdef TIMING
            auto start = zz::chrono::steady_clock::now();
#endif

            auto task = [&](const size_t begin, const size_t end) -> MinTravelInfo {

                const auto length = end - begin;
                const auto vectorCount = length - length % SimdSize;

                MinTravelInfo travel = vectorized_transform<SimdType, SimdSize>(begin, begin + vectorCount, dynamics,
                                                                                InfoType());

                if (vectorCount < length) { // Edge-case
                    travel = vectorized_transform<RealType, 1>(begin + vectorCount, end, dynamics, travel);
                }

                return travel;
            };

//            MinTravelInfo travel = (nThreads <= 1) ?
//                                   task(size_t(0), dimension) :
//                                   parallel_task_reduce(
//                                           size_t(0), dimension, MinTravelInfo(),
//                                           task,
//                                           [](MinTravelInfo lhs, MinTravelInfo rhs) {
//                                               return (lhs.time < rhs.time) ? lhs : rhs;
//                                           });

            MinTravelInfo travel;
            travel.time = 42.0;
            travel.index = static_cast<int>(seed);

#ifdef TIMING
            auto end = zz::chrono::steady_clock::now();
            duration["getNextBounceIrr"] += zz::chrono::duration_cast<chrono::TimingUnits>(end - start).count();
#endif

            return travel;
        }

        template<typename R>
        MinTravelInfo getNextBounce(const Dynamics<R> &dynamics) {

#ifdef TIMING
            auto start = zz::chrono::steady_clock::now();
#endif
            auto task = [&](const size_t begin, const size_t end) -> MinTravelInfo {

                const auto length = end - begin;
                const auto vectorCount = length - length % SimdSize;

                MinTravelInfo travel = vectorized_transform<SimdType, SimdSize>(begin, begin + vectorCount, dynamics,
                                                                                InfoType());

                if (vectorCount < length) { // Edge-case
                    travel = vectorized_transform<RealType, 1>(begin + vectorCount, end, dynamics, travel);
                }

                return travel;
            };

            MinTravelInfo travel = (nThreads <= 1) ?
                                   task(size_t(0), dimension) :
                                   parallel_task_reduce(
                                           size_t(0), dimension, MinTravelInfo(),
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

        template<typename T, typename F, typename G>
        inline T parallel_task_reduce(size_t begin, const size_t end, T sum, F transform, G reduce) {

#if 0
            auto block = (end - begin) / nThreads;
            for ( ; begin < (end - block); begin += block) {
                sum = reduce(sum, transform(begin, begin + block));
            }
            return reduce(sum, transform(begin, end));
#else
            return tbb::parallel_reduce(
                    tbb::blocked_range<size_t>(begin, end, (end - begin) / nThreads
                    ),
                    sum,
                    [transform, reduce](const tbb::blocked_range<size_t> &r, T sum) { // TODO Test &transform, &reduce
                        return reduce(sum, transform(r.begin(), r.end()));
                    },
                    reduce
            );
#endif
        }

        template<typename F>
        inline void parallel_task_for(size_t begin, const size_t end, F transform) {

            tbb::parallel_for(
                    tbb::blocked_range<size_t>(begin, end, (end - begin) / nThreads),
                    [transform](const tbb::blocked_range<size_t> &r) {
                        transform(r.begin(), r.end());
                    }
            );
        }

    protected:

    private:

        template<typename S, int SimdSize, typename R, typename I, typename Int>
        MinTravelInfo vectorized_transform(Int i, const Int end,
                                           const Dynamics<R> &dynamics, I result) {

            const auto *position = dynamics.position;
            const auto *velocity = dynamics.velocity;
            const auto *action = dynamics.action;
            const auto *gradient = dynamics.gradient;
            const auto *momentum = dynamics.momentum;
            const auto *observed = dynamics.observed;
            const auto *parameterSign = dynamics.parameterSign;

            for (; i < end; i += SimdSize) {

                const auto boundaryTime = findBoundaryTime(
                        SimdHelper<S, R>::get(position + i),
                        SimdHelper<S, R>::get(velocity + i),
                        SimdHelper<S, R>::get(observed + i),
                        SimdHelper<S, R>::get(parameterSign + i)
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

        template<typename T>
        void innerBounceImpl(Dynamics<T> &dynamics,
                             const T eventTime, const int eventIndex, const int eventType) {

            updatePosition<SimdType, SimdSize>(dynamics, eventTime);
            updateMomentum<SimdType, SimdSize>(dynamics, eventTime);

            if (eventType == BounceType::BOUNDARY) {

                reflectMomentum(dynamics, eventIndex);

            } else {

                setZeroMomentum(dynamics, eventIndex);

            }

            reflectVelocity(dynamics, eventIndex);
            updateGradient(dynamics, eventTime);
        }


        template<typename S, int Size, typename R>
        void updateDynamicsImpl(Dynamics<R> &dynamics,
                                const R time, const int index) {

            auto p = dynamics.position;
            auto v = dynamics.velocity;
            auto a = dynamics.action;
            auto g = dynamics.gradient;
            auto m = dynamics.momentum;
            const auto o = dynamics.observed;
            const auto ps = dynamics.parameterSign;
            const auto c = dynamics.column;

            const R halfTimeSquared = time * time / 2;
            const R twoV = 2 * v[index];

            auto scalar = [p, v, a, g, m, o, ps, c,
                    time, halfTimeSquared, twoV](size_t i) {
                const R gi = g[i];
                const R ai = a[i];

                p[i] = p[i] + time * v[i];
                m[i] = m[i] + time * gi - halfTimeSquared * ai;
                g[i] = gi - time * ai;
                a[i] = ai - twoV * c[i];
            };

            const S timeS = S(time);
            const S halfTimeSquaredS = S(halfTimeSquared);
            const S twoVS = S(twoV);

            auto simd = [p, v, a, g, m, o, ps, c,
                    timeS, halfTimeSquaredS, twoVS](size_t i) {
                const S gi = SimdHelper<S, R>::get(g + i);
                const S ai = SimdHelper<S, R>::get(a + i);

                SimdHelper<S, R>::put(
                        SimdHelper<S, R>::get(p + i) + timeS * SimdHelper<S, R>::get(v + i),
                        p + i);
                SimdHelper<S, R>::put(
                        SimdHelper<S, R>::get(m + i) + timeS * gi - halfTimeSquaredS * ai,
                        m + i);
                SimdHelper<S, R>::put(
                        gi - timeS * ai,
                        g + i);
                SimdHelper<S, R>::put(
                        ai - twoVS * SimdHelper<S, R>::get(c + i),
                        a + i);
            };

            if (nThreads <= 1) {
                simd_for_each<Size>(size_t(0), dimension, simd, scalar);
            } else {
                parallel_task_for(size_t(0), dimension,
                                  [simd, scalar](size_t begin, size_t end) { // TODO &task?
                                      simd_for_each<Size>(begin, end, simd, scalar);
                                  });
            }


        }


        template<typename R>
        BounceState doBounce(BounceState initialBounceState, MinTravelInfo firstBounce, Dynamics<R> &dynamics,
                             int dimension) {

            double remainingTime = initialBounceState.time;
            double eventTime = firstBounce.time;

            BounceState finalBounceState;
            if (remainingTime < eventTime) { // No event during remaining time
                updatePosition<SimdType, SimdSize>(dynamics, remainingTime);
                finalBounceState = BounceState(BounceType::NONE, -1, 0.0);

            } else {

                updatePosition<SimdType, SimdSize>(dynamics, eventTime);
                updateMomentum<SimdType, SimdSize>(dynamics, eventTime);

                const int eventType = firstBounce.type;
                const int eventIndex = firstBounce.index;

                DblSpan precisionColumn = fixedPrecision.subspan(eventIndex * dimension, dimension);

                if (eventType == BounceType::BOUNDARY) {

                    reflectMomentum(dynamics, eventIndex);

                } else {

                    setZeroMomentum(dynamics, eventIndex);

                }

                reflectVelocity(dynamics, eventIndex);
                updateGradient(dynamics, eventTime);
                updateAction(dynamics, eventIndex, precisionColumn);

                finalBounceState = BounceState(eventType, eventIndex, remainingTime - eventTime);
            }

            return finalBounceState;
        }

        template<typename S, int Size, typename R>
        inline void updatePosition(Dynamics<R> &dynamics, R time) {
            auto position = dynamics.position;
            const auto velocity = dynamics.velocity;

            auto scalar = [position, velocity, time](size_t i) {
                position[i] = position[i] + time * velocity[i];
            };

            auto simd = [position, velocity, time](size_t i) {
                SimdHelper<S, R>::put(
                        SimdHelper<S, R>::get(position + i)
                        + time * SimdHelper<S, R>::get(velocity + i),
                        position + i);
            };

            if (nThreads <= 1) {
                simd_for_each<Size>(size_t(0), dimension, simd, scalar);
            } else {
                parallel_task_for(size_t(0), dimension,
                                  [simd, scalar](size_t begin, size_t end) { // TODO &task?
                                      simd_for_each<Size>(begin, end, simd, scalar);
                                  });
            }
        }

        template<typename S, int Size, typename R>
        inline void updateMomentum(Dynamics<R> &dynamics, R time) {
            auto momentum = dynamics.momentum;
            const auto action = dynamics.action;
            const auto gradient = dynamics.gradient;
            const auto mk = mask.data(); // TODO Delegate

            const auto halfTimeSquared = time * time / 2;

            auto scalar = [momentum, action, gradient, mk, time, halfTimeSquared](size_t i) {
                momentum[i] = mk[i] * (momentum[i] + time * gradient[i] - halfTimeSquared * action[i]);
            };

            auto simd = [momentum, action, gradient, mk, time, halfTimeSquared](size_t i) {
                SimdHelper<S, R>::put(
                        SimdHelper<S, R>::get(mk + i) * (
                                SimdHelper<S, R>::get(momentum + i) +
                                time * SimdHelper<S, R>::get(gradient + i) -
                                halfTimeSquared * SimdHelper<S, R>::get(action + i)),
                        momentum + i);
            };

            if (mask.size() > 0) {
                if (nThreads <= 1) {
                    simd_for_each<Size>(size_t(0), dimension, simd, scalar);
                } else {
                    parallel_task_for(size_t(0), dimension,
                                      [simd, scalar](size_t begin, size_t end) { // TODO &task?
                                          simd_for_each<Size>(begin, end, simd, scalar);
                                      });
                }
            } else {
                exit(-1); // TODO Implement
            }
        }

        template<typename R>
        inline void updateAction(Dynamics<R> &dynamics, int index, const DblSpan precCol) {
            auto momentum = dynamics.momentum;
            const auto action = dynamics.action;
            const auto velocity = dynamics.velocity;
            const auto mk = mask.data(); // TODO Delegate

#ifdef TIMING
            auto start = zz::chrono::steady_clock::now();
#endif

            //const auto column = callback.getColumn(index);

#ifdef TIMING
            auto end = zz::chrono::steady_clock::now();
            duration["getColumn"] += zz::chrono::duration_cast<chrono::TimingUnits>(end - start).count();
#endif

            const auto twoV = 2 * velocity[index];

            if (mask.size() > 0) {
                vectorized_for_each(size_t(0), dimension,
                                    [action, precCol, mk, twoV](size_t i) {
                                        action[i] = mk[i] * (action[i] + twoV * precCol[i]);
                                    });
            } else {
                exit(-1); // TODO Implement
            }
        }

        template<typename R>
        inline void updateGradient(Dynamics<R> &dynamics, R time) {
            const auto action = dynamics.action;
            auto gradient = dynamics.gradient;

            auto task = [action, gradient, time](size_t i) {
                gradient[i] = gradient[i] - time * action[i];
            };

            if (nThreads <= 1) {
                vectorized_for_each(size_t(0), dimension, task);
            } else {
                parallel_task_for(size_t(0), dimension,
                                  [task](size_t begin, size_t end) { // TODO &task?
                                      vectorized_for_each(begin, end, task);
                                  });
            }
        }

        template<typename R>
        static inline void reflectMomentum(Dynamics<R> &dynamics, int index) {
            auto position = dynamics.position;
            auto momentum = dynamics.momentum;

            momentum[index] = -momentum[index];
            position[index] = R(0.0);

        }

        template<typename R>
        static inline void setZeroMomentum(Dynamics<R> &dynamics, int index) {
            auto momentum = dynamics.momentum;

            momentum[index] = R(0.0);
        }

        template<typename R>
        static inline void reflectVelocity(Dynamics<R> &dynamics, int index) {
            auto velocity = dynamics.velocity;

            velocity[index] = -velocity[index];
        }

        template<int Size, typename I, typename FV, typename FS>
        static inline void simd_for_each(I begin, const I end, FV vector, FS scalar) {

            if (Size > 1) { // TODO is this compile-time?
                const auto length = end - begin;
                const auto simdLength = length - length % SimdSize;
                const auto simdEnd = begin + simdLength;

                for (; begin < simdEnd; begin += SimdSize) {
                    vector(begin);
                }
            }

            for (; begin < end; ++begin) {
                scalar(begin);
            }
        }

        template<typename I, typename F>
        static inline void vectorized_for_each(I begin, const I end, F function) {
            for (; begin < end; ++begin) {
                function(begin);
            }
        }

        static inline void reduce_min(MinTravelInfo &result,
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

        static inline void reduce_min(DoubleSseMinTravelInfo &result,
                                      const D2 time, const int index, const int type) {
            const auto lessThan = time < result.time;
            if (xsimd::any(lessThan)) {
                result.time = select(lessThan, time, result.time);
                const auto mask = _mm_castpd_si128(lessThan);
                result.index = select(mask, makeSimdIndex<D2Index>(index),
                                      result.index); // TODO Merge into single register?
                result.type = select(mask, D2Index(type), result.type);
            }
        }

        static inline void reduce_min(DoubleAvxMinTravelInfo &result, // TODO Remove code-dup with above
                                      const D4 time, const int index, const int type) {
            const auto lessThan = time < result.time;
            if (xsimd::any(lessThan)) {
                result.time = select(lessThan, time, result.time);
                const auto mask = _mm256_castpd_si256(lessThan);
                result.index = select(mask, makeSimdIndex<D4Index>(index),
                                      result.index); // TODO Merge into single register?
                result.type = select(mask, D4Index(type), result.type);
            }
        }

        static inline MinTravelInfo horizontal_min(DoubleSseMinTravelInfo vector) {
            return (vector.time[0] < vector.time[1]) ?
                   MinTravelInfo(static_cast<int>(vector.type[0]), static_cast<int>(vector.index[0]), vector.time[0]) :
                   MinTravelInfo(static_cast<int>(vector.type[1]), static_cast<int>(vector.index[1]), vector.time[1]);
        }

        static inline MinTravelInfo horizontal_min(DoubleAvxMinTravelInfo vector) {

            auto const firstHalf = (vector.time[0] < vector.time[1]) ?
                                   MinTravelInfo(static_cast<int>(vector.type[0]), static_cast<int>(vector.index[0]),
                                                 vector.time[0]) :
                                   MinTravelInfo(static_cast<int>(vector.type[1]), static_cast<int>(vector.index[1]),
                                                 vector.time[1]);

            auto const secondHalf = (vector.time[2] < vector.time[3]) ?
                                    MinTravelInfo(static_cast<int>(vector.type[2]), static_cast<int>(vector.index[2]),
                                                  vector.time[2]) :
                                    MinTravelInfo(static_cast<int>(vector.type[3]), static_cast<int>(vector.index[3]),
                                                  vector.time[3]);

            return (firstHalf.time < secondHalf.time) ? firstHalf : secondHalf;
        }

        template<typename T>
        static inline T findBoundaryTime(const T position,
                                         const T velocity,
                                         const T observed,
                                         const T parameterSign) {
            return select(headingTowardsBoundary(parameterSign, velocity, observed),
                          fabs(position / velocity),
                          infinity<T>());
        }

        template<typename T>
        static inline auto headingTowardsBoundary(const T parameterSign,
                                                  const T velocity,
                                                  const T observed)
        -> decltype(T(1.0) > T(0.0)) {
            return parameterSign * velocity * observed < T(0.0);
        }

        template<typename T>
        static inline T minimumPositiveRoot(const T a, const T b, const T c) {

            const auto discriminant = b * b - 4 * a * c;
            const auto sqrtDiscriminant = select(c == T(0.0), b, sqrt(fabs(discriminant)));

            auto root1 = (-b - sqrtDiscriminant) / (2 * a);
            auto root2 = (-b + sqrtDiscriminant) / (2 * a);

            root1 = select(root1 > T(0.0), root1, infinity<T>());
            root2 = select(root2 > T(0.0), root2, infinity<T>());

            const auto root = select(root1 < root2, root1, root2);
            return select(discriminant < T(0.0), infinity<T>(), root);
        }

//        template <typename T>
//        static inline T sign(const T x) {
//            const auto zero = T(0.0);
//            return select(x == zero,
//                          x,
//                          select(x < zero,
//                                 T(-1.0),
//                                 T(+1.0))
//            );
//        }

        static mm::MemoryManager<MaskType> constructMask(double *raw, size_t length, bool zeroOneFlg) {

            mm::MemoryManager<MaskType> mask;
            mask.reserve(length);

            std::transform(raw, raw + length, std::back_inserter(mask),
                           [&zeroOneFlg](double x) {
                               if (zeroOneFlg) {
                                   return (x == 1.0) ? MaskType(1.0) : MaskType(0.0);
                               } else {
                                   return (x == 1.0) ? MaskType(1.0) : MaskType(-1.0);
                               }
                           });

            return mask;
        }

//        static double doMask(MaskType mask, double x) {
//            return mask * x;
//        }

        size_t dimension;
        mm::MemoryManager<MaskType> mask;
        mm::MemoryManager<MaskType> observed;
        mm::MemoryManager<MaskType> parameterSign;

        mm::MemoryManager<double> mmPosition;
        mm::MemoryManager<double> mmVelocity;
        mm::MemoryManager<double> mmAction;
        mm::MemoryManager<double> mmGradient;
        mm::MemoryManager<double> mmMomentum;

        DblSpan fixedPrecision;
        long flags;
        int nThreads;
        long seed;

        std::shared_ptr<tbb::global_control> control;

        std::vector<sitmo::threefry_20_64> rng;

#ifdef TIMING
        std::map<std::string, long long> duration;
#endif
    };

//    template <typename Integer, typename Transform, typename Reduce, typename Output>
//    inline Output transform_reduce(Integer begin, Integer end,
//                                   Output result, Transform transform, Reduce reduce) {
//        for (; begin != end; ++begin) {
//            result = reduce(result, transform(begin));
//        }
//        return result;
//    }

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

    std::unique_ptr<zz::AbstractZigZag> dispatch(
            int dimension,
            double *rawMask,
            double *rawObserved,
            double *rawParameterSign,
            long flags,
            int info,
            long seed,
            DblSpan precision=DblSpan()) {

        if (static_cast<unsigned long>(flags) & zz::Flags::AVX) {
            std::cerr << "Factory: AVX" << std::endl;
            return zz::make_unique<zz::ZigZag<zz::DoubleAvxTypeInfo>>(
                    dimension, rawMask, rawObserved, rawParameterSign, flags, info, seed, precision);
        } else if (static_cast<unsigned long>(flags) & zz::Flags::SSE) {
            std::cerr << "Factory: SSE" << std::endl;
            return zz::make_unique<zz::ZigZag<zz::DoubleSseTypeInfo>>(
                    dimension, rawMask, rawObserved, rawParameterSign, flags, info, seed, precision);
        } else {
            std::cerr << "Factory: No SIMD" << std::endl;
            return zz::make_unique<zz::ZigZag<zz::DoubleNoSimdTypeInfo>>(
                    dimension, rawMask, rawObserved, rawParameterSign, flags, info, seed, precision);
        }
    }
}

#pragma clang diagnostic pop

#endif //ZIG_ZAG_ZIGZAG_HPP
