//
// Created by Marc Suchard on 2019-12-03.
//

#ifndef NUTS_HPP
#define NUTS_HPP

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedMacroInspection" // Turn off warning for TBB_PREVIEW_GLOBAL_CONTROL

#include <vector>
#include <cmath>

#define TBB_PREVIEW_GLOBAL_CONTROL 1

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
#include "ZigZag.hpp"

namespace nuts {
    using DblSpan = tcb::span<double>;

    class TreeState {
        friend class NoUTurn;

    private:

        TreeState(DblSpan position, DblSpan momentum, DblSpan gradient,
                  int numNodes, bool flagContinue,
                  double cumAcceptProb, int numAcceptProbStates) : numNodes(numNodes),
                                                                   flagContinue(flagContinue),
                                                                   cumAcceptProb(cumAcceptProb),
                                                                   numAcceptProbStates(numAcceptProbStates),
                                                                   dim(position.size()) {
            std::vector<double> positionTri(dim * 3, 0);
            std::vector<double> momentumTri(dim * 3, 0);
            std::vector<double> gradientTri(dim * 3, 0);

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < dim; ++j) {
                    positionTri[i * dim + j] = position[j]
                    momentumTri[i * dim + j] = momentum[j]
                    gradientTri[i * dim + j] = gradient[j]
                }
            }
            positionTriSpan = DblSpan{positionTri};
            momentumTriSpan = DblSpan{momentumTri};
            gradientTriSpan = DblSpan{gradientTri};
        }

        DblSpan getPosition(int direction) {
            return positionTriSpan.subspan(getIndex(direction) * dim, dim)
        }

        DblSpan getMomentum(int direction) {
            return momentumTriSpan.subspan(getIndex(direction) * dim, dim)
        }

        DblSpan getGradient(int direction) {
            return gradientTriSpan.subspan(getIndex(direction) * dim, dim)
        }

        DblSpan getSample() {
            /*
            Returns a state chosen uniformly from the acceptable states along a hamiltonian dynamics trajectory tree.
            The sample is updated recursively while building trees.
            */
            return getPosition(0);
        }

        void setPosition(int direction, DblSpan position) {
            for (int j = 0; j < dim; ++j) {
                gradientTriSpan[getIndex(direction) * dim + j] = position[j]
            }
        }

        void setMomentum(int direction, DblSpan momentum) {
            for (int j = 0; j < dim; ++j) {
                momentumTriSpan[getIndex(direction) * dim + j] = momentum[j]
            }
        }

        void setGradient(int direction, DblSpan gradient) {
            for (int j = 0; j < dim; ++j) {
                gradientTriSpan[getIndex(direction) * dim + j] = gradient[j]
            }
        }

        void setSample(DblSpan position) { setPosition(0, position); }

        int getIndex(int direction) { // valid directions: -1, 0, +1
            assert (direction == -1 || direction == 1 || direction == 0);
            return direction + 1;
        }

        void mergeNextTree(TreeState nextTree, int direction) {

            setPosition(direction, nextTree.getPosition(direction));
            setMomentum(direction, nextTree.getMomentum(direction));
            setGradient(direction, nextTree.getGradient(direction));

            updateSample(nextTree);

            numNodes += nextTree.numNodes;
            flagContinue = computeStopCriterion(nextTree.flagContinue, this);

            cumAcceptProb += nextTree.cumAcceptProb;
            numAcceptProbStates += nextTree.numAcceptProbStates;
        }

        void updateSample(TreeState nextTree) {
            double uniform = getUniform();
            if (nextTree.numNodes > 0
                && uniform < ((double) nextTree.numNodes / (double) (numNodes + nextTree.numNodes))) {
                setSample(nextTree.getSample());
            }
        }

        DblSpan positionTriSpan;
        DblSpan momentumTriSpan;
        DblSpan gradientTriSpan;
        int dim;

        int numNodes;
        bool flagContinue;

        double cumAcceptProb;
        int numAcceptProbStates;
    };

    class NoUTurn {
    public:
        NoUTurn(double logProbErrorTol,
                int findMax,
                int maxHeight) : logProbErrorTol(logProbErrorTol),
                                 findMax(findMax),
                                 maxHeight(maxHeight) {
            std::cerr << "nuts constructed" << '\n' << std::endl;
        }

        ~NoUTurn() = default;

        static void takeOneStep(DblSpan initialPosition, DblSpan initialMomentum, DblSpan gradient) {

            DblSpan endPosition = initialPosition;

            const double initialJointDensity = zz::getJointProbability(initialPosition, initialMomentum); //to finish

            double logSliceU = log(getUniform()) + initialJointDensity;

            TreeState trajectoryTree = TreeState(initialPosition, initialMomentum, grandient, 1, true,
                                                 0, 0); //to finish
            int height = 0;

            while (trajectoryTree.flagContinue) {

                tmp = updateTrajectoryTree(trajectoryTree, height, logSliceU, initialJointDensity);
                if (!tmp.empty()) {
                    endPosition = tmp;
                }

                height++;

                if (height > options.maxHeight) {
                    trajectoryTree.flagContinue = false;
                }
            }
            return endPosition;
        }

        DblSpan updateTrajectoryTree(TreeState trajectoryTree,
                                     int depth,
                                     double logSliceU,
                                     double initialJointDensity
        ) {
            DblSpan endPosition;

            const double uniform1 = getUniform();
            int direction = (uniform1 < 0.5) ? -1 : 1;
            TreeState nextTrajectoryTree = buildTree(
                    trajectoryTree.getPosition(direction), trajectoryTree.getMomentum(direction),
                    trajectoryTree.getGradient(direction),
                    direction, logSliceU, depth, stepSizeInformation.getStepSize(), initialJointDensity);

            if (nextTrajectoryTree.flagContinue) {

                const double uniform = getUniform();
                const double acceptProb = (double) nextTrajectoryTree.numNodes / (double) trajectoryTree.numNodes;
                if (uniform < acceptProb) {
                    endPosition = nextTrajectoryTree.getSample();
                }
            }

            trajectoryTree.mergeNextTree(nextTrajectoryTree, direction);

            return endPosition;
        }

        TreeState buildTree(DblSpan position, DblSpan momentum, DblSpan gradient, int direction,
                            double logSliceU, int height, double stepSize, double initialJointDensity) {

            if (height == 0) {
                return buildBaseCase(position, momentum, gradient, direction, logSliceU, stepSize, initialJointDensity);
            } else {
                return buildRecursiveCase(position, momentum, gradient, direction, logSliceU, height, stepSize,
                                          initialJointDensity);
            }
        }

        TreeState buildBaseCase(DblSpan inPosition, DblSpan inMomentum, DblSpan inGradient, int direction,
                                double logSliceU, double stepSize, double initialJointDensity) {
            // Make deep copy of position and momentum
            std::vector<double> positionVec;
            std::vector<double> momentumVec;
            std::vector<double> gradientVec;

            position.assign(inPosition.begin(), inPosition.end());
            momentum.assign(inMomentum.begin(), inMomentum.end());
            gradient.assign(inGradient.begin(), inGradient.end());

            DblSpan position{positionVec};
            DblSpan momentum{momentumVec};
            DblSpan gradient{gradientVec};
            //hmcProvider.setParameter(position.getBuffer());
            // "one reversibleHMC integral
            ZigZag.reversiblePositionMomentumUpdate(position, momentum, gradient, direction, stepSize);

            double logJointProbAfter = ZigZag.getJointProbability(momentum);

            const int numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);

            const bool flagContinue = (logSliceU < options.logProbErrorTol + logJointProbAfter);

            // Values for dual-averaging
            const double acceptProb = Math.min(1.0, Math.exp(logJointProbAfter - initialJointDensity));
            const int numAcceptProbStates = 1;

            //hmcProvider.setParameter(inPosition);

            return new TreeState(position, momentum, gradient, numNodes,
                                 flagContinue, acceptProb,
                                 numAcceptProbStates);
        }
        
        TreeState buildRecursiveCase(DblSpan inPosition, DblSpan inMomentum, DblSpan gradient, int direction,
                                     double logSliceU, int height, double stepSize, double initialJointDensity) {

            TreeState subtree = buildTree(inPosition, inMomentum, gradient, direction, logSliceU,
                                          height - 1, // Recursion
                                          stepSize, initialJointDensity);

            if (subtree.flagContinue) {

                TreeState nextSubtree = buildTree(subtree.getPosition(direction), subtree.getMomentum(direction),
                                                  subtree.getGradient(direction), direction,
                                                  logSliceU, height - 1, stepSizeInformation.getStepSize(),
                                                  initialJointDensity);

                subtree.mergeNextTree(nextSubtree, direction);
            }
            return subtree;
        }

        static bool computeStopCriterion(bool flagContinue, TreeState state) {
            return computeStopCriterion(flagContinue,
                                        state.getPosition(1), state.getPosition(-1),
                                        state.getMomentum(1), state.getMomentum(-1));
        }

//        StepSize findReasonableStepSize(DblSpan initialPosition, DblSpan initialGradient,
//                                        double forcedInitialStepSize) {
//
//            if (forcedInitialStepSize != 0) {
//                return new StepSize(forcedInitialStepSize);
//            } else {
//                double stepSize = 0.1;
//
//                WrappedVector momentum = hmcProvider.drawMomentum();
//                int count = 1;
//                int dim = initialPosition.length;
//                WrappedVector
//                position = new WrappedVector.Raw(Arrays.copyOf(initialPosition, dim));
//                WrappedVector
//                gradient = new WrappedVector.Raw(Arrays.copyOf(initialGradient, dim));
//
//                double probBefore = hmcProvider.getJointProbability(momentum);
//
//                hmcProvider.reversiblePositionMomentumUpdate(position, momentum, gradient, 1, stepSize);
//
//                double probAfter = hmcProvider.getJointProbability(momentum);
//
//                double a = ((probAfter - probBefore) > Math.log(0.5) ? 1 : -1);
//
//                double probRatio = Math.exp(probAfter - probBefore);
//
//                while (Math.pow(probRatio, a) > Math.pow(2, -a)) {
//
//                    probBefore = probAfter;
//                    hmcProvider.reversiblePositionMomentumUpdate(position, momentum, gradient, 1, stepSize);
//
//                    probAfter = hmcProvider.getJointProbability(momentum);
//                    probRatio = Math.exp(probAfter - probBefore);
//
//                    stepSize = Math.pow(2, a) * stepSize;
//                    count++;
//
//                    if (count > options.findMax) {
//                        throw new RuntimeException("Cannot find a reasonable step-size in " + options.findMax + " " +
//                                                   "iterations");
//                    }
//                }
//                hmcProvider.setParameter(initialPosition);
//                return new StepSize(stepSize);
//            }
//        }

        static bool computeStopCriterion(bool flagContinue,
                                            DblSpan positionPlus, DblSpan positionMinus,
                                            DblSpan momentumPlus, DblSpan momentumMinus) {

            Dblspan positionDifference = subtractArray(positionPlus, positionMinus);


            return flagContinue &&
                   getDotProduct(positionDifference, momentumMinus) >= 0 &&
                   getDotProduct(positionDifference, momentumPlus) >= 0;
        }

        static double getDotProduct(DblSpan x, DblSpan y) {

            assert (x.size() == y.size());
            const int dim = x.size();

            double total = 0.0;
            for (int i = 0; i < dim; i++) {
                total += x[i] * y[i];
            }
            return total;
        }

        DblSpan subtractArray(DblSpan a, DblSpan b) {
            assert (a.size() == b.size());
            const int dim = a.length;

            double result[dim];
            for (int i = 0; i < dim; i++) {
                result[i] = a[i] - b[i];
            }

            DblSpan res{result};
            return res;
        }

        static double getUniform() {
            static std::default_random_engine e;
            static std::uniform_real_distribution<> dis(0, 1); // rage 0 - 1
            return dis(e);
        }

        double logProbErrorTol = 100.0;
        const int findMax = 100;
        const int maxHeight = 100;
    };
}

#pragma clang diagnostic pop

#endif //ZIG_ZAG_ZIGZAG_HPP
