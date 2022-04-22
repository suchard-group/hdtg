//
// Created by Zhenyu Zhang on 1/30/22.
//
#include <vector>
#include <cmath>

#define TBB_PREVIEW_GLOBAL_CONTROL 1

#define TIMING

#ifdef TIMING

#include <map>
#include <iomanip>
#include <random>
#include "Timing.h"

#endif // TIMING

#include "threefry.h"
//#include "dr_evomodel_operators_NativeZigZag.h"
#include "MemoryManagement.hpp"
#include "Simd.hpp"
#include "ZigZag.hpp"
#include "NoUTurn.hpp"
#include "UniformGenerator.h"

namespace nuts {
    using DblSpan = tcb::span<double>;

    class TreeState {
        friend class NoUTurn;

    public:

        TreeState(DblSpan position, DblSpan momentum, DblSpan gradient,
                  int numNodes, bool flagContinue,
                  double cumAcceptProb, int numAcceptProbStates, UniformGenerator &generator) :
                                                        dim(position.size()),
                                                        positionTri(position.size() * 3, 0),
                                                        momentumTri(position.size() * 3, 0),
                                                        gradientTri(position.size() * 3, 0), numNodes(numNodes),
                                                        flagContinue(flagContinue),
                                                        cumAcceptProb(cumAcceptProb),
                                                        numAcceptProbStates(numAcceptProbStates),
                                                        uniGenerator(generator) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < dim; ++j) {
                    positionTri[i * dim + j] = position[j];
                    momentumTri[i * dim + j] = momentum[j];
                    gradientTri[i * dim + j] = gradient[j];
                }
            }
        }

        DblSpan getPosition(int direction) {
            return DblSpan(&positionTri[getIndex(direction) * dim], dim);
        }

        DblSpan getMomentum(int direction) {
            return DblSpan(&momentumTri[getIndex(direction) * dim], dim);
        }

        DblSpan getGradient(int direction) {
            return DblSpan(&gradientTri[getIndex(direction) * dim], dim);
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
                positionTri[getIndex(direction) * dim + j] = position[j];
            }
        }

        void setMomentum(int direction, DblSpan momentum) {
            for (int j = 0; j < dim; ++j) {
                momentumTri[getIndex(direction) * dim + j] = momentum[j];
            }
        }

        void setGradient(int direction, DblSpan gradient) {
            for (int j = 0; j < dim; ++j) {
                gradientTri[getIndex(direction) * dim + j] = gradient[j];
            }
        }

        void setSample(DblSpan position) { setPosition(0, position); }

        int getIndex(int direction) { // valid directions: -1, 0, +1
            assert (direction == -1 || direction == 1 || direction == 0);
            return direction + 1;
        }

        bool computeStopCriterion() {
            DblSpan positionPlus = getPosition(1);
            DblSpan positionMinus = getPosition(-1);
            DblSpan momentumPlus = getMomentum(1);
            DblSpan momentumMinus = getMomentum(-1);

            double sumPlus = 0;
            double sumMinus = 0;
            double positionDiffI;
            for (int i = 0; i < positionPlus.size(); ++i) {
                positionDiffI = positionPlus[i] - positionMinus[i];
                sumPlus += positionDiffI * momentumPlus[i];
                sumMinus += positionDiffI * momentumMinus[i];
            }

            return (sumPlus > 0) && (sumMinus > 0);
        }

        void mergeNextTree(TreeState nextTree, int direction) {

            setPosition(direction, nextTree.getPosition(direction));
            setMomentum(direction, nextTree.getMomentum(direction));
            setGradient(direction, nextTree.getGradient(direction));

            updateSample(nextTree);

            numNodes += nextTree.numNodes;
            flagContinue = nextTree.flagContinue && computeStopCriterion();

            cumAcceptProb += nextTree.cumAcceptProb;
            numAcceptProbStates += nextTree.numAcceptProbStates;
        }

        void updateSample(TreeState nextTree) {
            if (nextTree.numNodes > 0
                && uniGenerator.getUniform() < ((double) nextTree.numNodes / (double) (numNodes + nextTree.numNodes))) {
                setSample(nextTree.getSample());
            }
        }

        int dim;
        std::vector<double> positionTri;
        std::vector<double> momentumTri;
        std::vector<double> gradientTri;

        int numNodes;
        bool flagContinue;

        double cumAcceptProb;
        int numAcceptProbStates;
        UniformGenerator &uniGenerator;
    };
}

