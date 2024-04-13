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
#include "MemoryManagement.h"
#include "Simd.h"
#include "ZigZag.h"
#include "NoUTurn.h"
#include "UniformGenerator.h"

namespace nuts {
    using DblSpan = tcb::span<double>;

    class TreeState {
        friend class NoUTurn;

    public:

        TreeState(DblSpan position, DblSpan momentum, DblSpan gradient,
                  int numAcceptableStates, bool flagContinue,
                  double cumAcceptProb, int numStates, UniformGenerator &generator) :
                                                        dim(position.size()),
                                                        positionTri(position.size() * 3, 0),
                                                        momentumTri(position.size() * 3, 0),
                                                        gradientTri(position.size() * 3, 0), 
                                                        numAcceptableStates(numAcceptableStates),
                                                        flagContinue(flagContinue),
                                                        cumAcceptProb(cumAcceptProb),
                                                        numStates(numStates),
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

        bool checkNoUturn() {
            DblSpan positionFront = getPosition(1);
            DblSpan positionRear = getPosition(-1);
            DblSpan momentumFront = getMomentum(1);
            DblSpan momentumRear = getMomentum(-1);

            // Derivatives of 0.5 * \| positionFront - positionRear \|^2, which are
            // inner products with the position difference and front/rear momenta.
            double distDerivFront = 0; 
            double distDerivRear = 0;
            for (int i = 0; i < positionFront.size(); ++i) {
                distDerivFront += (positionFront[i] - positionRear[i]) * momentumFront[i];
                distDerivRear += (positionFront[i] - positionRear[i]) * momentumRear[i];
            }

            return (distDerivFront > 0) && (distDerivRear > 0);
        }

        void mergeNextTree(TreeState nextTree, int direction, bool swapSampling) {

            setPosition(direction, nextTree.getPosition(direction));
            setMomentum(direction, nextTree.getMomentum(direction));
            setGradient(direction, nextTree.getGradient(direction));

            updateSample(nextTree, swapSampling);

            numAcceptableStates += nextTree.numAcceptableStates;
            flagContinue = nextTree.flagContinue && checkNoUturn();

            cumAcceptProb += nextTree.cumAcceptProb;
            numStates += nextTree.numStates;
        }

        void updateSample(TreeState nextTree, bool swapSampling) {
            double samplingWeightOnNext;
            if (swapSampling) {
              samplingWeightOnNext = (double) nextTree.numAcceptableStates / (double) numAcceptableStates;
            } else {
              samplingWeightOnNext = 
                (double) nextTree.numAcceptableStates / (double) (numAcceptableStates + nextTree.numAcceptableStates);
            }
            if (uniGenerator.getUniform() < samplingWeightOnNext) {
                setSample(nextTree.getSample());
            }
        }

        int dim;
        std::vector<double> positionTri;
        std::vector<double> momentumTri;
        std::vector<double> gradientTri;

        int numAcceptableStates;
        bool flagContinue;

        double cumAcceptProb;
        int numStates;
        UniformGenerator &uniGenerator;
    };
}

