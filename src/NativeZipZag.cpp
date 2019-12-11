#include <iostream>
#include <vector>
#include <memory>

#include "dr_evomodel_operators_NativeZigZag.h"
#include "ZigZag.hpp"

std::vector<std::unique_ptr<zz::AbstractZigZag>> implementation;

static jclass classMTL;
static jmethodID cid;

JNIEXPORT jint JNICALL JNI_OnLoad(JavaVM* vm, void* reserved) {

    JNIEnv *env;

    (void)reserved;

    if (vm->GetEnv(reinterpret_cast<void**>(&env), JNI_VERSION_1_6) != JNI_OK) {
        return JNI_ERR;
    }

    jclass tmpClassMTL = env->FindClass("dr/inference/operators/hmc/MinimumTravelInformation");
    classMTL =  (jclass) env->NewGlobalRef(tmpClassMTL);
    env->DeleteLocalRef(tmpClassMTL);

    cid = env->GetMethodID(classMTL, "<init>", "(DII)V");

    return JNI_VERSION_1_6;
}

JNIEXPORT void JNICALL JNI_OnUnload(JavaVM *vm, void *reserved) {
    JNIEnv *env;

    (void)reserved;

    if (vm->GetEnv(reinterpret_cast<void**>(&env), JNI_VERSION_1_6) != JNI_OK)
        return;

    env->DeleteGlobalRef(classMTL);
}


JNIEXPORT jint JNICALL Java_dr_evomodel_operators_NativeZigZag_create(
        JNIEnv *env,
        jobject obj,
        jint dimension,
        jobject provider,
        jdoubleArray mask,
        jdoubleArray observed) {

    double *rawMask = env->GetDoubleArrayElements(mask, NULL);
    double *rawObserved = env->GetDoubleArrayElements(observed, NULL);

    long flags = zz::Flags::TBB;
    int nThreads = 4;

    int instanceNumber = implementation.size();
    implementation.emplace_back(zz::make_unique<
//            zz::ZigZag<zz::DoubleNoSimdTypeInfo>
zz::ZigZag<zz::DoubleSseTypeInfo>
    >(dimension, rawMask, rawObserved, flags, nThreads));

    env->ReleaseDoubleArrayElements(mask, rawMask, JNI_ABORT);
    env->ReleaseDoubleArrayElements(observed, rawObserved, JNI_ABORT);

    return instanceNumber;
}


class JniCallback : public zz::PrecisionColumnCallback {
public:
    JniCallback(JNIEnv *env, jobject provider)
            : PrecisionColumnCallback(), env(env), provider(provider), array(nullptr), data(nullptr), column(-1) {
        jclass cls = env->GetObjectClass(provider);
        mid = env->GetMethodID(cls, "getColumn", "(I)[D");
        if (mid == 0) {
            throw;
        }
    }

    ~JniCallback() { releaseColumn(); }

    const double* getColumn(int index) {

        if (index != column) {

            releaseColumn();

            jobject object = env->CallObjectMethod(provider, mid, index);
            array = reinterpret_cast<jdoubleArray *>(&object);
            data = env->GetDoubleArrayElements(*array, &isCopy);
        }

        column = index;
        return data;
    }

    void releaseColumn() {

        if (data != nullptr) {
            if (isCopy == JNI_TRUE) {
                env->ReleaseDoubleArrayElements(*array, data, JNI_ABORT);
            }
            data = nullptr;
        }

        column = -1;
    }

private:
    JNIEnv *env;
    jobject provider;
    jmethodID mid;

    jdoubleArray *array;
    double *data;
    int column;
    jboolean isCopy;
};

JNIEXPORT jint JNICALL Java_dr_evomodel_operators_NativeZigZag_operate(
        JNIEnv *env,
        jobject obj,
        jint instanceNumber,
        jobject provider,
        jdoubleArray jPosition,
        jdoubleArray jVelocity,
        jdoubleArray jAction,
        jdoubleArray jGradient,
        jdoubleArray jMomentum,
        jdouble time) {

    JniCallback callback(env, provider);

    jboolean isPositionCopy, isVelocityCopy, isActionCopy, isGradientCopy, isMomentumCopy;

    double *position = env->GetDoubleArrayElements(jPosition, &isPositionCopy);
    double *velocity = env->GetDoubleArrayElements(jVelocity, &isVelocityCopy);
    double *action = env->GetDoubleArrayElements(jAction, &isActionCopy);
    double *gradient = env->GetDoubleArrayElements(jGradient, &isGradientCopy);
    double *momentum = env->GetDoubleArrayElements(jMomentum, &isMomentumCopy);

    // TODO Check all dimensions

    jsize dim = env->GetArrayLength(jPosition);

    auto returnValue = implementation[instanceNumber]->operate(
            std::span<double>(position, dim), std::span<double>(velocity, dim),
            std::span<double>(action, dim), std::span<double>(gradient, dim), std::span<double>(momentum, dim),
            time, callback);

    auto release = [&](jdoubleArray parent, double *child, jboolean isCopy, jint mode) {
        if (isCopy == JNI_TRUE) {
            env->ReleaseDoubleArrayElements(parent, child, mode);
        }
    };

    release(jPosition, position, isPositionCopy, 0); // or JNI_ABORT to avoid copy-back
    release(jVelocity, velocity, isVelocityCopy, 0);
    release(jAction, action, isActionCopy, 0);
    release(jGradient, gradient, isGradientCopy, 0);
    release(jMomentum, momentum, isMomentumCopy, 0);

    return returnValue;
}

#define CRIT

class JniCriticalHandler {
public:
    JniCriticalHandler(JNIEnv *env,
                       int dim,
                       const jdoubleArray jPosition,
                       const jdoubleArray jVelocity,
                       const jdoubleArray jAction,
                       const jdoubleArray jGradient,
                       const jdoubleArray jMomentum)
            : env(env), dim(dim), jArray({jPosition, jVelocity, jAction, jGradient, jMomentum}),
              array(jArray.size()), isCopy(jArray.size()) {

        for (int i = 0; i < jArray.size(); ++i) {
#ifdef CRIT
            array[i] = (double *) env->GetPrimitiveArrayCritical(jArray[i], &isCopy[i]);
#else
            array[i] = env->GetDoubleArrayElements(jArray[i], &isCopy[i]);
#endif
        }
    }

    ~JniCriticalHandler() {
        for (int i = 0; i < jArray.size(); ++i) {
#ifdef CRIT
            env->ReleasePrimitiveArrayCritical(jArray[i], (void *) array[i], JNI_ABORT);
#else
            if (isCopy[i] == JNI_TRUE) {
                env->ReleaseDoubleArrayElements(jArray[i], array[i], JNI_ABORT);
            }
#endif
        }
    }

    double* getArray(int i) { return array[i]; }

    std::span<double> getSpan(int i) { return std::span<double>(array[i], dim); }

    int getSize() const { return dim; }

private:

    JNIEnv *env;
    const int dim;

    std::vector<jdoubleArray> jArray;
    std::vector<double*> array;
    std::vector<jboolean> isCopy;
};

JNIEXPORT jobject JNICALL Java_dr_evomodel_operators_NativeZigZag_getNextEvent
        (JNIEnv *env, jobject obj,
                jint instanceNumber,
                jdoubleArray jPosition,
                jdoubleArray jVelocity,
                jdoubleArray jAction,
                jdoubleArray jGradient,
                jdoubleArray jMomentum) {

    JniCriticalHandler handler(env, env->GetArrayLength(jPosition),
            jPosition, jVelocity, jAction, jGradient, jMomentum);

    auto firstBounce = implementation[instanceNumber]->getNextBounce(
            handler.getSpan(0), handler.getSpan(1),
            handler.getSpan(2), handler.getSpan(3), handler.getSpan(4));

    return env->NewObject(classMTL, cid, firstBounce.time, firstBounce.index, firstBounce.type);
}

std::unique_ptr<JniCriticalHandler> handler;

JNIEXPORT jint JNICALL Java_dr_evomodel_operators_NativeZigZag_enterCriticalRegion
        (JNIEnv *env, jobject obj,
         jint instanceNumber,
         jdoubleArray jPosition,
         jdoubleArray jVelocity,
         jdoubleArray jAction,
         jdoubleArray jGradient,
         jdoubleArray jMomentum) {

    handler = zz::make_unique<JniCriticalHandler>(env,
                                                  env->GetArrayLength(jPosition),
                                                  jPosition, jVelocity, jAction, jGradient, jMomentum);

    return 0;
}

JNIEXPORT jint JNICALL Java_dr_evomodel_operators_NativeZigZag_exitCriticalRegion
        (JNIEnv *env, jobject obj, jint instanceNumber) {
    handler = nullptr;
    return 0;
}

JNIEXPORT jobject JNICALL Java_dr_evomodel_operators_NativeZigZag_getNextEventInCriticalRegion
        (JNIEnv *env, jobject obj, jint instanceNumber) {

    auto firstBounce = implementation[instanceNumber]->getNextBounce(
            handler->getSpan(0), handler->getSpan(1),
            handler->getSpan(2), handler->getSpan(3), handler->getSpan(4));

   return env->NewObject(classMTL, cid, firstBounce.time, firstBounce.index, firstBounce.type);
}
