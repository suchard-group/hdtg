#include <iostream>
#include <vector>
#include <memory>

#include "dr_evomodel_operators_NativeZigZag.h"
#include "ZigZag.hpp"

std::vector<std::unique_ptr<zz::AbstractZigZag>> implementation;

static jclass classMTL;
static jmethodID cid;
static jclass classNZZO;
static jfieldID flagsFid;
static jfieldID seedFid;
static jfieldID infoFid;

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

    jclass tmpClassNZZO = env->FindClass("dr/evomodel/operators/NativeZigZagOptions");
    classNZZO = (jclass) env->NewGlobalRef(tmpClassNZZO);
    env->DeleteLocalRef(tmpClassNZZO);

    flagsFid = env->GetFieldID(classNZZO, "flags", "J");
    seedFid = env->GetFieldID(classNZZO, "seed", "J");
    infoFid = env->GetFieldID(classNZZO, "info", "I");

    return JNI_VERSION_1_6;
}

JNIEXPORT void JNICALL JNI_OnUnload(JavaVM *vm, void *reserved) {
    JNIEnv *env;

    (void)reserved;

    if (vm->GetEnv(reinterpret_cast<void**>(&env), JNI_VERSION_1_6) != JNI_OK)
        return;

    env->DeleteGlobalRef(classMTL);
    env->DeleteGlobalRef(classNZZO);
}

// std::unique_ptr<zz::AbstractZigZag> dispatch(
//         int dimension,
//         double *rawMask,
//         double *rawObserved,
//         long flags,
//         int info,
//         long seed) {
// 
//     if (static_cast<unsigned long>(flags) & zz::Flags::AVX) {
//         std::cerr << "Factory: AVX" << std::endl;
//         return zz::make_unique<zz::ZigZag<zz::DoubleAvxTypeInfo>>(
//                 dimension, rawMask, rawObserved, flags, info, seed);
//     } else if (static_cast<unsigned long>(flags) & zz::Flags::SSE) {
//         std::cerr << "Factory: SSE" << std::endl;
//         return zz::make_unique<zz::ZigZag<zz::DoubleSseTypeInfo>>(
//                 dimension, rawMask, rawObserved, flags, info, seed);
//     } else {
//         std::cerr << "Factory: No SIMD" << std::endl;
//         return zz::make_unique<zz::ZigZag<zz::DoubleNoSimdTypeInfo>>(
//                 dimension, rawMask, rawObserved, flags, info, seed);
//     }
// }

JNIEXPORT jint JNICALL Java_dr_evomodel_operators_NativeZigZag_create(
        JNIEnv *env,
        jobject obj,
        jint dimension,
        jobject options,
        jdoubleArray mask,
        jdoubleArray observed,
        jdoubleArray parameterSign) {

    (void)obj;

    double *rawMask = env->GetDoubleArrayElements(mask, nullptr);
    double *rawObserved = env->GetDoubleArrayElements(observed, nullptr);
    double *rawParameterSign = env->GetDoubleArrayElements(parameterSign, nullptr);

    long flags = env->GetLongField(options, flagsFid);
    long seed = env->GetLongField(options, seedFid);
    int info = env->GetIntField(options, infoFid);

    int instanceNumber = static_cast<int>(implementation.size());
    implementation.emplace_back(zz::dispatch(dimension, rawMask, rawObserved, rawParameterSign, flags, info, seed));

    env->ReleaseDoubleArrayElements(mask, rawMask, JNI_ABORT);
    env->ReleaseDoubleArrayElements(observed, rawObserved, JNI_ABORT);
    env->ReleaseDoubleArrayElements(parameterSign, rawParameterSign, JNI_ABORT);

return instanceNumber;
}


class JniCallback : public zz::PrecisionColumnCallback {
public:
    JniCallback(JNIEnv *env, jobject provider)
            : PrecisionColumnCallback(), env(env), provider(provider), array(nullptr), data(nullptr),
              column(-1), isCopy(static_cast<jboolean>(false)) {
        jclass cls = env->GetObjectClass(provider);
        mid = env->GetMethodID(cls, "getColumn", "(I)[D");
        if (mid == nullptr) {
            throw;
        }
    }

    ~JniCallback() { releaseColumn(); }

    const double* getColumn(int index) override {

        if (index != column) {

            releaseColumn();

            jobject object = env->CallObjectMethod(provider, mid, index);
            array = reinterpret_cast<jdoubleArray *>(&object);
            data = env->GetDoubleArrayElements(*array, &isCopy);
        }

        column = index;
        return data;
    }

    void releaseColumn() override {

        if (data != nullptr) {
            env->ReleaseDoubleArrayElements(*array, data, JNI_ABORT);
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

    (void)obj;

    JniCallback callback(env, provider);

    jboolean isPositionCopy, isVelocityCopy, isActionCopy, isGradientCopy, isMomentumCopy;

    double *position = env->GetDoubleArrayElements(jPosition, &isPositionCopy);
    double *velocity = env->GetDoubleArrayElements(jVelocity, &isVelocityCopy);
    double *action = env->GetDoubleArrayElements(jAction, &isActionCopy);
    double *gradient = env->GetDoubleArrayElements(jGradient, &isGradientCopy);
    double *momentum = env->GetDoubleArrayElements(jMomentum, &isMomentumCopy);

    // TODO Check all dimensions

    const auto dim = static_cast<size_t>(env->GetArrayLength(jPosition));

    implementation[instanceNumber]->operate(
            zz::DblSpan(position, dim), zz::DblSpan(velocity, dim),
            zz::DblSpan(action, dim), zz::DblSpan(gradient, dim), zz::DblSpan(momentum, dim),
            time, callback);

    auto release = [&](jdoubleArray parent, double *child, jboolean isCopy, jint mode) {
//        if (isCopy == JNI_TRUE) {
            env->ReleaseDoubleArrayElements(parent, child, mode);
//        }
    };

    release(jPosition, position, isPositionCopy, 0); // or JNI_ABORT to avoid copy-back
    release(jVelocity, velocity, isVelocityCopy, 0);
    release(jAction, action, isActionCopy, 0);
    release(jGradient, gradient, isGradientCopy, 0);
    release(jMomentum, momentum, isMomentumCopy, 0);

    return 0;
}

#define CRITICAL

class JniCriticalHandler {
public:
    JniCriticalHandler(JNIEnv *env,
                       int dim,
                       int releaseMode,
                       jdoubleArray jPosition,
                       jdoubleArray jVelocity,
                       jdoubleArray jAction,
                       jdoubleArray jGradient,
                       jdoubleArray jMomentum)
            : env(env), dim(static_cast<size_t>(dim)), releaseMode(releaseMode),
              jArray({jPosition, jVelocity, jAction, jGradient, jMomentum}),
              array(jArray.size()), isCopy(jArray.size()) {

        for (int i = 0; i < jArray.size(); ++i) {
            if (jArray[i] != nullptr) {
#ifdef CRITICAL

                array[i] = (double *) env->GetPrimitiveArrayCritical(jArray[i], &isCopy[i]);
#else
                array[i] = env->GetDoubleArrayElements(jArray[i], &isCopy[i]);
#endif
            } else {
                array[i] = nullptr;
            }
        }
    }

    ~JniCriticalHandler() {
        for (int i = 0; i < jArray.size(); ++i) {
            if (jArray[i] != nullptr) {
#ifdef CRITICAL
                env->ReleasePrimitiveArrayCritical(jArray[i], (void *) array[i], releaseMode);
#else
                env->ReleaseDoubleArrayElements(jArray[i], array[i], releaseMode);
#endif
            }
        }
    }

//    double* getArray(int i) { return array[i]; }

    zz::DblSpan getSpan(int i) { return zz::DblSpan(array[i], dim); }

//    size_t getSize() const { return dim; }

private:

    JNIEnv *env;
    const size_t dim;
    const int releaseMode;

    std::vector<jdoubleArray> jArray;
    std::vector<double*> array;
    std::vector<jboolean> isCopy;
};

JNIEXPORT void JNICALL Java_dr_evomodel_operators_NativeZigZag_innerBounce
        (JNIEnv *env, jobject obj, jint instanceNumber,
                jdoubleArray position,
                jdoubleArray velocity,
                jdoubleArray action,
                jdoubleArray gradient,
                jdoubleArray momentum,
                jdouble time, jint index, jint type) {
    (void)obj;

    JniCriticalHandler handler(env, env->GetArrayLength(position), 0,
            position, velocity, action, gradient, momentum);

    implementation[instanceNumber]->innerBounce(
            handler.getSpan(0), handler.getSpan(1),
            handler.getSpan(2), handler.getSpan(3), handler.getSpan(4),
            time, index, type);
}

JNIEXPORT void JNICALL Java_dr_evomodel_operators_NativeZigZag_updateDynamics
        (JNIEnv *env, jobject, jint instanceNumber,
                jdoubleArray position,
                jdoubleArray velocity,
                jdoubleArray action,
                jdoubleArray gradient,
                jdoubleArray momentum,
                jdoubleArray jColumn,
                jdouble time, jint index, jint) {

    auto dim = env->GetArrayLength(position);

    JniCriticalHandler handler(env, dim, 0,
            position, velocity, action, gradient, momentum);

    jboolean isColumnCopy;

    auto *column = (double *) env->GetPrimitiveArrayCritical(jColumn, &isColumnCopy);

    implementation[instanceNumber]->updateDynamics(
            handler.getSpan(0), handler.getSpan(1),
            handler.getSpan(2), handler.getSpan(3), handler.getSpan(4),
            zz::DblSpan(column, static_cast<std::size_t>(dim)),
            time, index);

    env->ReleasePrimitiveArrayCritical(jColumn, (void *) column, JNI_ABORT);
}

JNIEXPORT jobject JNICALL Java_dr_evomodel_operators_NativeZigZag_getNextEvent
        (JNIEnv *env, jobject obj,
                jint instanceNumber,
                jdoubleArray jPosition,
                jdoubleArray jVelocity,
                jdoubleArray jAction,
                jdoubleArray jGradient,
                jdoubleArray jMomentum) {

    (void)obj;

    JniCriticalHandler handler(env, env->GetArrayLength(jPosition), JNI_ABORT,
            jPosition, jVelocity, jAction, jGradient, jMomentum);

    auto firstBounce = implementation[instanceNumber]->getNextBounce(
            handler.getSpan(0), handler.getSpan(1),
            handler.getSpan(2), handler.getSpan(3), handler.getSpan(4));

    return env->NewObject(classMTL, cid, firstBounce.time, firstBounce.index, firstBounce.type);
}

JNIEXPORT jobject JNICALL Java_dr_evomodel_operators_NativeZigZag_getNextEventIrreversible
        (JNIEnv *env, jobject obj,
                jint instanceNumber,
                jdoubleArray jPosition,
                jdoubleArray jVelocity,
                jdoubleArray jAction,
                jdoubleArray jGradient) {
    (void)obj;

    JniCriticalHandler handler(env, env->GetArrayLength(jPosition), JNI_ABORT,
            jPosition, jVelocity, jAction, jGradient, nullptr);

    auto firstBounce = implementation[instanceNumber]->getNextBounceIrreversible(
            handler.getSpan(0), handler.getSpan(1),
            handler.getSpan(2), handler.getSpan(3));

    return env->NewObject(classMTL, cid, firstBounce.time, firstBounce.index, firstBounce.type);
}


std::unique_ptr<JniCriticalHandler> handler;

JNIEXPORT void JNICALL Java_dr_evomodel_operators_NativeZigZag_innerBounceCriticalRegion
        (JNIEnv *env, jobject obj, jint instanceNumber,
         jdouble time, jint index, jint type) {

    (void)env;
    (void)obj;

    implementation[instanceNumber]->innerBounce(
            handler->getSpan(0), handler->getSpan(1),
            handler->getSpan(2), handler->getSpan(3), handler->getSpan(4),
            time, index, type);
}

JNIEXPORT jint JNICALL Java_dr_evomodel_operators_NativeZigZag_enterCriticalRegion
        (JNIEnv *env, jobject obj,
         jint instanceNumber,
         jdoubleArray jPosition,
         jdoubleArray jVelocity,
         jdoubleArray jAction,
         jdoubleArray jGradient,
         jdoubleArray jMomentum) {

    (void)obj;
    (void)instanceNumber;

    handler = zz::make_unique<JniCriticalHandler>(env,
                                                  env->GetArrayLength(jPosition), JNI_ABORT,
                                                  jPosition, jVelocity, jAction, jGradient, jMomentum);

    return 0;
}

JNIEXPORT jint JNICALL Java_dr_evomodel_operators_NativeZigZag_exitCriticalRegion
        (JNIEnv *env, jobject obj, jint instanceNumber) {

    (void)env;
    (void)obj;
    (void)instanceNumber;

    handler = nullptr;
    return 0;
}

JNIEXPORT jboolean JNICALL Java_dr_evomodel_operators_NativeZigZag_inCriticalRegion
        (JNIEnv *env, jobject obj, jint instanceNumber) {

    (void)env;
    (void)obj;
    (void)instanceNumber;

    return (handler != nullptr);
}



JNIEXPORT jobject JNICALL Java_dr_evomodel_operators_NativeZigZag_getNextEventInCriticalRegion
        (JNIEnv *env, jobject obj, jint instanceNumber) {

    (void)obj;

    auto firstBounce = implementation[instanceNumber]->getNextBounce(
            handler->getSpan(0), handler->getSpan(1),
            handler->getSpan(2), handler->getSpan(3), handler->getSpan(4));

   return env->NewObject(classMTL, cid, firstBounce.time, firstBounce.index, firstBounce.type);
}
