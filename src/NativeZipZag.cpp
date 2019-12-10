#include <iostream>
#include <vector>
#include <memory>

#include "dr_evomodel_operators_NativeZigZag.h"
#include "ZigZag.hpp"

//int loaded = 0; // Indicates is the initial library constructors have been run
// This patches a bug with JVM under Linux that calls the finalizer twice

//#ifdef __GNUC__
//void __attribute__ ((constructor)) beagle_gnu_init(void) {
////    beagle_library_initialize();
//    std::cerr << "gnuInit()" << std::endl;
//}
//void __attribute__ ((destructor)) beagle_gnu_finalize(void) {
////    beagle_library_finalize();
//    std::cerr << "gnuFinalize()" << std::endl;
//}
//#endif
//
//#ifdef _WIN32
//BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD fdwReason, LPVOID lpvReserved) {
//    switch (fdwReason) {
//    case DLL_PROCESS_ATTACH:
////        beagle_library_initialize();
//        break;
//    case DLL_PROCESS_DETACH:
////        beagle_library_finalize();
//        break;
//    }
//    return TRUE;
//}
//#endif

std::vector<std::unique_ptr<zz::AbstractZigZag>> implementation;

//JNIEXPORT jint JNICALL JNI_OnLoad(JavaVM* vm, void* reserved) {
//
//    JNIEnv *env;
//    jclass  cls;
//    jint    res;
//
//    (void)reserved;
//
//    if (vm->GetEnv((void **)&env, JNI_VERSION_1_6) != JNI_OK) {
//        return JNI_ERR;
//    }
//
////    cls = (*env)->FindClass(env, JNIT_CLASS);
////    if (cls == NULL)
////        return JNI_ERR;
////
////    res = (*env)->RegisterNatives(env, cls, funcs, sizeof(funcs)/sizeof(*funcs));
////    if (res != 0)
////        return -1;
//
//
//    std::cerr << "onLoad()" << std::endl;
//
//    return JNI_VERSION_1_6;
//}

//JNIEXPORT void JNICALL JNI_OnUnload(JavaVM *vm, void *reserved) {
//    JNIEnv *env;
//    jclass  cls;
//
//    (void)reserved;
//
//    std::cerr << "onUnload()" << std::endl;
//
//    if (vm->GetEnv((void **)&env, JNI_VERSION_1_6) != JNI_OK)
//        return;
//
////    cls = (*env)->FindClass(env, JNIT_CLASS);
////    if (cls == NULL)
////        return;
////
////    (*env)->UnregisterNatives(env, cls);
//
//    std::cerr << "onUnload()" << std::endl;
//}

//template <typename L, typename R>
//double bitwise_and(const L lhs, const R rhs) {
//    return static_cast<double>(
//            static_cast<uint64_t>(lhs) & static_cast<uint64_t>(rhs)
//    );
//}

JNIEXPORT jint JNICALL Java_dr_evomodel_operators_NativeZigZag_create(
        JNIEnv *env,
        jobject obj,
        jint dimension,
        jobject provider,
        jdoubleArray mask,
        jdoubleArray observed) {

    std::cerr << "In create()" << std::endl;

    jclass cls = env->GetObjectClass(provider);
    jmethodID mid = env->GetMethodID(cls, "getColumn", "(I)[D");
    if (mid == 0) {
        return JNI_ERR;
    }

    std::cerr << "Trying callback ..." << std::endl;

    jobject mvdata = env->CallObjectMethod(provider, mid, 0);
    jdoubleArray *arr = reinterpret_cast<jdoubleArray*>(&mvdata);
    double *data = env->GetDoubleArrayElements(*arr, NULL);

    for (int i = 0; i < dimension; ++i) {
        std::cerr << " " << data[i];
    }
    std::cerr << std::endl;

    env->ReleaseDoubleArrayElements(*arr, data, JNI_ABORT);

    double *rawMask = env->GetDoubleArrayElements(mask, NULL);
    double *rawObserved = env->GetDoubleArrayElements(observed, NULL);

    long flags = zz::Flags::TBB;
    int nThreads = 2;

    int instanceNumber = implementation.size();
    implementation.emplace_back(zz::make_unique<
//            zz::ZigZag<zz::DoubleNoSimdTypeInfo>
zz::ZigZag<zz::DoubleSseTypeInfo>
    >(dimension, rawMask, rawObserved, mid, flags, nThreads));

    env->ReleaseDoubleArrayElements(mask, rawMask, JNI_ABORT);
    env->ReleaseDoubleArrayElements(observed, rawObserved, JNI_ABORT);


    std::cerr << "Success" << std::endl;

    return instanceNumber;
}

JNIEXPORT jint JNICALL Java_dr_evomodel_operators_NativeZigZag_operate(
        JNIEnv *env,
        jobject obj,
        jint instanceNumber,
        jdoubleArray jPosition,
        jdoubleArray jVelocity,
        jdoubleArray jAction,
        jdoubleArray jGradient,
        jdoubleArray jMomentum,
        jdouble time) {

    jboolean isPositionCopy, isVelocityCopy, isActionCopy, isGradientCopy, isMomentumCopy;

    double *position = env->GetDoubleArrayElements(jPosition, &isPositionCopy);
    double *velocity = env->GetDoubleArrayElements(jVelocity, &isVelocityCopy);
    double *action = env->GetDoubleArrayElements(jAction, &isActionCopy);
    double *gradient = env->GetDoubleArrayElements(jGradient, &isGradientCopy);
    double *momentum = env->GetDoubleArrayElements(jMomentum, &isMomentumCopy);

    // TODO Check all dimensions

    jsize dim = env->GetArrayLength(jPosition);

    implementation[instanceNumber]->operate(
            std::span<double>(position, dim), std::span<double>(velocity, dim),
            std::span<double>(action, dim), std::span<double>(gradient, dim), std::span<double>(momentum, dim),
            time);

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

    return 0;
}

JNIEXPORT jobject JNICALL Java_dr_evomodel_operators_NativeZigZag_getNextEvent
        (JNIEnv *env, jobject obj,
                jint instanceNumber,
                jdoubleArray jPosition,
                jdoubleArray jVelocity,
                jdoubleArray jAction,
                jdoubleArray jGradient,
                jdoubleArray jMomentum){

    jboolean isPositionCopy, isVelocityCopy, isActionCopy, isGradientCopy, isMomentumCopy;

    double *position = env->GetDoubleArrayElements(jPosition, &isPositionCopy);
    double *velocity = env->GetDoubleArrayElements(jVelocity, &isVelocityCopy);
    double *action = env->GetDoubleArrayElements(jAction, &isActionCopy);
    double *gradient = env->GetDoubleArrayElements(jGradient, &isGradientCopy);
    double *momentum = env->GetDoubleArrayElements(jMomentum, &isMomentumCopy);

    jsize dim = env->GetArrayLength(jPosition);

    auto firstBounce = implementation[instanceNumber]->getNextBounce(
            std::span<double>(position, dim), std::span<double>(velocity, dim),
            std::span<double>(action, dim), std::span<double>(gradient, dim), std::span<double>(momentum, dim));

    auto release = [&](jdoubleArray parent, double *child, jboolean isCopy, jint mode) {
        if (isCopy == JNI_TRUE) {
            env->ReleaseDoubleArrayElements(parent, child, mode);
        }
    };

    release(jPosition, position, isPositionCopy, JNI_ABORT);
    release(jVelocity, velocity, isVelocityCopy, JNI_ABORT);
    release(jAction, action, isActionCopy, JNI_ABORT);
    release(jGradient, gradient, isGradientCopy, JNI_ABORT);
    release(jMomentum, momentum, isMomentumCopy, JNI_ABORT);

    jclass classMTL = env->FindClass("dr/inference/operators/hmc/MinimumTravelInformation");

    jmethodID cid = env->GetMethodID(classMTL, "<init>", "(DII)V");

    return env->NewObject(classMTL, cid, firstBounce.time, firstBounce.index, firstBounce.type);
}
