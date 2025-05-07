cmake_minimum_required(VERSION 3.0)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
project(zig_zag)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Base optimizations (always enabled)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

# Conditional SIMD (enabled only if supported)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx" COMPILER_SUPPORTS_AVX)
check_cxx_compiler_flag("-msse4.2" COMPILER_SUPPORTS_SSE42)
check_cxx_compiler_flag("-mfma" COMPILER_SUPPORTS_FMA)

if(COMPILER_SUPPORTS_AVX AND COMPILER_SUPPORTS_FMA)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -mfma")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mfma")
endif()
if(COMPILER_SUPPORTS_SSE42)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.2")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
endif()

# Debug flags (no SIMD)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")

find_package(Java REQUIRED)
find_package(JNI REQUIRED)
include(UseJava)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

OPTION(USE_SYSTEM_TBB "Use system TBB" ON)
IF(USE_SYSTEM_TBB)
    FIND_PACKAGE(TBB REQUIRED)
ELSE()
    MESSAGE(STATUS "Compiling local TBB")
    if (NO_RTM)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-rtm")
    endif()

    SET(TBB_BUILD_STATIC OFF CACHE BOOL         "Build TBB static library")
    SET(TBB_BUILD_TBBMALLOC OFF CACHE BOOL      "Build TBB malloc library")
    SET(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL "Build TBB malloc proxy library")
    SET(TBB_BUILD_TESTS OFF CACHE BOOL           "Build TBB tests and enable testing infrastructure")

    ADD_SUBDIRECTORY(external/tbb)
    SET(TBB_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/external/tbb/include")
    SET(TBB_LIBRARIES tbb)
    #INCLUDE("${CMAKE_MODULE_PATH}/build_tbb.cmake")
ENDIF()
INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
MESSAGE(STATUS "TBB directory: ${TBB_INCLUDE_DIRS}")

option(USE_SYSTEM_XSIMD "Use system XSIMD" OFF)
if(USE_SYSTEM_XSIMD)
    find_package(xsimd REQUIRED)
else()
    include("${CMAKE_MODULE_PATH}/build_xsimd.cmake")
endif()
message(STATUS "xsimd directory: ${xsimd_INCLUDE_DIRS}")
include_directories(${xsimd_INCLUDE_DIRS})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
message(STATUS, "Eigen 3 directory: ${}")

add_library (zig_zag SHARED
        src/NativeZigZag.cpp
        src/AbstractZigZag.h src/ZigZag.h src/NoUTurn.h src/Timing.h
        src/aligned_allocator.h src/MemoryManagement.h src/Simd.h src/PrecisionColumn.h
        src/threefry.h src/NutsTreeState.h)

target_link_libraries (zig_zag ${JAVA_JVM_LIBRARY} ${TBB_LIBRARIES})
target_include_directories(zig_zag PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${JNI_INCLUDE_DIRS}
)

target_link_libraries(zig_zag Eigen3::Eigen)

add_executable(main
        extra/main.cpp)

target_link_libraries(main
        zig_zag)

#add_jar(${PROJECT_NAME}
#        DemoFuncs.java
#        Main.java
#        Name.java
#    ENTRY_POINT Main
#)
