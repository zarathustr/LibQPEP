cmake_minimum_required(VERSION 3.17)
project(LibQPEP)

#set(CMAKE_CXX_COMPILER g++-10)

set(CMAKE_BUILD_TYPE "Debug")

if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math -Xclang -fopenmp -pthread")
elseif(LINUX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math -fopenmp -pthread")
endif()

set(CMAKE_CXX_STANDARD 11)
find_package(Eigen3 REQUIRED)
find_package(OpenCV)
find_package(LAPACK)

include_directories(${EIGEN3_INCLUDE_DIR} ./csdp/include)

if(OpenCV_FOUND)
    include_directories(${OpenCV_INCLUDE_DIRS})
    add_definitions(-D USE_OPENCV)
endif()

if(LAPACK_FOUND)
    add_definitions(-D EIGEN_USE_BLAS)
endif()

add_definitions(-D BIT64)
add_definitions(-D USEOPENMP)
add_definitions(-D SETNUMTHREADS)
add_definitions(-D USESIGTERM)
add_definitions(-D USEGETTIME)
add_definitions(-D CAPSLAPACK)

set(OMP_LIB "")
if(APPLE OR LINUX)
    if(EXISTS /usr/local/lib/libomp.a)
        set(OMP_LIB /usr/local/lib/libomp.a)
    elseif(EXISTS /usr/lib/libomp.a)
        set(OMP_LIB /usr/lib/libomp.a)
    else()
        message(FATAL_ERROR "No OpenMP Library Found!")
    endif()
endif()


add_executable(LibQPEP
        main.cpp
        generateProjectedPoints.cpp
        generateProjectedPoints.h
        pnp_WQD.cpp
        pnp_WQD.h
        solver_WQ_1_2_3_4_5_9_13_17_33_49_approx.cpp
        solver_WQ_1_2_3_4_5_9_13_17_33_49_approx.h
        solver_WQ_1_2_3_4_5_9_13_17_33_49_approx_helper.cpp
        solver_WQ_1_2_3_4_5_9_13_17_33_49_approx_helper.h
        QPEP_grobner.h
        QPEP_grobner.cpp
        utils.h
        misc_pnp_funcs.h
        misc_pnp_funcs.cpp
        QPEP_lm_single.cpp
        QPEP_lm_single.h
        pTop_WQD.cpp
        pTop_WQD.h
        misc_pTop_funcs.cpp
        misc_pTop_funcs.h
        solver_WQ_approx_helper.cpp
        solver_WQ_approx_helper.h
        solver_WQ_approx.cpp
        solver_WQ_approx.h
        QPEP.h
        solver_WQ.cpp
        solver_WQ.h
        solver_WQ_helper.cpp
        solver_WQ_helper.h
        CVLib.cpp
        CVLib.h
        ./csdp/lib/Fnorm.c
        ./csdp/lib/chol.c
        ./csdp/lib/linesearch.c
        ./csdp/lib/norms.c
        ./csdp/lib/qreig.c
        ./csdp/lib/sym_mat.c
        ./csdp/lib/zero_mat.c
        ./csdp/lib/add_mat.c
        ./csdp/lib/copy_mat.c
        ./csdp/lib/make_i.c
        ./csdp/lib/op_a.c
        ./csdp/lib/readprob.c
        ./csdp/lib/trace_prod.c
        ./csdp/lib/addscaledmat.c
        ./csdp/lib/easysdp.c
        ./csdp/lib/makefill.c
        ./csdp/lib/op_at.c
        ./csdp/lib/readsol.c
        ./csdp/lib/tweakgap.c
        ./csdp/lib/allocmat.c
        ./csdp/lib/freeprob.c
        ./csdp/lib/mat_mult.c
        ./csdp/lib/op_o.c
        ./csdp/lib/sdp.c
        ./csdp/lib/user_exit.c
        ./csdp/lib/calc_dobj.c
        ./csdp/lib/initparams.c
        ./csdp/lib/mat_multsp.c
        ./csdp/lib/packed.c
        ./csdp/lib/solvesys.c
        ./csdp/lib/writeprob.c
        ./csdp/lib/calc_pobj.c
        ./csdp/lib/initsoln.c
        ./csdp/lib/matvec.c
        ./csdp/lib/psd_feas.c
        ./csdp/lib/sortentries.c
        ./csdp/lib/writesol.c CovEstimation.cpp CovEstimation.h)

target_link_libraries(LibQPEP ${OMP_LIB} ${OpenCV_LIBS} ${LAPACK_LIBRARIES})

