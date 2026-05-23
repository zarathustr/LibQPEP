cmake_minimum_required(VERSION 3.6)
project(LibQPEP)

include(LibQPEP-basics.cmake)


set(srcs
        src/generateProjectedPoints.cpp
        src/pnp_WQD.cpp
        src/solver_WQ_1_2_3_4_5_9_13_17_33_49_approx.cpp
        src/QPEP_grobner.cpp
        src/misc_pnp_funcs.cpp
        src/QPEP_lm_single.cpp
        src/misc_pTop_funcs.cpp
        src/misc_hand_eye_funcs.cpp
        src/pTop_WQD.cpp
        src/hand_eye_WQD.cpp
        src/hand_eye_small_rotation.cpp
        src/pTop_small_rotation.cpp
        src/solver_WQ_approx_helper.cpp
        src/solver_WQ_1_2_3_4_5_9_13_17_33_49_approx_helper.cpp
        src/solver_WQ_approx.cpp
        src/solver_WQ.cpp
        src/solver_WQ_helper.cpp
        src/CVLib.cpp
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
        ./csdp/lib/writesol.c
        src/CovEstimation.cpp
        src/utils.cpp
        )

set(incs
        include/generateProjectedPoints.h
        include/pnp_WQD.h
        include/solver_WQ_1_2_3_4_5_9_13_17_33_49_approx.h
        include/solver_WQ_1_2_3_4_5_9_13_17_33_49_approx_helper.h
        include/QPEP_grobner.h
        include/utils.h
        include/misc_pnp_funcs.h
        include/QPEP_lm_single.h
        include/pTop_WQD.h
        include/hand_eye_WQD.h
        include/hand_eye_small_rotation.h
        include/pTop_small_rotation.h
        include/misc_pTop_funcs.h
        include/misc_hand_eye_funcs.h
        include/solver_WQ_approx_helper.h
        include/solver_WQ_approx.h
        include/QPEP.h
        include/solver_WQ.h
        include/solver_WQ_helper.h
        include/CVLib.h
        include/CovEstimation.h)
include_directories(include fsolve_wrapper)

add_executable(LibQPEP-test src/main.cpp)
add_library (LibQPEP SHARED ${srcs})
add_executable(LibQPEP-benchmark-linear-solvers src/benchmark_linear_solvers.cpp)

if(LINUX_OS)
    if(NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "armhf" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "armv7l" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "arm64")
        target_compile_options(LibQPEP PRIVATE -msse -msse2 -msse3 -msse4 -msse4.1 -msse4.2)
    elseif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "armhf" OR "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "armv7l")
        target_compile_options(LibQPEP PRIVATE -mfpu=neon -mfloat-abi=hard)
    elseif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" OR "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "arm64")
        target_compile_options(LibQPEP PRIVATE -march=native)
    endif()
endif()

if(NOT LINUX_OS)
    find_library(CORE_GRAPHICS CoreGraphics QUIET)
    target_link_libraries(LibQPEP ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${CORE_GRAPHICS})
    target_link_libraries(LibQPEP-test LibQPEP ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES})
    target_link_libraries(LibQPEP-benchmark-linear-solvers LibQPEP ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES})
else()
    if(NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "armhf" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "arm64")
    	target_link_libraries(LibQPEP ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${MKL_LIBS} -lX11)
    	target_link_libraries(LibQPEP-test LibQPEP ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${MKL_LIBS} -lX11)
        target_link_libraries(LibQPEP-benchmark-linear-solvers LibQPEP ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${MKL_LIBS} -lX11)
    else()
        target_link_libraries(LibQPEP ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${OPENCL_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} -lX11)
        target_link_libraries(LibQPEP-test LibQPEP ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${OPENCL_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} -lX11)
        target_link_libraries(LibQPEP-benchmark-linear-solvers LibQPEP ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${OPENCL_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} -lX11)
    endif()
endif()

include(GNUInstallDirs)
install(TARGETS LibQPEP
        EXPORT ${PROJECT_NAME}-config
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/../lib/
        )

if(LINUX_OS)
	INSTALL(CODE "execute_process( \
    		COMMAND ${CMAKE_COMMAND} -E create_symlink \
   		 ${CMAKE_INSTALL_PREFIX}/lib/libLibQPEP.so \
    		${CMAKE_INSTALL_PREFIX}/lib/libLibQPEP.a   \
    		)"
	)
endif()
install(FILES ${incs} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

export(
        TARGETS LibQPEP
        FILE ${PROJECT_NAME}-config.cmake
)
install(
        EXPORT ${PROJECT_NAME}-config
        DESTINATION share/${PROJECT_NAME}/cmake
)


# Create variable for the local build tree
get_property(CONF_INCLUDE_DIRS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)

# Create variable with the library location
set(QPEP_LIB ${PROJECT_NAME})

# Configure config file for local build tree
configure_file(${PROJECT_NAME}Config.cmake.in
        "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY)

# 2- installation build #

# Change the include location for the case of an install location
set(CONF_INCLUDE_DIRS ${INSTALL_INCLUDE_DIR} ${CONF_INCLUDE_DIRS} )

configure_file(${PROJECT_NAME}Config.cmake.in
        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}Config.cmake" @ONLY)

install (
        FILES
        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}Config.cmake"
        DESTINATION "share/${PROJECT_NAME}/cmake/"
)
