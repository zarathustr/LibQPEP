cmake_minimum_required(VERSION 3.6)
project(LibQPEP)

message("Compiler version is: ${CMAKE_SYSTEM_VERSION}")
message("System version is: ${CMAKE_SYSTEM_NAME}")
message("Architecture is: ${CMAKE_SYSTEM_PROCESSOR}")
set(NO_OPENMP FALSE)
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    message("No Build Type Specified")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    set(CMAKE_CXX_COMPILER clang++)
    set(CMAKE_C_COMPILER clang)
    add_definitions(-D COMPILER_CLANG)
    add_definitions(-D USE_DARWIN)

    if("${CMAKE_SYSTEM_VERSION}" GREATER "20.0.0")
        message("Using OS X Big Sur")
        set(OSX_BIGSUR True)
    elseif("${CMAKE_SYSTEM_VERSION}" GREATER "14.6.0" AND "${CMAKE_SYSTEM_VERSION}" LESS "20.0.0")
        message("Detected Mac OS X")
        set(OSX True)
    else()
        message("Detected old Mac OS X system, while 10.7 is the lowest version supported")
        set(OSX_10_9 True)
    endif()

    if(OSX_BIGSUR)
        add_definitions(-D OSX_BIG_SUR)
	    set(CMAKE_EXE_LINKER_FLAGS "-framework CoreGraphics")
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -Ofast -ffast-math -pthread -g")
    elseif(OSX_10_9)
        add_definitions(-D NO_OMP)
        add_definitions(-D OSX_10_9)
        set(CMAKE_EXE_LINKER_FLAGS "-framework CoreGraphics")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math -pthread -stdlib=libc++ -g")
        message("Operating system is Mac OS X 10.9 or below, disabling OpenMP")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -Ofast -ffast-math -pthread -g")
    endif()
endif()

if(NO_OMP)
    add_definitions(-DNO_OMP)
endif()

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-D COMPILER_GCC)
    set(USE_CLANG False)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_definitions(-D COMPILER_CLANG)
    set(USE_CLANG True)
    message("Using Clang compiler.")
endif()

if("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    set(LINUX_OS True)
endif()

if(LINUX_OS)
    if(NOT USE_CLANG)
    	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math -pthread -lpthread -lrt -lm -lz -g")
    else()
    	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -Ofast -ffast-math -pthread -lpthread -lrt -lm -lz -g")
    endif()

    if(NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "armhf" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "arm64")
    	include_directories(/opt/intel/mkl/include /opt/intel/oneapi/mkl/latest/include)
    	link_directories(
    		/opt/intel/mkl/lib/intel64 
    		/opt/intel/lib/intel64
    		/opt/intel/oneapi/mkl/latest/lib/intel64
    		/opt/intel/oneapi/intelpython/python3.7/lib/)
	set(MKL_LIBS "-lmkl_gf_lp64 -lmkl_intel_thread -lmkl_gnu_thread -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -liomp5")
    endif()
endif()

add_definitions(-DCURRENT_SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_CXX_STANDARD 11)
find_package(Eigen3 QUIET)
if(NOT EIGEN3_FOUND)
  # Fallback to cmake_modules
  #find_package(cmake_modules REQUIRED)
  find_package(Eigen REQUIRED)
  set(EIGEN3_INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS})
  set(EIGEN3_LIBRARIES ${EIGEN_LIBRARIES})  # Not strictly necessary as Eigen is head only
  # Possibly map additional variables to the EIGEN3_ prefix.
else()
  set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
endif()

if(LINUX_OS AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "armhf" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "arm64")
    find_package(CUDA QUIET)
    if(CUDA_FOUND)
        message("Found CUDA in ${CUDA_TOOLKIT_ROOT_DIR}!")
        message("Found OpenCL in ${CUDA_OpenCL_LIBRARY}")
        message("Found CUBLAS Libs in ${CUDA_CUBLAS_LIBRARIES}")
        add_definitions(-DUSE_CUDA)
        set(OPENCL_LIBRARY ${CUDA_OpenCL_LIBRARY})
        set(OPENCL_INCLUDE_DIR "${CUDA_TOOLKIT_ROOT_DIR}/include")
    endif()
endif()

if(NOT NO_OPENCL)
    find_package(OpenCL QUIET)
    if(OPENCL_FOUND)
        find_package(ViennaCL QUIET)
        if(VIENNACL_FOUND)
            message("Found ViennaCL!")
            include_directories(${VIENNACL_INCLUDE_DIR})
            add_definitions(-DUSE_OPENCL)
        else()
            message(WARNING "The ViennaCL Library is not installed! Will disable OpenCL features!")
        endif()
    endif()
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(X11)
if(NOT DEFINED OPENCV_VER)
    find_package(OpenCV)
else()
    find_package(OpenCV ${OPENCV_VER} EXACT REQUIRED PATHS ${OPENCV_PATH})
endif()

find_package(superlu QUIET)
if(SUPERLU_FOUND)
    add_definitions(-DUSE_SUPERLU)
    message("Using SuperLU for Acceleration of Linear System Solving")
endif()

include_directories(${EIGEN3_INCLUDE_DIRS} ./csdp/include)

if(OpenCV_FOUND)
    include_directories(${X11_INCLUDE_DIR})

    if("${OpenCV_VERSION}" GREATER_EQUAL "4")
        include_directories(${OpenCV_INCLUDE_DIRS})
        add_definitions(-D USE_OPENCV)
        add_definitions(-D USE_OPENCV4)
    elseif("${OpenCV_VERSION}" GREATER_EQUAL "3")
        include_directories(${OpenCV_INCLUDE_DIRS})
        add_definitions(-D USE_OPENCV)
        add_definitions(-D USE_OPENCV3)
    else()
        message(FATAL_ERROR "OpenCV Version: ${OPENCV_VERSION} is too old for this project!")
    endif()
endif()

if(LAPACK_FOUND AND NOT NO_MKL)
    add_definitions(-D EIGEN_USE_BLAS)
else()
    message("Disabling MKL and LAPACK")
endif()

add_definitions(-D BIT64)
if(NOT OSX_10_9 AND NOT LINUX_OS AND NOT NO_OPENMP)
    add_definitions(-D USEOPENMP)
    add_definitions(-D SETNUMTHREADS)
    add_definitions(-D USESIGTERM)
    add_definitions(-D USEGETTIME)
endif()

if(OSX OR OSX_BIGSUR)
    find_package(OpenMP QUIET)
    if(NOT NO_OPENMP OR NOT OPENMP_FOUND)
    	if(EXISTS /usr/local/lib/libomp.dylib)
        		set(OMP_LIB /usr/local/lib/libomp.dylib)
        elseif(EXISTS /usr/local/lib/libomp.a)
        		set(OMP_LIB /usr/local/lib/libomp.a)
	elseif(EXISTS /opt/local/lib/libomp.dylib)
        		set(OMP_LIB /opt/local/lib/libomp.dylib)
        		include_directories(/opt/local/include /opt/local/include/libomp)
    	elseif(EXISTS /opt/local/lib/libomp.a)
        		set(OMP_LIB /opt/local/lib/libomp.a)
        		include_directories(/opt/local/include /opt/local/include/libomp)
	elseif(EXISTS /usr/lib/libomp.dylib)
        		set(OMP_LIB /usr/lib/libomp.dylib)
    	elseif(EXISTS /usr/lib/libomp.a)
        		set(OMP_LIB /usr/lib/libomp.a)
    	endif()
         message("OpenMP: ${OMP_LIB}")
    endif()
elseif(LINUX_OS AND NOT USE_CLANG)
    find_package(OpenMP REQUIRED)
    set(OMP_LIB ${OpenMP_CXX_FLAGS})
endif()


set(srcs
        generateProjectedPoints.cpp
        pnp_WQD.cpp
        solver_WQ_1_2_3_4_5_9_13_17_33_49_approx.cpp
        QPEP_grobner.cpp
        misc_pnp_funcs.cpp
        QPEP_lm_single.cpp
        misc_pTop_funcs.cpp
        pTop_WQD.cpp
        solver_WQ_approx_helper.cpp
        solver_WQ_1_2_3_4_5_9_13_17_33_49_approx_helper.cpp
        solver_WQ_approx.cpp
        solver_WQ.cpp
        solver_WQ_helper.cpp
        CVLib.cpp
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
        CovEstimation.cpp
        utils.cpp
        )

set(incs
        generateProjectedPoints.h
        pnp_WQD.h
        solver_WQ_1_2_3_4_5_9_13_17_33_49_approx.h
        solver_WQ_1_2_3_4_5_9_13_17_33_49_approx_helper.h
        QPEP_grobner.h
        utils.h
        misc_pnp_funcs.h
        QPEP_lm_single.h
        pTop_WQD.h
        misc_pTop_funcs.h
        solver_WQ_approx_helper.h
        solver_WQ_approx.h
        QPEP.h
        solver_WQ.h
        solver_WQ_helper.h
        CVLib.h
        CovEstimation.h)
include_directories(fsolve_wrapper)

add_executable(LibQPEP-test main.cpp)
add_library (LibQPEP SHARED ${srcs})

if(NOT LINUX_OS)
    find_library(CORE_GRAPHICS CoreGraphics QUIET)
    target_link_libraries(LibQPEP ${OMP_LIB} ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${CORE_GRAPHICS})
    target_link_libraries(LibQPEP-test LibQPEP ${OMP_LIB} ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES})
else()
    if(NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "armhf" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "arm64")
    	target_link_libraries(LibQPEP ${OMP_LIB} ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${MKL_LIBS} -lX11)
    	target_link_libraries(LibQPEP-test LibQPEP ${OMP_LIB} ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${MKL_LIBS} -lX11)
    else()
        target_link_libraries(LibQPEP ${OMP_LIB} ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${OMP_LIB} ${OPENCL_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} -lX11)
        target_link_libraries(LibQPEP-test LibQPEP ${OMP_LIB} ${OpenCV_LIBS} ${LAPACK_LIBRARIES} ${X11_LIBRARIES} ${BLAS_LIBRARIES} ${OMP_LIB} ${OPENCL_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} -lX11)
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

