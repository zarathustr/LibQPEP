message("Compiler version is: ${CMAKE_SYSTEM_VERSION}")
message("System version is: ${CMAKE_SYSTEM_NAME}")
message("Architecture is: ${CMAKE_SYSTEM_PROCESSOR}")

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    message("No Build Type Specified")
    set(CMAKE_BUILD_TYPE "Release")
endif()
message("CMake Build Type: ${CMAKE_BUILD_TYPE}")

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
    elseif(OSX_10_9)
        add_definitions(-D OSX_10_9)
        set(CMAKE_EXE_LINKER_FLAGS "-framework CoreGraphics")
    endif()
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math")
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
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math -lrt -lm -lz -g")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -Ofast -ffast-math -lrt -lm -lz -g")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lstdc++")
    endif()

    if(NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "armhf" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" AND NOT "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "arm64")
        include_directories(/opt/intel/mkl/include /opt/intel/oneapi/mkl/latest/include)
        link_directories(
                /opt/intel/mkl/lib/intel64
                /opt/intel/lib/intel64
                /opt/intel/oneapi/mkl/latest/lib/intel64
                /opt/intel/oneapi/intelpython/python3.7/lib/)
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

message("Eigen: ${EIGEN3_INCLUDE_DIRS}")

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

add_definitions(-DBIT64 -DEIGEN_DONT_PARALLELIZE)