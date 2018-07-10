cmake_minimum_required(VERSION 3.5)
project(da3d)

# Find Matlab
find_package(Matlab COMPONENTS MEX_COMPILER MX_LIBRARY REQUIRED)

# Are we using gcc?
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # GCC on MacOs needs this option to use the clang assembler
  if (APPLE)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wa,-q")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-q")
  endif ()
  # Optimize to the current CPU and enable warnings
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -Wall -Wextra")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -Wextra")
endif ()

# Enable OpenMP
find_package (OpenMP)
if(OPENMP_FOUND)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

# Find fftw3
find_path (FFTW_INCLUDES fftw3.h)
find_library (FFTWF_LIBRARIES NAMES fftw3f)

set(SOURCE_FILES
    bm3d.h     lib_transforms.h      utilities.h
    bm3d.cpp   lib_transforms.cpp    utilities.cpp    mex.cpp
    )

if (NOT WIN32)
   set_property(SOURCE iio.c PROPERTY COMPILE_FLAGS "-Wno-unused-variable -Wno-unused-parameter -Wno-pointer-sign -Wno-parentheses -Wno-deprecated-declarations -Wno-unused-function")
endif ()

matlab_add_mex(NAME bm3d ${SOURCE_FILES} LINK_TO ${FFTWF_LIBRARIES})
target_include_directories(bm3d PUBLIC ${FFTW_INCLUDES})