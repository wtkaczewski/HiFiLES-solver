cmake_minimum_required(VERSION 3.10)
project(HiFiLES LANGUAGES CXX)

include_directories("${PROJECT_SOURCE_DIR}/include")
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

## Build type ##
if(NOT CMAKE_BUILD_TYPE) # force set build type in cache to debug
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb -D_DEBUG" CACHE STRING "Flags used by the linker during debug builds" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3" CACHE STRING "Flags used by the linker during release builds" FORCE)

## Building options ##
set(BUILDTYPE "CPU" CACHE STRING "Build for 'CPU' or 'GPU'")
set(PARALLEL OFF CACHE BOOL "Build with MPI support")

#source file list
set(SRCLIST
src/HiFiLES.cpp
src/bdy_inters.cpp
src/cubature_1d.cpp
src/cubature_hexa.cpp
src/cubature_quad.cpp
src/cubature_tet.cpp
src/cubature_tri.cpp
src/eles.cpp
src/eles_hexas.cpp
src/eles_pris.cpp
src/eles_quads.cpp
src/eles_tets.cpp
src/eles_tris.cpp
src/flux.cpp
src/funcs.cpp
src/geometry.cpp
src/global.cpp
src/input.cpp
src/inters.cpp
src/int_inters.cpp
src/linear_solvers_structure.cpp
src/matrix_structure.cpp
src/mesh.cpp
src/output.cpp
src/solver.cpp
src/source.cpp
src/vector_structure.cpp)

#cpu options
if(${BUILDTYPE} STREQUAL "CPU" )
    add_definitions(-D_CPU)
    #BLAS
    set(BLAS "NO" CACHE STRING "Build with external BLAS support, 'OPENBLAS', 'MKL', 'CBLAS', 'ATLAS','ACCELERATE' or 'NO'")
    if (NOT(${BLAS} STREQUAL "NO"))# if use blas
          set(BLAS_INCLUDE "/opt" CACHE PATH "path to BLAS include")
          set(BLAS_LD "/opt" CACHE PATH "path to BLAS library")
           if(${BLAS} STREQUAL "OPENBLAS")
                set(BLAS_LIB openblas pthread)
                add_definitions(-D_STANDARD_BLAS)
	  elseif (${BLAS} STREQUAL "MKL")
                set(BLAS_LIB mkl_intel_lp64 mkl_sequential mkl_core pthread m dl)
                add_definitions(-D_MKL_BLAS)
          elseif(${BLAS} STREQUAL "CBLAS")
                set(BLAS_LIB cblas blas gfortran)
                add_definitions(-D_STANDARD_BLAS)
          elseif(${BLAS} STREQUAL "ATLAS")
                set(BLAS_LIB f77blas cblas satlas)
                add_definitions(-D_STANDARD_BLAS)
          elseif(${BLAS} STREQUAL "ACCELERATE")
                set(BLAS_LIB -framework Accelerate)
                add_definitions(-flax-vector-conversions -D_ACCELERATE_BLAS)
          endif()
          include_directories(${BLAS_INCLUDE})
          set(CXX_LD ${CXX_LD} ${BLAS_LD})
          set(CXX_LIB ${CXX_LIB} ${BLAS_LIB})
    else()
        add_definitions(-D_NO_BLAS)
    endif()

#gpu options
elseif(${BUILDTYPE} STREQUAL "GPU")
    enable_language(CUDA)
    set(CUDA_LIB -lcudart -lcublas -lcusparse -lm)
    set(CUDA_PATH "/usr/local/cuda" CACHE PATH "path to CUDA toolkit")
    set(CUDA_ARCH "sm_35" CACHE STRING "GPU arch. FERMI sm_20,sm_21; KEPLER sm_30, sm_32,sm_35 etc")
    set(CUDA_INCLUDE ${CUDA_PATH}/include)
    set(CUDA_LD ${CUDA_PATH}/lib64)
    set(CMAKE_CUDA_COMPILER ${CUDA_PATH}/bin/nvcc)
    set(CMAKE_CUDA_FLAGS "-arch=${CUDA_ARCH} -Xcudafe '--diag_suppress=subscript_out_of_range'" )
    include_directories(${CUDA_INCLUDE})
    add_definitions(-D_GPU)
    set(CXX_LIB ${CXX_LIB} ${CUDA_LIB})
    set(CXX_LD ${CXX_LD} ${CUDA_LD}) 
    
    set(SRCLIST
    ${SRCLIST}
    ./src/cuda_kernels.cu)
else()
    message(SEND_ERROR  "Unknow BUILDTYPE. Please select CPU or GPU")
endif()

#MPI
if(${PARALLEL})
      find_package(MPI) 
      add_definitions(-D_MPI)
      set(PARMETIS_INCLUDE "/opt" CACHE PATH "path to parmetis include")
      set(PARMETIS_LD "/opt" CACHE PATH "path to parmetis library")
      if(${MPI_FOUND})
            include_directories(${MPI_CXX_INCLUDE_PATH} ${METIS_INCLUDE} ${PARMETIS_INCLUDE})
            set(CXX_LD ${CXX_LD} ${PARMETIS_LD} ${METIS_LD})
            set(CXX_LIB ${CXX_LIB} ${MPI_CXX_LIBRARIES} parmetis)
            set (SRCLIST ${SRCLIST} ./src/mpi_inters.cpp)
      else()
            message(SEND_ERROR  "Cannot find MPI library, please specify manually.")
      endif()
endif()

#build
LINK_DIRECTORIES(${CXX_LD})
if(${BUILDTYPE} STREQUAL "GPU")
    add_executable(HiFiLES_GPU ${SRCLIST})
    target_link_libraries(HiFiLES_GPU PRIVATE ${CXX_LIB})
endif()
if(${BUILDTYPE} STREQUAL "CPU")
    add_executable(HiFiLES_CPU ${SRCLIST})
    target_link_libraries(HiFiLES_CPU PRIVATE ${CXX_LIB})
endif()
