set(CMAKE_C_COMPILER nvc)
set(CMAKE_Fortran_COMPILER nvfortran)
set(CMAKE_CXX_COMPILER nvc++)

# Note: OpenMP_Fortran_FLAGS gets overwritten by the FindOpenMP module
# unless it is stored as a cache variable.
set(OpenMP_Fortran_FLAGS "-mp=bind,allcores,numa" CACHE STRING "" FORCE)

# Match the existing benchmark GPU setup for HPC2020.
set(OpenACC_Fortran_FLAGS "-acc=gpu -gpu=cc80" CACHE STRING "")

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 80)
endif()

# These flags reflect the NVHPC settings already used in the benchmark scripts.
set(ECBUILD_Fortran_FLAGS "-O3 -march=core-avx2" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=core-avx2" CACHE STRING "")
