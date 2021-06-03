cmake_minimum_required(VERSION 3.5)

# project name
project (CUPORTFOLIO)

 

add_definitions(-std=c++11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++11  -mavx -mfma")

include_directories(.)
INCLUDE_DIRECTORIES("/usr/include/eigen3")

find_package(Eigen3 REQUIRED)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")

set(CUDA_SEPARABLE_COMPILATION ON)
find_package(CUDA REQUIRED)

INCLUDE_DIRECTORIES(./${CUDA_INCLUDE_DIRS})
include_directories(.)
set(CUDA_NVCC_FLAGS  --fmad=true -Xptxas -dlcm=ca -O3
    -gencode=arch=compute_35,code=sm_35
   -gencode=arch=compute_50,code=sm_50
    -gencode=arch=compute_60,code=sm_60
   -gencode=arch=compute_70,code=sm_70
    -gencode=arch=compute_75,code=sm_75
   -gencode=arch=compute_75,code=compute_75 )

set(CUDA_NVCC_FLAGS
    --relocatable-device-code=true -Xcompiler
    -fPIC)
LINK_DIRECTORIES("/usr/local/cuda-10.0/lib64")
INCLUDE_DIRECTORIES ("/usr/local/cuda-10.0/include")
 
INCLUDE_DIRECTORIES (Common)
add_subdirectory(Common)
 

 


cuda_add_executable(CuPortfolio   main.cpp portfoliomodel.cpp  portfoliosolver.cu
    portfoliosolver.cpp kernel_function.cu dfp_cu.cu )
target_compile_options(CuPortfolio PUBLIC ${LEMON_COMPILE_OPTIONS})
target_include_directories(CuPortfolio PUBLIC ${LEMON_INCLUDE_DIRECTORY})
target_link_libraries(CuPortfolio
        ${LEMON_LIBRARY}
        Common
        -lcublas
        -fopenmp
        -lcudart -lcusolver
        )

cuda_add_executable(LogCuPortfolio   main_log.cpp portfoliomodel.cpp  portfoliosolver.cu
    portfoliosolver.cpp kernel_function.cu dfp_cu.cu )
target_compile_options(LogCuPortfolio PUBLIC ${LEMON_COMPILE_OPTIONS})
target_include_directories(LogCuPortfolio PUBLIC ${LEMON_INCLUDE_DIRECTORY})
target_link_libraries(LogCuPortfolio
        ${LEMON_LIBRARY}
        Common
        -lcublas
        -fopenmp
        -lcudart -lcusolver
        )
 



