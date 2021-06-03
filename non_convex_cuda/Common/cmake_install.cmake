# Install script for directory: /home/dongdong/matlab_cuda/Common

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/dongdong/matlab_cuda/../bin/libCommon.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/dongdong/matlab_cuda/../bin" TYPE STATIC_LIBRARY FILES "/home/dongdong/matlab_cuda/Common/libCommon.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/dongdong/matlab_cuda/../include/Common/eps_handler.h;/home/dongdong/matlab_cuda/../include/Common/family.h;/home/dongdong/matlab_cuda/../include/Common/ErrorCodes.h;/home/dongdong/matlab_cuda/../include/Common/GridTools.h;/home/dongdong/matlab_cuda/../include/Common/MultiScaleTools.h;/home/dongdong/matlab_cuda/../include/Common/PythonTypes.h;/home/dongdong/matlab_cuda/../include/Common/TCostFunctionProvider.h;/home/dongdong/matlab_cuda/../include/Common/TCouplingHandler.h;/home/dongdong/matlab_cuda/../include/Common/TEpsScaling.h;/home/dongdong/matlab_cuda/../include/Common/THierarchicalCostFunctionProvider.h;/home/dongdong/matlab_cuda/../include/Common/THierarchicalPartition.h;/home/dongdong/matlab_cuda/../include/Common/THierarchyBuilder.h;/home/dongdong/matlab_cuda/../include/Common/Tools.h;/home/dongdong/matlab_cuda/../include/Common/TVarListHandler.h;/home/dongdong/matlab_cuda/../include/Common/Verbose.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/dongdong/matlab_cuda/../include/Common" TYPE FILE FILES
    "/home/dongdong/matlab_cuda/Common/eps_handler.h"
    "/home/dongdong/matlab_cuda/Common/family.h"
    "/home/dongdong/matlab_cuda/Common/ErrorCodes.h"
    "/home/dongdong/matlab_cuda/Common/GridTools.h"
    "/home/dongdong/matlab_cuda/Common/MultiScaleTools.h"
    "/home/dongdong/matlab_cuda/Common/PythonTypes.h"
    "/home/dongdong/matlab_cuda/Common/TCostFunctionProvider.h"
    "/home/dongdong/matlab_cuda/Common/TCouplingHandler.h"
    "/home/dongdong/matlab_cuda/Common/TEpsScaling.h"
    "/home/dongdong/matlab_cuda/Common/THierarchicalCostFunctionProvider.h"
    "/home/dongdong/matlab_cuda/Common/THierarchicalPartition.h"
    "/home/dongdong/matlab_cuda/Common/THierarchyBuilder.h"
    "/home/dongdong/matlab_cuda/Common/Tools.h"
    "/home/dongdong/matlab_cuda/Common/TVarListHandler.h"
    "/home/dongdong/matlab_cuda/Common/Verbose.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/dongdong/matlab_cuda/Common/Models/cmake_install.cmake")

endif()

