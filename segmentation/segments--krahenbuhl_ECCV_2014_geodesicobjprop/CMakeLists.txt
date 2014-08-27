cmake_minimum_required(VERSION 2.8.9)
project( pseg )
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

# Checking for OpenMP
find_package( OpenMP )
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(MSVC)
	set(Boost_USE_STATIC_LIBS       ON)
	set(Boost_USE_MULTITHREADED     OFF)
	set(Boost_USE_STATIC_RUNTIME    OFF)
	add_definitions(-DBOOST_ALL_NO_LIB)
	add_definitions(-DBOOST_PYTHON_STATIC_LIB)
endif()

# Finding python and boost-python
find_package(Boost COMPONENTS python-py34)
if(NOT Boost_FOUND)
	find_package(Boost COMPONENTS python3)
endif()
if(NOT Boost_FOUND)
	find_package(Boost COMPONENTS python)
endif()
set(USE_PYTHON ${Boost_FOUND})
if(USE_PYTHON)
	find_package(PythonInterp 3.3 REQUIRED)
	find_package(PythonLibs 3.3 REQUIRED)
	find_package(NumPy REQUIRED)
	set(USE_PYTHON ${PYTHONLIBS_FOUND})
	if( USE_PYTHON )
		include_directories( ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} external/Boost.NumPy  external/liblbfgs-1.10/include )
		link_directories( ${Boost_LIBRARY_DIR} )
	endif()
endif()

# Finding Eigen
find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

set(CMAKE_POSITION_INDEPENDENT_CODE ON) 
# Setup some unix build flags and enable c++11
if (UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -Wno-sign-compare -pedantic") 
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DEIGEN_DONT_PARALLELIZE")
include_directories( external/CImg/ external/rapidxml-1.13/ external/libsvm-3.17/ external/liblbfgs-1.10/include external/maxflow-v3.03.src/ )
add_definitions(-DLBFGS_FLOAT=32)


# Build MATLAB
if(USE_MATLAB)
	add_subdirectory( matlab )
endif()


add_subdirectory( lib )
add_subdirectory( external )
add_subdirectory( examples )

