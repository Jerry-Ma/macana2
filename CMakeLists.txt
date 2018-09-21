## cmake config
cmake_minimum_required(VERSION 3.1)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "Do not build in-source. Please remove CMakeCache.txt and the CMakeFiles/ directory. Then build out-of-source.")
endif()
# set(CMAKE_MACOSX_RPATH 1)
set(cmake_modules ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
LIST(APPEND CMAKE_MODULE_PATH "${cmake_modules}")

## source files
project(macana2)
# set(src ${CMAKE_CURRENT_SOURCE_DIR}/src)
file(GLOB srcs Analysis/*.cpp Clean/*.cpp Mapmaking/*.cpp Observatory/*.cpp Simulate/*.cpp Sky/*.cpp Sky/Novas/*.c Utilities/*.cpp)
# list(REMOVE_ITEM srcs "Utilities/fitswriter.cpp")
set(include ${CMAKE_CURRENT_SOURCE_DIR}/include)

## cmake magics to resolve deps
# automatic git version
find_package(Git)
set(VERSION_UPDATE_FROM_GIT 1)
include(${cmake_modules}/cmake-git-semver/GetVersionFromGitTag.cmake)
message("==================================================")
message("PROJECT NAME: " ${CMAKE_PROJECT_NAME})
message("VERSION     : " ${${CMAKE_PROJECT_NAME}_VERSION_STRING_FULL})
message("==================================================")
# TODO openmp
# find_package(OpenMP REQUIRED)
set(OpenMP_C "${CMAKE_C_COMPILER}")
set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
set(OPENMP_LIBRARIES "/usr/local/opt/llvm/lib")
set(OPENMP_INCLUDES "/usr/local/opt/llvm/include")
include_directories("${OPENMP_INCLUDES}")
link_directories("${OPENMP_LIBRARIES}")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

find_package(GSL REQUIRED)
find_package(FFTW3 REQUIRED )
set (NETCDF_CXX "YES")
find_package (NetCDF REQUIRED)
find_package(CCFits REQUIRED)
find_package(CXSparse REQUIRED)

# build executables
include_directories(${include})
include_directories(Sky/Novas)

set(execs macanap beammap fitswriter)
foreach(exec ${execs})
    add_executable(${exec} ${exec}.cpp ${srcs})
    target_compile_options (${exec} PRIVATE -c -g -Wall -DHAVE_INLINE -O2 -fexceptions)
    target_link_libraries(${exec} GSL::gsl GSL::gslcblas ${FFTW3_LIBRARIES} ${NETCDF_CXX_LIBRARIES} ${CFITSIO_LIBRARY} ${CCFITS_LIBRARY} ${CXSPARSE_LIBRARIES})
endforeach(exec ${execs})

## docs
# https://vicrucann.github.io/tutorials/quick-cmake-doxygen/
#option(BUILD_DOC "Build documentation" ON)
#find_package(Doxygen)
#if (DOXYGEN_FOUND)
#    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
#    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
#    set(DOXYGEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs)
#    set(DOXYGEN_VER ${${CMAKE_PROJECT_NAME}_VERSION})
#    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
#    message("Doxygen build started")
#    add_custom_target(
#        doc ALL
#        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#        COMMENT "Generating API documentation with Doxygen"
#        VERBATIM)
#else()
#  message("Doxygen need to be installed to generate the doxygen documentation")
#endif()
