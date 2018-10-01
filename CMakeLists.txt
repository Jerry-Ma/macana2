cmake_minimum_required(VERSION 3.12)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "Do not build in-source. Please remove CMakeCache.txt and the CMakeFiles/ directory. Then build out-of-source.")
endif()
if (NOT APPLE)
    set(CMAKE_LINK_WHAT_YOU_USE TRUE)  # lwyu does not work with macos linker
endif()
set(CMAKE_CXX_STANDARD 11) # c++11
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # without compiler extensions like gnu++11
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to RelWithDebInfo")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# set(CMAKE_MACOSX_RPATH 1)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(macana2)

# load some cmake modules print out some useful info
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(cmake/PrintProperties.cmake)  # helper to print target attrs
find_package(Git)  # automatic git version
set(VERSION_UPDATE_FROM_GIT 1)
include(GetVersionFromGitTag)
set(_build_version ${${CMAKE_PROJECT_NAME}_VERSION_STRING_FULL})
string(TIMESTAMP _time_stamp)
message("==================================================")
message("PROJECT NAME: " ${CMAKE_PROJECT_NAME})
message("VERSION     : " ${_build_version})
message("TIME        : " ${_time_stamp})
message("==================================================")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/gitversion.h.in
    ${CMAKE_BINARY_DIR}/gitversion.h
    @ONLY
)
message("Update gitversion.h")

# dependencies
# TODO: better handle varieties of openmp on macos
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(FFTW3 REQUIRED COMPONENTS DOUBLE SERIAL)
# TODO: wrap these in imported targets
find_package(NetCDF REQUIRED COMPONENTS CXX)
find_package(CCFits REQUIRED)  # this one also find cfitsio
find_package(CXSparse REQUIRED)

print_target_properties(GSL::gsl)
print_target_properties(GSL::gslcblas)
print_target_properties(fftw3::double::serial)

# build executables
set(compile_options -Wall -Wextra -DHAVE_INLINE -fexceptions)  # -Werror
include_directories(${CMAKE_BINARY_DIR})  # version header

# source files
set(exec_list macanap beammap fitswriter)
set(
    macanap_sources
    Analysis/AnalParams.cpp Analysis/SimParams.cpp

    )
file(GLOB macanap_srcs Analysis/*.cpp Clean/*.cpp Mapmaking/*.cpp Observatory/*.cpp Simulate/*.cpp Sky/*.cpp Sky/Novas/*.c Utilities/*.cpp)
set(macanap_incs ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/Sky/Novas)
set(
    macanap_libs
    OpenMP::OpenMP_CXX
    GSL::gsl GSL::gslcblas
    fftw3::double::serial
    ${CXSPARSE_LIBRARIES}
    ${NETCDF_CXX_LIBRARIES}
    ${CFITSIO_LIBRARY} ${CCFITS_LIBRARY}
    )

file(GLOB beammap_srcs Analysis/*.cpp Clean/*.cpp Mapmaking/*.cpp Observatory/*.cpp Simulate/*.cpp Sky/*.cpp Sky/Novas/*.c Utilities/*.cpp)
set(beammap_incs ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/Sky/Novas)
set(
    beammap_libs
    OpenMP::OpenMP_CXX
    GSL::gsl GSL::gslcblas
    fftw3::double::serial
    ${CXSPARSE_LIBRARIES}
    ${NETCDF_CXX_LIBRARIES}
    ${CFITSIO_LIBRARY} ${CCFITS_LIBRARY}
    )

set(fitswriter_srcs) # dummy
set(fitswriter_incs ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(
    fitswriter_libs
    ${NETCDF_CXX_LIBRARIES}
    ${CFITSIO_LIBRARY} ${CCFITS_LIBRARY}
    )

foreach(exec ${exec_list})
    add_executable(${exec})
    target_compile_options(
    ${exec} PRIVATE
    ${compile_options})
    target_sources(
        ${exec} PRIVATE
        ${exec}.cpp
        ${${exec}_srcs}
        )
    target_include_directories(
        ${exec} PRIVATE
        ${${exec}_incs}
        )
    target_link_libraries(
        ${exec}
        PRIVATE
        ${${exec}_libs}
        )
endforeach(exec ${exec_list})

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
