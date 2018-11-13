cmake_minimum_required(VERSION 3.12)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "Do not build in-source. Please remove CMakeCache.txt and the CMakeFiles/ directory. Then build out-of-source.")
endif()
if (NOT APPLE)
    set(CMAKE_LINK_WHAT_YOU_USE TRUE)  # lwyu does not work with macos linker
endif()
set(CMAKE_CXX_STANDARD 17) # c++17
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # without compiler extensions like gnu++11
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to RelWithDebInfo")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()
# compiling options
OPTION(WITH_OPENMP "Enable OpenMP support?" ON)

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
find_package(OpenMP)
find_package(GSL REQUIRED)
find_package(FFTW3 REQUIRED COMPONENTS DOUBLE SERIAL)
# TODO: wrap these in imported targets
find_package(NetCDF REQUIRED COMPONENTS CXX)
find_package(CCFits REQUIRED)  # this one also find cfitsio
find_package(CXSparse REQUIRED)

print_target_properties(GSL::gsl)
print_target_properties(GSL::gslcblas)
print_target_properties(fftw3::double::serial)

# build options
set(compile_options -Wall -Wextra -DHAVE_INLINE -fexceptions -Wno-source-uses-openmp)  # -Werror
set(link_libraries
    GSL::gsl GSL::gslcblas
    fftw3::double::serial
    ${CXSPARSE_LIBRARIES}
    ${NETCDF_CXX_LIBRARIES}
    ${CFITSIO_LIBRARY} ${CCFITS_LIBRARY}
    )
if (OpenMP_CXX_FOUND AND ${WITH_OPENMP})
    list(APPEND link_libraries OpenMP::OpenMP_CXX)
elseif (OpenMP_CXX_FOUND AND (NOT ${WITH_OPENMP}))
    message("Found OpenMP, use -DWITH_OPENMP=1 to enable")
elseif (${WITH_OPENMP})
    message(FATAL_ERROR "Cannot find OpenMP, abort")
endif()
set(include_directories
    {CMAKE_BINARY_DIR}   # version header
    ${CXSPARSE_INCLUDE_DIRS}
)

# manana core
add_library(
    macana-core OBJECT
    Analysis/AnalParams.cpp
    Analysis/SimParams.cpp
    Clean/AzElTemplateCalculator.cpp
    Clean/Clean.cpp
    Clean/Clean2dStripe.cpp
    Clean/CleanBspline.cpp
    Clean/CleanHigh.cpp
    Clean/CleanPCA.cpp
    Clean/CleanSelector.cpp
    Mapmaking/Coaddition.cpp
    Mapmaking/CompletenessSim.cpp
    Mapmaking/Map.cpp
    Mapmaking/NoiseRealizations.cpp
    Mapmaking/Observation.cpp
    Mapmaking/PointSource.cpp
    Mapmaking/WienerFilter.cpp
    Observatory/Array.cpp
    Observatory/Detector.cpp
    Observatory/Telescope.cpp
    Observatory/TimePlace.cpp
    Simulate/MapNcFile.cpp
    Simulate/SimulationInserter.cpp
    Simulate/Subtractor.cpp
    Sky/Source.cpp
    Sky/astron_utilities.cpp
    Utilities/BinomialStats.cpp
    Utilities/GslRandom.cpp
    Utilities/SBSM.cpp
    Utilities/convolution.cpp
    Utilities/gaussFit.cpp
    Utilities/mpfit.cpp
    Utilities/sparseUtilities.cpp
    Utilities/tinyxml2.cpp
    Utilities/vector_utilities.cpp
    Sky/Novas/eph_manager.c
    Sky/Novas/novas.c
    Sky/Novas/novascon.c
    Sky/Novas/nutation.c
    Sky/Novas/readeph0.c
    Sky/Novas/solsys1.c
    )
target_include_directories(
    macana-core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/Sky/Novas
    ${include_directories}
    )
target_link_libraries(
    macana-core PUBLIC
    ${link_libraries}
    )
target_compile_options(
    macana-core PUBLIC
    ${compile_options})

# setup main executables
# TODO better separate the dependances for macana-core from the main execs
set(macanap_incs)
set(macanap_libs macana-core)

set(beammap_incs)
set(beammap_libs macana-core)

# setup utility executables
set(fitswriter_incs ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(fitswriter_libs ${NETCDF_CXX_LIBRARIES} ${CFITSIO_LIBRARY} ${CCFITS_LIBRARY})

# build
set(exec_list macanap beammap fitswriter)
foreach(exec ${exec_list})
    add_executable(${exec})
    target_sources(
        ${exec} PRIVATE
        ${exec}.cpp
        )
    target_include_directories(
        ${exec} PRIVATE
        ${${exec}_incs}
        )
    target_link_libraries(
        ${exec} PRIVATE
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
