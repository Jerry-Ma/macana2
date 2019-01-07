cmake_minimum_required(VERSION 3.12)
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()
# set(CMAKE_MACOSX_RPATH 1)
if (NOT APPLE)
    set(CMAKE_LINK_WHAT_YOU_USE TRUE)  # lwyu does not work with macos linker
endif()
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type specified, default to RelWithDebInfo")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(macana2 LANGUAGES CXX C)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -------
# options
# -------
option(BUILD_SHARED_LIBS "build shared libraries" OFF)
option(MACANA2_BUILD_TESTING "build tests for macana2" ON)
option(MACANA2_BUILD_PORT_TESTING "build tests for citlali-port" ON)
option(MACANA2_BUILD_BEAMMAP_GUI "build beammap GUI" ON)

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
    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/gitversion.h
    @ONLY
)
message("Update gitversion.h")

# add compiler prefix to the search path, useful for llvm compiler on macOS
include(GetCompilerPrefix)
GetCompilerPrefix(cxx_compiler_prefix CXX)
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${cxx_compiler_prefix}")

# ------------
# dependencies
# ------------

# dependencies installed system-wide
find_package(OpenMP)
find_package(GSL REQUIRED)
foreach(gsl_target IN ITEMS GSL::gsl GSL::gslcblas)
    set_target_properties(${gsl_target} PROPERTIES INTERFACE_COMPILE_DEFINITIONS "HAVE_INLINE")
endforeach(gsl_target)
find_package(FFTW3 REQUIRED COMPONENTS DOUBLE SERIAL)

# TODO: wrap these in imported targets
find_package(NetCDF REQUIRED COMPONENTS CXX)
find_package(CCFits REQUIRED)  # this one also find cfitsio
find_package(CXSparse REQUIRED)

# helper to print target attrs
include(PrintProperties)
print_target_properties(OpenMP::OpenMP_CXX)
print_target_properties(GSL::gsl)
print_target_properties(GSL::gslcblas)
print_target_properties(fftw3::double::serial)

# helper to handle dependencies pulled during configure-time
include(FetchContentHelper)

# setup testing/benchmarking framework
if (MACANA2_BUILD_TESTING OR MACANA2_BUILD_PORT_TESTING)
    # gtest
    FetchContentHelper(googletest https://github.com/google/googletest.git release-1.8.1)
    enable_testing()
    # https://gitlab.kitware.com/cmake/community/wikis/doc/tutorials/EmulateMakeCheck
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
    include(GoogleTest)
    # google benchmark
    FetchContentHelper(benchmark https://github.com/google/benchmark.git v1.4.1)
endif()

# -------------
# build targets
# -------------

# global build settings
# add_compile_options(-Weverything)  # -Werror
# include_directories({CMAKE_BINARY_DIR})  # version header

# citlali port lib
add_subdirectory(port)

# macana core lib
add_library(macana2 STATIC)
target_include_directories(macana2
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
        ${CXSPARSE_INCLUDE_DIRS}
        ${NETCDF_INCLUDE_DIRS}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/Sky/Novas
    )
target_link_libraries(macana2
    PUBLIC
        OpenMP::OpenMP_CXX
        GSL::gsl GSL::gslcblas
        fftw3::double::serial
        ${CXSPARSE_LIBRARIES}
        ${NETCDF_CXX_LIBRARIES}
        ${CFITSIO_LIBRARY} ${CCFITS_LIBRARY}
        port
    )
target_compile_options(macana2
    PUBLIC
        -fexceptions
    PRIVATE
        -Weverything
        )
target_sources(macana2
    PRIVATE
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
# macana2 executables
add_executable(macanap macanap.cpp)
target_link_libraries(macanap PRIVATE macana2)
add_custom_command(TARGET macanap
    POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/config
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/config/apDefault.xml ${CMAKE_CURRENT_BINARY_DIR}/config/apDefault.xml
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/config/apSimulate.xml ${CMAKE_CURRENT_BINARY_DIR}/config/apSimulate.xml
)
add_executable(beammap beammap.cpp)
target_link_libraries(beammap PRIVATE macana2)
add_custom_command(TARGET beammap
    POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/config
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/config/apbDefault.xml ${CMAKE_CURRENT_BINARY_DIR}/config/apbDefault.xml
)
add_executable(fitswriter fitswriter.cpp)
target_include_directories(fitswriter
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
        ${NETCDF_INCLUDE_DIRS}
    )
target_link_libraries(fitswriter
    PRIVATE
        ${NETCDF_CXX_LIBRARIES}
        ${CFITSIO_LIBRARY} ${CCFITS_LIBRARY}
    )
# macana2 tests
if (MACANA2_BUILD_TESTING)
    add_subdirectory(test)
endif()
# macana2 beammap gui
if (MACANA2_BUILD_BEAMMAP_GUI)
    add_subdirectory(beammap_gui)
endif()

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
