project(beammap_gui LANGUAGES CXX C)
if (NOT Qt5_PREFIX)
    # set some heuristic Qt5 paths
    # TODO detect Qt5 versions
    if (APPLE)
        set(Qt5_PREFIX "/Applications/Qt/5.11.2/clang_64")
        message("Use Qt5 from ${Qt5_PREFIX}")
    else()
        message(FATAL_ERROR "please set -DQt5_PREFIX=</path/to/Qt/prefix> to compile beammap_gui")
    endif()
endif()
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${Qt5_PREFIX}")
message("module path ${CMAKE_MODULE_PATH}")
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

# Find the QtWidgets library
find_package(Qt5 COMPONENTS Widgets PrintSupport Xml REQUIRED)
add_executable(beammap_gui)
target_compile_options(beammap_gui
    PRIVATE -Wall -Wextra)
target_include_directories(beammap_gui
    PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(beammap_gui
    PRIVATE
        macana2
        Qt5::Widgets
        Qt5::PrintSupport
        Qt5::Xml
        )
target_sources(beammap_gui
    PRIVATE
        beammap_gui.cpp
        domitem.cpp
        doublerangeslider.cpp
        mapviewer.cpp
        colormapcontrol.cpp
        dommodel.cpp
        imageview.cpp
        sciencespinbox.cpp
        qcustomplot/qcustomplot.cpp
        beammap_gui.ui
        imageview.ui
        colormapcontrol.ui
    )
