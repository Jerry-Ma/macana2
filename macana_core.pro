#-------------------------------------------------
#
# Project created by QtCreator 2018-10-09T09:51:00
#
#-------------------------------------------------
include(macana2.pri)
TARGET = $$MACANA_LIB_OUT
TEMPLATE = lib

CONFIG -= app_bundle
CONFIG -= qt
CONFIG += staticlib

SOURCES += \
    Analysis/AnalParams.cpp \
    Analysis/SimParams.cpp \
    Clean/AzElTemplateCalculator.cpp \
    Clean/Clean.cpp \
    Clean/Clean2dStripe.cpp \
    Clean/CleanBspline.cpp \
    Clean/CleanHigh.cpp \
    Clean/CleanPCA.cpp \
    Clean/CleanSelector.cpp \
    Mapmaking/Coaddition.cpp \
    Mapmaking/CompletenessSim.cpp \
    Mapmaking/Map.cpp \
    Mapmaking/NoiseRealizations.cpp \
    Mapmaking/Observation.cpp \
    Mapmaking/PointSource.cpp \
    Mapmaking/WienerFilter.cpp \
    Observatory/Array.cpp \
    Observatory/Detector.cpp \
    Observatory/Telescope.cpp \
    Observatory/TimePlace.cpp \
    Simulate/MapNcFile.cpp \
    Simulate/SimulationInserter.cpp \
    Simulate/Subtractor.cpp \
    Sky/Source.cpp \
    Sky/astron_utilities.cpp \
    Utilities/BinomialStats.cpp \
    Utilities/GslRandom.cpp \
    Utilities/SBSM.cpp \
    Utilities/convolution.cpp \
    Utilities/gaussFit.cpp \
    Utilities/mpfit.cpp \
    Utilities/sparseUtilities.cpp \
    Utilities/tinyxml2.cpp \
    Utilities/vector_utilities.cpp \
    Sky/Novas/eph_manager.c \
    Sky/Novas/novas.c \
    Sky/Novas/novascon.c \
    Sky/Novas/nutation.c \
    Sky/Novas/readeph0.c \
    Sky/Novas/solsys1.c
