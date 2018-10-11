#-------------------------------------------------
#
# Project created by QtCreator 2018-10-09T09:51:00
#
#-------------------------------------------------

QT       += core gui xml

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = macana2
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += \
        c++14 \
        sdk_no_version_check

SOURCES += \
        beammap_gui.cpp \
        qcustomplot/qcustomplot.cpp \
        beammap_gui/dommodel.cpp \
        beammap_gui/domitem.cpp \
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

HEADERS += \
        beammap_gui.h \
        qcustomplot/qcustomplot.h \
        beammap_gui/dommodel.h \
        beammap_gui/domitem.h \

FORMS += \
        beammap_gui.ui

LIBS += -L/usr/local/lib -lgsl -lgslcblas -lm -lnetcdf_c++ -lnetcdf -lfftw3 -lcxsparse -lCCfits -lcfitsio
INCLUDEPATH += \
        /usr/local/include \
        $$PWD/include \
        $$PWD/Sky/Novas

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
