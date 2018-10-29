include(../macana2.pri)
TARGET = macana_test
TEMPLATE = app

CONFIG -= app_bundle
CONFIG -= qt
CONFIG += thread

SOURCES += \
    test.cpp \
    AnalParamsTest.cpp \
    MapTest.cpp \
    GaussFitTest.cpp

LIBS += \
    -L /usr/local/lib -lgtest -lgmock \
    $$MACANA_LIB_DEPS -L../ -l$$MACANA_LIB_OUT
PRE_TARGETDEPS += ../lib$${MACANA_LIB_OUT}.a
testdata.path    = $${OUT_PWD}/data
testdata.files   = data/*
INSTALLS += testdata
