include(../macana2.pri)
TARGET = macana_test
TEMPLATE = app

CONFIG -= app_bundle
CONFIG -= qt
CONFIG += thread

SOURCES += \
    test.cpp \
    AnalParamsTest.cpp

LIBS += \
    -L /usr/local/lib -lgtest -lgmock \
    $$MACANA_LIB_DEPS -L../ -l$$MACANA_LIB_OUT

testdata.path    = $${OUT_PWD}/data
testdata.files   = data/*
INSTALLS += testdata
