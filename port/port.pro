DEFINES += QT_DEPRECATED_WARNINGS
CONFIG += \
    c++1z \
    sdk_no_version_check

TARGET = port_test
TEMPLATE = app

CONFIG -= app_bundle
CONFIG -= qt
CONFIG += thread

SOURCES += \
    port_generic_curvefit.cpp \
    test/test_port_generic_curvefit.cpp

HEADERS += \
    include/port_generic_curvefit.cpp

LIBS += \
    -L /usr/local/lib -lgtest -lgtest_main -lgmock
INCLUDEPATH += \
    /usr/local/include

port_test_data.path    = $${OUT_PWD}/data
port_test_data.files   = data/*
INSTALLS += port_test_data
