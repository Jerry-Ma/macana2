DEFINES += QT_DEPRECATED_WARNINGS
MACANA_LIB_DEPS = -L/usr/local/lib -lgsl -lgslcblas -lm -lnetcdf_c++ -lnetcdf -lfftw3 -lcxsparse -lCCfits -lcfitsio
MACANA_LIB_OUT = macana2
CONFIG += \
    c++1z \
    sdk_no_version_check
INCLUDEPATH += \
    /usr/local/include \
    $$IN_PWD/include \
    $$IN_PWD/Sky/Novas