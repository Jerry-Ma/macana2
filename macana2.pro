CONFIG += ordered
TEMPLATE = subdirs
SUBDIRS += \
    macana_core \
    beammap_gui \
    test \
    port
macana_core.file = macana_core.pro
test.depends = macana_core
beammap_gui.depends = macana_core
