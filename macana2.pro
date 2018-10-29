CONFIG += ordered
TEMPLATE = subdirs
SUBDIRS += \
    macana_core \
    test \
    beammap_gui
macana_core.file = macana_core.pro
test.depends = macana_core
beammap_gui.depends = macana_core
