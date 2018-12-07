
# Minimum Qt version
lessThan(QT_MAJOR_VERSION,5): error("Requires Qt 5")
lessThan(QT_MINOR_VERSION,7): error("Requires Qt 5.7")

# Default setting for GUI
isEmpty(GUI) { GUI = "yes" }

# Basic settings
TEMPLATE = subdirs

# Subdir projects
SUBDIRS += \
    core \
    cli \
    tests

# Dependencies
cli.depends = core
tests.depends = core

# This is if GUI is enabled
equals(GUI,"yes") {
    # Add GUI project and dependencies
    SUBDIRS += gui
    gui.depends = core
}
