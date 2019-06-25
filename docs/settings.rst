Settings
========

KINC has some settings which are global to all analytics and are stored in a file on disk. These settings can be changed by using the ``kinc settings`` command. Running this command alone will list all of the settings and their current values.

To change a setting, run ``kinc settings set <parameter> <value>``. For example, to disable CUDA, run ``kinc settings set cuda none``.

For more information, consult the help text by running ``kinc help settings``.
