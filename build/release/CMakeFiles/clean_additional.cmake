# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "gui/CMakeFiles/ligint-gui_autogen.dir/AutogenUsed.txt"
  "gui/CMakeFiles/ligint-gui_autogen.dir/ParseCache.txt"
  "gui/CMakeFiles/ligint_gui_lib_autogen.dir/AutogenUsed.txt"
  "gui/CMakeFiles/ligint_gui_lib_autogen.dir/ParseCache.txt"
  "gui/ligint-gui_autogen"
  "gui/ligint_gui_lib_autogen"
  )
endif()
