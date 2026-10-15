# Qt GUI: library with the widgets (also used by the smoke test) + ligint-gui executable
set(CMAKE_AUTOMOC ON)

add_library(ligint_gui_lib STATIC
    CommandLine.cpp CommandLine.h
    HitsModel.cpp HitsModel.h
    MainWindow.cpp MainWindow.h
    Style.cpp Style.h
    SvgView.cpp SvgView.h
    ViewerWidget.cpp ViewerWidget.h)
target_include_directories(ligint_gui_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ligint_gui_lib
    PUBLIC ligint::core Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Concurrent
    PRIVATE ligint_options)

# Optional: the SVG preview of the fingerprint needs the Qt Svg module (without it the fingerprint is shown as a table)
option(LIGINT_GUI_WITH_SVG "Use Qt Svg for the fingerprint figure in the GUI when it is available" ON)
if(LIGINT_GUI_WITH_SVG)
    find_package(Qt${QT_VERSION_MAJOR} QUIET COMPONENTS Svg)
endif()
if(LIGINT_GUI_WITH_SVG AND Qt${QT_VERSION_MAJOR}Svg_FOUND)
    message(STATUS "ligint: GUI fingerprint figure enabled (Qt Svg found)")
    target_link_libraries(ligint_gui_lib PUBLIC Qt${QT_VERSION_MAJOR}::Svg)
    target_compile_definitions(ligint_gui_lib PUBLIC LIGINT_GUI_HAS_SVG)
else()
    message(STATUS "ligint: Qt Svg not found: the GUI will show the fingerprint as a table only")
endif()

add_executable(ligint-gui main.cpp)
target_link_libraries(ligint-gui PRIVATE ligint_gui_lib ligint_options)
install(TARGETS ligint-gui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Tool to generate the documentation screenshots (optional)
option(LIGINT_BUILD_SCREENSHOT_TOOL "Build ligint-screenshot (images for the documentation)" OFF)
if(LIGINT_BUILD_SCREENSHOT_TOOL)
    add_executable(ligint-screenshot screenshot.cpp)
    target_link_libraries(ligint-screenshot PRIVATE ligint_gui_lib ligint_options)
endif()
