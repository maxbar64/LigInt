# Qt GUI: library with the widgets (also used by the smoke test) + ligint-gui executable
set(CMAKE_AUTOMOC ON)
find_package(Threads REQUIRED)    # the 3D rasterizer uses std::thread

add_library(ligint_gui_lib STATIC
    CommandLine.cpp CommandLine.h
    HitsModel.cpp HitsModel.h
    MainWindow.cpp MainWindow.h
    Style.cpp Style.h
    SvgView.cpp SvgView.h
    HelpWindow.cpp HelpWindow.h
    ${PROJECT_SOURCE_DIR}/docs/manual/html/manual.qrc
    PocketWorkspace.cpp PocketWorkspace.h
    ViewerWidget.cpp ViewerWidget.h
    ViewerOptionsWidget.cpp ViewerOptionsWidget.h
    Raster.cpp Raster.h
    SceneBuilder.cpp SceneBuilder.h ViewStyle.h)
# the manual (HTML pages and images) is embedded in the library as Qt resources
set_target_properties(ligint_gui_lib PROPERTIES AUTORCC ON)
target_include_directories(ligint_gui_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ligint_gui_lib
    PUBLIC ligint::core Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Concurrent Threads::Threads
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

set(CMAKE_AUTORCC ON)
add_executable(ligint-gui main.cpp resources.qrc)
target_link_libraries(ligint-gui PRIVATE ligint_gui_lib ligint_options)
ligint_add_windows_resources(ligint-gui ligint-gui "ligint - ligand-protein interactions")
if(WIN32)
    set_target_properties(ligint-gui PROPERTIES WIN32_EXECUTABLE ON)     # no console window
endif()
install(TARGETS ligint-gui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime)

# desktop integration (Linux): menu entry, PDB association and icon
if(UNIX AND NOT APPLE)
    install(FILES ${PROJECT_SOURCE_DIR}/packaging/ligint-gui.desktop DESTINATION ${CMAKE_INSTALL_DATADIR}/applications COMPONENT Runtime)
    install(FILES ${PROJECT_SOURCE_DIR}/packaging/ligint.svg DESTINATION ${CMAKE_INSTALL_DATADIR}/icons/hicolor/scalable/apps COMPONENT Runtime)
endif()

# Windows: copy the Qt libraries and plugins next to the program when installing (needs Qt 6.3 or later)
if(WIN32 AND QT_VERSION_MAJOR EQUAL 6)
    if(COMMAND qt_generate_deploy_app_script)
        set(LIGINT_DEPLOY_OPTIONS NO_UNSUPPORTED_PLATFORM_ERROR)
        if(Qt6_VERSION VERSION_GREATER_EQUAL 6.5)
            list(APPEND LIGINT_DEPLOY_OPTIONS NO_TRANSLATIONS)     # the GUI is in English: skip the Qt translation files
        endif()
        qt_generate_deploy_app_script(TARGET ligint-gui OUTPUT_SCRIPT LIGINT_QT_DEPLOY_SCRIPT ${LIGINT_DEPLOY_OPTIONS})
        install(SCRIPT ${LIGINT_QT_DEPLOY_SCRIPT} COMPONENT Runtime)
    else()
        message(WARNING "ligint: qt_generate_deploy_app_script not available (Qt older than 6.3): the Windows installer will not include the Qt libraries")
    endif()
endif()

# Tool to generate the documentation screenshots (optional)
option(LIGINT_BUILD_SCREENSHOT_TOOL "Build ligint-screenshot (images for the documentation)" OFF)
if(LIGINT_BUILD_SCREENSHOT_TOOL)
    add_executable(ligint-screenshot screenshot.cpp)
    target_link_libraries(ligint-screenshot PRIVATE ligint_gui_lib ligint_options)
endif()
