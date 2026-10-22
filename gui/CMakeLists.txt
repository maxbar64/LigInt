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
target_compile_definitions(ligint_gui_lib PRIVATE LIGINT_SOURCE_ROOT=\"${PROJECT_SOURCE_DIR}\")
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

# Copy the Qt libraries and plugins next to the program when installing (Windows: needs Qt 6.3 or later; Linux: Qt 6.5 or later).
# On by default on Windows; on Linux it is for the self-contained .tar.gz (the .deb uses the Qt of the distribution).
if(WIN32)
    set(_ligint_bundle_qt_default ON)
else()
    set(_ligint_bundle_qt_default OFF)
endif()
option(LIGINT_BUNDLE_QT "Install the Qt libraries and plugins with ligint-gui (Linux: also json-c and the launchers bin/*.sh)" ${_ligint_bundle_qt_default})
if(LIGINT_BUNDLE_QT AND QT_VERSION_MAJOR EQUAL 6)
    if(UNIX AND NOT APPLE)
        # with a Qt in a system prefix the deployment cannot tell Qt from the system libraries, and would bundle glibc too
        if(QT6_INSTALL_PREFIX MATCHES "^/usr(/|$)")
            message(FATAL_ERROR "LIGINT_BUNDLE_QT needs a Qt installed outside /usr (this one is in ${QT6_INSTALL_PREFIX})")
        endif()
        # the bundled Qt is found next to the program, wherever the package is unpacked
        set_target_properties(ligint-gui PROPERTIES INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
        # besides the default platform (xcb), the offscreen one: headless runs (QT_QPA_PLATFORM=offscreen), as in the CI
        qt_import_plugins(ligint-gui INCLUDE Qt6::QOffscreenIntegrationPlugin)
        # the list of plugins to deploy (platforms, image formats...) is written only for finalized targets (windeployqt finds them by itself)
        qt_finalize_target(ligint-gui)
    endif()
    if(COMMAND qt_generate_deploy_app_script)
        set(LIGINT_DEPLOY_OPTIONS NO_UNSUPPORTED_PLATFORM_ERROR)
        if(Qt6_VERSION VERSION_GREATER_EQUAL 6.5)
            list(APPEND LIGINT_DEPLOY_OPTIONS NO_TRANSLATIONS)     # the GUI is in English: skip the Qt translation files
        endif()
        qt_generate_deploy_app_script(TARGET ligint-gui OUTPUT_SCRIPT LIGINT_QT_DEPLOY_SCRIPT ${LIGINT_DEPLOY_OPTIONS})
        install(SCRIPT ${LIGINT_QT_DEPLOY_SCRIPT} COMPONENT Runtime)
    else()
        message(WARNING "ligint: qt_generate_deploy_app_script not available (Qt older than 6.3): the package will not include the Qt libraries")
    endif()
endif()

# Tool to generate the documentation screenshots (optional)
option(LIGINT_BUILD_SCREENSHOT_TOOL "Build ligint-screenshot (images for the documentation)" OFF)
if(LIGINT_BUILD_SCREENSHOT_TOOL)
    add_executable(ligint-screenshot screenshot.cpp)
    target_link_libraries(ligint-screenshot PRIVATE ligint_gui_lib ligint_options)
endif()
