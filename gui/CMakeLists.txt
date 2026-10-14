# Qt GUI: library with the widgets (also used by the smoke test) + ligint-gui executable
set(CMAKE_AUTOMOC ON)

add_library(ligint_gui_lib STATIC
    HitsModel.cpp HitsModel.h
    MainWindow.cpp MainWindow.h
    Style.cpp Style.h
    ViewerWidget.cpp ViewerWidget.h)
target_include_directories(ligint_gui_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ligint_gui_lib
    PUBLIC ligint::core Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Concurrent
    PRIVATE ligint_options)

add_executable(ligint-gui main.cpp)
target_link_libraries(ligint-gui PRIVATE ligint_gui_lib ligint_options)
install(TARGETS ligint-gui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Tool to generate the documentation screenshots (optional)
option(LIGINT_BUILD_SCREENSHOT_TOOL "Build ligint-screenshot (images for the documentation)" OFF)
if(LIGINT_BUILD_SCREENSHOT_TOOL)
    add_executable(ligint-screenshot screenshot.cpp)
    target_link_libraries(ligint-screenshot PRIVATE ligint_gui_lib ligint_options)
endif()
