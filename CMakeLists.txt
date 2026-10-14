cmake_minimum_required(VERSION 3.16)

project(ligint
    VERSION 0.3.0
    DESCRIPTION "Non-covalent ligand-protein interactions from PDB files"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(LIGINT_BUILD_CLI "Build the command-line executable (ligint)" ON)
option(LIGINT_BUILD_GUI "Build the Qt GUI (ligint-gui); requires Qt6 or Qt5 (Widgets, Concurrent)" ON)
option(LIGINT_BUILD_TESTS "Build the tests and register them in CTest" ON)
option(LIGINT_SANITIZE "Build with AddressSanitizer and UndefinedBehaviorSanitizer (GCC/Clang)" OFF)
option(LIGINT_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

include(GNUInstallDirs)

# Common compile options (warnings, sanitizers) applied to the project targets.
add_library(ligint_options INTERFACE)
if(MSVC)
    target_compile_options(ligint_options INTERFACE /utf-8 /W4 $<$<BOOL:${LIGINT_WARNINGS_AS_ERRORS}>:/WX>)
else()
    target_compile_options(ligint_options INTERFACE -Wall -Wextra
        $<$<BOOL:${LIGINT_WARNINGS_AS_ERRORS}>:-Werror>)
    if(LIGINT_SANITIZE)
        target_compile_options(ligint_options INTERFACE -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer -g)
        target_link_options(ligint_options INTERFACE -fsanitize=address,undefined)
    endif()
endif()

add_subdirectory(src)

if(LIGINT_BUILD_GUI)
    find_package(QT NAMES Qt6 Qt5 QUIET COMPONENTS Widgets Concurrent)
    if(QT_FOUND)
        find_package(Qt${QT_VERSION_MAJOR} QUIET COMPONENTS Widgets Concurrent)
    endif()
    if(QT_FOUND AND Qt${QT_VERSION_MAJOR}_FOUND)
        message(STATUS "ligint: GUI enabled with Qt ${QT_VERSION_MAJOR}")
        add_subdirectory(gui)
    else()
        message(WARNING "ligint: Qt (Widgets, Concurrent) not found: the GUI will not be built. "
                        "Install Qt6 or use -DLIGINT_BUILD_GUI=OFF to silence this warning.")
    endif()
endif()

if(LIGINT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ----------------------------------------------------------------------------------------------
# Installation (library, headers, programs) and CMake package for those who want to use the library
# ----------------------------------------------------------------------------------------------
include(CMakePackageConfigHelpers)
install(DIRECTORY include/ligint DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT ligintTargets NAMESPACE ligint:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ligint)
configure_package_config_file(cmake/ligintConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/ligintConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ligint)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/ligintConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ligintConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/ligintConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ligint)
install(FILES README.md LICENSE DESTINATION ${CMAKE_INSTALL_DOCDIR})
