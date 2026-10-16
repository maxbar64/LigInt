cmake_minimum_required(VERSION 3.16)

project(ligint
    VERSION 0.9.0
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

include(cmake/LigintWindowsResources.cmake)
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
set(LIGINT_PACKAGE_MAINTAINER "ligint authors <ligint@invalid.example>" CACHE STRING
    "Maintainer of the .deb/.rpm packages, as: Name <email> (replace the placeholder before publishing)")

# Two components: Runtime (programs, documentation, desktop integration) and Development (headers, library, CMake package)
include(CMakePackageConfigHelpers)
install(DIRECTORY include/ligint DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT Development)
install(EXPORT ligintTargets NAMESPACE ligint:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ligint COMPONENT Development)
configure_package_config_file(cmake/ligintConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/ligintConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ligint)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/ligintConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ligintConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/ligintConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ligint COMPONENT Development)

# documentation
install(FILES README.md LICENSE CHANGELOG.md DESTINATION ${CMAKE_INSTALL_DOCDIR} COMPONENT Runtime)
install(DIRECTORY docs/img docs/examples DESTINATION ${CMAKE_INSTALL_DOCDIR} COMPONENT Runtime)
# the HTML manual (it is also embedded in ligint-gui): open manual/index.html with a web browser
install(DIRECTORY docs/manual/html/ DESTINATION ${CMAKE_INSTALL_DOCDIR}/manual COMPONENT Runtime PATTERN "manual.qrc" EXCLUDE)
install(FILES docs/ARCHITECTURE.md docs/BENCHMARK.md docs/FINGERPRINT.md docs/INPUT_FORMATS.md docs/KLIFS.md docs/RCSB.md docs/MATRIX.md DESTINATION ${CMAKE_INSTALL_DOCDIR} COMPONENT Runtime)

# manual pages and copyright file (Unix); the pages are compressed when gzip is available, as Debian requires
if(UNIX AND NOT APPLE)
    set(LIGINT_MAN_DIR ${CMAKE_CURRENT_BINARY_DIR}/man)
    set(LIGINT_MAN_PAGES ligint.1)
    if(TARGET ligint-gui)                       # only when the GUI was really built (Qt found)
        list(APPEND LIGINT_MAN_PAGES ligint-gui.1)
    endif()
    find_program(LIGINT_GZIP gzip)
    set(LIGINT_MAN_FILES "")
    foreach(page ${LIGINT_MAN_PAGES})
        set(VERSION ${PROJECT_VERSION})
        configure_file(packaging/${page} ${LIGINT_MAN_DIR}/${page} @ONLY)
        if(LIGINT_GZIP)
            add_custom_command(OUTPUT ${LIGINT_MAN_DIR}/${page}.gz
                COMMAND ${LIGINT_GZIP} -9 -n -c ${LIGINT_MAN_DIR}/${page} > ${LIGINT_MAN_DIR}/${page}.gz
                DEPENDS ${LIGINT_MAN_DIR}/${page} VERBATIM)
            list(APPEND LIGINT_MAN_FILES ${LIGINT_MAN_DIR}/${page}.gz)
        else()
            list(APPEND LIGINT_MAN_FILES ${LIGINT_MAN_DIR}/${page})
        endif()
    endforeach()
    # changelog in Debian format (native package), compressed like the manual pages
    string(TIMESTAMP LIGINT_CHANGELOG_DATE "%a, %d %b %Y %H:%M:%S +0000" UTC)
    configure_file(packaging/changelog.in ${LIGINT_MAN_DIR}/changelog @ONLY)
    if(LIGINT_GZIP)
        add_custom_command(OUTPUT ${LIGINT_MAN_DIR}/changelog.gz
            COMMAND ${LIGINT_GZIP} -9 -n -c ${LIGINT_MAN_DIR}/changelog > ${LIGINT_MAN_DIR}/changelog.gz
            DEPENDS ${LIGINT_MAN_DIR}/changelog VERBATIM)
        set(LIGINT_CHANGELOG_FILE ${LIGINT_MAN_DIR}/changelog.gz)
    else()
        set(LIGINT_CHANGELOG_FILE ${LIGINT_MAN_DIR}/changelog)
    endif()
    add_custom_target(ligint_man ALL DEPENDS ${LIGINT_MAN_FILES} ${LIGINT_CHANGELOG_FILE})
    install(FILES ${LIGINT_MAN_FILES} DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 COMPONENT Runtime)
    install(FILES packaging/copyright ${LIGINT_CHANGELOG_FILE} DESTINATION ${CMAKE_INSTALL_DOCDIR} COMPONENT Runtime)
    # the development package has its own documentation directory (Debian policy: copyright and changelog)
    install(FILES packaging/copyright ${LIGINT_CHANGELOG_FILE} DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/doc/ligint-dev COMPONENT Development)
endif()

# installers and packages (cpack): .deb / .tar.gz on Linux, NSIS installer / .zip on Windows
include(cmake/Packaging.cmake)
