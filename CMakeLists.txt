cmake_minimum_required(VERSION 3.16)

project(ligint
    VERSION 0.1.0
    DESCRIPTION "Interazioni non covalenti ligando-proteina da file PDB"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Tipo di build" FORCE)
endif()

option(LIGINT_BUILD_CLI "Compila l'eseguibile a riga di comando (ligint)" ON)
option(LIGINT_BUILD_GUI "Compila la GUI Qt (ligint-gui); richiede Qt6 o Qt5 (Widgets, Concurrent)" ON)
option(LIGINT_BUILD_TESTS "Compila i test e li registra in CTest" ON)
option(LIGINT_SANITIZE "Compila con AddressSanitizer e UndefinedBehaviorSanitizer (GCC/Clang)" OFF)
option(LIGINT_WARNINGS_AS_ERRORS "Tratta gli avvisi del compilatore come errori" OFF)

include(GNUInstallDirs)

# Opzioni di compilazione comuni (avvisi, sanitizer) applicate ai target del progetto.
add_library(ligint_options INTERFACE)
if(MSVC)
    target_compile_options(ligint_options INTERFACE /W4 $<$<BOOL:${LIGINT_WARNINGS_AS_ERRORS}>:/WX>)
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
        message(STATUS "ligint: GUI abilitata con Qt ${QT_VERSION_MAJOR}")
        add_subdirectory(gui)
    else()
        message(WARNING "ligint: Qt (Widgets, Concurrent) non trovato: la GUI non verra' compilata. "
                        "Installa Qt6 oppure usa -DLIGINT_BUILD_GUI=OFF per silenziare l'avviso.")
    endif()
endif()

if(LIGINT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ----------------------------------------------------------------------------------------------
# Installazione (libreria, header, programmi) e pacchetto CMake per chi vuole usare la libreria
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
