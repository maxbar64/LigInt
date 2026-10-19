# Library with the analysis (no dependencies beyond the C++17 standard library)
add_library(ligint_core
    core/analysis.cpp
    core/secondary.cpp
    core/cartoon.cpp
    core/surface.cpp
    core/pocket.cpp
    core/chemistry.cpp
    core/mmcif.cpp
    core/pdb.cpp
    core/fingerprint.cpp
    core/matrix.cpp
    core/gpcrdb.cpp
    core/report.cpp
    core/svg.cpp
    core/templates.cpp
    klifs/klifs.cpp
    core/rcsb.cpp
    core/internal.hpp)
add_library(ligint::core ALIAS ligint_core)
target_include_directories(ligint_core
    PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE core)
target_compile_definitions(ligint_core PUBLIC LIGINT_VERSION="${PROJECT_VERSION}")
target_link_libraries(ligint_core PRIVATE ligint_options)
set_target_properties(ligint_core PROPERTIES POSITION_INDEPENDENT_CODE ON EXPORT_NAME core)   # exported as ligint::core
install(TARGETS ligint_core ligint_options EXPORT ligintTargets
    COMPONENT Development                                   # before the artifact kinds: it applies to all of them
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Command-line program
if(LIGINT_BUILD_CLI)
    add_executable(ligint cli/main.cpp)
    target_link_libraries(ligint PRIVATE ligint::core ligint_options)
    ligint_add_windows_resources(ligint ligint "ligint - ligand-protein interactions (command line)")
    install(TARGETS ligint RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime)
endif()

# Native KLIFS client (online HTTPS optional unless explicitly required).
# Use both pkg-config (Unix) and native CMake/vcpkg discovery (Windows).
option(LIGINT_KLIFS_REQUIRE_ONLINE "Fail configure if native KLIFS HTTPS support is unavailable" OFF)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(JSONC QUIET IMPORTED_TARGET json-c)
    pkg_check_modules(KLIFS_CURL QUIET IMPORTED_TARGET libcurl)
endif()

set(_klifs_curl_target "")
if(TARGET PkgConfig::KLIFS_CURL)
    set(_klifs_curl_target PkgConfig::KLIFS_CURL)
else()
    find_package(CURL QUIET)
    if(CURL_FOUND)
        set(_klifs_curl_target CURL::libcurl)
    endif()
endif()

set(_klifs_json_target "")
if(TARGET PkgConfig::JSONC)
    set(_klifs_json_target PkgConfig::JSONC)
else()
    # json-c may be installed without a CMake package config file.
    find_path(KLIFS_JSONC_INCLUDE_DIR NAMES json-c/json.h)
    find_library(KLIFS_JSONC_LIBRARY NAMES json-c json-c-static)
    if(KLIFS_JSONC_INCLUDE_DIR AND KLIFS_JSONC_LIBRARY)
        add_library(ligint_klifs_jsonc UNKNOWN IMPORTED)
        set_target_properties(ligint_klifs_jsonc PROPERTIES
            IMPORTED_LOCATION "${KLIFS_JSONC_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${KLIFS_JSONC_INCLUDE_DIR}")
        set(_klifs_json_target ligint_klifs_jsonc)
    endif()
endif()

if(_klifs_curl_target AND _klifs_json_target)
    target_compile_definitions(ligint_core PRIVATE LIGINT_KLIFS_ONLINE)
    # ligint_core is a static library: its dependencies reach the programs that link it, also through the installed CMake package.
    # The imported targets found here (PkgConfig::..., vcpkg) exist only in this build, so the package names the system libraries instead.
    target_link_libraries(ligint_core PRIVATE $<BUILD_INTERFACE:${_klifs_curl_target}> $<BUILD_INTERFACE:${_klifs_json_target}>
                                              $<INSTALL_INTERFACE:curl> $<INSTALL_INTERFACE:json-c>)
    set(LIGINT_NET_CURL ON CACHE INTERNAL "ligint is built with libcurl")
    set(LIGINT_NET_JSONC ON CACHE INTERNAL "ligint is built with json-c")
    message(STATUS "ligint: native KLIFS HTTPS enabled (curl=${_klifs_curl_target}, json=${_klifs_json_target})")
else()
    if(LIGINT_KLIFS_REQUIRE_ONLINE)
        message(FATAL_ERROR "Native KLIFS HTTPS requires both libcurl and json-c. Install their development packages or vcpkg ports curl and json-c and configure with the vcpkg toolchain file.")
    endif()
    message(WARNING "ligint: KLIFS online disabled (curl target='${_klifs_curl_target}', json-c target='${_klifs_json_target}'); local KLIFS mappings still work")
endif()

# RCSB downloads only need curl, independently of JSON-C / KLIFS availability.
option(LIGINT_RCSB_REQUIRE_ONLINE "Fail configure if native RCSB HTTPS is unavailable" OFF)
if(_klifs_curl_target)
    target_compile_definitions(ligint_core PRIVATE LIGINT_RCSB_ONLINE)
    target_link_libraries(ligint_core PRIVATE $<BUILD_INTERFACE:${_klifs_curl_target}> $<INSTALL_INTERFACE:curl>)
    set(LIGINT_NET_CURL ON CACHE INTERNAL "ligint is built with libcurl")
    message(STATUS "ligint: RCSB HTTPS enabled")
elseif(LIGINT_RCSB_REQUIRE_ONLINE)
    message(FATAL_ERROR "RCSB HTTPS requires libcurl. Install curl via vcpkg and configure with the vcpkg toolchain.")
else()
    message(WARNING "ligint: RCSB HTTPS disabled (libcurl not detected)")
endif()
