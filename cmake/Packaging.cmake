# Installers and packages, built with `cpack` after `cmake --build`:
#   Linux:   ligint_<version>_<arch>.deb (+ ligint-dev), ligint-<version>-linux-<arch>.tar.gz; `cpack -G RPM` for .rpm
#   Windows: ligint-<version>-win64.exe (NSIS installer), ligint-<version>-win64.zip (portable)
set(LIGINT_HOMEPAGE "https://github.com/maxbar64/LigInt" CACHE STRING "Project home page (package metadata)")

set(CPACK_PACKAGE_NAME "ligint")
set(CPACK_PACKAGE_VENDOR "ligint authors")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Non-covalent ligand-protein interactions from PDB files")
# wrapped: Debian wants lines of at most 80 characters in the extended description
set(CPACK_PACKAGE_DESCRIPTION
"ligint finds hydrogen bonds, salt bridges, hydrophobic contacts,
pi-stacking, cation-pi, halogen bonds, metal coordination and water
bridges between a ligand and a protein in a PDB file. It can summarize
them as an interaction fingerprint (CSV, JSON, SVG).
.
Includes a command-line program and a Qt graphical interface.")
set(CPACK_PACKAGE_HOMEPAGE_URL ${LIGINT_HOMEPAGE})
set(CPACK_PACKAGE_CONTACT ${LIGINT_PACKAGE_MAINTAINER})
set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/LICENSE)
set(CPACK_PACKAGE_INSTALL_DIRECTORY "ligint")
set(CPACK_STRIP_FILES ON)

if(WIN32)
    set(LIGINT_PACKAGE_SYSTEM "win64")
    set(CPACK_GENERATOR "NSIS;ZIP")
elseif(APPLE)
    set(LIGINT_PACKAGE_SYSTEM "macos-${CMAKE_SYSTEM_PROCESSOR}")
    set(CPACK_GENERATOR "TGZ")
else()
    set(LIGINT_PACKAGE_SYSTEM "linux-${CMAKE_SYSTEM_PROCESSOR}")
    set(CPACK_GENERATOR "DEB;TGZ")
endif()
set(CPACK_PACKAGE_FILE_NAME "ligint-${PROJECT_VERSION}-${LIGINT_PACKAGE_SYSTEM}")

# the component names that cpack needs, and the values that the per-generator file below reads
set(CPACK_COMPONENTS_ALL Runtime Development)
set(CPACK_LIGINT_MAINTAINER ${LIGINT_PACKAGE_MAINTAINER})
set(CPACK_LIGINT_BUILD_DIR ${CMAKE_BINARY_DIR})
# network features (RCSB downloads, KLIFS service) make the static library depend on libcurl and json-c
set(CPACK_LIGINT_DEV_DEBDEPS "")
set(CPACK_LIGINT_DEV_RPMDEPS "")
if(LIGINT_NET_CURL)
    set(CPACK_LIGINT_DEV_DEBDEPS "libcurl4-openssl-dev | libcurl-dev")
    set(CPACK_LIGINT_DEV_RPMDEPS "libcurl-devel")
endif()
if(LIGINT_NET_JSONC)
    set(CPACK_LIGINT_DEV_DEBDEPS "${CPACK_LIGINT_DEV_DEBDEPS}, libjson-c-dev")
    set(CPACK_LIGINT_DEV_RPMDEPS "${CPACK_LIGINT_DEV_RPMDEPS}, json-c-devel")
endif()
if(TARGET ligint-gui)
    set(CPACK_LIGINT_HAS_GUI ON)
else()
    set(CPACK_LIGINT_HAS_GUI OFF)
endif()
set(CPACK_LIGINT_ICON ${PROJECT_SOURCE_DIR}/packaging/ligint.ico)
set(CPACK_PROJECT_CONFIG_FILE ${PROJECT_SOURCE_DIR}/cmake/CPackOptions.cmake)

include(CPack)
cpack_add_component(Runtime DISPLAY_NAME "ligint" DESCRIPTION "Programs, documentation and desktop integration" REQUIRED)
cpack_add_component(Development DISPLAY_NAME "Development files" DESCRIPTION "Headers, static library and CMake package")
