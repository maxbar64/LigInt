# Per-generator settings; cpack reads this file once for each generator with CPACK_GENERATOR set.

if(CPACK_GENERATOR STREQUAL "DEB")
    set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
    set(CPACK_DEB_COMPONENT_INSTALL ON)
    set(CPACK_COMPONENTS_ALL Runtime Development)         # two packages: ligint and ligint-dev
    set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "${CPACK_LIGINT_MAINTAINER}")
    set(CPACK_DEBIAN_PACKAGE_SECTION "science")
    set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "${CPACK_PACKAGE_HOMEPAGE_URL}")
    set(CPACK_DEBIAN_COMPRESSION_TYPE "xz")
    set(CPACK_DEBIAN_RUNTIME_PACKAGE_NAME "ligint")
    set(CPACK_DEBIAN_RUNTIME_PACKAGE_SHLIBDEPS ON)        # Qt and C++ libraries found from the installed binaries
    if(CPACK_LIGINT_HAS_GUI)
        set(CPACK_DEBIAN_RUNTIME_PACKAGE_RECOMMENDS "qt6-qpa-plugins")   # platform plugins of Qt (xcb, wayland)
    endif()
    set(CPACK_DEBIAN_DEVELOPMENT_PACKAGE_NAME "ligint-dev")
    set(CPACK_DEBIAN_DEVELOPMENT_PACKAGE_SECTION "devel")
    if(CPACK_LIGINT_DEV_DEBDEPS)
        set(CPACK_DEBIAN_DEVELOPMENT_PACKAGE_DEPENDS "${CPACK_LIGINT_DEV_DEBDEPS}")
    endif()
    set(CPACK_COMPONENT_RUNTIME_DESCRIPTION "${CPACK_PACKAGE_DESCRIPTION}")
    set(CPACK_COMPONENT_DEVELOPMENT_DESCRIPTION
"Development files for ligint: headers, static library and CMake
package (find_package(ligint)).")

elseif(CPACK_GENERATOR STREQUAL "RPM")
    set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
    set(CPACK_RPM_COMPONENT_INSTALL ON)
    set(CPACK_COMPONENTS_ALL Runtime Development)
    set(CPACK_RPM_FILE_NAME RPM-DEFAULT)
    set(CPACK_RPM_PACKAGE_LICENSE "MIT")
    set(CPACK_RPM_PACKAGE_URL "${CPACK_PACKAGE_HOMEPAGE_URL}")
    set(CPACK_RPM_PACKAGE_GROUP "Applications/Science")
    set(CPACK_RPM_RUNTIME_PACKAGE_NAME "ligint")
    set(CPACK_RPM_DEVELOPMENT_PACKAGE_NAME "ligint-devel")
    if(CPACK_LIGINT_DEV_RPMDEPS)
        set(CPACK_RPM_DEVELOPMENT_PACKAGE_REQUIRES "${CPACK_LIGINT_DEV_RPMDEPS}")
    endif()
    set(CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST_ADDITION "/usr/share/man;/usr/share/man/man1;/usr/share/applications;/usr/share/icons;/usr/share/icons/hicolor;/usr/share/icons/hicolor/scalable;/usr/share/icons/hicolor/scalable/apps")

elseif(CPACK_GENERATOR STREQUAL "TGZ" OR CPACK_GENERATOR STREQUAL "ZIP")
    # portable archive: programs and documentation only (a single archive, without the development files)
    set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_LIGINT_BUILD_DIR};ligint;Runtime;/")
    set(CPACK_ARCHIVE_COMPONENT_INSTALL OFF)

elseif(CPACK_GENERATOR STREQUAL "NSIS")
    set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_LIGINT_BUILD_DIR};ligint;Runtime;/")   # end users do not need the development files
    set(CPACK_NSIS_COMPONENT_INSTALL OFF)
    string(REPLACE "/" "\\" _ico "${CPACK_LIGINT_ICON}")      # NSIS wants backslashes
    set(CPACK_NSIS_MUI_ICON "${_ico}")
    set(CPACK_NSIS_MUI_UNIICON "${_ico}")
    set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")
    set(CPACK_NSIS_DISPLAY_NAME "ligint ${CPACK_PACKAGE_VERSION}")
    set(CPACK_NSIS_PACKAGE_NAME "ligint")
    set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\ligint-gui.exe")
    set(CPACK_NSIS_URL_INFO_ABOUT "${CPACK_PACKAGE_HOMEPAGE_URL}")
    set(CPACK_NSIS_HELP_LINK "${CPACK_PACKAGE_HOMEPAGE_URL}")
    set(CPACK_NSIS_CONTACT "${CPACK_LIGINT_MAINTAINER}")
    set(CPACK_NSIS_MODIFY_PATH ON)                        # offers to add the programs to PATH, so `ligint` works in a terminal
    set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)    # an installed older version is removed first
    set(CPACK_NSIS_COMPRESSOR "/SOLID lzma")
    set(CPACK_PACKAGE_EXECUTABLES "ligint-gui" "ligint")  # Start menu shortcut
endif()
