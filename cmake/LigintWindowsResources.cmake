# Adds the Windows resource file (application icon and version information) to an executable target.
# Does nothing on other platforms.
function(ligint_add_windows_resources target name description)
    if(NOT WIN32)
        return()
    endif()
    set(LIGINT_RC_NAME ${name})
    set(LIGINT_RC_DESCRIPTION ${description})
    set(LIGINT_ICON_PATH ${PROJECT_SOURCE_DIR}/packaging/ligint.ico)
    set(LIGINT_RC_VERSION "${PROJECT_VERSION_MAJOR},${PROJECT_VERSION_MINOR},${PROJECT_VERSION_PATCH},0")
    set(rc ${CMAKE_CURRENT_BINARY_DIR}/${name}.rc)
    configure_file(${PROJECT_SOURCE_DIR}/packaging/windows/ligint.rc.in ${rc} @ONLY)
    target_sources(${target} PRIVATE ${rc})
endfunction()
