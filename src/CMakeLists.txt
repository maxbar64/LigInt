# Libreria con l'analisi (nessuna dipendenza oltre alla libreria standard C++17)
add_library(ligint_core
    core/analysis.cpp
    core/chemistry.cpp
    core/report.cpp
    core/templates.cpp
    core/internal.hpp)
add_library(ligint::core ALIAS ligint_core)
target_include_directories(ligint_core
    PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE core)
target_compile_definitions(ligint_core PUBLIC LIGINT_VERSION="${PROJECT_VERSION}")
target_link_libraries(ligint_core PRIVATE ligint_options)
set_target_properties(ligint_core PROPERTIES POSITION_INDEPENDENT_CODE ON EXPORT_NAME core)   # esportato come ligint::core
install(TARGETS ligint_core ligint_options EXPORT ligintTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Programma a riga di comando
if(LIGINT_BUILD_CLI)
    add_executable(ligint cli/main.cpp)
    target_link_libraries(ligint PRIVATE ligint::core ligint_options)
    install(TARGETS ligint RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
