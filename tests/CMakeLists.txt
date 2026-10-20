# ---- C++ tests of the library (no external dependencies)
add_executable(ligint_core_tests unit/test_core.cpp)
target_link_libraries(ligint_core_tests PRIVATE ligint::core ligint_options)
target_compile_definitions(ligint_core_tests PRIVATE LIGINT_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
add_test(NAME core_unit COMMAND ligint_core_tests)

# ---- Python suite (chemical rules, SDF/MOL2, parsing, malformed input, real complexes in tests/data)
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND AND TARGET ligint)
    add_test(NAME python_suite
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.py
                --bin $<TARGET_FILE:ligint> --real-dir ${CMAKE_CURRENT_SOURCE_DIR}/data)
    set_tests_properties(python_suite PROPERTIES TIMEOUT 600)
endif()

if(Python3_FOUND)
    add_test(NAME gpcrdb_shared_fingerprints_offline
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_gpcrdb_shared_fingerprints.py)
    set_tests_properties(gpcrdb_shared_fingerprints_offline PROPERTIES TIMEOUT 30)
    add_test(NAME gpcrdb_auto_mapping_offline
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_gpcrdb_auto_mapping.py)
    set_tests_properties(gpcrdb_auto_mapping_offline PROPERTIES TIMEOUT 30)
endif()

# ---- GUI smoke test (Qt "offscreen" platform: no display needed)
if(TARGET ligint_gui_lib)
    add_executable(ligint_gui_smoke gui/gui_smoke.cpp)
    target_link_libraries(ligint_gui_smoke PRIVATE ligint_gui_lib ligint_options)
    if(TARGET ligint)
        # the second argument lets the test run the generated command lines with the real CLI
        add_test(NAME gui_smoke COMMAND ligint_gui_smoke ${CMAKE_CURRENT_SOURCE_DIR}/data/1AZM.pdb $<TARGET_FILE:ligint>)
    else()
        add_test(NAME gui_smoke COMMAND ligint_gui_smoke ${CMAKE_CURRENT_SOURCE_DIR}/data/1AZM.pdb)
    endif()
    set_tests_properties(gui_smoke PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 120)
endif()
