# CMake generated Testfile for 
# Source directory: /home/massimo/devel/LigInt/tests
# Build directory: /home/massimo/devel/LigInt/build/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(core_unit "/home/massimo/devel/LigInt/build/tests/ligint_core_tests")
set_tests_properties(core_unit PROPERTIES  _BACKTRACE_TRIPLES "/home/massimo/devel/LigInt/tests/CMakeLists.txt;4;add_test;/home/massimo/devel/LigInt/tests/CMakeLists.txt;0;")
add_test(python_suite "/home/massimo/miniconda3/bin/python3" "/home/massimo/devel/LigInt/tests/run_tests.py" "--bin" "/home/massimo/devel/LigInt/build/src/ligint" "--real-dir" "/home/massimo/devel/LigInt/tests/data")
set_tests_properties(python_suite PROPERTIES  TIMEOUT "600" _BACKTRACE_TRIPLES "/home/massimo/devel/LigInt/tests/CMakeLists.txt;9;add_test;/home/massimo/devel/LigInt/tests/CMakeLists.txt;0;")
