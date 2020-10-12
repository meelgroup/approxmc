# Copyright (c) 2020, Mate Soos
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# -----------------------------------------------------------------------------
# Find GTest library which will be used to drive tests
# -----------------------------------------------------------------------------
# GoogleTest devs don't recommend using a pre-built GTest library
# ( https://code.google.com/p/googletest/wiki/FAQ#Why_is_it_not_recommended_to_install_a_pre-compiled_copy_of_Goog ).
# Because of this, distros like Ubuntu don't provide a pre-built GTest library
# so ``find_package(GTest REQUIRED)`` fails.
#
# Instead it is recommended that projects build their own copy of GTest. Detecting
# the location of GTest source code is probably error prone so using a copy in the
# repository seems like the easiest thing to do. This also has the added benefit that
# everyone uses the same version of GTest.
set(GTEST_PREFIX ${PROJECT_SOURCE_DIR}/utils/gtest)
message(STATUS "NOTE: if adding the 'gtest' subdirectory fails, you need to issue 'git submodule init' and 'git submodule update'")
add_subdirectory(${GTEST_PREFIX} gtest)
set(GTEST_BOTH_LIBRARIES gtest gtest_main)

include_directories(${GTEST_PREFIX}/include)
add_definitions( -DCMS_TESTING_ENABLED )

# Add handy macros/functions
include(AddSTPGTest)
include(AddGTestSuite)

# -----------------------------------------------------------------------------


include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/)
include_directories(${PROJECT_BINARY_DIR}/apx-src)
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS})

# README test
# add_executable(readme_test
#     readme_test.cpp
# )
# target_link_libraries(readme_test
#     cryptominisat5
# )
# add_test (
#     NAME readme_test
#     COMMAND readme_test
#     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
# )

# unit tests
set (MY_TESTS
    simpletest
)

foreach(F ${MY_TESTS})
    add_executable(${F}
        ${F}.cpp
    )
    target_link_libraries(${F}
        ${GMP_LIBRARY}
        ${GTEST_BOTH_LIBRARIES}
        approxmc
        ${CRYPTOMINISAT5_LIBRARIES}
    )
    add_test (
        NAME ${F}
        COMMAND ${F}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endforeach()
