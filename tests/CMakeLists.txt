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

if (NOT WIN32)
    add_cxx_flag_if_supported("-Wno-sign-compare")
endif()

set(GTEST_PREFIX ${PROJECT_SOURCE_DIR}/tests/gtest)
message(STATUS "NOTE: if adding the 'gtest' subdirectory fails, you need to issue 'git submodule init' and 'git submodule update'")
set(GTEST_BOTH_LIBRARIES gtest gtest_main)

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_definitions( -DCMS_TESTING_ENABLED )

# Add handy macros/functions
include(AddSTPGTest)
include(AddGTestSuite)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/)
include_directories(${PROJECT_BINARY_DIR}/apx-src)
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS})
include_directories(${ARJUN_INCLUDE_DIRS})

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
