# Copyright (c) 2017, Mate Soos
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

include_directories(
    ${PROJECT_SOURCE_DIR}
)

if (NOT WIN32)
    add_cxx_flag_if_supported("-Wno-bitfield-constant-conversion")
    #add_cxx_flag_if_supported("-Wduplicated-cond")
    #add_cxx_flag_if_supported("-Wduplicated-branches")
    add_cxx_flag_if_supported("-Wlogical-op")
    add_cxx_flag_if_supported("-Wrestrict")
    add_cxx_flag_if_supported("-Wnull-dereference")
    add_cxx_flag_if_supported("-Wjump-misses-init")
    add_cxx_flag_if_supported("-Wdouble-promotion")
    add_cxx_flag_if_supported("-Wshadow")
    add_cxx_flag_if_supported("-Wformat=2")
    add_cxx_flag_if_supported("-Wextra-semi")
    add_cxx_flag_if_supported("-pedantic")
    #add_cxx_flag_if_supported("-Wdeprecated")
endif()

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS})

add_executable(scalmc
    signalcode.cpp
    scalmc.cpp
)

set(scalmc_exec_link_libs
    ${Boost_LIBRARIES}
)

IF (ZLIB_FOUND)
    SET(scalmc_exec_link_libs ${scalmc_exec_link_libs} ${ZLIB_LIBRARY})
ENDIF()

target_link_libraries(scalmc
    ${scalmc_exec_link_libs}
)

if (STATICCOMPILE)
    SET_TARGET_PROPERTIES(scalmc PROPERTIES LINK_SEARCH_START_STATIC 1)
    SET(cms_libs ${CRYPTOMINISAT5_STATIC_LIBRARIES})
    message(STATUS "bimm: ${CRYPTOMINISAT5_STATIC_LIBRARIES}")
else()
    SET(cms_libs ${CRYPTOMINISAT5_LIBRARIES})
endif()
set_target_properties(scalmc PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE)
if (STATICCOMPILE)
    SET_TARGET_PROPERTIES(scalmc PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

target_link_libraries(scalmc
    ${cms_libs}
)

install(TARGETS scalmc
    EXPORT ${SCALMC_EXPORT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
SET(CPACK_PACKAGE_EXECUTABLES "scalmc")

# to compile statically:
# /usr/lib/ccache/bin/c++   -Wall -Wextra -Wunused -Wsign-compare -fno-omit-frame-pointer -Wtype-limits -Wuninitialized -Wno-deprecated -Wstrict-aliasing -Wpointer-arith -Wpointer-arith -Wformat-nonliteral -Winit-self -Wparentheses -Wunreachable-code -ggdb3 -Wlogical-op -Wrestrict -Wnull-dereference -Wdouble-promotion -Wshadow -Wformat=2 -pedantic    -Wl,--discard-all -Wl,--build-id=sha1 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -static    CMakeFiles/scalmc.dir/signalcode.cpp.o CMakeFiles/scalmc.dir/scalmc.cpp.o  -o ../scalmc /usr/lib/libboost_program_options.a /usr/lib/libz.a /usr/local/lib/libcryptominisat5.a /usr/local/lib/libm4ri.a -pthread

