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
include_directories(${GMP_INCLUDE_DIRS})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

if (NOT EMSCRIPTEN)
    add_executable(approxmc
        approxmc.cpp
        main.cpp
        constants.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
    )
else()
    add_executable(approxmc
        approxmc.cpp
        main_min.cpp
        constants.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
    )

    set_target_properties(
      approxmc PROPERTIES LINK_FLAGS "-s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\", \"cwrap\"]'  -s LINKABLE=1 -s EXPORT_ALL=1"
    )


    SET_SOURCE_FILES_PROPERTIES(main_min.cpp
        PROPERTIES COMPILE_FLAGS "-s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\", \"cwrap\"]'  -s LINKABLE=1 -s EXPORT_ALL=1"
    )
endif()


set(approxmc_exec_link_libs
    ${Boost_LIBRARIES}
    ${GMP_LIBRARY}
)

IF (ZLIB_FOUND)
    SET(approxmc_exec_link_libs ${approxmc_exec_link_libs} ${ZLIB_LIBRARY})
ENDIF()

target_link_libraries(approxmc
    ${approxmc_exec_link_libs}
)

if (STATICCOMPILE)
    SET_TARGET_PROPERTIES(approxmc PROPERTIES LINK_SEARCH_START_STATIC 1)
    SET(cms_libs ${CRYPTOMINISAT5_STATIC_LIBRARIES})
    message(STATUS "bimm: ${CRYPTOMINISAT5_STATIC_LIBRARIES}")
else()
    SET(cms_libs ${CRYPTOMINISAT5_LIBRARIES})
endif()
set_target_properties(approxmc PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE)
if (STATICCOMPILE)
    SET_TARGET_PROPERTIES(approxmc PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

target_link_libraries(approxmc
    ${cms_libs}
)


SET(CPACK_PACKAGE_EXECUTABLES "approxmc")



approxmc_add_public_header(approxmc ${CMAKE_CURRENT_SOURCE_DIR}/approxmc.h )



set(HEADER_DEST "${PROJECT_BINARY_DIR}/include/approxmc")

add_custom_target(CopyPublicHeaders ALL)

get_target_property(approxmc_public_headers approxmc PUBLIC_HEADER)

foreach(public_header ${approxmc_public_headers})

    get_filename_component(HEADER_NAME ${public_header} NAME)

    add_custom_command(TARGET CopyPublicHeaders PRE_BUILD

                       COMMAND ${CMAKE_COMMAND} -E make_directory

                               "${HEADER_DEST}"

                       COMMAND ${CMAKE_COMMAND} -E echo

                       "Copying ${HEADER_NAME} to ${HEADER_DEST}"

                       COMMAND ${CMAKE_COMMAND} -E

                           copy_if_different

                           ${public_header}

                           "${HEADER_DEST}"

                      )

endforeach()

install(TARGETS approxmc
    EXPORT ${APPROXMC_EXPORT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/approxmc"

)
