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

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# -----------------------------------------------------------------------------
# Make RelWithDebInfo the default build type if otherwise not set
# -----------------------------------------------------------------------------
set(build_types Debug Release RelWithDebInfo MinSizeRel)
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "You can choose the type of build, options are:${build_types}")
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Options are ${build_types}"
        FORCE
    )

    # Provide drop down menu options in cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${build_types})
endif()
message(STATUS "Doing a ${CMAKE_BUILD_TYPE} build")

# -----------------------------------------------------------------------------
# Option to enable/disable assertions
# -----------------------------------------------------------------------------

# Filter out definition of NDEBUG from the default build configuration flags.
# We will add this ourselves if we want to disable assertions
foreach(build_config ${build_types})
    string(TOUPPER ${build_config} upper_case_build_config)
    foreach(language CXX C)
        set(VAR_TO_MODIFY "CMAKE_${language}_FLAGS_${upper_case_build_config}")
        string(REGEX REPLACE "(^| )[/-]D *NDEBUG($| )"
                             " "
                             replacement
                             "${${VAR_TO_MODIFY}}")
        set(${VAR_TO_MODIFY} "${replacement}" CACHE STRING "Default flags for ${build_config} configuration" FORCE)
    endforeach()
endforeach()

set(CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY ON)
set(CMAKE_EXPORT_PACKAGE_REGISTRY OFF)
project(approxmc)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(GenerateExportHeader)
include(GNUInstallDirs)

message(STATUS "LIB directory is '${CMAKE_INSTALL_LIBDIR}'")
message(STATUS "BIN directory is '${CMAKE_INSTALL_BINDIR}'")

# contains some library search cmake scripts
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# generate JSON file of compile commands -- useful for code extension
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# static compilation
option(BUILD_SHARED_LIBS "Build the shared library" ON)
# STATIC_BINARY controls whether the approxmc executable is linked fully
# statically.  Defaults to ON when BUILD_SHARED_LIBS=OFF (the historical
# behaviour), but can be set to OFF independently — e.g. when building the
# Python wheel: BUILD_SHARED_LIBS=OFF gives static dep .a files embedded in
# the .so, but we don't need (or want) a static approxmc binary.
if(NOT BUILD_SHARED_LIBS)
    option(STATIC_BINARY "Link the approxmc binary fully statically" ON)
else()
    option(STATIC_BINARY "Link the approxmc binary fully statically" OFF)
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if((${CMAKE_SYSTEM_NAME} MATCHES "Linux") OR (${CMAKE_SYSTEM_NAME} MATCHES "Darwin"))
    if(STATIC_BINARY AND ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        message(STATUS "Linking approxmc binary statically")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif()
endif()

option(SANITIZE "Use Clang sanitizers. You MUST use clang++ as the compiler for this to work" OFF)
if(SANITIZE)
    message(WARNING " --Using clang sanitizers -- you MUST use clang++ or the compile WILL fail")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lasan -lubsan")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lasan -lubsan")
    add_compile_options("-fsanitize=address")
    # below warns on overflows EVEN when that's OK, because it's well-defined, disabling
    # add_compile_options("-fsanitize=integer")
    add_compile_options("-fsanitize=undefined")

    add_compile_options("-fsanitize=null")
    add_compile_options("-fsanitize=alignment")
    #add_compile_options("-fno-sanitize-recover")

    add_compile_options("-fsanitize=return")
    add_compile_options("-fsanitize=bounds")
    add_compile_options("-fsanitize=float-divide-by-zero")
    add_compile_options("-fsanitize=integer-divide-by-zero")
    #add_compile_options("-fsanitize=unsigned-integer-overflow")
    add_compile_options("-fsanitize=signed-integer-overflow")
    add_compile_options("-fsanitize=bool")
    add_compile_options("-fsanitize=enum")
    add_compile_options("-fsanitize=float-cast-overflow")

    #add_compile_options("-Weverything")
    #add_compile_options("-Wshorten-64-to-32")
    #add_compile_options("-Wweak-vtables")
    #add_compile_options("-Wsign-conversion")
    #add_compile_options("-Wconversion")
endif()

include(CheckCXXCompilerFlag)
macro(add_cxx_flag_if_supported flagname)
  check_cxx_compiler_flag("${flagname}" HAVE_FLAG_${flagname})

  if(HAVE_FLAG_${flagname})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flagname}")
  endif()
endmacro()

option(ENABLE_ASSERTIONS "Build with assertions enabled" ON)
message(STATUS "build type is ${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(ENABLE_ASSERTIONS OFF)
endif()

if(ENABLE_ASSERTIONS)
    # NDEBUG was already removed.
else()
    # Note this definition doesn't appear in the cache variables.
    add_compile_definitions(NDEBUG)
    add_cxx_flag_if_supported("-fno-stack-protector")
    add_compile_definitions(_FORTIFY_SOURCE=0)
endif()


# Note: O3 gives slight speed increase, 1 more solved from SAT Comp'14 @ 3600s
if(NOT MSVC)
    add_compile_options(-g)

    add_compile_options("$<$<CONFIG:RelWithDebInfo>:-O3>")

    add_compile_options("$<$<CONFIG:Release>:-O3>")
    add_compile_options("$<$<CONFIG:Release>:-g0>")

    add_compile_options("$<$<CONFIG:Debug>:-O0>")

else()
    # see https://msdn.microsoft.com/en-us/library/fwkeyyhe.aspx for details
    # /ZI = include debug info
    # /Wall = all warnings

    add_compile_options("$<$<CONFIG:RELWITHDEBINFO>:/Ox>")
    add_compile_options("$<$<CONFIG:RELWITHDEBINFO>:/ZI>")

    add_compile_options("$<$<CONFIG:RELEASE>:/Ox>")
    add_compile_options("$<$<CONFIG:RELEASE>:/DNDEBUG>")
    add_compile_options("$<$<CONFIG:RELEASE>:/ZI>")

    add_compile_options("$<$<CONFIG:DEBUG>:/Od>")

    if(NOT BUILD_SHARED_LIBS)
        # We statically link to reduce dependencies
        foreach(flag_var CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
            # /MD -- Causes the application to use the multithread-specific and DLL-specific version of the run-time library.
            #        Defines _MT and _DLL and causes the compiler to place the library name MSVCRT.lib into the .obj file.
            if(${flag_var} MATCHES "/MD")
                string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
            endif()

            # /MDd	-- Defines _DEBUG, _MT, and _DLL and causes the application to use the debug multithread-specific and DLL-specific version of the run-time library.
            #          It also causes the compiler to place the library name MSVCRTD.lib into the .obj file.
            if(${flag_var} MATCHES "/MDd")
                string(REGEX REPLACE "/MDd" "/MTd" ${flag_var} "${${flag_var}}")
            endif()
        endforeach()

        # Creates a multithreaded executable (static) file using LIBCMT.lib.
        add_compile_options(/MT)
    endif()

    # buffers security check
    add_compile_options(/GS)

    # Proper warning level
    add_compile_options(/W1)

    # Disable STL used in DLL-boundary warning
    add_compile_options(/wd4251)
    add_compile_options(/D_CRT_SECURE_NO_WARNINGS)

    # Wall is MSVC's Weverything, so annoying unless used from the start
    # and with judiciously used warning disables
    # add_compile_options(/Wall)

    # /Za = only ansi C98 & C++11
    # /Za is not recommended for use, not tested, etc.
    # see: http://stackoverflow.com/questions/5489326/za-compiler-directive-does-not-compile-system-headers-in-vs2010
    # add_compile_options(/Za)

    add_compile_options(/fp:precise)

    # exception handling. s = The exception-handling model that catches C++ exceptions only and tells the compiler to assume that functions declared as extern "C" may throw an exception.
    # exception handling. c = If used with s (/EHsc), catches C++ exceptions only and tells the compiler to assume that functions declared as extern "C" never throw a C++ exception.
    add_compile_options(/EHsc)

    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /INCREMENTAL:NO")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /PDBCOMPRESS")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:1572864")

    set(DEF_INSTALL_CMAKE_DIR CMake)
endif()

option(ENABLE_TESTING "Enable testing" OFF)

if(NOT WIN32)
    add_cxx_flag_if_supported("-Wall")
    add_cxx_flag_if_supported("-Wextra")
    add_cxx_flag_if_supported("-Wunused")
    add_cxx_flag_if_supported("-Wsign-compare")
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
        add_cxx_flag_if_supported("-fno-omit-frame-pointer")
        add_cxx_flag_if_supported("-g")
    endif()
    add_cxx_flag_if_supported("-Wtype-limits")
    add_cxx_flag_if_supported("-Wuninitialized")
    add_cxx_flag_if_supported("-Wno-deprecated")
    add_cxx_flag_if_supported("-Wstrict-aliasing")
    add_cxx_flag_if_supported("-Wpointer-arith")
    add_cxx_flag_if_supported("-Wheader-guard")
    add_cxx_flag_if_supported("-msse4.2")
    if(NOT ENABLE_TESTING AND ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        add_cxx_flag_if_supported("-fvisibility=hidden")
    endif()
    add_cxx_flag_if_supported("-Wformat-nonliteral")
    add_cxx_flag_if_supported("-Winit-self")
    add_cxx_flag_if_supported("-Wparentheses")
    add_cxx_flag_if_supported("-Wunreachable-code")
endif()

if((NOT WIN32) AND (NOT ENABLE_TESTING))
    add_cxx_flag_if_supported("-fvisibility=hidden")
endif()

set(DEF_INSTALL_CMAKE_DIR lib/cmake/approxmc)
set(APPROXMC_INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
    "Installation directory for approxmc CMake files")

# -----------------------------------------------------------------------------
# Add GIT version
# -----------------------------------------------------------------------------
function(SetVersionNumber PREFIX VERSION_MAJOR VERSION_MINOR VERSION_PATCH)
  set(${PREFIX}_VERSION_MAJOR ${VERSION_MAJOR} PARENT_SCOPE)
  set(${PREFIX}_VERSION_MINOR ${VERSION_MINOR} PARENT_SCOPE)
  set(${PREFIX}_VERSION_PATCH ${VERSION_PATCH} PARENT_SCOPE)
  set(${PREFIX}_VERSION
        "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}"
        PARENT_SCOPE)
endfunction()

find_program(GIT_EXECUTABLE git)
if(GIT_EXECUTABLE)
  include(GetGitRevisionDescription)
  get_git_head_revision(GIT_REFSPEC GIT_SHA1)
  message(STATUS "GIT hash found: ${GIT_SHA1}")
else()
  set(GIT_SHA1 "GIT-hash-notfound")
endif()
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml" PYPROJECT_TOML)
string(REGEX MATCH "version = \"([0-9]+\\.[0-9]+\\.[0-9]+)\"" _ "${PYPROJECT_TOML}")
set(APPROXMC_FULL_VERSION "${CMAKE_MATCH_1}")
if(NOT APPROXMC_FULL_VERSION)
    message(FATAL_ERROR "Could not parse version from pyproject.toml")
endif()

string(REPLACE "." ";" APPROXMC_FULL_VERSION_LIST ${APPROXMC_FULL_VERSION})
SetVersionNumber("PROJECT" ${APPROXMC_FULL_VERSION_LIST})
message(STATUS "PROJECT_VERSION: ${PROJECT_VERSION}")
message(STATUS "PROJECT_VERSION_MAJOR: ${PROJECT_VERSION_MAJOR}")
message(STATUS "PROJECT_VERSION_MINOR: ${PROJECT_VERSION_MINOR}")
message(STATUS "PROJECT_VERSION_PATCH: ${PROJECT_VERSION_PATCH}")

# ----------
# manpage
# ----------
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND NOT EMSCRIPTEN)
    find_program(HELP2MAN_FOUND help2man)
    if(HELP2MAN_FOUND)
        add_custom_target(man_approxmc
            ALL
            DEPENDS approxmc-bin
        )

        add_custom_command(
            TARGET man_approxmc
            POST_BUILD
            COMMAND help2man
            ARGS --version-string=${APPROXMC_FULL_VERSION} --help-option="-h" $<TARGET_FILE:approxmc-bin> -o ${CMAKE_CURRENT_BINARY_DIR}/approxmc.1
        )

        install(
            FILES ${CMAKE_CURRENT_BINARY_DIR}/approxmc.1
            DESTINATION ${CMAKE_INSTALL_PREFIX}/man/man1)

        message(STATUS "Manpage will be created and installed")
    else()
        message(STATUS "Cannot find help2man, not creating manpage")
    endif()
else()
    message(STATUS "Not on Linux, not creating manpage")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Collect all compile definitions set so far into a string for embedding
# via @COMPILE_DEFINES@ in GitSHA1.cpp.in
get_directory_property(DirDefs COMPILE_DEFINITIONS)
set(COMPILE_DEFINES "")
foreach(d ${DirDefs})
    set(COMPILE_DEFINES "${COMPILE_DEFINES} -D${d}")
endforeach()
message(STATUS "All defines at startup: ${COMPILE_DEFINES}")

macro(approxmc_add_public_header LIBTARGET HEADER)
    get_target_property(EXISTING_PUBLIC_HEADERS ${LIBTARGET} PUBLIC_HEADER)
    if(EXISTING_PUBLIC_HEADERS)
        list(APPEND EXISTING_PUBLIC_HEADERS "${HEADER}")
    else()
        # Do not append to empty list
        set(EXISTING_PUBLIC_HEADERS "${HEADER}")
    endif()
    set_target_properties(
        ${LIBTARGET}
        PROPERTIES
        PUBLIC_HEADER "${EXISTING_PUBLIC_HEADERS}"
     )
endmacro()

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
# When building statically, prefer .a over .so for all find_library/find_package calls
if(NOT BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()
find_package(PkgConfig REQUIRED)
if(APPLE)
    execute_process(COMMAND brew --prefix OUTPUT_VARIABLE _homebrew_prefix OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(_homebrew_prefix)
        set(ENV{PKG_CONFIG_PATH} "${_homebrew_prefix}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
endif()
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)
find_library(GMPXX_LIBRARY NAMES gmpxx HINTS ${GMP_LIBRARY_DIRS})
if(GMPXX_LIBRARY)
    get_property(_gmp_iface_libs TARGET PkgConfig::GMP PROPERTY INTERFACE_LINK_LIBRARIES)
    set_property(TARGET PkgConfig::GMP PROPERTY INTERFACE_LINK_LIBRARIES "${GMPXX_LIBRARY}" ${_gmp_iface_libs})
endif()
pkg_check_modules(MPFR REQUIRED IMPORTED_TARGET mpfr)

include(FetchContent)
# FetchContent_Populate with declared details is deprecated in CMake 3.30+ (CMP0169).
# We call it intentionally to clone-only without add_subdirectory, so silence the warning.
if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
endif()


# ---- CryptoMiniSat5 and its deps (cadical, cadiback) ----
# cadical and cadiback are only needed when auto-building CryptoMiniSat5.
# cadiback's configure looks for ../cadical, so both must be siblings under deps/.
set(cadical_DIR "" CACHE PATH "Directory where CaDiCaL is built (contains libcadical.a). Auto-fetched if empty.")
set(cadiback_DIR "" CACHE PATH "Directory where CaDiBaCk is built (contains libcadiback.a). Auto-fetched if empty.")
set(cryptominisat5_DIR "" CACHE PATH "Directory where CryptoMiniSat5 is built. Auto-fetched if empty.")
if(NOT cryptominisat5_DIR)
    if(NOT cadical_DIR)
        set(_cadical_src "${CMAKE_BINARY_DIR}/deps/cadical")
        set(_cadical_lib "${_cadical_src}/build/libcadical.a")
        FetchContent_Declare(cadical_src
            GIT_REPOSITORY https://github.com/meelgroup/cadical
            GIT_TAG        master
            SOURCE_DIR     "${_cadical_src}"
            GIT_SHALLOW    TRUE
        )
        FetchContent_GetProperties(cadical_src)
        if(NOT cadical_src_POPULATED)
            message(STATUS "CaDiCaL: fetching from meelgroup/cadical master")
            FetchContent_Populate(cadical_src)
        endif()
        if(NOT EXISTS "${_cadical_lib}")
            message(STATUS "CaDiCaL: building in ${_cadical_src}")
            execute_process(
                COMMAND bash -c "CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=\"-fPIC\" ./configure --no-contracts --no-tracing --fp"
                WORKING_DIRECTORY "${_cadical_src}"
                RESULT_VARIABLE _cadical_cfg_ret
            )
            if(NOT _cadical_cfg_ret EQUAL 0)
                message(FATAL_ERROR "CaDiCaL configure failed (exit ${_cadical_cfg_ret})")
            endif()
            execute_process(
                COMMAND make -j8
                WORKING_DIRECTORY "${_cadical_src}"
                RESULT_VARIABLE _cadical_make_ret
            )
            if(NOT _cadical_make_ret EQUAL 0)
                message(FATAL_ERROR "CaDiCaL build failed (exit ${_cadical_make_ret})")
            endif()
        endif()
        set(cadical_DIR "${_cadical_src}/build" CACHE PATH "Directory where CaDiCaL is built" FORCE)
        set(cadical_SRC_DIR "${_cadical_src}/src" CACHE PATH "CaDiCaL source dir (for cadical.hpp)" FORCE)
        message(STATUS "CaDiCaL: using auto-built ${cadical_DIR}")
    endif()
    if(NOT cadiback_DIR)
        set(_cadiback_src "${CMAKE_BINARY_DIR}/deps/cadiback")
        set(_cadiback_lib "${_cadiback_src}/libcadiback.a")
        FetchContent_Declare(cadiback_src
            GIT_REPOSITORY https://github.com/meelgroup/cadiback
            GIT_TAG        main
            SOURCE_DIR     "${_cadiback_src}"
            GIT_SHALLOW    TRUE
        )
        FetchContent_GetProperties(cadiback_src)
        if(NOT cadiback_src_POPULATED)
            message(STATUS "CaDiBaCk: fetching from meelgroup/cadiback main")
            FetchContent_Populate(cadiback_src)
        endif()
        if(NOT EXISTS "${_cadiback_lib}")
            message(STATUS "CaDiBaCk: building in ${_cadiback_src}")
            execute_process(
                COMMAND bash -c "CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=\"-fPIC\" ./configure --release"
                WORKING_DIRECTORY "${_cadiback_src}"
                RESULT_VARIABLE _cadiback_cfg_ret
            )
            if(NOT _cadiback_cfg_ret EQUAL 0)
                message(FATAL_ERROR "CaDiBaCk configure failed (exit ${_cadiback_cfg_ret})")
            endif()
            execute_process(
                COMMAND make -j8 libcadiback.a
                WORKING_DIRECTORY "${_cadiback_src}"
                RESULT_VARIABLE _cadiback_make_ret
            )
            if(NOT _cadiback_make_ret EQUAL 0)
                message(FATAL_ERROR "CaDiBaCk build failed (exit ${_cadiback_make_ret})")
            endif()
        endif()
        set(cadiback_DIR "${_cadiback_src}" CACHE PATH "Directory where CaDiBaCk is built" FORCE)
        message(STATUS "CaDiBaCk: using auto-built ${cadiback_DIR}")
    endif()
    set(_cms_src "${CMAKE_BINARY_DIR}/deps/cryptominisat")
    set(_cms_build "${_cms_src}/build")
    FetchContent_Declare(cryptominisat5_src
        GIT_REPOSITORY https://github.com/msoos/cryptominisat
        GIT_TAG        master
        SOURCE_DIR     "${_cms_src}"
        GIT_SHALLOW    TRUE
    )
    FetchContent_GetProperties(cryptominisat5_src)
    if(NOT cryptominisat5_src_POPULATED)
        message(STATUS "CryptoMiniSat5: fetching from msoos/cryptominisat master")
        FetchContent_Populate(cryptominisat5_src)
    endif()
    if(NOT EXISTS "${_cms_build}/cryptominisat5Config.cmake")
        message(STATUS "CryptoMiniSat5: configuring and building in ${_cms_build}")
        execute_process(
            COMMAND ${CMAKE_COMMAND} -S "${_cms_src}" -B "${_cms_build}"
                -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                -Dcadical_DIR=${cadical_DIR}
                -Dcadiback_DIR=${cadiback_DIR}
            RESULT_VARIABLE _cms_cfg_ret
        )
        if(NOT _cms_cfg_ret EQUAL 0)
            message(FATAL_ERROR "CryptoMiniSat5 cmake configure failed (exit ${_cms_cfg_ret})")
        endif()
        execute_process(
            COMMAND ${CMAKE_COMMAND} --build "${_cms_build}" -j8
            RESULT_VARIABLE _cms_build_ret
        )
        if(NOT _cms_build_ret EQUAL 0)
            message(FATAL_ERROR "CryptoMiniSat5 build failed (exit ${_cms_build_ret})")
        endif()
    endif()
    set(cryptominisat5_DIR "${_cms_build}" CACHE PATH "Directory where CryptoMiniSat5 is built" FORCE)
    message(STATUS "CryptoMiniSat5: using auto-built ${cryptominisat5_DIR}")
endif()
find_package(cryptominisat5 CONFIG REQUIRED)
message(STATUS "CryptoMiniSat5 dynamic lib: ${CRYPTOMINISAT5_LIBRARIES}")
message(STATUS "CryptoMiniSat5 static lib:  ${CRYPTOMINISAT5_STATIC_LIBRARIES}")
message(STATUS "CryptoMiniSat5 static lib deps: ${CRYPTOMINISAT5_STATIC_LIBRARIES_DEPS}")
message(STATUS "CryptoMiniSat5 include dirs: ${CRYPTOMINISAT5_INCLUDE_DIRS}")

# ---- Arjun and its deps (sbva, treedecomp) ----
# sbva and treedecomp are only needed when auto-building Arjun.
set(sbva_DIR "" CACHE PATH "Directory where sbva is built. Auto-fetched if empty.")
set(treedecomp_DIR "" CACHE PATH "Directory where treedecomp is built. Auto-fetched if empty.")
set(arjun_DIR "" CACHE PATH "Directory where Arjun is built. Auto-fetched if empty.")
if(NOT arjun_DIR)
    if(NOT sbva_DIR)
        set(_sbva_src "${CMAKE_BINARY_DIR}/deps/sbva")
        set(_sbva_build "${_sbva_src}/build")
        FetchContent_Declare(sbva_src
            GIT_REPOSITORY https://github.com/meelgroup/SBVA.git
            GIT_TAG        master
            SOURCE_DIR     "${_sbva_src}"
            GIT_SHALLOW    TRUE
        )
        FetchContent_GetProperties(sbva_src)
        if(NOT sbva_src_POPULATED)
            message(STATUS "SBVA: fetching from meelgroup/SBVA master")
            FetchContent_Populate(sbva_src)
        endif()
        if(NOT EXISTS "${_sbva_build}/sbvaConfig.cmake")
            message(STATUS "SBVA: configuring and building in ${_sbva_build}")
            execute_process(
                COMMAND ${CMAKE_COMMAND} -S "${_sbva_src}" -B "${_sbva_build}"
                    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                RESULT_VARIABLE _sbva_cfg_ret
            )
            if(NOT _sbva_cfg_ret EQUAL 0)
                message(FATAL_ERROR "SBVA cmake configure failed (exit ${_sbva_cfg_ret})")
            endif()
            execute_process(
                COMMAND ${CMAKE_COMMAND} --build "${_sbva_build}" -j8
                RESULT_VARIABLE _sbva_build_ret
            )
            if(NOT _sbva_build_ret EQUAL 0)
                message(FATAL_ERROR "SBVA build failed (exit ${_sbva_build_ret})")
            endif()
        endif()
        set(sbva_DIR "${_sbva_build}" CACHE PATH "Directory where sbva is built" FORCE)
        message(STATUS "SBVA: using auto-built ${sbva_DIR}")
    endif()
    if(NOT treedecomp_DIR)
        set(_tdec_src "${CMAKE_BINARY_DIR}/deps/treedecomp")
        set(_tdec_build "${_tdec_src}/build")
        FetchContent_Declare(treedecomp_src
            GIT_REPOSITORY https://github.com/meelgroup/treedecomp.git
            GIT_TAG        main
            SOURCE_DIR     "${_tdec_src}"
            GIT_SHALLOW    TRUE
        )
        FetchContent_GetProperties(treedecomp_src)
        if(NOT treedecomp_src_POPULATED)
            message(STATUS "treedecomp: fetching from meelgroup/treedecomp main")
            FetchContent_Populate(treedecomp_src)
        endif()
        if(NOT EXISTS "${_tdec_build}/treedecompConfig.cmake")
            message(STATUS "treedecomp: configuring and building in ${_tdec_build}")
            execute_process(
                COMMAND ${CMAKE_COMMAND} -S "${_tdec_src}" -B "${_tdec_build}"
                    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                RESULT_VARIABLE _tdec_cfg_ret
            )
            if(NOT _tdec_cfg_ret EQUAL 0)
                message(FATAL_ERROR "treedecomp cmake configure failed (exit ${_tdec_cfg_ret})")
            endif()
            execute_process(
                COMMAND ${CMAKE_COMMAND} --build "${_tdec_build}" -j8
                RESULT_VARIABLE _tdec_build_ret
            )
            if(NOT _tdec_build_ret EQUAL 0)
                message(FATAL_ERROR "treedecomp build failed (exit ${_tdec_build_ret})")
            endif()
        endif()
        set(treedecomp_DIR "${_tdec_build}" CACHE PATH "Directory where treedecomp is built" FORCE)
        message(STATUS "treedecomp: using auto-built ${treedecomp_DIR}")
    endif()
    set(_arjun_src "${CMAKE_BINARY_DIR}/deps/arjun")
    set(_arjun_build "${_arjun_src}/build")
    FetchContent_Declare(arjun_src
        GIT_REPOSITORY https://github.com/meelgroup/arjun.git
        GIT_TAG        master
        SOURCE_DIR     "${_arjun_src}"
        GIT_SHALLOW    TRUE
    )
    FetchContent_GetProperties(arjun_src)
    if(NOT arjun_src_POPULATED)
        message(STATUS "Arjun: fetching from meelgroup/arjun master")
        FetchContent_Populate(arjun_src)
    endif()
    if(NOT EXISTS "${_arjun_build}/arjunConfig.cmake")
        message(STATUS "Arjun: configuring and building in ${_arjun_build}")
        execute_process(
            COMMAND ${CMAKE_COMMAND} -S "${_arjun_src}" -B "${_arjun_build}"
                -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                -Dcadical_DIR=${cadical_DIR}
                -Dcryptominisat5_DIR=${cryptominisat5_DIR}
                -Dsbva_DIR=${sbva_DIR}
                -Dtreedecomp_DIR=${treedecomp_DIR}
            RESULT_VARIABLE _arjun_cfg_ret
        )
        if(NOT _arjun_cfg_ret EQUAL 0)
            message(FATAL_ERROR "Arjun cmake configure failed (exit ${_arjun_cfg_ret})")
        endif()
        execute_process(
            COMMAND ${CMAKE_COMMAND} --build "${_arjun_build}" -j8
            RESULT_VARIABLE _arjun_build_ret
        )
        if(NOT _arjun_build_ret EQUAL 0)
            message(FATAL_ERROR "Arjun build failed (exit ${_arjun_build_ret})")
        endif()
    endif()
    set(arjun_DIR "${_arjun_build}" CACHE PATH "Directory where Arjun is built" FORCE)
    message(STATUS "Arjun: using auto-built ${arjun_DIR}")
endif()
find_package(arjun CONFIG REQUIRED)
message(STATUS "Arjun dynamic lib : ${ARJUN_LIBRARIES}")
message(STATUS "Arjun static lib  : ${ARJUN_STATIC_LIBRARIES}")
message(STATUS "Arjun include dirs: ${ARJUN_INCLUDE_DIRS}")

# ---- EvalMaxSAT (cmake-based, optional) ----
option(USE_EVALMAXSAT "Build with EvalMaxSAT support" OFF)
if(USE_EVALMAXSAT)
    set(EvalMaxSAT_DIR "" CACHE PATH "Directory where EvalMaxSAT is built. Auto-fetched if empty.")
    if(NOT EvalMaxSAT_DIR)
        set(_evalmaxsat_src "${CMAKE_BINARY_DIR}/deps/EvalMaxSAT")
        set(_evalmaxsat_build "${_evalmaxsat_src}/build")
        FetchContent_Declare(evalmaxsat_src
            GIT_REPOSITORY https://github.com/meelgroup/EvalMaxSAT.git
            GIT_TAG        master
            SOURCE_DIR     "${_evalmaxsat_src}"
            GIT_SHALLOW    TRUE
        )
        FetchContent_GetProperties(evalmaxsat_src)
        if(NOT evalmaxsat_src_POPULATED)
            message(STATUS "EvalMaxSAT: fetching from meelgroup/EvalMaxSAT master")
            FetchContent_Populate(evalmaxsat_src)
        endif()
        if(NOT EXISTS "${_evalmaxsat_build}/EvalMaxSATConfig.cmake")
            message(STATUS "EvalMaxSAT: configuring and building in ${_evalmaxsat_build}")
            execute_process(
                COMMAND ${CMAKE_COMMAND} -S "${_evalmaxsat_src}" -B "${_evalmaxsat_build}"
                    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                RESULT_VARIABLE _evalmaxsat_cfg_ret
            )
            if(NOT _evalmaxsat_cfg_ret EQUAL 0)
                message(FATAL_ERROR "EvalMaxSAT cmake configure failed (exit ${_evalmaxsat_cfg_ret})")
            endif()
            execute_process(
                COMMAND ${CMAKE_COMMAND} --build "${_evalmaxsat_build}" -j8
                RESULT_VARIABLE _evalmaxsat_build_ret
            )
            if(NOT _evalmaxsat_build_ret EQUAL 0)
                message(FATAL_ERROR "EvalMaxSAT build failed (exit ${_evalmaxsat_build_ret})")
            endif()
        endif()
        set(EvalMaxSAT_DIR "${_evalmaxsat_build}" CACHE PATH "Directory where EvalMaxSAT is built" FORCE)
        message(STATUS "EvalMaxSAT: using auto-built ${EvalMaxSAT_DIR}")
    endif()
    find_package(EvalMaxSAT REQUIRED)
endif()

option(NOZLIB "Don't use zlib" OFF)
if(NOT NOZLIB AND NOT (NOT BUILD_SHARED_LIBS AND WIN32))
    find_package(ZLIB)
    if(ZLIB_FOUND)
        message(STATUS "OK, Found ZLIB!")
        add_compile_definitions(USE_ZLIB)
    else()
        message(STATUS "WARNING: Did not find ZLIB, gzipped file support will be disabled")
    endif()
endif()

# -----------------------------------------------------------------------------
# Provide an export name to be used by targets that wish to export themselves.
# -----------------------------------------------------------------------------
set(APPROXMC_EXPORT_NAME "approxmcTargets")

if(ENABLE_TESTING)
    add_subdirectory(tests/gtest)
endif()
add_subdirectory(src apx-src)

# =============================================================================
# Python extension module (pyapproxmc)
# Activated by scikit-build-core via -DBUILD_PYTHON_EXTENSION=ON.
# Requires CMake >= 3.18 (python_add_library).
# =============================================================================
option(BUILD_PYTHON_EXTENSION "Build the pyapproxmc Python extension module" OFF)
if(BUILD_PYTHON_EXTENSION)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

    # GitSHA1 symbols are already provided by the approxmc static library.
    python_add_library(pyapproxmc MODULE WITH_SOABI
        python/src/pyapproxmc.cpp
    )

    target_compile_definitions(pyapproxmc PRIVATE
        APPMC_FULL_VERSION="${APPROXMC_FULL_VERSION}"
    )

    target_include_directories(pyapproxmc PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${CRYPTOMINISAT5_INCLUDE_DIRS}
        ${ARJUN_INCLUDE_DIRS}
    )

    target_link_libraries(pyapproxmc PRIVATE
        approxmc
        ${ARJUN_STATIC_LIBRARIES}
        PkgConfig::GMP
        PkgConfig::MPFR
    )
    target_compile_features(pyapproxmc PRIVATE cxx_std_20)

    # Only install the extension into the wheel (not the library/headers/binary).
    install(TARGETS pyapproxmc DESTINATION . COMPONENT python)
endif()

# -----------------------------------------------------------------------------
# Add uninstall target for makefiles
# -----------------------------------------------------------------------------
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)

# -----------------------------------------------------------------------------
# Testing
# -----------------------------------------------------------------------------
if(ENABLE_TESTING)
    enable_testing()

    message(STATUS "Testing is enabled")
    set(UNIT_TEST_EXE_SUFFIX "Tests" CACHE STRING "Suffix for Unit test executable")
    add_subdirectory(tests)

else()
    message(WARNING "Testing is disabled")
endif()

# -----------------------------------------------------------------------------
# Export our targets so that other CMake based projects can interface with
# the build of approxmc in the build-tree
# -----------------------------------------------------------------------------
set(APPROXMC_TARGETS_FILENAME "approxmcTargets.cmake")
set(APPROXMC_CONFIG_FILENAME "approxmcConfig.cmake")

# Export targets
export(
    TARGETS approxmc
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${APPROXMC_TARGETS_FILENAME}"
)

# Build-tree config (absolute paths, for use without installation)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/approxmcConfig-buildtree.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${APPROXMC_CONFIG_FILENAME}"
    @ONLY
)

# Export this package to the CMake user package registry
# Now the user can just use find_package(approxmc) on their system
export(PACKAGE approxmc)

# Install-tree config and version files (relocatable via CMakePackageConfigHelpers)
set(APPROXMC_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}")
set(APPROXMC_VERSION_FILENAME "approxmcConfigVersion.cmake")

include(CMakePackageConfigHelpers)
configure_package_config_file(approxmcConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${APPROXMC_CONFIG_FILENAME}"
    INSTALL_DESTINATION "${APPROXMC_INSTALL_CMAKE_DIR}"
    PATH_VARS APPROXMC_INSTALL_INCLUDEDIR
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${APPROXMC_VERSION_FILENAME}"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${APPROXMC_CONFIG_FILENAME}"
    "${CMAKE_CURRENT_BINARY_DIR}/${APPROXMC_VERSION_FILENAME}"
    DESTINATION "${APPROXMC_INSTALL_CMAKE_DIR}"
)

# Install the export set for use with the install-tree
install(
    EXPORT ${APPROXMC_EXPORT_NAME}
    DESTINATION "${APPROXMC_INSTALL_CMAKE_DIR}"
)
