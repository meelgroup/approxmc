find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_MPFR QUIET mpfr)
endif()

if(APPLE)
    execute_process(
        COMMAND brew --prefix mpfr
        OUTPUT_VARIABLE HOMEBREW_MPFR_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
endif()

find_path(
    MPFR_INCLUDE_DIR
    NAMES mpfr.h
    HINTS ${PC_MPFR_INCLUDEDIR} ${HOMEBREW_MPFR_PREFIX}/include
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /opt/homebrew/include
          /usr/local/include
          /usr/include
)

set(_mpfr_old_suffixes "${CMAKE_FIND_LIBRARY_SUFFIXES}")
if(NOT BUILD_SHARED_LIBS)
    # Prefer static MPFR, but fall back to shared — mpfr-devel on RHEL/AlmaLinux
    # only ships libmpfr.so (no libmpfr.a).
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".dylib")
else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".dylib" ".a")
endif()

find_library(
    MPFR_LIBRARY
    NAMES mpfr
    HINTS ${PC_MPFR_LIBDIR} ${HOMEBREW_MPFR_PREFIX}/lib
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /opt/homebrew/lib
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
)

set(CMAKE_FIND_LIBRARY_SUFFIXES "${_mpfr_old_suffixes}")
set(MPFR_INCLUDE_DIRS ${MPFR_INCLUDE_DIR})
set(MPFR_LIBRARIES ${MPFR_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR
    DEFAULT_MSG
    MPFR_LIBRARY
    MPFR_INCLUDE_DIR)
mark_as_advanced(MPFR_LIBRARY MPFR_INCLUDE_DIR)
