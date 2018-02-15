# - this module looks for bioplib
# Defines:
#   BIOPLIB_LIBRARIES
#   BIOPLIB_FOUND
#   BIOPLIB_INCLUDE
#   BIOPLIB_VERSION
#   BIOPLIB_MAJOR_VERSION
#   BIOPLIB_MINOR_VERSION

find_package(PackageHandleStandardArgs)

find_library(BIOPLIB_LIBRARIES biops HINTS $ENV{HOME}/lib/)

get_filename_component(BIOPLIB_REAL_PATH ${BIOPLIB_LIBRARIES} REALPATH)

string(REGEX MATCH "[0-9]+" BIOPLIB_MAJOR_VERSION ${BIOPLIB_REAL_PATH})
string(REGEX MATCH "[0-9]+.[0-9]+$" BIOPLIB_MINOR_VERSION ${BIOPLIB_REAL_PATH})

find_path(BIOPLIB_INCLUDE 00bioplib.c HINTS $ENV{HOME}/bioplib/src)

set(BIOPLIB_VERSION ${BIOPLIB_MAJOR_VERSION}.${BIOPLIB_MINOR_VERSION})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(bioplib
        REQUIRED_VARS BIOPLIB_LIBRARIES BIOPLIB_INCLUDE
        VERSION_VAR BIOPLIB_VERSION)

MARK_AS_ADVANCED (
        BIOPLIB_LIBRARIES
        BIOPLIB_INCLUDE
        BIOPLIB_VERSION
        BIOPLIB_MAJOR_VERSION
        BIOPLIB_MINOR_VERSION
)