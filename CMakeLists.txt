cmake_minimum_required(VERSION 3.8)
project(bioptools)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Bioplib)
find_package(Genlib)
find_package(LibXml2 REQUIRED)

set(BINDIR $ENV{HOME}/bin CACHE PATH "Installation directory")
set(DATADIR $ENV{HOME}/data CACHE PATH "Data storage directory")


###############################################################################
#                        COMPILER AND LINKER FLAGS
###############################################################################

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ansi -Wall -pedantic")

###############################################################################
#                           DEFINE SOURCE FILES
###############################################################################

set(DATA_FILES
        ${bioptools_SOURCE_DIR}/data/radii.dat)

set(SOURCE_FILES
        src/naccess2bval.c
        src/pdbaddhet.c
        src/pdbatomcount.c
        src/pdbatomsel.c
        src/pdbavbr.c
        src/pdbcalcrms.c
        src/pdbcentralres.c
        src/pdbchain.c
        src/pdbcheckforres.c
        src/pdbcount.c
        src/pdbdummystrip.c
        src/pdbfindresrange.c
        src/pdbflip.c
        src/pdbgetchain.c
        src/pdbgetresidues.c
        src/pdbgetzone.c
        src/pdbhetstrip.c
        src/pdbmakepatch.c
        src/pdbpatchbval.c
        src/pdbpatchnumbering.c
        src/pdbrenum.c
        src/pdbrotate.c
        src/pdbsolv.c
        src/pdbsplitchains.c
        src/pdbsumbval.c
        src/pdbtorsions.c
        src/pdbtranslate.c
        src/setpdbnumbering.c)

###############################################################################
#                           DEFINE BINARY NAMES
###############################################################################

# final name of binaries
set(BINARY_NAMES
        as2bval
        addhet
        atomcount
        atomsel
        avbr
        rmspdb
        centralres
        chainpdb
        checkforres
        countpdb
        dummystrip
        findresrange
        flip
        getchain
        getresidues
        getpdb
        hetstrip
        makepatch
        patchbval
        patchpdbnum
        renumpdb
        rotate
        solv
        splitchains
        sumbval
        torsions
        transpdb
        numpdb
        )

include_directories(${BIOPLIB_INCLUDE} ${LIBXML2_INCLUDE_DIR})

###############################################################################
#                           BUILD RULES
###############################################################################

# get length of SOURCE_FILES list to iterate over
# SOURCE_FILES and BINARY_NAMES in parallel
list(LENGTH SOURCE_FILES len)
math(EXPR len "${len} - 1")

foreach(val RANGE ${len})
    list(GET SOURCE_FILES ${val} SOURCE_FILE)
    list(GET BINARY_NAMES ${val} BIN_NAME)
    add_executable(${BIN_NAME} ${SOURCE_FILE} ${BIOPLIB_INCLUDE})
    # link bioplib, genlib, lxml2 and math standard library (m)
    target_link_libraries(${BIN_NAME} ${BIOPLIB_LIBRARIES} ${GENLIB_LIBRARIES} ${LIBXML2_LIBRARIES} m)
endforeach()

###############################################################################
#                           MAKE RULES
###############################################################################

# message to print at the end of build -> `make`
add_custom_target( FinalMessage ALL
        ${CMAKE_COMMAND} -E cmake_echo_color --cyan
        " "
        "                   --- BUILD COMPLETE ---"
        " "
        "To use some programs, you will need to set the environment variable"
        " "
        "DATADIR to ${DATADIR}"
        " "
        "                        --- END ---"
        COMMENT "")
add_dependencies( FinalMessage ${BINARY_NAMES} )

# make install
install(TARGETS ${BINARY_NAMES} DESTINATION ${BINDIR} EXPORT myproj-targets)
install(FILES   ${DATA_FILES}  DESTINATION ${DATADIR})
install(EXPORT  myproj-targets DESTINATION ${BINDIR})
# message at the end of install
# next block comes from:
# https://stackoverflow.com/questions/18968979/how-to-get-colorized-output-with-cmake
if(NOT WIN32)
    string(ASCII 27 Esc)
    set(Cyan "${Esc}[36m")
endif()

install(CODE "MESSAGE(\"${Cyan}
                    --- INSTALL COMPLETE ---

To use some programs, you will need to set the environment variable

DATADIR to ${DATADIR}

                          --- END ---
\")")
