if(NOT VIF_FOUND)
    # Find vif headers
    find_path(VIF_INCLUDE_DIR vif.hpp
        HINTS ${VIF_ROOT_DIR} PATH_SUFFIXES include)

    find_path(VIF_COMPILER_DIR cvif
        HINTS ${VIF_ROOT_DIR} PATH_SUFFIXES bin)

    find_path(VIF_REFGEN_DIR vif-refgen
        HINTS ${VIF_ROOT_DIR} PATH_SUFFIXES bin)

    set(VIF_COMPILER ${VIF_COMPILER_DIR}/cvif)
    set(VIF_REFGEN ${VIF_REFGEN_DIR}/vif-refgen)

    mark_as_advanced(VIF_INCLUDE_DIR VIF_REFGEN VIF_COMPILER)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(VIF DEFAULT_MSG VIF_INCLUDE_DIR)

    set(VIF_INCLUDE_DIRS ${VIF_INCLUDE_DIR})

    # Configure compilers
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.3))
            message(STATUS "clang version >= 3.3 (${CMAKE_CXX_COMPILER_VERSION})")
        else()
            message(FATAL_ERROR "vif requires advanced features from the C++11 norm that are only available with clang 3.3 or higher (your version: ${CMAKE_CXX_COMPILER_VERSION}). Please upgrade your compiler.")
        endif()

        add_definitions(-Weverything)
        add_definitions(-Wno-c++98-compat-pedantic)
        add_definitions(-Wno-c++98-compat)
        add_definitions(-Wno-unused-parameter)
        add_definitions(-Wno-sign-conversion)
        add_definitions(-Wno-conversion)
        add_definitions(-Wno-missing-variable-declarations)
        add_definitions(-Wno-missing-prototypes)
        add_definitions(-Wno-padded)
        add_definitions(-Wno-float-equal)
        add_definitions(-Wno-unused-variable)
        add_definitions(-Wno-global-constructors)
        add_definitions(-Wno-exit-time-destructors)
        add_definitions(-Wno-weak-vtables)
        add_definitions(-Wno-covered-switch-default)
        add_definitions(-Wno-documentation-unknown-command)
        add_definitions(-Wno-unneeded-internal-declaration)
        add_definitions(-Wno-unused-function)
        add_definitions(-Wno-unused-macros)

        if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5))
            add_definitions(-Wno-old-style-cast)
        endif()
        if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8))
            add_definitions(-Wno-double-promotion)
        endif()

        add_definitions(-std=c++11)
        add_definitions(-ftemplate-backtrace-limit=0)
        add_definitions(-ferror-limit=5)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        if ("${CMAKE_CXX_COMPILER_VERSION}" STREQUAL "")
            message(WARNING "could not figure out the version of gcc, let's hope it is >= 4.7")
        else()
            if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7))
                message(STATUS "gcc version >= 4.7 (${CMAKE_CXX_COMPILER_VERSION})")
            else()
                message(FATAL_ERROR "vif requires advanced features from the C++11 norm that are only available with gcc 4.7 or higher (your version: ${CMAKE_CXX_COMPILER_VERSION}). Please upgrade your compiler.")
            endif()
        endif()

        add_definitions(-Wall)
        add_definitions(-std=c++11)
        add_definitions(-fmax-errors=5)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        if ("${CMAKE_CXX_COMPILER_VERSION}" STREQUAL "")
            message(WARNING "could not figure out the version of Intel compiler, let's hope it is >= 14")
        else()
            if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14))
                message(STATUS "Intel version >= 14 (${CMAKE_CXX_COMPILER_VERSION})")
            else()
                message(FATAL_ERROR "vif requires advanced features from the C++11 norm that are only available with Intel compiler 14 or higher (your version: ${CMAKE_CXX_COMPILER_VERSION}). Please upgrade your compiler.")
            endif()
        endif()

        add_definitions(-Wall)
        add_definitions(-std=c++11)
        add_definitions(-diag-error-limit=5)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
      message(ERROR "Microsoft Visual C++ compiler is not supported")
    endif()

    # find required libraries
    find_package(Threads REQUIRED)

    set(VIF_LIBRARIES ${VIF_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    # find optional libraries
    if (NOT NO_CFITSIO)
        find_package(CFITSIO)
    endif()
    if (NOT NO_LIBUNWIND)
        find_package(LibUnwind)
    endif()
    if (NOT NO_LIBDWARF)
        find_package(LibElf)
        find_package(LibDwarf)
    endif()
    if (NOT NO_FFTW)
        find_package(FFTW 3)
    endif()
    if (NOT NO_PROFILER OR NOT NO_TCMALLOC)
        find_package(GooglePerfTools)
    endif()
    if (NOT NO_LAPACK)
        find_package(LAPACK)
    endif()
    if (NOT NO_GSL AND NOT NO_LAPACK)
        find_package(GSL)
    endif()
    if (NOT NO_WCSLIB)
        find_package(WCSLib)
    endif()

    # Handle conditional reflection support
    # WIP: just disabled for now
    set(NO_REFLECTION 1)
    add_definitions(-DNO_REFLECTION)

    # Handle conditional CFITSIO support
    if (NOT CFITSIO_FOUND OR NO_CFITSIO)
        add_definitions(-DNO_CFITSIO)
    else()
        set(VIF_INCLUDE_DIRS ${VIF_INCLUDE_DIRS} ${CFITSIO_INCLUDES})
        set(VIF_LIBRARIES ${VIF_LIBRARIES} ${CFITSIO_LIBRARIES})
    endif()

    # Handle conditional LAPACK support
    if (NOT LAPACK_FOUND OR NO_LAPACK)
        add_definitions(-DNO_LAPACK)
    else()
        set(VIF_LIBRARIES ${VIF_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()

    # Handle conditional GSL support
    if (NOT GSL_FOUND OR NO_GSL OR NO_LAPACK OR NOT LAPACK_FOUND)
        add_definitions(-DNO_GSL)
    else()
        set(VIF_INCLUDE_DIRS ${VIF_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})
        set(VIF_LIBRARIES ${VIF_LIBRARIES} ${GSL_LIBRARIES})
    endif()

    # Handle conditional WCSLib support
    if (NOT WCSLIB_FOUND OR NO_WCSLIB OR NO_CFITSIO)
        add_definitions(-DNO_WCSLIB)
    else()
        set(VIF_INCLUDE_DIRS ${VIF_INCLUDE_DIRS} ${WCSLIB_INCLUDE_DIRS})
        set(VIF_LIBRARIES ${VIF_LIBRARIES} ${WCSLIB_LIBRARIES})

        if (WCSLIB_VERSION_STRING VERSION_LESS 5.0)
            add_definitions(-DWCSLIB_NO_DIS)
        endif()
    endif()

    # Handle conditional FFTW support
    if (NOT FFTW_FOUND OR NO_FFTW)
        add_definitions(-DNO_FFTW)
    else()
        set(VIF_INCLUDE_DIRS ${VIF_INCLUDE_DIRS} ${FFTW_INCLUDES})
        set(VIF_LIBRARIES ${VIF_LIBRARIES} ${FFTW_LIBRARIES})
    endif()

    # Handle conditional LibUnwind support
    if (NOT LIBUNWIND_FOUND OR NO_LIBUNWIND)
        set(NO_UNWIND 1)
        add_definitions(-DNO_LIBUNWIND)
    else()
        set(NO_UNWIND 0)
        set(VIF_INCLUDE_DIRS ${VIF_INCLUDE_DIRS} ${LIBUNWIND_INCLUDE_DIR})
        set(VIF_LIBRARIES ${VIF_LIBRARIES} ${LIBUNWIND_LIBRARIES})
    endif()

    # Handle conditional LibDwarf support
    if (NO_UNWIND OR NOT LIBELF_FOUND OR NOT LIBDWARF_FOUND OR NO_LIBDWARF)
        add_definitions(-DNO_LIBDWARF)
    else()
        set(VIF_INCLUDE_DIRS ${VIF_INCLUDE_DIRS} ${LIBDWARF_INCLUDE_DIRS})
        set(VIF_LIBRARIES ${VIF_LIBRARIES} ${LIBDWARF_LIBRARIES} ${LIBELF_LIBRARIES})
    endif()

    # Handle conditional Google perftools support
    if (TCMALLOC_LIBRARY)
        set(VIF_LIBRARIES ${VIF_LIBRARIES} ${TCMALLOC_LIBRARY})
    endif()
    if (PROFILER_LIBRARY)
        set(VIF_LIBRARIES ${VIF_LIBRARIES} ${PROFILER_LIBRARY})
    endif()
endif()

# Function to compile a vif program in CMake
# This feature does not support including/linking other external libraries, as well as
# "#define" commands, and is therefore not very powerful. But it is sufficient for basic
# needs. Supports reflection.
function(add_phypp_target CPP_FILE_NAME)
    # Generate binary name from c++ file
    get_filename_component(FILE_BASE ${CPP_FILE_NAME} NAME_WE)

    # Define the command to generate the binary file
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_custom_command(OUTPUT "${FILE_BASE}-make" VERBATIM COMMAND
            ${VIF_COMPILER} debug "${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME}" -o "${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make"
            DEPENDS ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME})
    else()
        add_custom_command(OUTPUT "${FILE_BASE}-make" VERBATIM COMMAND
            ${VIF_COMPILER} optimize "${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME}" -o "${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make"
            DEPENDS ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME})
    endif()

    # Create a target that will call this command
    add_custom_target(${FILE_BASE} ALL DEPENDS
        ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME}
        ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make)

    # Specify installation of the binary
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make
        DESTINATION bin RENAME ${FILE_BASE})
endfunction()
