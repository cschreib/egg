if(NOT PHYPP_FOUND)
    find_path(PHYPP_INCLUDE_DIR phypp.hpp
        HINTS ${PHYPP_ROOT_DIR} PATH_SUFFIXES include)

    file(GLOB PHYPP_HEADERS ${PHYPP_INCLUDE_DIR}/phypp/*.hpp)

    find_path(PHYPP_COMPILER_DIR cphy++
        HINTS ${PHYPP_ROOT_DIR} PATH_SUFFIXES bin)

    find_path(PHYPP_REFGEN_DIR phy++-refgen
        HINTS ${PHYPP_ROOT_DIR} PATH_SUFFIXES bin)

    set(PHYPP_COMPILER ${PHYPP_COMPILER_DIR}/cphy++)
    set(PHYPP_REFGEN ${PHYPP_REFGEN_DIR}/phy++-refgen)

    mark_as_advanced(PHYPP_INCLUDE_DIR PHYPP_HEADERS PHYPP_REFGEN PHYPP_COMPILER)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(PHYPP DEFAULT_MSG
        PHYPP_INCLUDE_DIR PHYPP_HEADERS PHYPP_REFGEN PHYPP_COMPILER)

    set(PHYPP_INCLUDE_DIRS ${PHYPP_INCLUDE_DIR})
endif(NOT PHYPP_FOUND)

# Function to compile a phy++ program in CMake
function(add_phypp_target CPP_FILE_NAME)
    # Generate binary name from c++ file
    get_filename_component(FILE_BASE ${CPP_FILE_NAME} NAME_WE)

    # Define the command to generate the binary file
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_custom_command(OUTPUT ${FILE_BASE}-make VERBATIM COMMAND
            ${PHYPP_COMPILER} debug ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME} -o ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make
            DEPENDS ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME} ${PHYPP_HEADERS})
    else()
        add_custom_command(OUTPUT ${FILE_BASE}-make VERBATIM COMMAND
            ${PHYPP_COMPILER} optimize ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME} -o ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make
            DEPENDS ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME} ${PHYPP_HEADERS})
    endif()

    # Create a target that will call this command
    add_custom_target(${FILE_BASE} ALL DEPENDS
        ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME}
        ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make
        ${PHYPP_HEADERS})

    # Specify installation of the binary
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make
        DESTINATION bin RENAME ${FILE_BASE})
endfunction()
