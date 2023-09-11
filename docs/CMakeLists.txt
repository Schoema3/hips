#################### Doxygen documentation


find_package(Doxygen REQUIRED)

if(DOXYGEN_FOUND)
    # Set the input and output directories for Doxygen
    set(DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}/src")
    set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/doxygen")

    configure_file(Doxyfile.in Doxyfile @ONLY)

    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating Doxygen documentation"
        VERBATIM
    )
else()
    message(WARNING "Doxygen not found, cannot generate documentation")
endif()

