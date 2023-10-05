
# ---- Create the hips library target ----
add_library(hips_lib STATIC
    cvodeDriver.cc
    hips.cc
    hips.h
    cvodeDriver.h
)

# ---- Set C++ standard for the library ----
target_compile_features(hips_lib PUBLIC cxx_std_17)

# ---- Include directories for the library ----
target_include_directories(hips_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

find_package(Threads REQUIRED)
target_link_libraries(hips_lib PRIVATE Threads::Threads)


# ---- the option and variable for custom integrator ----
option(USE_CUSTOM_INTEGRATOR "Use a custom integrator" OFF)
set(CUSTOM_INTEGRATOR_PATH "" CACHE STRING "Path to custom integrator library")

# ---- the source files for the library based on the option ----
if(USE_CUSTOM_INTEGRATOR)
    # Include the headers and source files of the custom integrator
    target_sources(hips_lib PRIVATE cvodeDriver.cc hips.cc hips.h cvodeDriver.h ${CUSTOM_INTEGRATOR_PATH}/customIntegrator.cc)
    target_include_directories(hips_lib PRIVATE ${CUSTOM_INTEGRATOR_PATH})
else()
    #  the headers and source files of cvodeDriver (default)
    target_sources(hips_lib PRIVATE cvodeDriver.cc hips.cc hips.h cvodeDriver.h)
endif()

# ---- Link the necessary libraries to hips_lib (update as needed) ----
target_link_libraries(hips_lib PRIVATE cantera sundials_cvodes python3.10 Threads::Threads)

# ---- Optionally install the library ----
install(TARGETS hips_lib
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

# ---- Optionally install the header files ----
install(FILES hips.h cvodeDriver.h DESTINATION include)
