##################### hips library #####################

add_library(hips)

target_sources(hips PRIVATE 
    hips.cc                 
    hips.h 
)

if(REACTIONS_ENABLED)
    target_sources(hips PRIVATE
        batchreactor_cvode.cc   
        batchreactor_cvode.h   
        integrator_cvode.h
        batchreactor_cantera.cc 
        batchreactor_cantera.h
        batchreactor.h
    )
endif()


# Option to enable reactions


# Setting target properties
set_target_properties(hips PROPERTIES 
    PUBLIC_HEADER "hips.h"
)

# Find necessary libraries
find_package(Threads REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(LAPACK REQUIRED)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(hips ${Python_LIBRARIES} ${Python} ${ACCELERATE})


# Include boost if necessary
if(reactions_enabled)
    find_library(CVODES sundials_cvodes)
    find_library(CANTERA NAMES cantera_shared cantera)
    get_filename_component(cantera_path ${CANTERA} DIRECTORY)
    set_property(TARGET hips APPEND PROPERTY INSTALL_RPATH ${cantera_path})
    target_link_libraries(hips ${CANTERA})
    target_link_libraries(hips ${CANTERA} ${YAML} ${CVODES} ${Python_LIBRARIES} ${Python} ${ACCELERATE})

endif()

# Link libraries
include_directories("/usr/local/include")
target_link_directories(hips PRIVATE "/usr/local/lib")


############# Installation ################

include(GNUInstallDirs)

install(TARGETS hips
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" 
)

