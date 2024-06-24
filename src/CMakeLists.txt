##################### hips library

add_library(hips SHARED)
add_library(hips::hips ALIAS hips)

target_sources(hips PRIVATE 
    hips.cc                 
    PUBLIC
    FILE_SET hips_headers
        TYPE HEADERS
        BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" 
        FILES hips.h
)

if(REACTIONS_ENABLED)
    target_sources(hips PRIVATE
        batchReactor_cvode.cc
        batchReactor_cantera.cc
    )
endif()

##################### build library

find_library(YAML   yaml-cpp)

target_include_directories(hips PUBLIC "$ENV{CONDA_PREFIX}/include")
target_include_directories(hips PUBLIC "/usr/local/include")

target_link_directories(hips PUBLIC "$ENV{CONDA_PREFIX}/lib")
target_link_directories(hips PUBLIC "/usr/local/lib")

target_link_libraries(hips ${YAML} )

if(REACTIONS_ENABLED)
    find_library(CVODES sundials_cvodes)
    find_library(CANTERA NAMES cantera_shared cantera)
    target_link_libraries(hips ${CANTERA} ${CVODES})
endif()

############# Installation ################

include(GNUInstallDirs)

install(
    TARGETS hips                                        # install library
    EXPORT  hips_exp                                    # define export var for next comand
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    FILE_SET hips_headers DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

install(
    EXPORT      hips_exp                                # hips export file
    FILE        hips_exp.cmake  
    NAMESPACE   hips::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hips"
)

install(
    FILES       "hipsConfig.cmake"                      # install config files
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hips"    # edit file hipsConfig.cmake
)
