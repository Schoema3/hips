#################### PREAMBLE #####################

cmake_minimum_required(VERSION 3.12)

project(HiPS VERSION 1.0.0
        LANGUAGES CXX C)

############### PROJECT-WIDE SETUP ################

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)      # Release or Debug
endif()

set(CMAKE_CXX_STANDARD 11)             # Set C++ language standard
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../run")

#------------ Executable targets --------------

add_executable(hips-run)

#------------ Target properties ---------------

target_compile_features(hips-run PUBLIC cxx_std_11)

#################### preprocessor directives

if(PARALLEL STREQUAL "yes")
    message("MPI is on")
    add_definitions(-DDOMPI)
endif()

if(NOT CHEMISTRY STREQUAL "CANTERARR" )
    add_definitions(-DPROBLEMSPECIFICRR)
endif()

if(SILENT STREQUAL "yes" )
    message("Suppressing some runtime output")
    add_definitions(-DSILENT)
endif()

add_definitions(-D_USESTDVECTOR_)

#################### Link: libraries and flags


#------- Yaml: look in source/yaml for installation, if not found, then install it there

find_library(YAML yaml-cpp          HINTS ${CMAKE_CURRENT_SOURCE_DIR}/yaml/lib)
find_path(YAML_INCLUDE_DIR yaml-cpp HINTS ${CMAKE_CURRENT_SOURCE_DIR}/yaml/include)
if(NOT YAML OR NOT YAML_INCLUDE_DIR)
    message("\n*************************************************")
    message("Yaml not found, building it in hips/source/yaml")
    message("*************************************************\n")
    find_package(Git)
    if(GIT_FOUND)
        execute_process(COMMAND rm -rf ../source/yaml/include ../source/yaml/lib ../source/yaml/yaml-cpp)
        execute_process(COMMAND git clone https://github.com/jbeder/yaml-cpp.git ../source/yaml/yaml-cpp)
        execute_process(COMMAND mkdir ../source/yaml/yaml-cpp/build)
        execute_process(COMMAND cmake .. -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_SOURCE_DIR}/yaml 
                        WORKING_DIRECTORY ../source/yaml/yaml-cpp/build)
        execute_process(COMMAND make WORKING_DIRECTORY ../source/yaml/yaml-cpp/build)
        execute_process(COMMAND make install WORKING_DIRECTORY ../source/yaml/yaml-cpp/build)
    else()
        message("\n*************************************************")
        message("Cannot get and build yaml because git not found")
        message("*************************************************\n")
    endif()
    find_library(YAML yaml-cpp HINTS ${CMAKE_CURRENT_SOURCE_DIR}/yaml/lib)
    find_path(YAML_INCLUDE_DIR yaml-cpp HINTS ${CMAKE_CURRENT_SOURCE_DIR}/yaml/include)
endif()

target_link_libraries(hips-run ${YAML})
target_include_directories(hips-run PRIVATE ${YAML_INCLUDE_DIR})



#------- Cantera needs to come after YAML in the library list to avoid duplicate symbol link errors

find_library(CANTERA cantera HINTS ${CANTERA_LIB_DIR})
if(CANTERA)
    message("**********CANTERA WAS FOUND*************")
    target_include_directories(hips-run PRIVATE ${CANTERA_INCLUDE_DIR}) 
    target_link_libraries(hips-run ${CANTERA})
    add_definitions(-DDOCANTERA)
    find_package(Boost REQUIRED COMPONENTS system)
    target_include_directories(hips-run PRIVATE ${Boost_INCLUDE_DIRS})
else()
    message("\n*************************************************")
    message("Building the code without Cantera")
    message("*************************************************\n")
endif()


#-------

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    target_link_libraries(hips-run "-framework Accelerate")
endif()

target_link_libraries(hips-run "pthread")

#################### Local source files

target_include_directories(hips-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(hips-run PRIVATE main.cc
        inputoutput.cc              inputoutput.h
        domain.cc                   domain.h
        streams.cc                  streams.h
        param.cc                    param.h
        micromixer.cc               micromixer.h
        solver.cc                   solver.h
        cvodeDriver.cc              cvodeDriver.h
        cantera_shell_functions.cc  cantera_shell_functions.h
                                    interp_linear.h
        processor.cc                processor.h
)

#################### Subdirectories

add_subdirectory(domainvariables)
add_subdirectory(cvode)
add_subdirectory(domaincases)
add_subdirectory(user_chemical_mechanisms)
#add_subdirectory(yaml)

#################### Doxygen documentation

find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_custom_target(doxygen 
        COMMAND rm -rf ../doc/doxygen/html
        COMMAND mkdir  ../doc/doxygen/html
        COMMAND cp    ../doc/doxygen_build_pages/montage_odt.png ../doc/doxygen/html
        COMMAND doxygen -s ../doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building doxygen documentation"
        VERBATIM
    )
else()
    message("\n********************************************************")
    message("Cannot build documentation because doxygen was not found")
    message("********************************************************\n")
endif()

