set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib") # to include hips (built in this project)

add_executable(test.x)
target_sources(test.x PRIVATE test.cc)

#--------------------

set_property(TARGET test.x PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)  # adds LC_RPATH to executable

target_include_directories(test.x PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_include_directories(test.x PRIVATE "$ENV{CONDA_PREFIX}/include")
target_link_directories(   test.x PRIVATE "$ENV{CONDA_PREFIX}/lib")

#--------------------

find_package(Catch2 3 REQUIRED)
target_link_libraries(test.x PRIVATE Catch2::Catch2WithMain hips)

#--------------------

install(TARGETS test.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
