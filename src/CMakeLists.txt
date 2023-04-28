
add_executable(hips.x)
target_sources(hips.x PRIVATE driver.cc cvodeDriver.cc cvodeDriver.h hips.cc hips.h)

if(APPLE)
    find_library(ACCELERATE Accelerate)
endif()
find_package(Threads REQUIRED)
target_include_directories(hips.x PRIVATE "/usr/local/include" "/usr/include")
target_link_directories(hips.x PRIVATE "/usr/local/lib" "/usr/lib/x86_64-linux-gnu" "/usr/lib/python3.10/config-3.10-x86_64-linux-gnu")
target_link_libraries(hips.x cantera sundials_cvodes python3.10 Threads::Threads)

install(TARGETS hips.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
