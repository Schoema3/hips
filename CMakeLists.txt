
# ---- PREAMBLE ----
cmake_minimum_required(VERSION 3.15)
project(hips VERSION 0.1.0 LANGUAGES CXX)

# ---- PROJECT SETUP ----
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)       # Release or Debug
endif()

# ---- Set C++ language standard ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

# ---- Installation location ----
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/installation")

# ---- Build options ----
option(HIPS_BUILD_EXAMPLES   "Build hips examples"        ON)
option(HIPS_BUILD_DOCS       "Build hips documentation"   OFF)

# ---- MAIN PROJECT TARGETS ----
add_subdirectory(src)

# ---- OPTIONAL PROJECT TARGETS ----
# -------- EXAMPLES --------
if (HIPS_BUILD_EXAMPLES)
    set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/examples) # Optional, to set output path for examples
    add_subdirectory(examples)
endif()

# -------- DOCUMENTATION ----------

if (HIPS_BUILD_DOCS )
    add_subdirectory(docs)
endif()
