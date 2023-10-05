
# Example 1
add_executable(ex_1)
target_sources(ex_1 PRIVATE ex_1.cc)

# Link the hips_lib library to example1
target_link_libraries(ex_1 PRIVATE hips_lib)

# Example 2
add_executable(ex_2)
target_sources(ex_2 PRIVATE ex_2.cc)

# Link the hips_lib library to example2
target_link_libraries(ex_2 PRIVATE hips_lib)

# Example 3
add_executable(ex_3)
target_sources(ex_3 PRIVATE ex_3.cc)

# Link the hips_lib library to example2
target_link_libraries(ex_3 PRIVATE hips_lib)

# Add more examples as needed...
