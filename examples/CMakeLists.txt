
# Example 1
add_executable(ex_1)
target_sources(ex_1 PRIVATE ex_1.cc)

# Link the hips_lib library to example1
target_link_libraries(ex_1 PRIVATE hips_lib)

# Example 2
#add_executable(example2)
#target_sources(example2 PRIVATE example2.cc)

# Link the hips_lib library to example2
#target_link_libraries(example2 PRIVATE hips_lib)

# Add more examples as needed...
