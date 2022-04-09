target_include_directories(hips-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(hips-run PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase.h
        )
