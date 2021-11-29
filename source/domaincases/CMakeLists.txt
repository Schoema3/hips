target_include_directories(sec.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(sec.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase.cc                      ${CMAKE_CURRENT_SOURCE_DIR}/domaincase.h
        )
