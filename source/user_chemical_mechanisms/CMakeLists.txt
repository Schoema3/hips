target_include_directories(sec.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(    CHEMISTRY STREQUAL "ONESTEPC2H4")
    target_sources(sec.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/onestep_c2h4.cc)
elseif(CHEMISTRY STREQUAL "ONESTEPCH4")
    target_sources(sec.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/onestep_ch4.cc)
endif()

