target_include_directories(sec.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(sec.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/dv.cc                 ${CMAKE_CURRENT_SOURCE_DIR}/dv.h 
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_rho.cc             ${CMAKE_CURRENT_SOURCE_DIR}/dv_rho.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_dvisc.cc           ${CMAKE_CURRENT_SOURCE_DIR}/dv_dvisc.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_temp.cc            ${CMAKE_CURRENT_SOURCE_DIR}/dv_temp.h
       )
