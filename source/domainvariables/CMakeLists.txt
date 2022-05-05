target_include_directories(hips-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(hips-run PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/dv.cc                 ${CMAKE_CURRENT_SOURCE_DIR}/dv.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_chi_hips.cc        ${CMAKE_CURRENT_SOURCE_DIR}/dv_chi_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_enth_hips.cc       ${CMAKE_CURRENT_SOURCE_DIR}/dv_enth_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_mixf_hips.cc       ${CMAKE_CURRENT_SOURCE_DIR}/dv_mixf_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_rho.cc             ${CMAKE_CURRENT_SOURCE_DIR}/dv_rho.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_dvisc.cc           ${CMAKE_CURRENT_SOURCE_DIR}/dv_dvisc.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_temp.cc            ${CMAKE_CURRENT_SOURCE_DIR}/dv_temp.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_cold_hips.cc  ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_cold_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_hips.cc       ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_hips.h
)
