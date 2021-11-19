target_include_directories(sec.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(sec.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_enth_hips.cc              ${CMAKE_CURRENT_SOURCE_DIR}/dv_enth_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_mixf_hips.cc              ${CMAKE_CURRENT_SOURCE_DIR}/dv_mixf_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_chi_hips.cc               ${CMAKE_CURRENT_SOURCE_DIR}/dv_chi_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_hips.cc              ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_cold_hips.cc         ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_cold_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_hips.cc              ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_hips.h
                                                                 ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_hips_CQMOM.h
                                                                 ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_hips_LOGN.h
                                                                 ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_hips_MOMIC.h
                                                                 ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_hips_MONO.h
                                                                 ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_hips_SECT.h
                                                                 ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_hips_QMOM.h
        ${CMAKE_CURRENT_SOURCE_DIR}/micromixer_hips.cc           ${CMAKE_CURRENT_SOURCE_DIR}/micromixer_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/solver_hips.cc               ${CMAKE_CURRENT_SOURCE_DIR}/solver_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips.cc           ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips_comb.cc      ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips_comb.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips_simpleRxn.cc ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips_simpleRxn.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips_soot.cc      ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips_soot.h
)
