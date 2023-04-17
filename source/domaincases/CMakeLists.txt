target_include_directories(hips PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(hips PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase.cc                   ${CMAKE_CURRENT_SOURCE_DIR}/domaincase.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips.cc              ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips_comb.cc         ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips_comb.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips_simpleRxn.cc    ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_hips_simpleRxn.h
               )
