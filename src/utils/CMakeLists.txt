add_library (utils INTERFACE)

target_include_directories (utils
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
) 

target_sources(utils INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}/UtilsOxygenationMap.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UtilsOxygenationMap.hpp
)

target_link_libraries(utils
    INTERFACE Eigen3::Eigen
)
