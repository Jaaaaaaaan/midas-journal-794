add_library (classes INTERFACE)

target_include_directories (classes
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
) 

target_sources(classes INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}/NodeTable.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NodeTable.hpp
    ${CMAKE_CURRENT_LIST_DIR}/OBB.cpp
    ${CMAKE_CURRENT_LIST_DIR}/OBB.hpp
    ${CMAKE_CURRENT_LIST_DIR}/OxygenationMap.cpp
    ${CMAKE_CURRENT_LIST_DIR}/OxygenationMap.hpp
    ${CMAKE_CURRENT_LIST_DIR}/SupplyMap.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SupplyMap.hpp
    ${CMAKE_CURRENT_LIST_DIR}/VascularTree.cpp
    ${CMAKE_CURRENT_LIST_DIR}/VascularTree.hpp
)

target_link_libraries(classes
    INTERFACE Eigen3::Eigen
)
