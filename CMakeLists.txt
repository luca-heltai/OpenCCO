cmake_minimum_required(VERSION 3.16...3.31)
project(CCO_Implementation)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(CTest)


SET(SOURCE_BASE_HEADER_FILES
  ${PROJECT_SOURCE_DIR}/src/helpers/GeomHelpers.h
  ${PROJECT_SOURCE_DIR}/src/helpers/ExpandTreeHelpers.h
  ${PROJECT_SOURCE_DIR}/src/CoronaryArteryTree.h
  ${PROJECT_SOURCE_DIR}/src/helpers/XmlHelpers.h
  ${PROJECT_SOURCE_DIR}/src/DomainController.h
  ${PROJECT_SOURCE_DIR}/tests/TreeTest.h        
)
SET(SOURCE_BASE_FILES
  ${PROJECT_SOURCE_DIR}/src/CoronaryArteryTree.ih
)

find_package(Ceres REQUIRED)



# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/ext/" "${PROJECT_SOURCE_DIR}/src/helpers/" "${PROJECT_SOURCE_DIR}/src/" )
FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS} )

option(OPENCCO_WITH_VTK "Enable original-STL wall certificates with VTK" OFF)
if(OPENCCO_WITH_VTK)
  find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel FiltersCore IOGeometry)
endif()

add_subdirectory(tests)
add_subdirectory(bin)
add_subdirectory(tools)

