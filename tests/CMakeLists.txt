SET(TESTS_SRC
  testAnatomicalTree
  testAnatomicalGrowth
  testFrustumGeometry
  testDomainCoordinates
  testSamplingBudget
  testProjection
  testImplicitDomain
  testCompCCO
  testGenerateTree
  testGeom
  testKamiya
  testSchreiner
)


FOREACH(FILE ${TESTS_SRC})
  add_executable(${FILE} ${FILE}.cpp ${SOURCE_BASE_FILES} ${SOURCE_BASE_HEADER_FILES} )
  target_link_libraries (${FILE} ${CERES_LIBRARIES} ${DGTAL_LIBRARIES} )
ENDFOREACH(FILE)
target_compile_definitions(testGenerateTree PUBLIC SAMPLE_DIR="${PROJECT_SOURCE_DIR}/Samples/")
target_compile_definitions(testGeom PUBLIC SAMPLE_DIR="${PROJECT_SOURCE_DIR}/Samples/")

if(OPENCCO_WITH_VTK)
  add_executable(testSurfaceWall testSurfaceWall.cpp)
  target_compile_definitions(testSurfaceWall PRIVATE OPENCCO_WITH_VTK=1
    OPENCCO_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(testSurfaceWall ${DGTAL_LIBRARIES} ${VTK_LIBRARIES})
endif()

if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  add_test(NAME testLiverInputBoundary
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_liver_inputs.py
      $<TARGET_FILE:inspectAnatomy>)
  set_tests_properties(testLiverInputBoundary PROPERTIES TIMEOUT 30)
  add_test(NAME testGenerationCLI
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_generation_cli.py
      $<TARGET_FILE:generateTree2D> $<TARGET_FILE:generateTree3D>)
  set_tests_properties(testGenerationCLI PROPERTIES TIMEOUT 120)
  if(OPENCCO_WITH_VTK)
    set(surface_test_flag --with-surface)
    add_test(NAME testSurfaceWall COMMAND testSurfaceWall)
  endif()
  add_test(NAME testAnatomicalGrowthCLI
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_anatomical_growth_cli.py
      $<TARGET_FILE:growAnatomicalTree> ${surface_test_flag})
  set_tests_properties(testAnatomicalGrowthCLI PROPERTIES TIMEOUT 120)
  # testCompCCO requires external reference data not distributed with OpenCCO.
  foreach(test_name testAnatomicalTree testAnatomicalGrowth testFrustumGeometry testDomainCoordinates testSamplingBudget testProjection testImplicitDomain testGenerateTree testGeom testKamiya testSchreiner)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${test_name}-output")
    add_test(NAME ${test_name} COMMAND ${test_name})
    set_tests_properties(${test_name} PROPERTIES
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${test_name}-output"
      TIMEOUT 180
      FAIL_REGULAR_EXPRESSION "TEST ERRORS; KO;:KO")
  endforeach()
endif()
