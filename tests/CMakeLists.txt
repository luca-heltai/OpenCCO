SET(TESTS_SRC
  testAnatomicalTree
  testDomainCoordinates
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

if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  add_test(NAME testLiverInputBoundary
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_liver_inputs.py
      $<TARGET_FILE:inspectAnatomy>)
  set_tests_properties(testLiverInputBoundary PROPERTIES TIMEOUT 30)
  # testCompCCO requires external reference data not distributed with OpenCCO.
  foreach(test_name testAnatomicalTree testDomainCoordinates testGenerateTree testGeom testKamiya testSchreiner)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${test_name}-output")
    add_test(NAME ${test_name} COMMAND ${test_name})
    set_tests_properties(${test_name} PROPERTIES
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${test_name}-output"
      TIMEOUT 180
      FAIL_REGULAR_EXPRESSION "TEST ERRORS; KO;:KO")
  endforeach()
endif()
