SET(BIN_SRC
  growAnatomicalTree
  generateTree2D
  generateTree3D
)

FOREACH(FILE ${BIN_SRC})
  add_executable(${FILE} ${FILE}.cpp ${SOURCE_BASE_FILES} ${SOURCE_BASE_HEADER_FILES} )
  target_link_libraries (${FILE} ${CERES_LIBRARIES} ${DGTAL_LIBRARIES} )
ENDFOREACH(FILE)
if(OPENCCO_WITH_VTK)
  target_compile_definitions(growAnatomicalTree PRIVATE OPENCCO_WITH_VTK=1)
  target_link_libraries(growAnatomicalTree ${VTK_LIBRARIES})
endif()
