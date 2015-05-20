start_demo_dependencies()
find_package(Doxygen REQUIRED)
end_demo_dependencies()

configure_file(${cpp_tech_demo_SOURCE_DIR}/docs/Doxyfile.in ${cpp_tech_demo_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
  ${DOXYGEN_EXECUTABLE} ${cpp_tech_demo_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${cpp_tech_demo_BINARY_DIR}
  COMMENT "Generating API documentation with Doxygen" VERBATIM)
set(GENERATE_DOCS ON)
