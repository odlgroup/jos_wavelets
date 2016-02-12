set(Headers
    test.h)

set(Sources
    test.c
    compimages3.c)

include_directories ("${PROJECT_SOURCE_DIR}/lib")
link_directories ("${PROJECT_BINARY_DIR}/lib")

add_executable (jos_wavelets_test ${Headers} ${Sources})
target_link_libraries(jos_wavelets_test jos_wavelets)

add_executable (jos_wavelets_testadjoint testadjoint.c)
target_link_libraries(jos_wavelets_testadjoint jos_wavelets)

# Link against math lib
IF(UNIX)
  TARGET_LINK_LIBRARIES(jos_wavelets_test m)
  TARGET_LINK_LIBRARIES(jos_wavelets_testadjoint m)
ENDIF(UNIX)

