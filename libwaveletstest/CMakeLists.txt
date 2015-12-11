set(Headers
	test.h)

set(Sources
	test.c
	compimages3.c)

include_directories ("${PROJECT_SOURCE_DIR}/libwavelets")
link_directories ("${PROJECT_BINARY_DIR}/libwavelets")

add_executable (libwaveletstest ${Headers} ${Sources})
target_link_libraries(libwaveletstest libwavelets)

add_executable (libwaveletstestadjoint testadjoint.c)
target_link_libraries(libwaveletstestadjoint libwavelets)

# Link against math lib
IF(UNIX)
  TARGET_LINK_LIBRARIES(libwaveletstest m)
  TARGET_LINK_LIBRARIES(libwaveletstestadjoint m)
ENDIF(UNIX)

