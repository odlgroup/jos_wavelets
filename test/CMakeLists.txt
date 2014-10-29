set(Headers
	test.h)

set(Sources
	test.c
	compimages3.c)

include_directories ("${PROJECT_SOURCE_DIR}/libwavelets")
link_directories ("${PROJECT_BINARY_DIR}/libwavelets")

add_executable (test ${Headers} ${Sources})

target_link_libraries(test libwavelets)