set(CPP_FILES 
	hello
	world)

#Generate projects for each file
foreach(CPP_FILE ${CPP_FILES})

	# set up MATLAB libraries
	set(CPP_FILE_NAME ${CPP_FILE}.cpp)
	INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIR})
	add_library(${CPP_FILE} SHARED ${CPP_FILE_NAME} ${CMAKE_SOURCE_DIR}/Matlabdef.def)
	target_link_libraries(${CPP_FILE} ${MATLAB_LIBRARIES})

	## Change naming for the correct MATLAB wrapper
	if(WIN32)
		if (CMAKE_CL_64)
			#Windows64
			SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexw64 LINK_FLAGS /export:mexFunction)
		else()
			#Windows32
			SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexw32 LINK_FLAGS /export:mexFunction)
		endif()
	else(WIN32)
		if(APPLE)
			if (CMAKE_SIZEOF_VOID_P MATCHES "8")
				#Linux64
				SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexmaci64)
			else()
				#Linux32
				SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexmaci32)
			endif (CMAKE_SIZEOF_VOID_P MATCHES "8")
		else(APPLE)
			if (CMAKE_SIZEOF_VOID_P MATCHES "8")
				#Linux64
				SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexa64)
			else()
				#Linux32
				SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexglx)
			endif (CMAKE_SIZEOF_VOID_P MATCHES "8")
		endif(APPLE)
	endif(WIN32)

	# install to /bin by default
	install(TARGETS ${CPP_FILE} DESTINATION ../bin)
endforeach()