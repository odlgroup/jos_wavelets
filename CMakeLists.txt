cmake_minimum_required (VERSION 2.6) 
project (libwavelets)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake) # add FindMatlab module


# Choosing order of wavelet coefficients

#Default order:     LLL_N, LLH_N ...   .... HHL_1, HHH_1 

#for order:   HHH_1,HHL_1 ..                ..,LLH_N,LLL_N   
#uncomment following line
# 	add_definitions(-DCOARSEEND)

# Enable or disable high precision
option(HIGH_PRECISION "Use double precision" TRUE)

if(HIGH_PRECISION)
	add_definitions(-DHIGH_PRECISION)
endif(HIGH_PRECISION)

# Enable or disable high precision
option(MEX_ENABLED "Use double precision" FALSE)

if(MEX_ENABLED)
	find_package(Matlab REQUIRED)
endif()

# Option to enable or disable debugging and optimization
option(DEBUGGING "Enable debugging (in make)" FALSE)
option(OPTIMIZATION "Enable optimization (in make)" TRUE)

if(OPTIMIZATION)
	if(DEBUGGING)
		SET(CMAKE_BUILD_TYPE RelWithDebInfo)
	else(DEBUGGING)
		SET(CMAKE_BUILD_TYPE Release)
	endif(DEBUGGING)
else(OPTIMIZATION)
	if(DEBUGGING)
		SET(CMAKE_BUILD_TYPE Debug)
	else(DEBUGGING)
		SET(CMAKE_BUILD_TYPE None)
	endif(DEBUGGING)
endif(OPTIMIZATION)

# Add sub directories as needed
add_subdirectory (libwavelets)
add_subdirectory (libwaveletstest)
if(MEX_ENABLED)
	add_subdirectory (libwaveletsmex)
endif(MEX_ENABLED)
