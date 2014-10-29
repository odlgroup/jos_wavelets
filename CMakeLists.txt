cmake_minimum_required (VERSION 2.6) 
project (libwavelets)

# Enable or disable high precision
add_definitions(-DHIGH_PRECISION)

# Optimized 
# set(CMAKE_BUILD_TYPE Release)
# Debug
# set(CMAKE_BUILD_TYPE Debug)

# Add sub directories as needed
add_subdirectory (libwavelets)
add_subdirectory (test)