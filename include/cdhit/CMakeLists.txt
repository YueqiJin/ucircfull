message("compiling cd-hit library --- ${PROJECT_SOURCE_DIR}")

# Compiler flags
set(CMAKE_CXX_FLAGS "-Wall -ggdb -pg -fopenmp")

# Support debugging
option(DEBUG "Enable debugging" OFF)
if(DEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_ZLIB")
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})
set(LINK_LIBRARIES ${LINK_LIBRARIES} ${ZLIB_LIBRARIES})

add_library(cdhit SHARED 
	cdhit-common.c++
	cdhit-utility.c++
	cdhit.c++
)