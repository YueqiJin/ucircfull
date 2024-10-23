# set the project name 
cmake_minimum_required(VERSION 3.4)

project(ucircfull VERSION 1.0.0 LANGUAGES CXX) 

# specify the C++ standard

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++2a -std=c++20" COMPILER_SUPPORTS_CXX20)
if (COMPILER_SUPPORTS_CXX20)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a -std=c++20")
else ()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++20 support. Please use a different C++ compiler.")
endif ()


if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(TARGET_DIR "debug")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -g")
else ()
    set(TARGET_DIR "release")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC")
endif ()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True) 
set(CMAKE_BUILD_DIR ${PROJECT_SOURCE_DIR}/build/${TARGET_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BUILD_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BUILD_DIR}/bin)

# generate version header.

configure_file(
	${PROJECT_SOURCE_DIR}/include/version.hpp.in
	${PROJECT_SOURCE_DIR}/include/version.hpp
	@ONLY
)

find_package(ZLIB REQUIRED)

include_directories($(ZLIB_INCLUDE_DIRS))

# import cd-hit library
add_subdirectory(include/cdhit)
add_subdirectory(include/cccs)
add_subdirectory(include/spoa)
add_subdirectory(include/minimap2)

#include dependence
include_directories (
	${PROJECT_SOURCE_DIR}/include
	${PROJECT_SOURCE_DIR}/include/spoa/include
	${PROJECT_SOURCE_DIR}/include/minimap2/
)

add_library(circfullcirccall SHARED
	src/circCall.cpp
	src/circFindccs.cpp
	src/circFindBsj.cpp
	src/circRecords.cpp
	src/circFullStruct.cpp
)

add_dependencies(circfullcirccall cccs)

add_library(circfullalign SHARED
	src/degenerateBase.cpp
	src/align.cpp
)

add_library(circfullclust SHARED
	src/clust.cpp
)

add_library(circfullUtils SHARED
	src/circfullUtils.cpp
)

add_library(circfull SHARED
	src/ucircfull.cpp
	src/extractUMI.cpp
	src/clustUMI.cpp
	src/filterCircCand.cpp
)

find_package(Threads REQUIRED)
find_package(SeqAn REQUIRED)
find_package(seqan3 3.2 REQUIRED)
find_package(Boost 1.66.0 REQUIRED)
include_directories (${SEQAN3_INCLUDE_DIRS})
add_definitions (${SEQAN3_DEFINITIONS})
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -fPIC -DWITH_ZLIB")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")

# add the executable 
add_executable(circfull_satcurve src/satcurve.cpp)
add_executable(ucircfull src/ucircfull.cpp)
#add_executable(circfull_prepGTF src/resultTransfer.cpp)

#add_dependencies(circfullalign MINIMAP2)
target_link_libraries(circfullalign
	${CMAKE_THREAD_LIBS_INIT}
	minimap2
	seqan3::seqan3
)

target_link_libraries(circfullclust
	${CMAKE_THREAD_LIBS_INIT}
	cdhit
	spoa::spoa
	ZLIB::ZLIB
	circfullUtils
)

target_link_libraries(circfullUtils
	${CMAKE_THREAD_LIBS_INIT}
	circfullalign
)

target_link_libraries(circfull_satcurve 
	${CMAKE_THREAD_LIBS_INIT}
	circfullUtils
)

target_link_libraries(circfullcirccall 
	${CMAKE_THREAD_LIBS_INIT}
	seqan3::seqan3
	${SEQAN_LIBRARIES}
	circfullUtils
	cccs
	circfullclust
)

target_link_libraries(circfull 
	${CMAKE_THREAD_LIBS_INIT}
	circfullUtils
	circfullcirccall
	circfullalign
	circfullclust
	seqan3::seqan3
)

target_link_libraries(ucircfull 
	${CMAKE_THREAD_LIBS_INIT}
	circfullUtils
	circfull
)

# install to destination
install(
	TARGETS ucircfull circfull_satcurve circfullalign circfullclust circfullUtils circfullcirccall circfull
	EXPORT ucircfullTargets
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
)