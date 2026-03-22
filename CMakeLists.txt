cmake_minimum_required(VERSION 3.16)

project(ucircfull VERSION 1.1.0 LANGUAGES CXX)

include(CheckCXXCompilerFlag)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(NOT MSVC)
	set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

check_cxx_compiler_flag("-std=c++23" COMPILER_SUPPORTS_CXX23)
if(NOT COMPILER_SUPPORTS_CXX23)
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support C++23.")
endif()

set(UCIRCFULL_GENERATED_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated/include)
file(MAKE_DIRECTORY ${UCIRCFULL_GENERATED_INCLUDE_DIR})

configure_file(
	${PROJECT_SOURCE_DIR}/include/version.hpp.in
	${UCIRCFULL_GENERATED_INCLUDE_DIR}/version.hpp
	@ONLY
)

if(UNIX AND NOT APPLE)
	set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
	set(CMAKE_BUILD_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
	set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

option(UCIRCFULL_BUNDLE_GNU_CXX_RUNTIME "Install GNU C++ runtime libraries (libstdc++, libgcc_s) for portable Linux deployment" ON)

find_package(ZLIB REQUIRED)

include_directories(${ZLIB_INCLUDE_DIRS})

# import cd-hit library
add_subdirectory(include/cdhit)
add_subdirectory(include/cccs)
add_subdirectory(include/spoa)
add_subdirectory(include/minimap2)

#include dependence
include_directories (
	${PROJECT_SOURCE_DIR}/include
	${UCIRCFULL_GENERATED_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}/include/spoa/include
	${PROJECT_SOURCE_DIR}/include/minimap2/
)

add_library(circfullcirccall SHARED
	src/circCall.cpp
	src/circFindccs.cpp
	src/circFindBsj.cpp
	src/circFusion.cpp
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
	PRIVATE
		Threads::Threads
		minimap2_dev
		seqan3::seqan3
)

target_link_libraries(circfullclust
	PRIVATE
		Threads::Threads
		cdhit
		spoa::spoa
		ZLIB::ZLIB
		circfullUtils
)

target_link_libraries(circfullUtils
	PRIVATE
		Threads::Threads
		circfullalign
)

target_link_libraries(circfull_satcurve 
	PRIVATE
		Threads::Threads
		circfullUtils
)

target_link_libraries(circfullcirccall 
	PRIVATE
		Threads::Threads
		seqan3::seqan3
		${SEQAN_LIBRARIES}
		circfullUtils
		cccs
		circfullclust
)

target_link_libraries(circfull 
	PRIVATE
		Threads::Threads
		circfullUtils
		circfullcirccall
		circfullalign
		circfullclust
		seqan3::seqan3
)

target_link_libraries(ucircfull 
	PRIVATE
		Threads::Threads
		circfullUtils
		circfull
		circfullcirccall
)

set(UCIRCFULL_EXPORT_NAME ucircfullTargets)
set(UCIRCFULL_CMAKE_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/ucircfull)

install(
	TARGETS
		ucircfull
		circfull_satcurve
		circfullalign
		circfullclust
		circfullUtils
		circfullcirccall
		circfull
		cdhit
		spoa
		minimap2
	EXPORT ${UCIRCFULL_EXPORT_NAME}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ucircfull
)

install(
	FILES
		${PROJECT_SOURCE_DIR}/include/align.hpp
		${PROJECT_SOURCE_DIR}/include/circCall.hpp
		${PROJECT_SOURCE_DIR}/include/circfullUtils.hpp
		${PROJECT_SOURCE_DIR}/include/clust.hpp
		${PROJECT_SOURCE_DIR}/include/degenerateBase.hpp
		${PROJECT_SOURCE_DIR}/include/ucircfull.hpp
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ucircfull
)

install(
	FILES ${PROJECT_SOURCE_DIR}/include/argparse/argparse.hpp
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ucircfull/argparse
)

install(
	FILES ${UCIRCFULL_GENERATED_INCLUDE_DIR}/version.hpp
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ucircfull
)

write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/ucircfullConfigVersion.cmake
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/ucircfullConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/ucircfullConfig.cmake
	INSTALL_DESTINATION ${UCIRCFULL_CMAKE_INSTALL_DIR}
)

install(
	EXPORT ${UCIRCFULL_EXPORT_NAME}
	NAMESPACE ucircfull::
	DESTINATION ${UCIRCFULL_CMAKE_INSTALL_DIR}
)

install(
	FILES
		${CMAKE_CURRENT_BINARY_DIR}/ucircfullConfig.cmake
		${CMAKE_CURRENT_BINARY_DIR}/ucircfullConfigVersion.cmake
	DESTINATION ${UCIRCFULL_CMAKE_INSTALL_DIR}
)

if(UNIX AND NOT APPLE AND UCIRCFULL_BUNDLE_GNU_CXX_RUNTIME)
	execute_process(
		COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.so.6
		OUTPUT_VARIABLE UCIRCFULL_LIBSTDCXX
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	execute_process(
		COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libgcc_s.so.1
		OUTPUT_VARIABLE UCIRCFULL_LIBGCC
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)

	if(EXISTS "${UCIRCFULL_LIBSTDCXX}")
		file(REAL_PATH "${UCIRCFULL_LIBSTDCXX}" UCIRCFULL_LIBSTDCXX_REAL)
		install(
			FILES "${UCIRCFULL_LIBSTDCXX_REAL}"
			DESTINATION ${CMAKE_INSTALL_LIBDIR}
			RENAME libstdc++.so.6
		)
	else()
		message(WARNING "UCIRCFULL_BUNDLE_GNU_CXX_RUNTIME is ON, but libstdc++.so.6 was not found via compiler path lookup.")
	endif()

	if(EXISTS "${UCIRCFULL_LIBGCC}")
		file(REAL_PATH "${UCIRCFULL_LIBGCC}" UCIRCFULL_LIBGCC_REAL)
		install(
			FILES "${UCIRCFULL_LIBGCC_REAL}"
			DESTINATION ${CMAKE_INSTALL_LIBDIR}
			RENAME libgcc_s.so.1
		)
	else()
		message(WARNING "UCIRCFULL_BUNDLE_GNU_CXX_RUNTIME is ON, but libgcc_s.so.1 was not found via compiler path lookup.")
	endif()
endif()

set(CPACK_PACKAGE_NAME "ucircfull")
set(CPACK_PACKAGE_VENDOR "yangence")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_SOURCE_IGNORE_FILES "/build/;/.git/;\\\\.vscode/;\\\\.idea/")
if(WIN32)
	set(CPACK_GENERATOR "ZIP")
else()
	set(CPACK_GENERATOR "TGZ")
endif()
include(CPack)