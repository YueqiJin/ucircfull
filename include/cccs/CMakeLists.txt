if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(CARGO_CMD cargo build)
	set(TARGET_DIR "debug")
else ()
	set(CARGO_CMD cargo build --release)
	set(TARGET_DIR "release")
endif ()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(RS_SO "${CMAKE_CURRENT_SOURCE_DIR}/target/${TARGET_DIR}/libcccs.a")

# build rs module
add_custom_target( CCCS_BUILD ALL
	COMMENT "Compiling rs module"
	COMMAND ${CARGO_CMD}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(NAME rs_test
COMMAND cargo test
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_library(cccs STATIC IMPORTED GLOBAL)
add_dependencies(cccs CCCS_BUILD)
set_target_properties(cccs PROPERTIES IMPORTED_LOCATION ${RS_SO})

install(
	FILES ${RS_SO}
	TYPE LIB
)