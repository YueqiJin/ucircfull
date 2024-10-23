cmake_minimum_required(VERSION 3.0)

project(minimap2)

# Define compiler flags and options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -O2 -Wc++-compat -fPIC")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D_DEBUG")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -DNDEBUG")

# Additional compiler flags
add_definitions(-DHAVE_KALLOC)

# Include directories
include_directories(${CMAKE_SOURCE_DIR})

# Source files
set(SOURCES
    main.c
    kthread.c
    kalloc.c
    misc.c
    bseq.c
    sketch.c
    sdust.c
    options.c
    index.c
    lchain.c
    align.c
    hit.c
    seed.c
    map.c
    format.c
    pe.c
    esterr.c
    splitidx.c
    ksw2_ll_sse.c
)

# SSE-specific targets on x86/x86_64
if(NOT arm_neon AND NOT sse2only)
    list(APPEND SOURCES
        ksw2_extz2_sse.c
        ksw2_extd2_sse.c
        ksw2_exts2_sse.c
        ksw2_dispatch.c
    )
    # Set compilation options for these source files
    set_source_files_properties(
        ksw2_extz2_sse.c
        ksw2_extd2_sse.c
        ksw2_exts2_sse.c
        ksw2_dispatch.c
        PROPERTIES COMPILE_FLAGS "-msse4.1"
    )
endif()

# NEON-specific targets on ARM
if(arm_neon)
    list(APPEND SOURCES
        ksw2_extz2_sse.c
        ksw2_extd2_sse.c
        ksw2_exts2_sse.c
    )
    include_directories(sse2neon)
    if(aarch64)
        add_definitions(-D_FILE_OFFSET_BITS=64 -fsigned-char)
    else()
        add_definitions(-D_FILE_OFFSET_BITS=64 -mfpu=neon -fsigned-char)
    endif()
endif()


# Create the library target
add_library(minimap2 STATIC ${SOURCES})

# Find and link ZLIB library
find_package(ZLIB REQUIRED)
# Link 'pthread' automatically
find_package(Threads REQUIRED)
target_link_libraries(minimap2
	ZLIB::ZLIB
	Threads::Threads
	m
)