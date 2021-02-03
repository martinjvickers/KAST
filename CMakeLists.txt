# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (kast CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)

SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})
add_definitions ( -DSEQAN_DISABLE_VERSION_CHECK )

# Add CXX flags found by find_package (SeqAn).
#set (CMAKE_CXX_FLAGS  "-static ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -I/usr/local/include/jellyfish-2.3.0 -Wl,--whole-archive -ljellyfish-2.0 -lpthread -Wl,-no-whole-archive")
set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (kast distance.cpp kast.cpp)
target_link_libraries (kast ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_executable (kast_tests distance.cpp kast_tests.cpp)
target_link_libraries (kast_tests ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
