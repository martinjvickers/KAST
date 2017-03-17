# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (alfsc CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS  "-static ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -Wl,--whole-archive -lpthread -Wl,-no-whole-archive")

# Add executable and link against SeqAn dependencies.
add_executable (alfsc alfsc.cpp)
target_link_libraries (alfsc ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_executable (unit_tests unit_tests.cpp)
target_link_libraries (unit_tests ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
