
ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)

INCLUDE_DIRECTORIES(
	${ZLIB_INCLUDE}
	${BAMTOOLS_INCLUDE}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

add_executable(detect
	detect.cpp
)

target_link_libraries(detect
	${CORE_LIB}
)
