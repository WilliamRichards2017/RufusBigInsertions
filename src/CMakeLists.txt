INCLUDE_DIRECTORIES(
	${BAMTOOLS_INCLUDE}
	${ZLIB_INCLUDE}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

add_library(rufus_big_insertions_core STATIC
	util.cpp
	contigs.cpp
	insertion.cpp
	parentGT.cpp
)

set(CORE_LIB rufus_big_insertions_core CACHE STRING "")

TARGET_LINK_LIBRARIES(rufus_big_insertions_core
	${BAMTOOLS_LIB}
	${ZLIB_LIBRARY}
)

add_dependencies(${CORE_LIB} ${RUFUS_BIG_INSERTIONS_EXTERNAL_PROJECT})
