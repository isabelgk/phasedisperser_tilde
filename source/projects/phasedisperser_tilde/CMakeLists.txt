cmake_minimum_required(VERSION 3.0)

set(C74_MIN_API_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../min-api)
include(${C74_MIN_API_DIR}/script/min-pretarget.cmake)


include_directories( 
	"${C74_INCLUDES}"
    "${CMAKE_CURRENT_LIST_DIR}"
)

add_library( 
	${PROJECT_NAME} 
	MODULE
        phasedisperser_tilde.cpp
)

include(${C74_MIN_API_DIR}/script/min-posttarget.cmake)

