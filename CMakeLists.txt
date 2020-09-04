IF(WIN32)
cmake_minimum_required(VERSION 2.6)
ELSE(WIN32)
cmake_minimum_required(VERSION 2.4)
ENDIF(WIN32)

IF(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

# Set name of project
SET(PROJECT_NAME "DetectLineModule")

# Start a project
PROJECT(${PROJECT_NAME})

##########################################################################################
# Section to include your own FindX cmake files and options
##########################################################################################

## Commented example for the EIGEN library. Replace this with what your module needs, 
## for example Boost libraries, OpenCV, inference libraries, etc.
#
# OPTION(BUILD_EIGEN OFF "Build module with Eigen support")
# CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_SOURCE_DIR}/inc/config.h)
# IF (BUILD_EIGEN) 
# 	MESSAGE("[*] Eigen support, this will enable all (sparse) matrix related code")
# 	FIND_PACKAGE(Eigen3 REQUIRED)
# 	INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
# ENDIF()

##########################################################################################

# First set empty
SET(PARENT_INCLUDE_DIRS "")

# Get include directories
FOREACH(INCLDIR ${PARENT_INCLUDE_DIRS})
  MESSAGE("[*] Add ${INCLDIR}")
  INCLUDE_DIRECTORIES(${INCLDIR})
ENDFOREACH()

# Manually add some CXX flag (disabled by default because most cross-compilers choke on it)
SET(CMAKE_CXX_FLAGS -std=c++11)

# Some debug information
MESSAGE("[*] \"${PROJECT_NAME}\" uses the following CXX flags: \"${CMAKE_CXX_FLAGS}\" (can be empty)")
MESSAGE("[*] Linked libraries: ${LIBS}")

# Source directory
SET(SOURCE_DIR "src")
SET(INCLUDE_DIR "include")

# Search for functional source code
FILE(GLOB FOLDER_SOURCE ${SOURCE_DIR}/*.cpp ${SOURCE_DIR}/*.cc ${SOURCE_DIR}/*.c)
FILE(GLOB FOLDER_HEADER ${INCLUDE_DIR}/*.h ${INCLUDE_DIR}/*.hh ${INCLUDE_DIR}/*.hpp)

# For Visual Studio group the files
SOURCE_GROUP("Source files" FILES ${FOLDER_SOURCE})
SOURCE_GROUP("Header files" FILES ${FOLDER_HEADER})

# Automatically add include directories if needed.
FOREACH(HEADER_FILE ${FOLDER_HEADER})
  GET_FILENAME_COMPONENT(HEADER_DIR ${HEADER_FILE} PATH)
  INCLUDE_DIRECTORIES(${HEADER_DIR})
ENDFOREACH(HEADER_FILE ${FOLDER_HEADER})

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "INCLUDE_DIRECTORY='${dir}'")
endforeach()

#message(STATUS "Enable testing framework")
#enable_testing()
#ADD_SUBDIRECTORY(test)

SET(COMMON_INCLUDES ${FOLDER_SOURCE} ${FOLDER_HEADER})
MESSAGE(STATUS "Files included: ${COMMON_INCLUDES}")

set(MAINFILENAME src/${PROJECT_NAME}Main.cpp)
FILE(GLOB MAINFILE "${MAINFILENAME}")

# For testing we have to include everything too, except for the main file
SET(TEST_INCLUDES ${FOLDER_SOURCE} ${FOLDER_HEADER} ${AIM_FILES})
list(REMOVE_ITEM TEST_INCLUDES ${MAINFILE})
SET(PROJECT_TESTLIB "${PROJECT_NAME}Test")

# Set up our main executable.
IF(FOLDER_SOURCE STREQUAL "")
  MESSAGE(FATAL_ERROR "No source code files found. Please add something")
ELSE ()
  ADD_EXECUTABLE(${PROJECT_NAME} ${COMMON_INCLUDES})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${LIBS})
  INSTALL(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
  
  # Create the 
  ADD_LIBRARY(${PROJECT_TESTLIB} ${TEST_INCLUDES})
ENDIF ()

