cmake_minimum_required (VERSION 2.8)
SET(CMAKE_CONFIGURATION_TYPES "Release")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
SET(CMAKE_LINKER_FLAGS "-lboost_regex -lboost_program_options")
project (MFE)

add_executable(MFE main.cpp fsp.cpp)

find_package(OpenCV REQUIRED)
if(OpenCV_FOUND)
  include_directories(${OpenCV_INCLUDE_DIRS})
  target_link_libraries(MFE ${OpenCV_LIBS})
endif(OpenCV_FOUND)

find_package(Boost 1.40.0 COMPONENTS regex program_options REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR})
  target_link_libraries(MFE ${Boost_LIBRARIES})
endif(Boost_FOUND)
