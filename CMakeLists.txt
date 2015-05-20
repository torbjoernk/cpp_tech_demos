cmake_minimum_required(VERSION 3.1)
project(cpp_tech_demo
  VERSION 0.0.1
  LANGUAGES C CXX
)

include(FeatureSummary)
include(ExternalProject)

option(USE_GIT_PROTOCOL "use 'git@host:repo' for 3rd-party libraries instead of 'https://host/repo'" ON)

option(DEMO_STD      "enable C++(03/11/14/17) STD demo" ON )
add_feature_info(DEMO_STD
  DEMO_STD
  "Demos utilizing exclusively the C++ standard libraries."
)
option(DEMO_CONCEPTS "enable C++ concepts demo"         OFF)
add_feature_info(DEMO_CONCEPTS
  DEMO_CONCEPTS
  "Demos for concepts and traits introspection provided by external libraries."
)
option(DEMO_QT5      "enable Qt5 demo"                  OFF)
add_feature_info(DEMO_QT5
  DEMO_QT5
  "Demos for Qt5 based libraries and applications."
)

message(STATUS "-------------------------------------------------------------------------------")
set(3RDPATY_ROOT_PATH ${CMAKE_SOURCE_DIR}/3rdparty)
set_directory_properties(PROPERTIES EP_BASE ${3RDPATY_ROOT_PATH})
add_subdirectory(3rdparty)
message(STATUS "")

if(DEMO_STD)
  message(STATUS "-------------------------------------------------------------------------------")
  add_subdirectory(std)
  message(STATUS "")
endif()

if(DEMO_CONCEPTS)
  message(STATUS "-------------------------------------------------------------------------------")
  add_subdirectory(concepts)
  message(STATUS "")
endif()

if(DEMO_QT5)
  message(STATUS "-------------------------------------------------------------------------------")
  add_subdirectory(qt5)
  message(STATUS "")
endif()

message(STATUS "-------------------------------------------------------------------------------")
feature_summary(WHAT ALL)
