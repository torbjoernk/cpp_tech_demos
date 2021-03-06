cmake_minimum_required(VERSION 3.1)
project(cpp_tech_demo
  VERSION 0.0.1
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH ${cpp_tech_demo_SOURCE_DIR}/cmake)
include(FeatureSummary)
include(ExternalProject)
include(utility_functions)

option(USE_GIT_PROTOCOL "use 'git@host:repo' for 3rd-party libraries instead of 'https://host/repo'" ON)

option(DEMO_STD      "enable C++(03/11/14/17) STD demo"   ON )
add_feature_info(DEMO_STD
  DEMO_STD
  "Demos utilizing exclusively the C++ standard libraries."
)
option(DEMO_CONCEPTS "enable C++ concepts demo"           OFF)
add_feature_info(DEMO_CONCEPTS
  DEMO_CONCEPTS
  "Demos for concepts and traits introspection provided by external libraries."
)
option(DEMO_MPI "enable MPI demos"                        OFF)
add_feature_info(DEMO_MPI
  DEMO_MPI
  "Demos for MPI communication with C++ either through native C interface or 3rd-party C++ bindings."
)
option(DEMO_QT5      "enable Qt5 demo"                    OFF)
add_feature_info(DEMO_QT5
  DEMO_QT5
  "Demos for Qt5 based libraries and applications."
)

option(DOCS          "enable generation of documentation" OFF)
add_feature_info(DOCS
  DOCS
  "enable generation of documentation."
)

message(STATUS "-------------------------------------------------------------------------------")
set(3RDPARTY_ROOT_PATH ${CMAKE_SOURCE_DIR}/3rdparty)
set_directory_properties(PROPERTIES EP_BASE ${3RDPARTY_ROOT_PATH})
add_subdirectory(3rdparty)

add_demo("Standard C++(03/11/14/17)" "std" DEMO_STD)
add_demo("C++ Concepts and Traits Introspection" "concepts" DEMO_CONCEPTS)
add_demo("MPI Examples of multiple scale" "mpi" DEMO_MPI)
add_demo("Qt5" "qt5" DEMO_QT5)

add_demo("Documentation" "docs" DOCS)

message(STATUS "-------------------------------------------------------------------------------")
feature_summary(WHAT ALL)
