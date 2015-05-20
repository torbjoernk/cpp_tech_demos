message(STATUS "Demo: C++ Concepts and Traits Introspection")

# Tick library
set(Tick_GITHUB_REPO "pfultz2/Tick")
message(STATUS "  retreiving Github repository: ${Tick_GITHUB_REPO}")

ExternalProject_Add(
  Tick
  LIST_SEPARATOR " "
  GIT_REPOSITORY "${GITHUB_PREFIX}${Tick_GITHUB_REPO}"
  TIMEOUT ${3RDPARTY_DOWNLOAD_TIMEOUT}
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  BUILD_IN_SOURCE ON
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  TEST_COMMAND ""
  INSTALL_DIR ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
  LOG_CONFIGURE OFF
  LOG_BUILD OFF
  EXCLUDE_FROM_ALL 1
)

set(Tick_INCLUDE "${3RDPARTY_ROOT_PATH}/Source/Tick"
    CACHE STRING "Tick trait introspection and concepts library")
