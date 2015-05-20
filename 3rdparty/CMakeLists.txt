message(STATUS "Preparing potential 3rd-party libraries")

set(3RDPARTY_DOWNLOAD_TIMEOUT 60 PARENT_SCOPE)

if(${USE_GIT_PROTOCOL})
  message(STATUS "  using git protocol for Github")
  set(GITHUB_PREFIX "git@github.com:" PARENT_SCOPE)
else()
  message(STATUS "  using HTTPS protocol for Github")
  set(GITHUB_PREFIX "https://github.com/" PARENT_SCOPE)
endif()
message(STATUS "    you can change that by setting USE_GIT_PROTOCOL=[ON,OFF]")
