start_demo_dependencies()
find_package(MPI QUIET)
if(MPI_FOUND)
  message(STATUS "  C++")
  message(STATUS "    Compiler Wrapper: ${MPI_CXX_COMPILER} (this is not used)")
  message(STATUS "    Include Path:     ${MPI_CXX_INCLUDE_PATH}")
  message(STATUS "    Libraries:        ${MPI_CXX_LIBRARIES}")
  message(STATUS "    Link Flags:       ${MPI_CXX_LINK_FLAGS}")
  message(STATUS "  Executable:         ${MPIEXEC}")
  end_demo_dependencies()
else()
  message(STATUS "MPI not found. Skipping MPI demo.")
  end_demo_dependencies()
  return()
endif()

start_demo_targets()
set(executables
  mpi_demo
  provoke_open_requests
  split_comms
)

include_directories(${MPI_CXX_INCLUDE_PATH})

foreach(exe ${executables})
  add_executable(${exe} ${exe}.cpp)
  set_target_properties(${exe} PROPERTIES
    CXX_STANDARD "11"
    CXX_STANDARD_REQUIRED TRUE
    LINK_FLAGS ${MPI_CXX_LINK_FLAGS}
  )
  target_link_libraries(${exe}
    ${MPI_CXX_LIBRARIES}
  )
  message_demo_target(${exe})
endforeach(exe)
