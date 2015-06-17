start_demo_dependencies()
message(STATUS "MPI")
find_package(MPI QUIET)
if(MPI_FOUND)
  message(STATUS "  C++")
  message(STATUS "    Compiler Wrapper: ${MPI_CXX_COMPILER} (this is not used)")
  message(STATUS "    Include Path:     ${MPI_CXX_INCLUDE_PATH}")
  message(STATUS "    Libraries:        ${MPI_CXX_LIBRARIES}")
  message(STATUS "    Link Flags:       ${MPI_CXX_LINK_FLAGS}")
  message(STATUS "  Executable:         ${MPIEXEC}")
else()
  message(STATUS "MPI not found. Skipping MPI demo.")
  end_demo_dependencies()
  return()
endif()

find_package(Boost)
if(NOT ${Boost_FOUND})
  message(STATUS "  simple_comm requires Boost.Format")
endif()
end_demo_dependencies()

start_demo_targets()
set(executables
  mpi_demo
  one_sided_fenced
  one_sided_mutex
  provoke_open_requests
  split_comms
)
if(${Boost_FOUND})
  list(APPEND executables simple_comm simple_comm_rma)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

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
  add_demo_target(${exe})
endforeach(exe)
