if(${Boost_FOUND})
  set(executables
    simple_comm_isend
    simple_comm_rma_get
    simple_comm_rma_get_onefile
  )

  include_directories(${Boost_INCLUDE_DIRS})

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
endif()
