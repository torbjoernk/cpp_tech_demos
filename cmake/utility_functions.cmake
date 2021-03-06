macro(ADD_DEMO _title _dir _option)
  if(${_option})
    message(STATUS "-------------------------------------------------------------------------------")
    message(STATUS "Demo: ${_title}")
    message(STATUS "")
    add_subdirectory(${_dir})
  endif()
endmacro()

macro(START_DEMO_DEPENDENCIES)
  message(STATUS "Configuring Dependencies")
endmacro()

macro(END_DEMO_DEPENDENCIES)
  message(STATUS "")
endmacro()

macro(START_DEMO_TARGETS)
  message(STATUS "Configuring Targets")
endmacro()

macro(ADD_DEMO_TARGET _name)
  message(STATUS "  ${_name} - run with 'run_${_name}'")
  add_custom_target(run_${_name}
    COMMAND ${_name}
    DEPENDS ${_name}
  )
endmacro()
