
start_demo_targets()

set(EXAMPLES
  exceptions
  inheritance
  mixins
  template_inherit
)

foreach(example ${EXAMPLES})
  add_executable(${example}
    ${example}.cpp
  )

  set_target_properties(${example} PROPERTIES
    CXX_STANDARD "11"
    CXX_STANDARD_REQUIRED TRUE
  )

  add_demo_target(${example})
endforeach(example)
