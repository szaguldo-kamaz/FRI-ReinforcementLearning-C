set(EXAMPLE_acrobot TRUE CACHE BOOL "acrobot example")
if (EXAMPLE_acrobot)
    add_subdirectory(acrobot)
endif()

set(EXAMPLE_cartpole TRUE CACHE BOOL "cartpole example")
if (EXAMPLE_cartpole)
    add_subdirectory(cartpole)
endif()

set(EXAMPLE_mountaincar TRUE CACHE BOOL "mountaincar example")
if (EXAMPLE_mountaincar)
    add_subdirectory(mountaincar)
endif()
