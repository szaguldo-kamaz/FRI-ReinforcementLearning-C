add_executable(cartpole
    cartpole.c
)

if(BUILD_VISUALIZATION)
    target_link_libraries(cartpole LINK_PUBLIC frirl frirlgui glut GL GLU)
else()
    target_link_libraries(cartpole LINK_PUBLIC frirl)
endif()
