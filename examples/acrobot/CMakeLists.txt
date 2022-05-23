add_executable(acrobot
    acrobot.c
)

if(BUILD_VISUALIZATION)
    target_link_libraries(acrobot LINK_PUBLIC frirl frirlgui glut GL GLU)
else()
    target_link_libraries(acrobot LINK_PUBLIC frirl)
endif()
