add_executable(mountaincar
    mountaincar.c
)

if(BUILD_VISUALIZATION)
    target_link_libraries(mountaincar LINK_PUBLIC frirl frirlgui glut GL GLU)
else()
    target_link_libraries(mountaincar LINK_PUBLIC frirl)
endif()
