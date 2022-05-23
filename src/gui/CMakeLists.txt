set(SRC
    color.c
    color.h
    draw.c
    draw.h
    gui.h
    init.c
    init.h
)

add_library(frirlgui STATIC ${SRC})
target_include_directories(frirlgui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(frirlgui PUBLIC ${PROJECT_BINARY_DIR})
