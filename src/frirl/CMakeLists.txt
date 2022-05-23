set(SRC
    frirl.h frirl_types.h frirl_types_def.h
    frirl_init.c frirl_deinit.c frirl_init_ve.c frirl_init_rb.c
    frirl_check_possible_states.c frirl_e_greedy_selection.c frirl_episode.c frirl_get_best_action.c frirl_update_sarsa.c
        frirl_sequential_run.c frirl_agent.c

    frirl_utils.h frirl_utils.c

    frirl_app_helpers.h frirl_app_helpers.c

    frirl_test.h frirl_test_run.c

    frirl_imitation.h frirl_imitation.c
)

add_library(frirl STATIC ${SRC})

target_include_directories(frirl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(frirl PUBLIC ${PROJECT_BINARY_DIR})

target_link_libraries(frirl PUBLIC five)

if(BUILD_MPI)
message(${MPI_C_INCLUDE_PATH})
    target_include_directories(frirl PUBLIC ${MPI_C_INCLUDE_PATH})
    target_link_libraries(frirl PUBLIC ${MPI_C_LIBRARIES})
endif()

if(BUILD_VISUALIZATION)
    target_link_libraries(frirl PUBLIC frirlgui)
endif()
