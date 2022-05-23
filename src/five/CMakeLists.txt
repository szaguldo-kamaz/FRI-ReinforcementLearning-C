set(SRC
    FIVE.h
    FIVEGScFunc.c
    FIVEGVagEnv.c
    FIVEInit.c
    FIVEVagConcl.c
    FIVEVagConcl_FRIRL_BestAct.c
    FIVEVagConclWeight.c
    FIVEValVag.c
    five_vague_distance.c
    five_vague_distance_parallel.c
    five_rule_distance.c
    five_add_rule.c
    five_deinit.c
    five_remove_rule.c
)

add_library(five STATIC ${SRC})
target_include_directories(five PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(five PUBLIC ${PROJECT_BINARY_DIR})
target_link_libraries(five PUBLIC m)

#TODO: resolve include mpi.h conflict
if(BUILD_MPI)
message(${MPI_C_INCLUDE_PATH})
        target_include_directories(five PUBLIC ${MPI_C_INCLUDE_PATH})
        #target_link_libraries(five PUBLIC ${MPI_C_LIBRARIES})
endif()
