include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
)
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory("${googletest_SOURCE_DIR}" "${googletest_BINARY_DIR}")
endif()

include(GoogleTest)

add_executable(simple_test simple_test.cpp)
target_link_libraries(simple_test BiLAD gtest_main)
#gtest_add_tests(simple_test "" AUTO)
gtest_discover_tests(simple_test)

add_executable(random_test random_test.cpp)
target_link_libraries(random_test BiLAD gtest_main)
gtest_discover_tests(random_test)

add_executable(m_larac PSQSR_test.cpp)
target_link_libraries(m_larac BiLAD gtest_main)
gtest_discover_tests(m_larac)

add_executable(new_bilad new_bilad_test.cpp)
target_link_libraries(new_bilad BiLAD gtest_main)
gtest_discover_tests(new_bilad)

message(STATUS "copying test data files")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/data" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

