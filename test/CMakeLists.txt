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
