list(APPEND TEST_sources ${CMAKE_CURRENT_SOURCE_DIR}/GTest.cpp)

list(APPEND TEST_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(TEST_sources ${TEST_sources} PARENT_SCOPE)
set(TEST_includes ${TEST_includes} PARENT_SCOPE)
