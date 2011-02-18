project(module-basicadapt)

add_library(hermes_basicadapt SHARED basicadapt.cpp)
target_link_libraries(hermes_basicadapt ${HERMES_REAL_BIN})

add_executable(${PROJECT_NAME} main.cpp disc.cpp ../basic/basic.cpp)
target_link_libraries(${PROJECT_NAME} hermes_basicadapt)
include (../CMake.common)

install(TARGETS hermes_basicadapt
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
