CYTHON_ADD_MODULE_PYX(basicadapt)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/modules/basic)
include_directories(${PROJECT_BINARY_DIR}/modules/basicadapt)
add_python_library(basicadapt basicadapt.cpp)
target_link_libraries(basicadapt hermes_basic hermes_basicadapt)

set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes2d/modules/basicadapt)
install(TARGETS basicadapt
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
