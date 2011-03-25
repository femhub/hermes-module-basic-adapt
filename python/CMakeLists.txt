include_directories(${PROJECT_BINARY_DIR}/src)

set(CYTHON_INCLUDE_DIRECTORIES -I.)
set(CYTHON_INCLUDE_DIRECTORIES ${CYTHON_INCLUDE_DIRECTORIES} -I${HERMES_PYTHON_WRAPPERS_DIR})

CYTHON_ADD_MODULE(basicadapt)
target_link_libraries(basicadapt hermes_basicadapt)

set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes_modules)
install(TARGETS basicadapt
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
install(FILES __init__.py basicadapt_defs.pxd
    DESTINATION ${PY_PATH}
    )
