# Links a Python extension module.
#
# The exact link flags are platform dependent and this macro makes it possible
# to write platform independent cmakefiles. All you have to do is to change
# this:
#
# add_library(simple_wrapper SHARED ${SRC})  # Linux only
# set_target_properties(simple_wrapper PROPERTIES PREFIX "")
#
# to this:
#
# add_python_library(simple_wrapper ${SRC})  # Platform independent
#
# Full example:
#
# set(SRC
#     iso_c_utilities.f90
#     pde_pointers.f90
#     example1.f90
#     example2.f90
#     example_eigen.f90
#     simple.f90
#     simple_wrapper.c
# )
# add_python_library(simple_wrapper ${SRC})

macro(ADD_PYTHON_LIBRARY name)
    # When linking Python extension modules, a special care must be taken about
    # the link flags, which are platform dependent:
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # on Mac, we need to use the "-bundle" gcc flag, which is what MODULE
        # does:
        add_library(${name} MODULE ${ARGN})
        # and "-flat_namespace -undefined suppress" link flags, that we need
        # to add by hand:
        set_target_properties(${name} PROPERTIES
            LINK_FLAGS "-flat_namespace -undefined suppress")
    ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # on Linux, we need to use the "-shared" gcc flag, which is what SHARED
        # does:
        add_library(${name} SHARED ${ARGN})
    ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${name} PROPERTIES PREFIX "")
    IF(${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
            target_link_libraries(${name} ${PYTHON_LIBRARY})
    ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
endmacro(ADD_PYTHON_LIBRARY)
