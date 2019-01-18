include(FetchContent)
function(FetchContentHelper name git tag)
FetchContent_Declare(
    ${name}
    GIT_REPOSITORY ${git}
    GIT_TAG        ${tag}
)
FetchContent_GetProperties(${name})
if(NOT ${name}_POPULATED)
    FetchContent_Populate(${name})
    set(${name}_SOURCE_DIR ${${name}_SOURCE_DIR} PARENT_SCOPE)
    set(${name}_BINARY_DIR ${${name}_BINARY_DIR} PARENT_SCOPE)
    message(argv3 ${ARGV3})
    if (ARGV3 STREQUAL "NO_ADD_SUBDIR")
    else()
        add_subdirectory(${${name}_SOURCE_DIR} ${${name}_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()
endif()
endfunction(FetchContentHelper)
