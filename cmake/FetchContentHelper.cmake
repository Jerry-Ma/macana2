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
    add_subdirectory(${${name}_SOURCE_DIR} ${${name}_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
endfunction(FetchContentHelper)
