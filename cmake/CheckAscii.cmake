cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED VELA_SOURCE_DIR)
    message(FATAL_ERROR "VELA_SOURCE_DIR must be set")
endif()

set(ASCII_CHECK_ROOTS
    tests
    src
    include
    examples
)

set(ASCII_CHECK_EXTENSIONS
    c
    cc
    cpp
    cxx
    h
    hh
    hpp
    hxx
    inl
    ipp
    tpp
    json
    txt
    md
    py
    cmake
)

set(files_to_check)
foreach(root IN LISTS ASCII_CHECK_ROOTS)
    if(NOT IS_DIRECTORY "${VELA_SOURCE_DIR}/${root}")
        continue()
    endif()

    foreach(extension IN LISTS ASCII_CHECK_EXTENSIONS)
        file(GLOB_RECURSE matched_files
            LIST_DIRECTORIES false
            "${VELA_SOURCE_DIR}/${root}/*.${extension}"
        )
        list(APPEND files_to_check ${matched_files})
    endforeach()
endforeach()

list(REMOVE_DUPLICATES files_to_check)
list(LENGTH files_to_check file_count)

set(offending_files)
foreach(file_path IN LISTS files_to_check)
    file(READ "${file_path}" file_hex HEX)
    string(REGEX REPLACE "(..)" "\\1;" file_bytes "${file_hex}")

    foreach(file_byte IN LISTS file_bytes)
        if(file_byte MATCHES "^[89A-Fa-f][0-9A-Fa-f]$")
            file(RELATIVE_PATH relative_path "${VELA_SOURCE_DIR}" "${file_path}")
            list(APPEND offending_files "${relative_path}")
            break()
        endif()
    endforeach()
endforeach()

if(offending_files)
    list(REMOVE_DUPLICATES offending_files)
    list(SORT offending_files)
    list(JOIN offending_files "\n  " offending_report)
    message(FATAL_ERROR
        "Non-ASCII bytes found in the configured ASCII check scope (tests/, src/, include/, examples/ and selected text/source extensions):\n"
        "  ${offending_report}\n\n"
        "Keep files under tests/, src/, include/, and examples/ ASCII-only for reliable MSYS2 UCRT64 builds and test selection."
    )
endif()

message(STATUS "ASCII check passed for ${file_count} files")
