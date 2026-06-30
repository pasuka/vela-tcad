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
    # Allow an optional leading UTF-8 BOM while rejecting all other non-ASCII bytes.
    string(REGEX REPLACE "^[Ee][Ff][Bb][Bb][Bb][Ff]" "" file_hex "${file_hex}")
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

set(SCALING_SCHEMA_DOCS
    README.md
    docs/config_schema.md
    docs/examples.md
)

set(missing_scaling_docs)
foreach(doc_path IN LISTS SCALING_SCHEMA_DOCS)
    set(full_doc_path "${VELA_SOURCE_DIR}/${doc_path}")
    if(NOT EXISTS "${full_doc_path}")
        list(APPEND missing_scaling_docs "${doc_path} (missing file)")
        continue()
    endif()

    file(READ "${full_doc_path}" doc_text)
    string(FIND "${doc_text}" "unit_scaling" has_unit_scaling)
    string(FIND "${doc_text}" "legacy SI" has_legacy_si)
    string(FIND "${doc_text}" "No `scaling` field" has_no_scaling)

    if(has_unit_scaling EQUAL -1
       OR has_legacy_si EQUAL -1
       OR has_no_scaling EQUAL -1)
        list(APPEND missing_scaling_docs "${doc_path}")
    endif()
endforeach()

if(missing_scaling_docs)
    list(JOIN missing_scaling_docs "\n  " missing_scaling_report)
    message(FATAL_ERROR
        "Scaling mode documentation is incomplete. Each public schema doc must mention "
        "`unit_scaling`, `legacy SI`, and `No `scaling` field`:\n"
        "  ${missing_scaling_report}"
    )
endif()

foreach(doc_path IN LISTS SCALING_SCHEMA_DOCS)
    file(READ "${VELA_SOURCE_DIR}/${doc_path}" doc_text)
    string(TOLOWER "${doc_text}" doc_text_lower)
    string(FIND "${doc_text_lower}" "sentaurus" has_forbidden_sentaurus)
    string(FIND "${doc_text_lower}" "\"system\": \"si\"" has_forbidden_system_si)
    string(FIND "${doc_text_lower}" "\"system\":\"si\"" has_forbidden_system_si_compact)

    if(NOT has_forbidden_sentaurus EQUAL -1
       OR NOT has_forbidden_system_si EQUAL -1
       OR NOT has_forbidden_system_si_compact EQUAL -1)
        message(FATAL_ERROR
            "Forbidden scaling schema terminology found in ${doc_path}. "
            "Do not document commercial software names or a `system: si` mode."
        )
    endif()
endforeach()

message(STATUS "Scaling mode documentation check passed")
