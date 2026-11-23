cmake_minimum_required(VERSION 3.20)
project(cadcore LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(CADCORE_SANITIZE "Build with AddressSanitizer + UndefinedBehaviorSanitizer" OFF)
option(CADCORE_TESTS "Build tests and tools" ON)


# --- library -------------------------------------------------------------------------------------
add_library(cadcore STATIC
  src/core/core.cpp
  src/core/stats.cpp
  src/core/compare.cpp
  src/core/text.cpp
  src/core/parallel.cpp
  src/codec/inflate.cpp
  src/codec/xml.cpp
  src/codec/json.cpp
  src/codec/meshopt.cpp
  src/codec/draco.cpp
  src/codec/draco/bitstream.cpp
  src/codec/draco/connectivity.cpp
  src/codec/draco/attributes.cpp
  src/codec/draco/kd_tree.cpp
  src/formats/registry.cpp
  src/formats/zip.cpp
  src/formats/stl.cpp
  src/formats/obj.cpp
  src/formats/ply.cpp
  src/formats/off.cpp
  src/formats/gltf.cpp
  src/formats/threemf.cpp
  src/formats/catpart.cpp
  src/formats/catproduct.cpp
  src/formats/threedxml.cpp
  src/formats/cgr.cpp
  src/formats/catia/osm.cpp
  src/formats/catia/spec.cpp
  src/formats/catia/document.cpp
  src/formats/catia/tessellation.cpp
  src/formats/catia/store.cpp
  src/formats/drc.cpp
  src/formats/dxf.cpp
  src/formats/dxf/groups.cpp
  src/formats/dxf/document.cpp
  src/formats/dxf/geometry.cpp
  src/formats/dxf/draw.cpp
  src/formats/dxf/colors.cpp
  src/formats/dxf/proxy.cpp
  src/formats/dxf/clip.cpp)
target_include_directories(cadcore PUBLIC include PRIVATE src)
find_package(Threads REQUIRED)
target_link_libraries(cadcore PUBLIC Threads::Threads)

set(CADCORE_WARNINGS -Wall -Wextra -Wpedantic -Wshadow -Wno-unused-parameter)
target_compile_options(cadcore PRIVATE ${CADCORE_WARNINGS})

if(CADCORE_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
  add_link_options(-fsanitize=address,undefined)
  target_compile_options(cadcore PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
endif()

# --- test-part generator, tools, tests ----------------------------------------------------------
if(CADCORE_TESTS AND NOT EMSCRIPTEN)
  add_library(cadcore_testgen STATIC testgen/models.cpp testgen/writers.cpp testgen/zipwriter.cpp)
  target_include_directories(cadcore_testgen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} PRIVATE src)
  target_link_libraries(cadcore_testgen PUBLIC cadcore)
  target_compile_options(cadcore_testgen PRIVATE ${CADCORE_WARNINGS})

  add_executable(cadcli tools/cadcli.cpp)
  target_link_libraries(cadcli PRIVATE cadcore cadcore_testgen)
  target_include_directories(cadcli PRIVATE src)  # debug commands read internal CATIA structures

  add_executable(cadcore_tests tests/test_main.cpp tests/test_generated.cpp tests/test_robustness.cpp tests/test_codecs.cpp tests/test_draco.cpp tests/test_catproduct.cpp tests/test_threedxml.cpp tests/test_cgr.cpp tests/test_catia_grammar.cpp tests/test_dxf.cpp tests/assembly_gt.cpp)
  target_include_directories(cadcore_tests PRIVATE src)
  set(CADCORE_CATIA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../samples/catia" CACHE PATH "Folder with the CATIA ground-truth projects (CATProduct + STEP exports)")
  target_compile_definitions(cadcore_tests PRIVATE CADCORE_SAMPLES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../samples"
                             CADCORE_FIXTURES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures" CADCORE_CATIA_DIR="${CADCORE_CATIA_DIR}")
  target_link_libraries(cadcore_tests PRIVATE cadcore cadcore_testgen)
  target_compile_options(cadcore_tests PRIVATE ${CADCORE_WARNINGS})

  # Speed / complexity (separate: timing-based, run on a quiet machine).
  add_executable(cadcore_perf tests/test_main.cpp tests/perf.cpp)
  target_include_directories(cadcore_perf PRIVATE src)
  target_link_libraries(cadcore_perf PRIVATE cadcore cadcore_testgen)
  target_compile_options(cadcore_perf PRIVATE ${CADCORE_WARNINGS})

  enable_testing()
  add_test(NAME unit COMMAND cadcore_tests)
  add_test(NAME perf COMMAND cadcore_perf)
endif()

# --- WebAssembly module ------------------------------------------------------------------------
if(EMSCRIPTEN)
  add_executable(cadcore_wasm wasm/bindings.cpp)
  target_link_libraries(cadcore_wasm PRIVATE cadcore)
  target_link_options(cadcore_wasm PRIVATE
    -lembind -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createCadcore -sALLOW_MEMORY_GROWTH=1
    -sENVIRONMENT=web,worker,node -sFILESYSTEM=0 -fexceptions)
  target_compile_options(cadcore_wasm PRIVATE -fexceptions)
  set_target_properties(cadcore_wasm PROPERTIES OUTPUT_NAME cadcore)
endif()
