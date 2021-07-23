# CMakeList.txt: projekt CMake dla elementu CTBN, dołącz źródło i zdefiniuj
# tutaj logikę specyficzną dla projektu.
#
cmake_minimum_required (VERSION 3.8)
add_library(Parameters "utils/parameters.cpp" "utils/parameters.h")
target_link_libraries(Parameters)
# Dodaj źródło do pliku wykonywalnego tego projektu.
add_executable (CTBN "CTBN.cpp" "CTBN.h")
target_link_libraries(CTBN PUBLIC Parameters)
# TODO: Dodaj testy i zainstaluj elementy docelowe w razie potrzeby.
