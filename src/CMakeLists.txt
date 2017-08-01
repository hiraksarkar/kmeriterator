set(kmeriterator_lib_srcs FileOps.cpp FastxParser.cpp )
add_library(hfasta STATIC ${kmeriterator_lib_srcs})

add_executable(testkmerit MakeFastq.cpp)
target_link_libraries(testkmerit z hfasta)
