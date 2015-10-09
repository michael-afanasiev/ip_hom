INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "/usr/local/Cellar/ipopt/3.12.1_3/include/coin"
  "/usr/local/include"
  )
LINK_DIRECTORIES(
  "/usr/local/lib"
  )

SET(cxxbind_SRCS
  src/main.cpp
  src/exodusFile.cpp
  src/exoMod1D_nlp.cpp
  src/model.cpp
  )

ADD_EXECUTABLE(HomTest ${cxxbind_SRCS})
TARGET_LINK_LIBRARIES(HomTest ipopt)
TARGET_LINK_LIBRARIES(HomTest exodus)
