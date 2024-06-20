add_input_file_links()

dune_add_test(NAME columneicpbox
              SOURCES eicp.cc
              COMPILE_DEFINITIONS TYPETAG=EICPColumnBoxTypeTag)

dune_add_test(NAME columnleobox
              SOURCES leo.cc
              COMPILE_DEFINITIONS TYPETAG=LEOColumnBoxTypeTag)
  
#install headers
install(FILES
eicp.cc
eicpcolumnproblem.hh
leo.cc
leocolumnproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/appl/icp/eicp)

