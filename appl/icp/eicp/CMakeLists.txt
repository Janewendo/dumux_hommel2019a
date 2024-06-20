add_input_file_links()

dune_add_test(NAME columneicpbox
              SOURCES eicp.cc
              COMPILE_DEFINITIONS TYPETAG=EICPColumnBoxTypeTag)

#install headers
install(FILES
eicp.cc
eicpcolumnproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/appl/icp/eicp)
