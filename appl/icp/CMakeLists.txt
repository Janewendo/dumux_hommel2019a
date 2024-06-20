add_subdirectory(eicp)

add_input_file_links()
dune_symlink_to_source_files(FILES injections)

#install headers
install(FILES
icpspatialparams.hh
co2tableslaboratoryhightemp.hh
co2valueslaboratoryhightemp.inc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/appl/icp)