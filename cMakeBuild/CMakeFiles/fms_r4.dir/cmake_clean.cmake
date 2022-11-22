file(REMOVE_RECURSE
  "libfms_r4.a"
  "libfms_r4.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C Fortran)
  include(CMakeFiles/fms_r4.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
