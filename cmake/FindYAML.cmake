# Finds libyaml and provides the imported target PkgConfig::YAML
# Tries pkg-config (yaml-0.1) first, then falls back to find_path/find_library
# so a plain -lyaml is still used when no .pc file is available.

if(TARGET PkgConfig::YAML)
  set(YAML_FOUND TRUE)
  return()
endif()

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(YAML QUIET IMPORTED_TARGET "yaml-0.1")
endif()

if(NOT YAML_FOUND)
  find_path(YAML_INCLUDE_DIR NAMES yaml.h)
  find_library(YAML_LIBRARY NAMES yaml)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(YAML DEFAULT_MSG YAML_LIBRARY YAML_INCLUDE_DIR)
  if(YAML_FOUND)
    add_library(PkgConfig::YAML UNKNOWN IMPORTED)
    set_target_properties(PkgConfig::YAML PROPERTIES
      IMPORTED_LOCATION "${YAML_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${YAML_INCLUDE_DIR}")
  endif()
endif()
