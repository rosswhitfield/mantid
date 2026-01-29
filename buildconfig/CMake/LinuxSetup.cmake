# ######################################################################################################################
# Who are we
# ######################################################################################################################
include(DetermineLinuxDistro)

# ######################################################################################################################
# If required, find the address sanitizer library: libasan
# ######################################################################################################################
string(TOLOWER "${USE_SANITIZER}" USE_SANITIZERS_LOWER)
if(${USE_SANITIZERS_LOWER} MATCHES "address")
  find_package(AsanLib REQUIRED)
  if(ASANLIB_FOUND)
    set(ASAN_LIBRARY ${ASAN_LIBRARIES})
  endif(ASANLIB_FOUND)
endif()

# Tag used by dynamic loader to identify directory of loading library
set(DL_ORIGIN_TAG \$ORIGIN)

# ######################################################################################################################
# Set up package scripts for this distro
# ######################################################################################################################
include(LinuxPackageScripts)
