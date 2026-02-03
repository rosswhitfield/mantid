# ######################################################################################################################
# Who are we
# ######################################################################################################################
include(DetermineLinuxDistro)

# ######################################################################################################################
# If required, find tcmalloc
# ######################################################################################################################
option(USE_TCMALLOC "If true, use LD_PRELOAD=libtcmalloc_minimal.so in startup scripts" ON)
# If not wanted, just carry on without it
if(USE_TCMALLOC)
  find_library(TCMALLOC_LIBRARIES tcmalloc_minimal)
  if(NOT TCMALLOC_LIBRARIES)
    message(STATUS "tcmalloc_minimal not found, skipping tcmalloc setup.")
    set(USE_TCMALLOC OFF)
  else()
    set(TCMALLOC_FOUND TRUE)
    message(STATUS "Found tcmalloc library: ${TCMALLOC_LIBRARIES}")
    set(TCMALLOC_RUNTIME_LIB "${TCMALLOC_LIBRARIES}")
  endif()
  # if it can't be found still carry on as the build will work. The package depenendencies will install it for the end
  # users
else(USE_TCMALLOC)
  message(STATUS "tcmalloc will not be included in startup scripts")
endif()

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
