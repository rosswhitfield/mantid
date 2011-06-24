###############################################################################
# Putting all the common CPack stuff in one place
###############################################################################

# Common description stuff
set ( CPACK_PACKAGE_DESCRIPTION_SUMMARY "Neutron Scattering Data Analysis" )
set ( CPACK_PACKAGE_VENDOR "ISIS Rutherford Appleton Laboratory and NScD Oak Ridge National Laboratory" )
set ( CPACK_PACKAGE_URL http://www.mantidproject.org/ )
set ( CPACK_PACKAGE_CONTACT mantid-help@mantidproject.org )
set ( CPACK_PACKAGE_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${MtdVersion_WC_LAST_CHANGED_REV} )
set ( CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR} )
set ( CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR} )
set ( CPACK_PACKAGE_VERSION_PATCH ${MtdVersion_WC_LAST_CHANGED_REV} )

# RPM information - only used if generating a rpm
set ( CPACK_RPM_PACKAGE_LICENSE GPLv3+ )
set ( CPACK_RPM_PACKAGE_RELEASE 1 )
set ( CPACK_RPM_PACKAGE_GROUP Applications/Engineering )

# DEB information - only used if generating a deb
set ( CPACK_DEBIAN_PACKAGE_RELEASE 1 )
