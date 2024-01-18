// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once
#include "../DllConfig.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "QPair"
#include "QString"
#include "QStringList"
#include <boost/optional.hpp>
#include <string>

namespace MantidQt {
namespace CustomInterfaces {
namespace WorkspaceManipulationUtils {

MANTIDQT_INELASTIC_DLL std::string getWorkspaceSuffix(const std::string &wsName);
MANTIDQT_INELASTIC_DLL QString getWorkspaceBasename(const QString &wsName);
MANTIDQT_INELASTIC_DLL std::unordered_map<std::string, size_t>
extractAxisLabels(const Mantid::API::MatrixWorkspace_const_sptr &workspace, const size_t &axisIndex);

MANTIDQT_INELASTIC_DLL std::string getEMode(const Mantid::API::MatrixWorkspace_sptr &ws);
MANTIDQT_INELASTIC_DLL double getEFixed(const Mantid::API::MatrixWorkspace_sptr &ws);

MANTIDQT_INELASTIC_DLL bool getResolutionRangeFromWs(const QString &workspace, QPair<double, double> &res);
MANTIDQT_INELASTIC_DLL bool getResolutionRangeFromWs(const Mantid::API::MatrixWorkspace_const_sptr &workspacee,
                                                     QPair<double, double> &res);

MANTIDQT_INELASTIC_DLL QPair<double, double>
getXRangeFromWorkspace(const Mantid::API::MatrixWorkspace_const_sptr &workspace, double precision = 0.00001);
MANTIDQT_INELASTIC_DLL QPair<double, double> getXRangeFromWorkspace(std::string const &workspaceName,
                                                                    double precision = 0.000001);

} // namespace WorkspaceManipulationUtils
} // namespace CustomInterfaces
} // namespace MantidQt
