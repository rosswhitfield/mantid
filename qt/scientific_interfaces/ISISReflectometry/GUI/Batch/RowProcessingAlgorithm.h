// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "Common/DllConfig.h"
#include "MantidQtWidgets/Common/BatchAlgorithmRunner.h"
#include <boost/optional.hpp>
#include <map>
#include <string>
#include <vector>

namespace MantidQt {
namespace CustomInterfaces {
namespace ISISReflectometry {
class Batch;
class Row;

MANTIDQT_ISISREFLECTOMETRY_DLL MantidQt::API::IConfiguredAlgorithm_sptr createConfiguredAlgorithm(Batch const &model,
                                                                                                  Row &row);

MANTIDQT_ISISREFLECTOMETRY_DLL API::IConfiguredAlgorithm::AlgorithmRuntimeProps
createAlgorithmRuntimeProps(Batch const &model, Row const &row);
MANTIDQT_ISISREFLECTOMETRY_DLL API::IConfiguredAlgorithm::AlgorithmRuntimeProps
createAlgorithmRuntimeProps(Batch const &model);
} // namespace ISISReflectometry
} // namespace CustomInterfaces
} // namespace MantidQt
