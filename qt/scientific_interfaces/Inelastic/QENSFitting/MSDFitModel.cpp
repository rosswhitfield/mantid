// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MSDFitModel.h"

using namespace Mantid::API;

namespace MantidQt::CustomInterfaces::Inelastic {

MSDFitModel::MSDFitModel() { m_fitType = MSDFIT_STRING; }

std::string MSDFitModel::getResultXAxisUnit() const { return "Temperature"; }

} // namespace MantidQt::CustomInterfaces::Inelastic
