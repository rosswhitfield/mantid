// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidDataObjects/PeaksWorkspace.h"
#include "MantidPythonInterface/api/WorkspacePropertyExporter.h"
#include "MantidPythonInterface/core/GetPointer.h"

using Mantid::API::WorkspaceProperty; // NOLINT
using Mantid::DataObjects::PeaksWorkspace;

GET_POINTER_SPECIALIZATION(WorkspaceProperty<PeaksWorkspace>)

// Explicit template instantiation to ensure visibility across shared libraries
namespace Mantid {
namespace API {
extern template class WorkspaceProperty<Mantid::DataObjects::PeaksWorkspace>;
} // namespace API
} // namespace Mantid

void export_PeaksWorkspaceProperty() {
  using Mantid::PythonInterface::WorkspacePropertyExporter;
  WorkspacePropertyExporter<PeaksWorkspace>::define("PeaksWorkspaceProperty");
}
