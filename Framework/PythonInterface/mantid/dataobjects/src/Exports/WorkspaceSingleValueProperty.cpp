// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2024 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidDataObjects/WorkspaceSingleValue.h"
#include "MantidPythonInterface/api/WorkspacePropertyExporter.h"
#include "MantidPythonInterface/core/GetPointer.h"

using Mantid::DataObjects::WorkspaceSingleValue;

GET_POINTER_SPECIALIZATION(Mantid::API::WorkspaceProperty<WorkspaceSingleValue>)

// Explicit template instantiation to ensure visibility across shared libraries
namespace Mantid {
namespace API {
extern template class WorkspaceProperty<Mantid::DataObjects::WorkspaceSingleValue>;
} // namespace API
} // namespace Mantid

void export_WorkspaceSingleValueProperty() {
  using Mantid::PythonInterface::WorkspacePropertyExporter;
  WorkspacePropertyExporter<WorkspaceSingleValue>::define("WorkspaceSingleValueProperty");
}
