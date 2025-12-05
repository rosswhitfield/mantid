// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidAPI/IPeaksWorkspace.h"
#include "MantidKernel/PropertyWithValue.hxx"
#include "MantidPythonInterface/api/WorkspacePropertyExporter.h"
#include "MantidPythonInterface/core/GetPointer.h"

using Mantid::API::IPeaksWorkspace;
using Mantid::API::WorkspaceProperty; // NOLINT

GET_POINTER_SPECIALIZATION(WorkspaceProperty<IPeaksWorkspace>)

// Explicit template instantiation to ensure visibility across shared libraries
namespace Mantid {
namespace API {
extern template class WorkspaceProperty<IPeaksWorkspace>;
} // namespace API
} // namespace Mantid

void export_IPeaksWorkspaceProperty() {
  using Mantid::PythonInterface::WorkspacePropertyExporter;
  WorkspacePropertyExporter<IPeaksWorkspace>::define("IPeaksWorkspaceProperty");
}

namespace Mantid::Kernel {
template class PropertyWithValue<std::shared_ptr<Mantid::API::IPeaksWorkspace>>;
}
