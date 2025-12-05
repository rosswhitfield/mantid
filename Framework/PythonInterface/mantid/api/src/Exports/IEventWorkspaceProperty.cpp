// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidAPI/IEventWorkspace.h"
#include "MantidKernel/PropertyWithValue.hxx"
#include "MantidPythonInterface/api/WorkspacePropertyExporter.h"
#include "MantidPythonInterface/core/GetPointer.h"

using Mantid::API::IEventWorkspace;
using Mantid::API::WorkspaceProperty; // NOLINT

GET_POINTER_SPECIALIZATION(WorkspaceProperty<IEventWorkspace>)

// Explicit template instantiation to ensure visibility across shared libraries
namespace Mantid {
namespace API {
extern template class WorkspaceProperty<IEventWorkspace>;
} // namespace API
} // namespace Mantid

void export_IEventWorkspaceProperty() {
  using Mantid::PythonInterface::WorkspacePropertyExporter;

  WorkspacePropertyExporter<IEventWorkspace>::define("IEventWorkspaceProperty");
}

// Explicit instantiation for PropertyWithValue<std::shared_ptr<IEventWorkspace>> for Python module
namespace Mantid::Kernel {
template class PropertyWithValue<std::shared_ptr<Mantid::API::IEventWorkspace>>;
} // namespace Mantid::Kernel
