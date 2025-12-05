// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidAPI/IEventWorkspace.h"
#include "MantidAPI/IMDEventWorkspace.h"
#include "MantidAPI/IMDHistoWorkspace.h"
#include "MantidAPI/IMDWorkspace.h"
#include "MantidAPI/IPeaksWorkspace.h"
#include "MantidAPI/ISplittersWorkspace.h"
#include "MantidAPI/ITableWorkspace.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/Workspace.h"
#include "MantidAPI/WorkspaceGroup.h"

// WorkspaceProperty implementation
#include "MantidAPI/WorkspaceProperty.hxx"
#include "MantidKernel/PropertyWithValue.hxx"

namespace Mantid::API {
// Explicit template instantiation definitions
template class MANTID_API_DLL WorkspaceProperty<Workspace>;
template class MANTID_API_DLL WorkspaceProperty<MatrixWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<IEventWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<IMDWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<IMDEventWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<IMDHistoWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<IPeaksWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<ITableWorkspace>;
// ISplittersWorkspace excluded - doesn't have getName() and other required methods
template class MANTID_API_DLL WorkspaceProperty<WorkspaceGroup>;
} // namespace Mantid::API

// Explicit instantiations for PropertyWithValue<std::shared_ptr<WorkspaceType>>
namespace Mantid::Kernel {
template class PropertyWithValue<std::shared_ptr<Mantid::API::Workspace>>;
template class PropertyWithValue<std::shared_ptr<Mantid::API::MatrixWorkspace>>;
template class PropertyWithValue<std::shared_ptr<Mantid::API::IEventWorkspace>>;
template class PropertyWithValue<std::shared_ptr<Mantid::API::IMDWorkspace>>;
template class PropertyWithValue<std::shared_ptr<Mantid::API::IMDEventWorkspace>>;
template class PropertyWithValue<std::shared_ptr<Mantid::API::IMDHistoWorkspace>>;
template class PropertyWithValue<std::shared_ptr<Mantid::API::IPeaksWorkspace>>;
template class PropertyWithValue<std::shared_ptr<Mantid::API::ITableWorkspace>>;
// ISplittersWorkspace excluded - doesn't work with WorkspaceProperty
template class PropertyWithValue<std::shared_ptr<Mantid::API::WorkspaceGroup>>;
} // namespace Mantid::Kernel
