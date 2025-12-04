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
#include "MantidAPI/Workspace.h"
#include "MantidAPI/WorkspaceGroup.h"

// WorkspaceProperty implementation
#include "MantidAPI/WorkspaceProperty.hxx"

namespace Mantid::API {
// Explicit template instantiation definitions
template class MANTID_API_DLL WorkspaceProperty<Workspace>;
template class MANTID_API_DLL WorkspaceProperty<IEventWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<IMDEventWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<IMDHistoWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<IMDWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<MatrixWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<IPeaksWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<ITableWorkspace>;
template class MANTID_API_DLL WorkspaceProperty<WorkspaceGroup>;
} // namespace Mantid::API
