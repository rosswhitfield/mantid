// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidDataObjects/MDHistoWorkspace.h"
#include "MantidKernel/PropertyWithValue.hxx"
#include "MantidPythonInterface/api/RegisterWorkspacePtrToPython.h"
#include "MantidPythonInterface/core/GetPointer.h"

#include <boost/python/class.hpp>

using Mantid::API::IMDHistoWorkspace;
using Mantid::DataObjects::MDHistoWorkspace;
using namespace Mantid::PythonInterface::Registry;
using namespace boost::python;

GET_POINTER_SPECIALIZATION(MDHistoWorkspace)

void export_MDHistoWorkspace() {
  class_<MDHistoWorkspace, bases<IMDHistoWorkspace>, boost::noncopyable>("MDHistoWorkspace", no_init);

  // register pointers
  RegisterWorkspacePtrToPython<MDHistoWorkspace>();
}

template class Mantid::Kernel::PropertyWithValue<std::shared_ptr<Mantid::DataObjects::MDHistoWorkspace>>;
