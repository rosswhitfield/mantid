// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidKernel/PropertyWithValue.h"
#include "MantidDataObjects/DllConfig.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidDataObjects/GroupingWorkspace.h"
#include "MantidDataObjects/LeanElasticPeaksWorkspace.h"
#include "MantidDataObjects/MDEvent.h"
#include "MantidDataObjects/MDEventWorkspace.h"
#include "MantidDataObjects/MDHistoWorkspace.h"
#include "MantidDataObjects/MaskWorkspace.h"
#include "MantidDataObjects/OffsetsWorkspace.h"
#include "MantidDataObjects/PeaksWorkspace.h"
#include "MantidDataObjects/RebinnedOutput.h"
#include "MantidDataObjects/SpecialWorkspace2D.h"
#include "MantidDataObjects/SplittersWorkspace.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidDataObjects/WorkspaceSingleValue.h"

// PropertyWithValue implementation
#include "MantidKernel/PropertyWithValue.hxx"

namespace Mantid {
namespace DataObjects {
template <size_t nd> using MDEventWS = MDEventWorkspace<MDEvent<nd>, nd>;
template <size_t nd> using MDLeanEventWS = MDEventWorkspace<MDLeanEvent<nd>, nd>;
} // namespace DataObjects
namespace Kernel {

// Explicit static logger instantiations
#define PROPERTYWITHVALUE_LOGGER(type) template <> Logger PropertyWithValue<type>::g_logger("PropertyWithValue");

PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::EventWorkspace>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::GroupingWorkspace>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MaskWorkspace>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDEventWS<1>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDEventWS<2>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDEventWS<3>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDEventWS<4>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDEventWS<5>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDEventWS<6>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDEventWS<7>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDEventWS<8>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDEventWS<9>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDLeanEventWS<1>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDLeanEventWS<2>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDLeanEventWS<3>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDLeanEventWS<4>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDLeanEventWS<5>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDLeanEventWS<6>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDLeanEventWS<7>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDLeanEventWS<8>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDLeanEventWS<9>>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::MDHistoWorkspace>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::OffsetsWorkspace>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::PeaksWorkspace>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::LeanElasticPeaksWorkspace>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::RebinnedOutput>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::SpecialWorkspace2D>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::SplittersWorkspace>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::TableWorkspace>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::Workspace2D>)
PROPERTYWITHVALUE_LOGGER(std::shared_ptr<DataObjects::WorkspaceSingleValue>)

/// @cond
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::EventWorkspace>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::GroupingWorkspace>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MaskWorkspace>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDEventWS<1>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDEventWS<2>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDEventWS<3>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDEventWS<4>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDEventWS<5>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDEventWS<6>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDEventWS<7>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDEventWS<8>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDEventWS<9>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDLeanEventWS<1>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDLeanEventWS<2>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDLeanEventWS<3>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDLeanEventWS<4>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDLeanEventWS<5>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDLeanEventWS<6>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDLeanEventWS<7>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDLeanEventWS<8>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDLeanEventWS<9>>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::MDHistoWorkspace>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::OffsetsWorkspace>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::PeaksWorkspace>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::LeanElasticPeaksWorkspace>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::RebinnedOutput>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::SpecialWorkspace2D>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::SplittersWorkspace>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::TableWorkspace>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::Workspace2D>>;
template class MANTID_DATAOBJECTS_DLL PropertyWithValue<std::shared_ptr<DataObjects::WorkspaceSingleValue>>;
/// @endcond

} // namespace Kernel
} // namespace Mantid
