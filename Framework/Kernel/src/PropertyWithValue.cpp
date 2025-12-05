// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidKernel/PropertyWithValue.h"
#include "MantidKernel/Matrix.h"

// PropertyWithValue implementation
#include "MantidKernel/OptionalBool.h"
#include "MantidKernel/PropertyManager.h"
#include "MantidKernel/PropertyWithValue.hxx"
#include "MantidNexus/NexusFile.h"

#include <vector>

namespace Mantid::Kernel {

// Explicit instantiations for basic types
template class PropertyWithValue<bool>;
template class PropertyWithValue<int>;
template class PropertyWithValue<long long>;
template class PropertyWithValue<int64_t>;
template class PropertyWithValue<uint16_t>;
template class PropertyWithValue<uint32_t>;
template class PropertyWithValue<unsigned long long>;
template class PropertyWithValue<uint64_t>;
template class PropertyWithValue<float>;
template class PropertyWithValue<double>;
template class PropertyWithValue<std::string>;

// Explicit instantiations for vector types used by ArrayProperty
template class PropertyWithValue<std::vector<bool>>;
template class PropertyWithValue<std::vector<int>>;
template class PropertyWithValue<std::vector<long long>>;
template class PropertyWithValue<std::vector<int64_t>>;
template class PropertyWithValue<std::vector<uint16_t>>;
template class PropertyWithValue<std::vector<uint32_t>>;
template class PropertyWithValue<std::vector<unsigned long long>>;
template class PropertyWithValue<std::vector<uint64_t>>;
template class PropertyWithValue<std::vector<float>>;
template class PropertyWithValue<std::vector<double>>;
template class PropertyWithValue<std::vector<std::string>>;

// OptionalBool type used in Python exports
template class PropertyWithValue<OptionalBool>;

} // namespace Mantid::Kernel
