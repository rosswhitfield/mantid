// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once
#include "Common/First.h"
#include "GetInstrumentParameter.h"
#include <optional>
#include <vector>

namespace MantidQt {
namespace CustomInterfaces {
namespace ISISReflectometry {

template <typename T>
std::optional<T> firstFromParameterFile(Mantid::Geometry::Instrument_const_sptr instrument,
                                        std::string const &parameterName) {
  return first(getInstrumentParameter<T>(instrument, parameterName));
}

template <typename... Ts>
std::optional<boost::variant<Ts...>> firstFromParameterFileVariant(Mantid::Geometry::Instrument_const_sptr instrument,
                                                                   std::string const &parameterName) {
  auto values = getInstrumentParameter<boost::variant<Ts...>>(instrument, parameterName);
  return boost::apply_visitor(FirstVisitor<Ts...>(), values);
}

class MissingInstrumentParameterValue {
public:
  explicit MissingInstrumentParameterValue(std::string const &parameterName) : m_parameterName(parameterName) {}

  std::string const &parameterName() const;

private:
  std::string m_parameterName;
};

class InstrumentParameters {
public:
  explicit InstrumentParameters(Mantid::Geometry::Instrument_const_sptr instrument);

  template <typename T> std::optional<T> optional(std::string const &parameterName) {
    return fromFile<T>(parameterName);
  }

  template <typename Default, typename T>
  T handleMandatoryIfMissing(std::optional<T> const &value, std::string const &parameterName) {
    if (value) {
      return value.value();
    }
    m_missingValueErrors.emplace_back(parameterName);
    return Default();
  }

  template <typename T> T mandatory(std::string const &parameterName) {
    try {
      return handleMandatoryIfMissing<T>(firstFromParameterFile<T>(m_instrument, parameterName), parameterName);
    } catch (InstrumentParameterTypeMissmatch const &ex) {
      m_typeErrors.emplace_back(ex);
      return T();
    }
  }

  /**
   * Tries to get the value of a property which may hold a value of any of the
   * specified types.
   *
   * Will try the types in the order specified from left to right.
   * Returns a variant containing a default constructed value of the first type
   * and records the parameter as missing in the event that the property is not
   * in the file.
   *
   * If the property is in the file but is not one of the specified types a type
   * missmatch error will be recorded.
   */
  template <typename T, typename... Ts> boost::variant<T, Ts...> mandatoryVariant(std::string const &parameterName) {
    try {
      return handleMandatoryIfMissing<T>(firstFromParameterFileVariant<T, Ts...>(m_instrument, parameterName),
                                         parameterName);
    } catch (InstrumentParameterTypeMissmatch const &ex) {
      m_typeErrors.emplace_back(ex);
      return T();
    }
  }

  std::vector<InstrumentParameterTypeMissmatch> const &typeErrors() const;
  bool hasTypeErrors() const;
  std::vector<MissingInstrumentParameterValue> const &missingValues() const;
  bool hasMissingValues() const;

private:
  template <typename T> std::optional<T> fromFile(std::string const &parameterName) {
    try {
      return firstFromParameterFile<T>(m_instrument, parameterName);
    } catch (InstrumentParameterTypeMissmatch const &ex) {
      m_typeErrors.emplace_back(ex);
      return std::nullopt;
    }
  }

  Mantid::Geometry::Instrument_const_sptr m_instrument;
  std::vector<InstrumentParameterTypeMissmatch> m_typeErrors;
  std::vector<MissingInstrumentParameterValue> m_missingValueErrors;
};
} // namespace ISISReflectometry
} // namespace CustomInterfaces
} // namespace MantidQt
