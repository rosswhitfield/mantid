// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL-3.0 +
#include "MantidPythonInterface/core/PythonObjectProperty.h"
#include "MantidJson/Json.h"
#include "MantidKernel/Exception.h"

// Forward declare the specialization before including .hxx to avoid "after instantiation" errors
namespace Mantid::Kernel {
template <>
PropertyWithValue<boost::python::object>::PropertyWithValue(std::string name, boost::python::object defaultValue,
                                                            IValidator_sptr validator, unsigned int direction);

template <>
PropertyWithValue<boost::python::object>::PropertyWithValue(const std::string &name,
                                                            const boost::python::object &defaultValue,
                                                            const std::string &defaultValueStr,
                                                            IValidator_sptr validator, unsigned int direction);

template <>
PropertyWithValue<boost::python::object> &
PropertyWithValue<boost::python::object>::operator=(const boost::python::object &value);
template <> const boost::python::object &PropertyWithValue<boost::python::object>::operator()() const;

template <> PropertyWithValue<boost::python::object>::operator const boost::python::object &() const;

template <> void PropertyWithValue<boost::python::object>::replaceValidator(IValidator_sptr newValidator);
} // namespace Mantid::Kernel

#include "MantidKernel/PropertyWithValue.hxx"
#include "MantidKernel/PropertyWithValueJSON.h"
#include "MantidPythonInterface/core/GlobalInterpreterLock.h"
#include "MantidPythonInterface/core/IsNone.h"

#include <boost/python/dict.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/import.hpp>
#include <boost/python/list.hpp>

namespace {
namespace bp = boost::python;

/** Atomic types which can be serialized by python's json.dumps */
std::set<std::string> const jsonAllowedTypes{"int", "float", "str", "NoneType", "bool"};

inline bool isJsonAtomic(bp::object const &obj) {
  std::string const objname = bp::extract<std::string>(obj.attr("__class__").attr("__name__"));
  return jsonAllowedTypes.count(objname) > 0;
}

/** Recurisvely JSONify the object, including replacing custom class objects with a dictionary representation
 * At each stage in the process, will check if the value is one of the above types, or a list-like object.
 * If the object is one of the atomic types, it remains.  If it is a list-like, the below is run on each individual
 * element. If not, the object is replaced with a dictionary, run recursively on each individual value
 * @param obj the object to be recursively written
 * @param depth the current recursion depth, which is limited to ten; beyond this, ellipses will be printed
 * @return a python object which is a dictionary corresponding to `obj`
 */
bp::object recursiveDictDump(bp::object const &obj, unsigned char depth = 0) {
  static unsigned char constexpr max_depth(10);
  bp::object ret;
  // limit recursion depth, to avoid infinity loops or possible segfaults
  if (depth >= max_depth) {
    ret = bp::str("...");
  }
  // if the object can be json-ified already, return this object
  else if (isJsonAtomic(obj)) {
    ret = obj;
  }
  // if the object is a list, json-ify each element of the list
  else if (PyList_Check(obj.ptr()) || PyTuple_Check(obj.ptr())) {
    bp::list ls;
    for (bp::ssize_t i = 0; i < bp::len(obj); i++) {
      ls.append(recursiveDictDump(obj[i], depth + 1));
    }
    ret = ls;
  }
  // if the object is a dictionary, json-ify all of the values
  else if (PyDict_Check(obj.ptr())) {
    bp::dict d;
    bp::list keyvals = bp::extract<bp::dict>(obj)().items();
    for (bp::ssize_t i = 0; i < bp::len(keyvals); i++) {
      bp::object key = keyvals[i][0];
      bp::object val = keyvals[i][1];
      d[key] = recursiveDictDump(val, depth + 1);
    }
    ret = d;
  }
  // if the object is not one of the above types, then extract its __dict__ method and repeat the above
  else if (PyObject_HasAttrString(obj.ptr(), "__dict__")) {
    bp::dict d = bp::extract<bp::dict>(obj.attr("__dict__"));
    ret = recursiveDictDump(d, depth); // NOTE re-run at same depth
  } else {
    ret = bp::str(obj);
  }
  return ret;
}
} // namespace

namespace Mantid::Kernel {

// Constructor specializations
template <>
MANTID_PYTHONINTERFACE_CORE_DLL
PropertyWithValue<boost::python::object>::PropertyWithValue(std::string name, boost::python::object defaultValue,
                                                            IValidator_sptr validator, unsigned int direction)
    : Property(std::move(name), typeid(boost::python::object), direction), m_value(defaultValue),
      m_initialValue(std::move(defaultValue)), m_validator(std::move(validator)) {}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL PropertyWithValue<boost::python::object>::PropertyWithValue(
    const std::string &name, const boost::python::object &defaultValue, const std::string &defaultValueStr,
    IValidator_sptr validator, unsigned int direction)
    : Property(name, typeid(boost::python::object), direction), m_value(defaultValue), m_initialValue(defaultValue),
      m_validator(std::move(validator)) {
  // Use the defaultValue parameter instead of trying to parse defaultValueStr
  // since boost::python::object can't use extractToValueVector
  UNUSED_ARG(defaultValueStr);
}

} // namespace Mantid::Kernel

namespace Mantid::PythonInterface {

// PythonObjectProperty constructor implementations
PythonObjectProperty::PythonObjectProperty(std::string const &name, PythonObject const &defaultValue,
                                           IValidator_sptr const &validator, unsigned int const direction)
    : BaseClass(name, defaultValue, validator, direction) {}

PythonObjectProperty::PythonObjectProperty(std::string const &name, PythonObject const &defaultValue,
                                           unsigned int const direction)
    : BaseClass(name, defaultValue, std::make_shared<NullValidator>(), direction) {}

PythonObjectProperty::PythonObjectProperty(std::string const &name, IValidator_sptr const &validator,
                                           unsigned int const direction)
    : BaseClass(name, PythonObject(), validator, direction) {}

PythonObjectProperty::PythonObjectProperty(std::string const &name, unsigned int const direction)
    : BaseClass(name, PythonObject(), std::make_shared<NullValidator>(), direction) {}

PythonObjectProperty::PythonObjectProperty(std::string const &name, std::string const &strvalue,
                                           IValidator_sptr const &validator, unsigned int const direction)
    : BaseClass(name, PythonObject(), strvalue, validator, direction) {}

} // namespace Mantid::PythonInterface

namespace Mantid::Kernel {

/**
 * Creates a string representation of the object
 */
template <> std::string toString(PythonObject const &obj) {
  Mantid::PythonInterface::GlobalInterpreterLock gil;

  // std::string ret;
  boost::python::object rep;
  // if the object is None type, return an empty string
  if (Mantid::PythonInterface::isNone(obj)) {
    rep = boost::python::str("");
  }
  // if the object can be read as a string, then return the object itself
  else if (boost::python::extract<std::string>(obj).check()) {
    rep = obj;
  }
  // otherwise, use json to return a string representation of the class
  else {
    // try loading as a json -- will work for most 'built-in' types
    boost::python::object json = boost::python::import("json");
    try {
      rep = json.attr("dumps")(obj);
    }
    // if json doesn't work, then build a json-like dictionary representation of the object and dump that
    catch (boost::python::error_already_set const &) {
      PyErr_Clear(); // NOTE must clear error registry, or bizarre errors will occur at unexpected lines
      boost::python::object dict = recursiveDictDump(obj);
      try {
        rep = json.attr("dumps")(dict);
      } catch (boost::python::error_already_set const &) {
        PyErr_Clear();
        rep = boost::python::str("<unrepresentable object>");
      }
    }
  }
  return boost::python::extract<std::string>(rep);
}

/**
 * Creates a pretty string representation of the object. In this case it matches
 * the simple string case.
 */
template <> std::string toPrettyString(PythonObject const &value, size_t /*maxLength*/, bool /*collapseLists*/) {
  return toString(value);
}

/**
 * Creates a Json representation of the object
 */
template <> Json::Value encodeAsJson(PythonObject const &) {
  throw Exception::NotImplementedError("encodeAsJson(const boost::python::object &)");
}

} // namespace Mantid::Kernel

namespace Mantid::PythonInterface {

using Kernel::PropertyWithValue;
using Kernel::Exception::NotImplementedError;

/**
 *  @return A string representati on of the default value
 */
std::string PythonObjectProperty::getDefault() const { return Mantid::Kernel::toString(m_initialValue); }

/** Set the property value.
 *  @param value :: The object definition string.
 *  @return Error message, or "" on success.
 */
std::string PythonObjectProperty::setValue(PythonObject const &value) {
  std::string ret;
  try {
    *this = value;
  } catch (std::invalid_argument const &except) {
    ret = except.what();
  }
  return ret;
}

/** Set the property value.
 *  @param value :: The value of the property as a string.
 *  @return Error message, or "" on success.
 */
std::string PythonObjectProperty::setValue(std::string const &value) {
  Mantid::PythonInterface::GlobalInterpreterLock gil;

  std::string ret;
  boost::python::object newVal;
  // try to load as a json object
  try {
    boost::python::object json = boost::python::import("json");
    newVal = json.attr("loads")(value);
  }
  // if it cannot be loaded as a json then it is probably a string
  catch (boost::python::error_already_set const &) {
    PyErr_Clear(); // NOTE must clear error registry, or bizarre errors will occur
    try {
      newVal = boost::python::str(value.c_str());
    } catch (boost::python::error_already_set const &) {
      PyErr_Clear(); // NOTE must clear error registry, or bizarre errors will occur
      return "Failed to interpret string as JSON or string property: " + value;
    }
  }

  // use the assignment operator, which also calls the validator
  try {
    *this = newVal;
  } catch (std::invalid_argument const &except) {
    ret = except.what();
  }
  return ret;
}

/**
 * Assumes the Json object is a string and parses it to create the object
 * @param value A Json::Value containing a string
 * @return An empty string indicating success otherwise the string will contain
 * the value of the error.
 */
std::string PythonObjectProperty::setValueFromJson(const Json::Value &json) {
  std::string jsonstr = JsonHelpers::jsonToString(json);
  return setValue(jsonstr);
}

std::string PythonObjectProperty::setDataItem(const std::shared_ptr<Kernel::DataItem> &) {
  throw NotImplementedError("PythonObjectProperty::setDataItem(const std::shared_ptr<Kernel::DataItem> &)");
}

/** Indicates if the value matches the value None
 *  @return true if the value is None
 *  NOTE: For reasons (surely good ones), trying to compare m_value to m_initialValue raises a segfault
 */
bool PythonObjectProperty::isDefault() const { return m_value.is_none(); }

} // namespace Mantid::PythonInterface

// Explicit specializations for PropertyWithValue<boost::python::object>
// We can't use the full template instantiation because it uses lexical_cast
// which doesn't work with boost::python::object. Instead, we provide explicit
// specializations for the virtual methods that boost::python needs.
namespace Mantid::Kernel {

template <>
MANTID_PYTHONINTERFACE_CORE_DLL PropertyWithValue<boost::python::object> *
PropertyWithValue<boost::python::object>::clone() const {
  return new PropertyWithValue<boost::python::object>(*this);
}

template <> MANTID_PYTHONINTERFACE_CORE_DLL std::string PropertyWithValue<boost::python::object>::isValid() const {
  return m_validator->isValid(m_value);
}

template <> MANTID_PYTHONINTERFACE_CORE_DLL bool PropertyWithValue<boost::python::object>::isDefault() const {
  // For python objects, we consider it default if it's None
  return m_value.is_none();
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL void PropertyWithValue<boost::python::object>::saveProperty(Nexus::File * /*file*/) {
  throw std::invalid_argument("PropertyWithValue::saveProperty - Cannot save '" + this->name() +
                              "', property type boost::python::object not implemented.");
}

template <> MANTID_PYTHONINTERFACE_CORE_DLL std::string PropertyWithValue<boost::python::object>::value() const {
  return toString(m_value);
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL std::string
PropertyWithValue<boost::python::object>::valueAsPrettyStr(const size_t maxLength, const bool collapseLists) const {
  return toPrettyString(m_value, maxLength, collapseLists);
}

template <> MANTID_PYTHONINTERFACE_CORE_DLL Json::Value PropertyWithValue<boost::python::object>::valueAsJson() const {
  return encodeAsJson((*this)());
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL std::string
PropertyWithValue<boost::python::object>::setValue(const std::string &value) {
  // This shouldn't be called directly - PythonObjectProperty overrides it
  // But we need to provide it for vtable completeness
  (void)value; // Suppress unused parameter warning
  throw std::runtime_error("PropertyWithValue<boost::python::object>::setValue should not be called directly. "
                           "Use PythonObjectProperty::setValue instead.");
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL std::string
PropertyWithValue<boost::python::object>::setValueFromJson(const Json::Value &value) {
  // This shouldn't be called directly - PythonObjectProperty overrides it
  // But we need to provide it for vtable completeness
  (void)value; // Suppress unused parameter warning
  throw std::runtime_error("PropertyWithValue<boost::python::object>::setValueFromJson should not be called directly. "
                           "Use PythonObjectProperty::setValueFromJson instead.");
}
template <>
MANTID_PYTHONINTERFACE_CORE_DLL std::string
PropertyWithValue<boost::python::object>::setValueFromProperty(const Property &right) {
  // Try to cast to the same type and copy the value
  if (auto prop = dynamic_cast<const PropertyWithValue<boost::python::object> *>(&right)) {
    m_value = prop->m_value;
    return "";
  } else {
    // Fall back to string conversion
    return setValue(right.value());
  }
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL std::string
PropertyWithValue<boost::python::object>::setDataItem(const std::shared_ptr<DataItem> & /*data*/) {
  throw Exception::NotImplementedError(
      "PropertyWithValue<boost::python::object>::setDataItem - Not implemented for python objects.");
}

template <> MANTID_PYTHONINTERFACE_CORE_DLL std::string PropertyWithValue<boost::python::object>::getDefault() const {
  return toString(m_initialValue);
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL bool PropertyWithValue<boost::python::object>::isMultipleSelectionAllowed() {
  return m_validator->isMultipleSelectionAllowed();
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL std::vector<std::string>
PropertyWithValue<boost::python::object>::allowedValues() const {
  return determineAllowedValues(m_value, *m_validator);
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL PropertyWithValue<boost::python::object> &
PropertyWithValue<boost::python::object>::operator+=(Property const *right) {
  // Adding python objects doesn't make much sense, but we provide this for vtable completeness
  // Just ignore it
  (void)right;
  return *this;
}

template <> MANTID_PYTHONINTERFACE_CORE_DLL int PropertyWithValue<boost::python::object>::size() const {
  // Python object is a single value
  return 1;
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL PropertyWithValue<boost::python::object> &
PropertyWithValue<boost::python::object>::operator=(const boost::python::object &value) {
  boost::python::object oldValue = m_value;
  m_value = value;
  std::string problem = this->isValid();
  if (problem.empty()) {
    return *this;
  } else {
    m_value = oldValue;
    throw std::invalid_argument("When setting value of property \"" + this->name() + "\": " + problem);
  }
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL const boost::python::object &
PropertyWithValue<boost::python::object>::operator()() const {
  return m_value;
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL PropertyWithValue<boost::python::object>::operator const boost::python::object &()
    const {
  return m_value;
}

template <>
MANTID_PYTHONINTERFACE_CORE_DLL void
PropertyWithValue<boost::python::object>::replaceValidator(IValidator_sptr newValidator) {
  m_validator = newValidator;
}

// These template specializations are already provided earlier in this file:
// - toString<PythonObject>
// - toPrettyString<PythonObject>
// - encodeAsJson<PythonObject>

} // namespace Mantid::Kernel
