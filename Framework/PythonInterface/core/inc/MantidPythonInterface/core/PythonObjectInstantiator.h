// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2011 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include "MantidKernel/Instantiator.h"
#include "MantidPythonInterface/core/GlobalInterpreterLock.h"

#include <boost/python/extract.hpp>
#include <boost/python/object.hpp>
// older versions of boost::python seem to require this last, after the object
// include
#include <memory>

namespace Mantid {
namespace PythonInterface {
/**
 * Special shared_ptr::deleter object that locks the GIL while deleting
 * the underlying Python object
 */
struct GILSharedPtrDeleter {
  GILSharedPtrDeleter(const boost::python::converter::shared_ptr_deleter &deleter)
      : m_deleter(boost::python::converter::shared_ptr_deleter(deleter.owner)) {}
  /**
   * Called when the shared_ptr reference count is zero
   * @param data A pointer to the data to be deleted
   */
  void operator()(void const *data) {
    GlobalInterpreterLock gil;
    m_deleter(data);
  }
  /// Main deleter object
  boost::python::converter::shared_ptr_deleter m_deleter;
};

/**
 * @tparam Base A pointer to the C++ type that is stored in within the created
 *              python object
 */
template <typename Base> class PythonObjectInstantiator : public Kernel::AbstractInstantiator<Base> {
public:
  /// Constructor taking a Python class object wrapped as a boost::python:object
  PythonObjectInstantiator(const boost::python::object &classObject) : m_classObject(classObject) {}

  /// Creates an instance of the object as shared_ptr to the Base type
  std::shared_ptr<Base> createInstance() const override;

  /// Creates an instance of the object as raw pointer to the Base type
  Base *createUnwrappedInstance() const override;

private:
  /// The class name
  boost::python::object m_classObject;
};

/**
 * Creates an instance of the object as shared_ptr to the Base type
 * @returns A shared_ptr to Base.
 */
template <typename Base> std::shared_ptr<Base> PythonObjectInstantiator<Base>::createInstance() const {
  using namespace boost::python;
  GlobalInterpreterLock gil;

  object instance = m_classObject();
  // The instantiator assumes that the exported type uses a
  // HeldType=std::shared_ptr<Adapter>,
  // where Adapter inherits from Base,
  // see
  // http://www.boost.org/doc/libs/1_42_0/libs/python/doc/v2/class.html#class_-spec.
  // The advantage is that the memory management is very simple as it is all
  // handled within the
  // boost layer.

  // Swap the deleter for one that acquires the GIL before actually deallocating
  // the
  // Python object or we get a segfault when deleting the objects on later
  // versions of Python 2.7
  auto instancePtr = extract<std::shared_ptr<Base>>(instance)();
  auto *deleter = std::get_deleter<converter::shared_ptr_deleter, Base>(instancePtr);
  instancePtr.reset(instancePtr.get(), GILSharedPtrDeleter(*deleter));
  return instancePtr;
}

/**
 * @throws std::runtime_error as we're unable to extract a non-shared ptr from
 * the wrapped object
 */
template <typename Base> Base *PythonObjectInstantiator<Base>::createUnwrappedInstance() const {
  throw std::runtime_error("Unable to create unwrapped instance of Python object");
}
} // namespace PythonInterface
} // namespace Mantid
