#ifndef MANTID_API_ICONSTRAINT_H_
#define MANTID_API_ICONSTRAINT_H_

//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAPI/IFunction.h"

namespace Mantid
{
namespace API
{
/** An interface to a constraint.

    @author Anders Markvardsen, ISIS, RAL
    @date 12/11/2009

    Copyright &copy; 2009 STFC Rutherford Appleton Laboratory

    This file is part of Mantid.

    Mantid is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Mantid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    File change history is stored at: <https://svn.mantidproject.org/mantid/trunk/Code/Mantid>.
    Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
class DLLExport IConstraint
{
public:
  /// Default constructor
  IConstraint() {}
  /// Virtual destructor
  virtual ~IConstraint() {}

  /// Returns a penalty number which is bigger than or equal to zero
  /// If zero it means that the constraint is not penalized. If larger
  /// than zero the constraint is penalized where the larger this number
  /// is the larger the penalty
  virtual double check(IFunction* fn) = 0;

  /// Returns the derivative of the penalty for each parameter
  virtual boost::shared_ptr<std::vector<double> > checkDeriv(IFunction* fn) = 0;

};


} // namespace API
} // namespace Mantid


#endif /*MANTID_API_ICONSTRAINT_H_*/
