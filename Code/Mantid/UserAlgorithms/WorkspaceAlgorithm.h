#ifndef USER_ALGORITHMS_WORKSPACEALGORITHM_H_
#define USER_ALGORITHMS_WORKSPACEALGORITHM_H_

#include <iostream>

#include "MantidAPI/Algorithm.h"

namespace Mantid
{
namespace Algorithms
{
	/** @class WorkspaceAlgorithm WorkspaceAlgorithm.h
   	
 	    An example algorithm iterating over a workspace.
   			    	
      @author Roman Tolchenov, ISIS, RAL
      @date 02/05/2008
   	    
      Copyright &copy; 2008 STFC Rutherford Appleton Laboratories
   	
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
class WorkspaceAlgorithm : public API::Algorithm
{
public:
  ///no arg constructor
  WorkspaceAlgorithm() : API::Algorithm() {}
  ///virtual destructor
  virtual ~WorkspaceAlgorithm() {}
  /// Algorithm's name for identification overriding a virtual method
  virtual const std::string name() const { return "WorkspaceAlgorithm";}
  /// Algorithm's version for identification overriding a virtual method
  virtual const int version()const { return (1);}

  ///Initialisation code
  void init();
  ///Execution code
  void exec();

  /// Static reference to the logger class
  static Mantid::Kernel::Logger& g_log;

};

} // namespace Algorithm
} // namespace Mantid

#endif /*USER_ALGORITHMS_WORKSPACEALGORITHM_H_*/

