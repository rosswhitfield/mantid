#ifndef MANTID_API_MEMORYMANAGER_H_
#define MANTID_API_MEMORYMANAGER_H_

//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
//#include "boost/shared_ptr.hpp"
#include <vector>
#include <string>
#include "MantidAPI/DllExport.h"
#include "MantidKernel/Logger.h"
#include "MantidKernel/SingletonHolder.h"

namespace Mantid
{
	namespace API
	{

        /// Information about the memory 
        struct MemoryInfo
        {
            int totalMemory;  ///< total physical memory in KB
            int availMemory;  ///< available physical memory in KB
            int freeRatio;    ///< percentage of the available memory ( 0 - 100 )
        };

		/** @class MemoryManagerImpl MemoryManager.h API/MemoryManager.h

		The MemoryManagerImpl class is responsible for memory management.

		@author Roman Tolchenov, Tessella Support Services plc
		@date 01/08/2008

		Copyright &copy; 2007 STFC Rutherford Appleton Laboratories

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

		class EXPORT_OPT_MANTID_API MemoryManagerImpl 
		{
		public:

            ///Returns available physical memory in the system in KB.
            MemoryInfo getMemoryInfo();

		private:
            friend struct Mantid::Kernel::CreateUsingNew<MemoryManagerImpl>;

			///Class cannot be instantiated by normal means
			MemoryManagerImpl();
            ///destructor
			~MemoryManagerImpl();
            ///Copy contrstrutor
			MemoryManagerImpl(const MemoryManagerImpl&);
    
			/// Standard Assignment operator    
			MemoryManagerImpl& operator = (const MemoryManagerImpl&);

			/// Static reference to the logger class
			Kernel::Logger& g_log;

		};

		///Forward declaration of a specialisation of SingletonHolder for AlgorithmManagerImpl (needed for dllexport/dllimport) and a typedef for it.
#ifdef _WIN32
		// this breaks new namespace declaraion rules; need to find a better fix
		template class EXPORT_OPT_MANTID_API Mantid::Kernel::SingletonHolder<MemoryManagerImpl>;
#endif /* _WIN32 */
		typedef Mantid::Kernel::SingletonHolder<MemoryManagerImpl> MemoryManager;

	} // namespace API
}  //Namespace Mantid

#endif /* MANTID_API_MEMORYMANAGER_H_ */
