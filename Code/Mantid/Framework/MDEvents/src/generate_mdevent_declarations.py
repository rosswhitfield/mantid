""" Simple script that generates references to all
needed MDLeanEvent<X> instantiations. """
import sys
import os
import time
import datetime

# List of the dimensions to instantiate
dimensions = range(1,10)

# List of every possible MDEvent or MDLeanEvent types.
mdevent_types = ["MDEvent", "MDLeanEvent"]

topheader = """/* Auto-generated by '%s' 
 *     on %s
 *
 * DO NOT EDIT!
 */ 
 
#include <boost/shared_ptr.hpp>
#include "MantidMDEvents/MDEventFactory.h"
#include "MantidAPI/IMDEventWorkspace.h"
#include "MantidMDEvents/MDEventWorkspace.h"
 
""" % (__file__, datetime.datetime.now())
 
header = """
namespace Mantid
{
namespace MDEvents
{
""" 


footer = """
} // namespace Mantid
} // namespace MDEvents 

/* THIS FILE WAS AUTO-GENERATED BY %s - DO NOT EDIT! */ 
""" % (__file__)




#============================================================================================================
#============================================================================================================
#============================================================================================================
#============================================================================================================
factory_top = """

/** Create a MDEventWorkspace of the given type
@param nd :: number of dimensions
@param eventType :: string describing the event type (MDEvent or MDLeanEvent)
@return shared pointer to the MDEventWorkspace created (as a IMDEventWorkspace).
*/
API::IMDEventWorkspace_sptr MDEventFactory::CreateMDWorkspace(size_t nd, const std::string & eventType)
{

"""
factory_if_event_type = """  if (eventType == "%s")
  {
    switch(nd)
    {
"""

factory_lines = """    case (%d):
      return boost::shared_ptr<MDEventWorkspace<%s,%d> >(new MDEventWorkspace<%s,%d>);
"""
factory_bottom = """    default:
      throw std::invalid_argument("Invalid number of dimensions passed to CreateMDWorkspace.");
    } // end switch
  } // end if eventType
"""

def write_factory(f):
    """ Write out a factory method """
    f.write(factory_top)
    for mdevent_type in mdevent_types:
        f.write(factory_if_event_type % mdevent_type)
        for nd in dimensions:
            eventType = "%s<%d>" % (mdevent_type, nd)
            f.write(factory_lines % (nd,eventType,nd,eventType,nd) )
        f.write(factory_bottom)
    # Throw for an invalid event type
    f.write('  // Unknown event type\n')
    f.write('  throw std::invalid_argument("Unknown event type passed to CreateMDWorkspace.");\n')
    f.write("}\n")
        


        
#============================================================================================================
#============================================================================================================
#============================================================================================================
#============================================================================================================

# for the calling function macro
macro_top = """
/** Macro that makes it possible to call a templated method for
 * a MDEventWorkspace using a IMDEventWorkspace_sptr as the input.
 *
 * @param funcname :: name of the function that will be called.
 * @param workspace :: IMDEventWorkspace_sptr input workspace.
 */
 
#define %sCALL_MDEVENT_FUNCTION%s(funcname, workspace) \\
{ \\
"""


macro = """%sMDEventWorkspace<%s, %d>::sptr %s = boost::dynamic_pointer_cast<%sMDEventWorkspace<%s, %d> >(workspace); \\
if (%s) funcname<%s, %d>(%s); \\
"""

def get_macro(min_dimension=1, const=""):
    """ Return the macro code CALL_MDEVENT_FUNCTION
    Parameter:
        min_dimension :: to avoid compiler warnings, limit to dimensions higher than this
        const :: set to "const " to make a const equivalent
    """
    suffix = ""
    prefix = ""
    if (min_dimension > 1): suffix = "%d" % min_dimension
    if const != "": prefix = "CONST_"
    s = macro_top % (prefix, suffix);
    
    for mdevent_type in mdevent_types:
        for nd in dimensions:
            if (nd >= min_dimension):
                eventType = "%s<%d>" % (mdevent_type, nd)
                varname = "MDEW_%s_%d" % (mdevent_type.upper(),nd)
                if const != "":
                    varname = "CONST_" + varname
                s += macro % (const, eventType,nd, varname, const, eventType,nd, varname, eventType,nd, varname) 
    s +=  "} \n"
    return s



#======================================================================================
def get_padding(line):
    """Return a string with the spaces padding the start of the given line."""
    out = ""
    for c in line:
        if c == " ":
            out += " "
        else:
            break
    return out

#======================================================================================
def find_line_number(lines, searchfor, startat=0):
    """Look line-by-line in lines[] for a line that starts with searchfor. Return
    the line number in source where the line was found, and the padding (in spaces) before it"""
    count = 0
    done = False
    padding = ""
    for n in xrange(startat, len(lines)):
        line = lines[n]
        s = line.strip()
        if s.startswith(searchfor):
            # How much padding?
            padding = get_padding(line)
            return (n, padding)
    return (None, None)


#======================================================================================
def insert_lines(lines, insert_lineno, extra, padding):
    """Insert a text, split by lines, inside a list of 'lines', at index 'insert_lineno'
    Adds 'padding' to each line."""
    # split
    extra_lines = extra.split("\n");
    #Pad 
    for n in xrange(len(extra_lines)):
        extra_lines[n] = padding+extra_lines[n]
    return lines[0:insert_lineno] + extra_lines + lines[insert_lineno:]

#============================================================================================================
#============================================================================================================
#============================================================================================================

def generate():
    print "Generating MDEventFactory.cpp"

    # Classes that have a .cpp file (and will get an Include line)
    classes_cpp = ["IMDBox", "MDBox", "MDEventWorkspace", "MDGridBox", "MDBin", "MDBoxIterator"]
    # All of the classes to instantiate
    classes = classes_cpp + mdevent_types
    
    # First, open the header and read all the lines
    f = open("../inc/MantidMDEvents/MDEventFactory.h", 'r')
    s = f.read();
    lines = s.split("\n")
    (n1, padding) = find_line_number(lines, "//### BEGIN AUTO-GENERATED CODE ###", startat=0)
    (n2, padding_ignored) = find_line_number(lines, "//### END AUTO-GENERATED CODE ###", startat=0)
    print n1
    print n2
    if n1 is None or n2 is None:
        raise Exception("Could not find the marker in the MDEventFactory.h file.")
    print n2
    header_before = lines[:n1+1]
    header_after = lines[n2:]
    f.close()
    
    # =========== Do the Source File ===========
    f = open("MDEventFactory.cpp", 'w');
    f.write(topheader)
    
    for c in classes:
        f.write('#include "MantidMDEvents/%s.h"\n' % c)
    
    f.write("\n")
    f.write("// We need to include the .cpp files so that the declarations are picked up correctly. Weird, I know. \n") 
    f.write("// See http://www.parashift.com/c++-faq-lite/templates.html#faq-35.13 \n") 
    for c in classes_cpp:
        f.write('#include "%s.cpp"\n' % c)
    f.write("\n")
    
    f.write(header)

    # MDEvent and MDLeanEvent type (just one template arg)
    classes = mdevent_types
    for c in classes:
        f.write("// Instantiations for %s\n" % c )
        for nd in dimensions:
            f.write("template DLLExport class %s<%d>;\n" % (c, nd) )
        f.write("\n\n")

    # Classes with MDLeanEvent<x>,x
    classes = classes_cpp
    for c in classes:
        f.write("// Instantiations for %s\n" % c )
        for mdevent_type in mdevent_types:
            for nd in dimensions:
                f.write("template DLLExport class %s<%s<%d>, %d>;\n" % (c, mdevent_type, nd, nd) )
        f.write("\n\n")
            
    write_factory(f)
    
    f.write(footer)
    f.close()
    
    
    # ========== Start the header file =============
    lines = header_before

    # Make the macro then pad it into the list of lines
    macro = get_macro()
    lines = insert_lines(lines, len(lines), macro, padding)
    # Again for 3+ dimensions 
    macro = get_macro(3)
    lines = insert_lines(lines, len(lines), macro, padding)
    # Again for const workspace
    macro = get_macro(1, "const ")
    lines = insert_lines(lines, len(lines), macro, padding)

    # Typedefs for MDEventWorkspace
    lines.append("\n");
    classes = ["MDBox", "IMDBox", "MDGridBox", "MDEventWorkspace", "MDBin"]
    for c in classes:
        lines.append("\n%s// ------------- Typedefs for %s ------------------\n" % (padding, c));
        mdevent_type = "MDEvent"
        for nd in dimensions:
            lines.append(padding + "/// Typedef for a %s with %d dimension%s " % (c, nd, ['','s'][nd>1]) )
            lines.append(padding + "typedef %s<%s<%d>, %d> %s%d;" % (c, mdevent_type, nd, nd, c, nd) )
        mdevent_type = "MDLeanEvent"
        for nd in dimensions:
            lines.append(padding + "/// Typedef for a %s with %d dimension%s " % (c, nd, ['','s'][nd>1]) )
            lines.append(padding + "typedef %s<%s<%d>, %d> %s%dLean;" % (c, mdevent_type, nd, nd, c, nd) )
            
        lines.append("\n");
            

    lines += header_after
    
    f = open("../inc/MantidMDEvents/MDEventFactory.h", 'w')
    for line in lines:
        f.write(line + "\n")
    f.close()


if __name__=="__main__":
    generate()
