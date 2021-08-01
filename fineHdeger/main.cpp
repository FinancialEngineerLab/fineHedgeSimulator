
#include <ql/quantlib.hpp>

#ifdef BOOST_MSVC
#   include <ql/auto_link.hpp>
#endif

#include <iostream>
#include <iomanip>

using namespace QuantLib;

#if defined(QL_ENBALE_SESSIONS)
namespace QuantLib
{
    Threadey sessionId() { return {};}
}
#endif



int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}
