
#include <ql/quantlib.hpp>

#ifdef BOOST_MSVC
#   include <ql/auto_link.hpp>
#endif

#include "ReplicationError.hpp"

#include <iostream>
#include <iomanip>

using namespace QuantLib;

#if defined(QL_ENBALE_SESSIONS)
namespace QuantLib
{
    Threadey sessionId() { return {};}
}
#endif



int main(int argc, const char * argv[])
{
    try
    {
        std::cout << std::endl;
        Time maturity = 1.0 / 12.0;
        Real strike = 100.0;
        Real underlying = 100.0;
        Volatility volatility = 0.20;
        Rate riskFreeRate = 0.01;
        Rate dividendRate = 0.01;
        ReplicationError rp(Option::Call, maturity, strike, underlying, volatility, riskFreeRate, dividendRate, 0);
        Size scenarios = 20000;
        Size hedgesNum;
        
        hedgesNum = 21;
        std::cout << " Hedge : " << hedgesNum << "hedging number " <<std::endl;
        rp.compute(hedgesNum, scenarios);
        
        hedgesNum = 84;
        std::cout << " Hedge : " << hedgesNum << "hedging number " <<std::endl;
        rp.compute(hedgesNum, scenarios);
        
        hedgesNum = 84;
        
        return 0;
    }
    catch (std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr <<"Unknown Error" << std::endl;
        return 1;
    }
}
