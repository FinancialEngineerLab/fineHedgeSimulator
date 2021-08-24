
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
        Rate riskFreeRate = 0.02;
        Rate dividendRate = 0.01;
        Real transactionCost = 0.00015;
        Size TimeStep = 21;

		
        ReplicationError rp(Option::Call, maturity, strike, underlying, volatility, riskFreeRate, dividendRate, transactionCost);
        
		rp.nSamples = 2000;
		rp.nTimeSteps = TimeStep;
		rp.PLMean = 0;
		rp.PLStddev = 0;
		rp.PLSkew = 0;
		rp.PLKurt = 0;
		rp.theorStD = 0;

		std::cout << std::endl;

        std::cout << std::endl;

        std::cout << " Hedge : " << rp.nTimeSteps << "hedging number " <<std::endl;
        rp.compute(rp.nTimeSteps, rp.nSamples);
		rp.printResult();

		std::cout << std::endl;


        rp.nTimeSteps = 84;
		std::cout << " Hedge : " << rp.nTimeSteps << "hedging number " << std::endl;
		rp.compute(rp.nTimeSteps, rp.nSamples);
		rp.printResult();

		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "*********** Mismatch of Drift *************" << std::endl;
		rp.nTimeSteps = 84;
		rp.u_ = 0.04;
		std::cout << " Hedge : " << rp.nTimeSteps << "hedging number " << std::endl;
		rp.compute(rp.nTimeSteps, rp.nSamples);
		rp.printResult();
		std::cout << std::endl;


		std::cout << std::endl;


		std::cout << "*********** Mismatch of Volatility *************" << std::endl;
		rp.simulSigma_ = 0.18;
		rp.nTimeSteps = 84;
        std::cout << " Realized Volatilty " <<  volatility << std::endl;
        std::cout << " Inputted Volatility " << rp.simulSigma_ << std::endl;
		std::cout << " Hedge : " << rp.nTimeSteps << "hedging number " << std::endl;
		rp.compute(rp.nTimeSteps, rp.nSamples);
		rp.printResult();
		std::cout << std::endl;


		std::cout << std::endl;


		std::cout << "******* Mismatch of Volatility & Drift Term *******" << std::endl;
		rp.u_ = 0.01;
		rp.simulSigma_ = 0.1;
		rp.nTimeSteps = 84;
		std::cout << " Hedge : " << rp.nTimeSteps << "hedging number " << std::endl;
		rp.compute(rp.nTimeSteps, rp.nSamples);
		rp.printResult();

		std::cout << std::endl;

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "***** Optimal Hedging with Transaction cost *****" << std::endl;
		rp.optimalHedging(100);
		std::cout << std::endl;
        rp.nTimeSteps = rp.optimalTimes;
        rp.compute(rp.optimalTimes, rp.nSamples);
        std::cout << "Optimal Hedge : " << rp.optimalTimes << "Optimal hedging number " << std::endl;
        rp.printResult();
        std::cout << std::endl;
        std::cout << std::endl;
        
        
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "***** Hedging with Real Volatility *****" << std::endl;
        
        rp.nTimeSteps = 42;
        vector < vector<Real> > TimePaths(rp.nTimeSteps, vector<Real>(rp.nSamples, 0));
        
        for(Size i = 1.0; i<=rp.nTimeSteps;i++)
        {
            Time n = rp.nTimeSteps;
            Time m = i;
            rp.maturity_ = m / n;
            rp.nTimeSteps = i;
            rp.compute(rp.nTimeSteps, rp.nSamples);
            
            for(Size j = 0; j<rp.nSamples;j++)
            {
                TimePaths[i-1][j] = rp.tempPaths[j];
            }
        }
        
        for(Size i = 0; i< TimePaths.size();i++)
        {
            for(Size j = 0; j < rp.nSamples;j++)
            {
                std::cout << "i = "<<i <<" j = " <<j << ", Result : " << TimePaths[i][j];
                std::cout << std::endl;
            }
        }
        
        
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
