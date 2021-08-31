#include <ql/quantlib.hpp>

#ifdef BOOST_MSVC
#   include <ql/auto_link.hpp>
#endif

#include "ReplicationError.hpp"
#include <iostream>
#include <iomanip>

using namespace QuantLib;
using namespace std;

#if defined(QL_ENBALE_SESSIONS)
namespace QuantLib
{
    Threadey sessionId() { return {}; }
}
#endif

Real normalCDF(Real value);
Real normalCDF(Real value)
{
    return 0.5*erfc(-value * M_SQRT1_2);
}

int main(int argc, const char * argv[])
{
    try
    {
        std::cout << std::endl;
        Time maturity = 3.0 / 12.0;
        Real strike = 100.0;
        Real underlying = 100.0;
        Volatility volatility = 0.20;
        Rate riskFreeRate = 0.02;
        Rate dividendRate = 0.01;
        Real transactionCost = 0.0;
        Size TimeStep = 21;


        ReplicationError rp(Option::Call, maturity, strike, underlying, volatility, riskFreeRate, dividendRate, transactionCost);
        rp.nSamples = 100;
        rp.nTimeSteps = TimeStep;
        rp.PLMean = 0;
        rp.PLStddev = 0;
        rp.PLSkew = 0;
        rp.PLKurt = 0;
        rp.theorStD = 0;

        std::cout << std::endl;
        std::cout << std::endl;

        std::cout << " Hedge : " << rp.nTimeSteps << "hedging number " << std::endl;
        rp.compute(rp.nTimeSteps, rp.nSamples);
        rp.printResult();
        //Real OptionRealVolValue;
        //OptionRealVolValue = rp.OptionValue_;

        std::cout << std::endl;
        std::cout << std::endl;
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
        std::cout << "*********** Hedging,Real Vol < Input, Implied Vol *************" << std::endl;
        rp.HedgeVol_ = 0.15; // for hedging
        rp.nTimeSteps = 84;

        std::cout << " Input Volatilty " << volatility << std::endl;
        std::cout << " Hedging Volatility " << rp.HedgeVol_ << std::endl;
        std::cout << " Hedge : " << rp.nTimeSteps << "hedging number " << std::endl;
        rp.compute(rp.nTimeSteps, rp.nSamples);
        rp.printResult();
        std::cout << std::endl;



        std::cout << "*********** Mismatch of Volatility *************" << std::endl;
        std::cout << "*********** Hedging,Real Vol > Input, Implied Vol *************" << std::endl;
        rp.HedgeVol_ = 0.25;
        rp.nTimeSteps = 84;

        std::cout << " Input Volatilty " << volatility << std::endl;
        std::cout << " Hedging Volatility " << rp.HedgeVol_ << std::endl;
        std::cout << " Hedge : " << rp.nTimeSteps << "hedging number " << std::endl;
        rp.compute(rp.nTimeSteps, rp.nSamples);
        rp.printResult();
        std::cout << std::endl;

        std::cout << std::endl;


        std::cout << "******* Mismatch of Volatility & Drift Term *******" << std::endl;
        std::cout << "******* drift : u > r *******" << std::endl;
        std::cout << "******* Vol : Hedging,Real Vol > Input, Implied Vol  *******" << std::endl;
        rp.u_ = 0.05;
        rp.HedgeVol_ = 0.25;
        rp.nTimeSteps = 84;

        std::cout << " Pricing Volatilty Small" << volatility << std::endl;
        std::cout << " Hedging Volatility Big!" << rp.HedgeVol_ << std::endl;
        std::cout << " Hedge : " << rp.nTimeSteps << "hedging number " << std::endl;
        rp.compute(rp.nTimeSteps, rp.nSamples);
        rp.printResult();

        std::cout << "******* Therefore, Hedge Vol bigger, PnL Bigger  *******" << std::endl;
        std::cout << "******* Therefore, drift bigger, PnL Loser  *******" << std::endl;
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << std::endl;

        std::cout << std::endl;


        std::cout << "***** Optimal Hedging with Transaction cost *****" << std::endl;

        ReplicationError rpOpt(Option::Call, maturity, strike, underlying, volatility, riskFreeRate, dividendRate, transactionCost);
        rpOpt.nSamples = 100;
        rpOpt.nTimeSteps = TimeStep;
        //rpOpt.PLMean = 0;
        //rpOpt.PLStddev = 0;
        //rpOpt.PLSkew = 0;
        //rpOpt.PLKurt = 0;
        //rpOpt.theorStD = 0;
        rpOpt.optimalHedging(100);
        std::cout << std::endl;
        rpOpt.nTimeSteps = rpOpt.optimalTimes;
        rpOpt.compute(rpOpt.optimalTimes, rpOpt.nSamples);
        std::cout << "Optimal Hedge : " << rpOpt.optimalTimes << "Optimal hedging number " << std::endl;
        rpOpt.printResult();

        std::cout << std::endl;
        std::cout << std::endl;


        std::cout << std::endl;
        std::cout << std::endl;

        std::cout << std::endl;



        std::cout << "***** Hedging with Underestimated Real Volatility *****" << std::endl;
        Size nSamples = 20;
        Size nTimeSteps = 1000;
        Volatility impVol = 0.2;
        Volatility HedgeVol = 0.3;
        ReplicationError RealVolHedge(Option::Call, maturity, strike, underlying, impVol, riskFreeRate, 0, 0);
        RealVolHedge.u_ = 0.05;
        RealVolHedge.HedgeVol_ = HedgeVol;
        RealVolHedge.nSamples = 20;
        RealVolHedge.nTimeSteps = 1000;
        RealVolHedge.maturity_ = 1.0;
        cout << RealVolHedge.nSamples << endl;
        RealVolHedge.computePnL(nTimeSteps, nSamples);
        std::cout << std::endl;
        RealVolHedge.compute(nTimeSteps, nSamples);
        RealVolHedge.printResult();
        std::cout << std::endl;


        /*

       rp.nTimeSteps = 1000;
       rp.nSamples = 20;

       Volatility HedgeVol; // for hedigng
       Volatility InputVol;  // for pricing
       HedgeVol = 0.3;
       InputVol = 0.2;

       ReplicationError rpInput(Option::Call, maturity, strike, underlying, HedgeVol, riskFreeRate, dividendRate, 0);
       rpInput.nTimeSteps = rp.nTimeSteps;
       rpInput.nSamples = rp.nSamples;
       rpInput.HedgeVol_ = HedgeVol;
       rpInput.compute(rp.nTimeSteps, rp.nSamples, rp.nTimeSteps);

       ReplicationError rpReal(Option::Call, maturity, strike, underlying, InputVol, riskFreeRate, dividendRate, 0);
       rpReal.nTimeSteps = rp.nTimeSteps;
       rpReal.nSamples = rp.nSamples;
       rpReal.HedgeVol_ = InputVol;
       rpReal.compute(rp.nTimeSteps, rp.nSamples, rp.nTimeSteps);

       Real OptionHedgeVolValue = rpInput.OptionValue_;
       Real OptionRealVolValue = rpReal.OptionValue_;


       ReplicationError rpMismatch(Option::Call, maturity, strike, underlying, InputVol, riskFreeRate, dividendRate, 0);
       rpMismatch.nTimeSteps = rp.nTimeSteps;
       rpMismatch.nSamples = rp.nSamples;
       rpMismatch.HedgeVol_ = HedgeVol; // for hedging
       rpMismatch.compute(rp.nTimeSteps, rp.nSamples, rp.nTimeSteps);

       vector < vector<Real> > TimePaths(rpMismatch.nTimeSteps, vector<Real>(rpMismatch.nSamples, 0));
       vector < Real> Lowerbounds;
       Real Upperbound = OptionHedgeVolValue - OptionRealVolValue;
       Real Lowerbound;
       Lowerbounds.clear();

       for(Size i = 1.0; i<= rp.nTimeSteps;i++)
       {
           rpMismatch.HedgeVol_ = HedgeVol;
           Time n = rp.nTimeSteps;
           Time m = i;
           rpMismatch.maturity_ = m / n * maturity;

           //rp.nTimeSteps = i;
           //rp.compute(i, rp.nSamples);
           rpMismatch.compute(rp.nTimeSteps, rp.nSamples, rp.nTimeSteps);

           for(Size j = 0; j< rp.nSamples;j++)
           {
               TimePaths[i-1][j] = rpMismatch.tempPaths[j];
           }
           Lowerbound = Upperbound - 2 * strike * std::exp(-(riskFreeRate-0) *(rpMismatch.maturity_)) *
               (normalCDF(0.5*(HedgeVol - InputVol)*std::sqrt(rpMismatch.maturity_)) - 0.5);
           Lowerbounds.push_back(Lowerbound);
       }

       for(Size i = 0; i< TimePaths.size();i++)
       {
           for(Size j = 0; j < rp.nSamples;j++)
           {
               std::cout << "i = " << i + 1 << " j = " << j + 1 << ", Result : " << TimePaths[i][j]
                   << ", Upper Bound : " << Upperbound << ", Lower Bound : " << Lowerbounds[i] << std::endl;
               std::cout << std::endl;
           }
       }
       */

        return 0;
    }
    catch (std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "Unknown Error" << std::endl;
        return 1;
    }
}
