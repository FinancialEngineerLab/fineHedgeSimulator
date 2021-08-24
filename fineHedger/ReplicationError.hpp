
#include <ql/quantlib.hpp>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#ifndef ReplicationError_h
#define ReplicationError_h

using namespace QuantLib;
using namespace std;

class ReplicationError
{
public:


    ReplicationError(Option::Type type, Time maturity, Real strike, Real s0,Rate sigma, Rate r, Rate q, Real transactionCost) :
    maturity_(maturity), payoff_(type, strike), s0_(s0), sigma_(sigma), r_(r), q_(q), transactionCost_(transactionCost)
    {
        // for option BSM Price
        DiscountFactor rDiscount = std::exp(-r_*maturity_);
        DiscountFactor qDiscount = std::exp(-q*maturity_);;
        Real forward = s0_*qDiscount/rDiscount;
        Real stdDev = std::sqrt(sigma_*sigma_*maturity_);
        ext::shared_ptr<StrikedTypePayoff> payoff(new PlainVanillaPayoff(payoff_));
        BlackCalculator black(payoff, forward, stdDev, rDiscount);
        std::cout << "Option Value: " << black.value() << std::endl;
        
        //for option BSM Vega
        vega_ = black.vega(maturity_);
        
        std::cout << std::endl;
        std::cout << std::setw(8) << " " << " | "
                  << std::setw(8) << " " << " | "
                  << std::setw(8) << "PnL" << " | "
                  << std::setw(8) << "PnL" << " | "
                  << std::setw(12) << "Derman Kamal" << " | "
                  << std::setw(8) << "PnL" << " | "
                  << std::setw(8) << "PnL" << std::endl;
        std::cout << std::setw(8) << "Samples" << " | "
                  << std::setw(8) << "Trades" << " | "
                  << std::setw(8) << "Mean" << " | "
                  << std::setw(8) << "Stddev" << " | "
                  << std::setw(12) << "formaula" << " | "
                  << std::setw(8) << "Skewness" << " | "
        << std::setw(8) << "Kurtosis" << std::endl;
        
        std::cout << std::string(78, '-') << std::endl;
		initialSimul();
    }
    
    void compute(Size nTimeSteps, Size nSamples);
	void optimalHedging(Size maxDt);
	void printResult();

    
    //vector<pair<Size, Size> getPnL(Time time, Size nTimeSteps, Size nSamples);
    
	void initialSimul()
	{
		u_ = r_;
		simulSigma_ = sigma_;
	}

	void initialHedge()
	{
		PLMean = 0;
		PLStddev = 0;
		PLSkew = 0;
		PLKurt = 0;
		theorStD = 0;
        tempPaths.clear();
        optimalTimes = nTimeSteps;
	}

	Rate u_;
	Rate simulSigma_;

	Size scenarios;
	Size nTimeSteps;
	//static Size hedgesNum;
    Volatility inputVol;
    Time time;
    
	Size nSamples =0;
	Real PLMean;
	Real PLStddev;
	Real PLSkew;
	Real PLKurt;
	Real theorStD;
    Real optimalTimes;
    
    std::vector <pair<Real, Real> > pnlPaths;
    std::vector <Real> tempPaths;
    Matrix pnlPaths2(Size i, Size j);
    
    
    Time maturity_;
private:
    PlainVanillaPayoff payoff_;
    Real s0_;
    Volatility sigma_;
    Rate r_;
    Rate q_;
    Real vega_;
    Real transactionCost_;

};


class ReplicationPathPricer : public PathPricer<Path>
{
public:
    ReplicationPathPricer(Option::Type type, Real strike, Real r, Real q, Time maturity, Volatility sigma, Real transactionCost):
    type_(type), strike_(strike), r_(r), q_(q), maturity_(maturity), sigma_(sigma), transactionCost_(transactionCost)
    {
        QL_REQUIRE(strike_ > 0.0, "strike must be positive");
        QL_REQUIRE(maturity_ > 0.0, "maturity must be positive");
        QL_REQUIRE(sigma_ >= 0.0, "volatility must be equal or larger than zero");
    }
    Real operator()(const Path& path) const override;
    
    
private:
    Option::Type type_;
    Real strike_;
    Rate r_;
    Rate q_;
    Time maturity_;
    Volatility sigma_;
    Real transactionCost_;
    
};

#endif /* ReplicationError_h */
