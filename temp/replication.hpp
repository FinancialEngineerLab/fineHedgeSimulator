
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
        strike_ = payoff_.strike();
        type_ = payoff_.optionType();
        BlackCalculator black(payoff, forward, stdDev, rDiscount);
        std::cout << "Option Value: " << black.value() << std::endl;
        
        initialSimul();
        //for option BSM Vega
        vega_ = black.vega(maturity_);
        OptionValue_ = black.value();

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
        //initialSimul();
    }
    
    //void compute(Size nTimeSteps, Size nSamples);

    void compute(Size nTimeSteps, Size nSamples);
    void optimalHedging(Size maxDt);
    void printResult();


    void computePnL(Size nTimeSteps, Size nSamples);
  
    //vector<pair<Size, Size> getPnL(Time time, Size nTimeSteps, Size nSamples);
    
    void initialSimul()
    {
        u_ = r_;
        HedgeVol_ = sigma_;
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
        //OptionValue_ = 0;
    }

    Rate u_;
    Rate simulSigma_;
    Rate HedgeVol_;

    Size scenarios;
    Size nTimeSteps;
    //static Size hedgesNum;
    Volatility inputVol;
    Time time;
    Real OptionValue_ =0;

    Size nSamples =0;
    Real PLMean;
    Real PLStddev;
    Real PLSkew;
    Real PLKurt;
    Real theorStD;
    Real optimalTimes;
    
    std::vector <pair<Real, Real> > pnlPaths;
    std::vector <Real> tempPaths;
    
    Real getd1(Volatility vol, Real s, Time t)
    {
        Real d1;
        d1 = (std::log(s / strike_) + (r_ - q_ + 0.5* vol * vol) * t) / (vol * std::sqrt(t));
        return d1;
    }
    Real getd2(Volatility vol, Real s, Time t)
    {
        return getd1(vol, s, t) - vol * std::sqrt(t);
    }
    Real stdNomCumDist(Real x, double mu, double sigma)
    {
        double v = (x - mu) / sigma;
        return 0.5*std::erfc(-v * M_SQRT1_2);
    }
    
    Real BSMDelta(Volatility vol, Real s, Time t)
    {
        Real delta(10);
        double d1 = getd1(vol, s, t);
        double nd1 = stdNomCumDist(d1, 0 ,1);
        if (type_ == Option::Put)
        {
            delta = std::exp(-q_ * t)*(nd1+ (-1- 1) / 2);
        }
        else
        {
            delta = std::exp(-q_ * t)*(nd1 + (1 - 1) / 2);
        }
        return delta;
    }

    Real BSMGamma(Volatility vol, Real s, Time t)
    {
        Real gamma;
        Real d1;
        d1 = getd1(vol, s, t);
        Real npd1 = stdNomCumDist(d1, 0, 1);
        gamma = (exp(-d1 * d1 / 2) * (1 / std::sqrt(2*3.1415926))) / (s * vol * std::sqrt(t));
        return gamma;
    }

    Real BSMPrice(Volatility vol, Real s, Time t)
    {
        Real price;
        double d1 = getd1(vol, s, t);
        double d2 = getd2(vol, s, t);
        if (type_ == Option::Put)
        {
            double nd1 = stdNomCumDist(-1 * d1,0,1);
            double nd2 = stdNomCumDist(-1 * d2, 0, 1);
            double price = -1 * (s * std::exp(-q_ * t)*nd1 - strike_ * std::exp(-r_ - t) * nd2);
            return price;
        }
        else
        {
            double nd1 = stdNomCumDist(1 * d1, 0, 1);
            double nd2 = stdNomCumDist(1 * d2, 0, 1);
            double price = 1 * (s * std::exp(-q_ * t)*nd1 - strike_ * std::exp(-r_ - t) * nd2);
            return price;
        }
    }

    //std::vector<vector<Real>> pnlProcess;
    //std::vector<vector<Real>> pnlResult;
    Time maturity_;
private:
    PlainVanillaPayoff payoff_;
    Real s0_;
    Real strike_;
    Option::Type type_;
    Volatility sigma_;
    Rate r_;
    Rate q_;
    Real vega_;
    Real transactionCost_;
    Real gamma_;

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
    Volatility simulSigma_;
    Real transactionCost_;
    
};

#endif /* ReplicationError_h */

