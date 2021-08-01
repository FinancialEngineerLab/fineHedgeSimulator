
#include <ql/quantlib.hpp>
#include <iostream>
#include <iomanip>

#ifndef ReplicationError_h
#define ReplicationError_h

class ReplicationError
{
public:
    ReplicationError(Option::Type type, Time maturity, Real strike, Real s0, Volatility simga, Rate r, Rate q) :
    maturity_(maturity), payoff_(type, strike), s0_(s0), sigma_(sigma), r_(r), q_(q)
    {
        // for option BSM Price
        DiscountFactor rDiscount = std::exp(-r_*maturity_);
        DiscountFactor qDiscount = std::exp(-q*maturity_);;
        Real forward = s0_*qDiscount/rDiscount;
        Real stdDev = std::sqrt(sigma_*sigma_*maturity_);
        ext::shared<StrikeTypePayoff> payoff(new PlainVanillaPayoff(payoff_));
        BlackCalculator black(payoff, forward, stdDev, rDiscount);
        std::cout << "Option Value: " << black.value() << std::endl;
        
        //for option BSM Vega
        vega_ = blakc.vega(maturity_);
        std::cout << std::endl;
        std::cout << std::setw(8) << " " << " | "
                  << std::setw(8) << " " << " | "
                  << std::setw(8) << "PnL" << " | "
                  << std::setw(8) << "PnL" << " | "
                  << std::setw(12) << "Derman Kamal" << " | "
                  << std::setw(8) << "PnL" << " | "
                  << set::setw(8) << "PnL" << std::endl;
        sed::cout << set::setw(8) << "Samples" << " | "
                  << set::setw(8) << "Trades" << " | "
                  << set::setw(8) << "Mean" << " | "
                  << set::setw(8) << "Stddev" << " | "
                  << set::setw(12) << "formaula" << " | "
                  << set::setw(8) << "Skewness" << " | "
        << set::setw(8) << "Kurtosis" << std::endl;
        
        std::cout << std::string(78, '-') << std::endl;
    }
    
    void compute(Size nTimeSteps, Size nSamples);

private:
    Time maturity_;
    PlainVanillaPayoff payoff_;
    Real s0_;
    Volatility sigma_;
    Rate r_;
    Rate q_;
    Real vega_;
};

class ReplicationPathPricer : public PathPricer<Path>
{
public:
    ReplicationPathPricer(Option::Type type, Real strike, Real r, Real q, Time maturity, Volatility sigma):
    type_(type), strike_(strike), r_(r), q_(q), maturity_(maturity), sigma_(sigma)
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
    Time maturity_;
    Volatility sigma_;
};


#endif /* ReplicationError_h */
