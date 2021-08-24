
#include "ReplicationError.hpp"
#include <stdio.h>
#include <vector>

using namespace QuantLib;
using namespace std;


void ReplicationError::compute(Size nTimeSteps, Size nSamples)
{
    QL_REQUIRE(nTimeSteps>0, "the number of steps must be > 0");
	initialHedge();
    // if interval set,
    // Time tau = maturity_ / nTimeSteps;
    
    Calendar calendar = SouthKorea();
    Date today = Date::todaysDate();
    DayCounter dayCount = Actual365Fixed();
    Handle<Quote> stateVariable(ext::shared_ptr<Quote>(new SimpleQuote(s0_)));
    Handle<YieldTermStructure> riskFreeRate(ext::shared_ptr<YieldTermStructure>(new FlatForward(today, u_, dayCount)));
    Handle<YieldTermStructure> dividendYield(ext::shared_ptr<YieldTermStructure>(new FlatForward(today, 0.0, dayCount)));
    Handle<BlackVolTermStructure> volatility(ext::shared_ptr<BlackVolTermStructure>(new BlackConstantVol(today, calendar, simulSigma_, dayCount)));
    ext::shared_ptr<StochasticProcess1D> diffusion(new BlackScholesMertonProcess(stateVariable, dividendYield, riskFreeRate, volatility));
    
    //
    PseudoRandom::rsg_type rsg = PseudoRandom::make_sequence_generator(nTimeSteps, 1);
    bool brownianBridge = false;
    typedef SingleVariate<PseudoRandom>::path_generator_type generator_type;
    ext::shared_ptr<generator_type> myPathGenerator(new generator_type(diffusion, maturity_, nTimeSteps, rsg, brownianBridge));
    
    ext::shared_ptr<PathPricer<Path> > myPathPricer(new ReplicationPathPricer(payoff_.optionType(), payoff_.strike(),
                                                                              r_, q_, maturity_, sigma_, transactionCost_));
    Statistics statisticsAccumulator;
    MonteCarloModel<SingleVariate, PseudoRandom> MCSimulation(myPathGenerator, myPathPricer, statisticsAccumulator, false);
    MCSimulation.addSamples(nSamples);
    
    PLMean = MCSimulation.sampleAccumulator().mean();
    PLStddev = MCSimulation.sampleAccumulator().standardDeviation();
    PLSkew = MCSimulation.sampleAccumulator().skewness();
    PLKurt = MCSimulation.sampleAccumulator().kurtosis();
    
    // Derman and Kamil //
    theorStD = std::sqrt(M_PI/4/nTimeSteps)*vega_*sigma_;
    
    // // //
    std::vector<pair<Real, Real> > temp1;
    temp1= MCSimulation.sampleAccumulator().data();
    
    for (Size i =0; i<temp1.size(); i++)
    {
        tempPaths.push_back(temp1[i].first);
    }
    /*
    for (Size i = 0; i < nSamples; i++)
    {
        for (Size j = 0; j < 1/nTimeSteps;j++)
        {
            pnlPaths.push_back(MCSimulation.sampleAccumulator().data());
        }
    }
   */
}

void ReplicationError::printResult()
{
	std::cout << std::fixed << std::setprecision(3)
		<< std::setw(8) << nSamples << " | "
		<< std::setw(8) << nTimeSteps << " | "
		<< std::setw(8) << PLMean << " | "
		<< std::setw(8) << PLStddev << " | "
		<< std::setw(12) << theorStD << " | "
		<< std::setw(8) << PLSkew << " | "
		<< std::setw(8) << PLKurt << std::endl;
}


Real ReplicationPathPricer::operator()(const Path& path) const
{
    Size n = path.length() - 1;
    QL_REQUIRE(n>0, "the path cannot be empty");
    
    Time dt = maturity_/n;
    //Rate stockDividendYield = q_;
    
    // simulation
    Time t =0;
    Real stock= path.front();
    Real money_account = 0.0;
    
    DiscountFactor rDiscount = std::exp(-r_*maturity_);
    DiscountFactor qDiscount = std::exp(-q_*maturity_);
    Real forward = stock * qDiscount / rDiscount;
    Real stdDev = std::sqrt(sigma_*sigma_*maturity_);
    ext::shared_ptr<StrikedTypePayoff> payoff(new PlainVanillaPayoff(type_, strike_));
    
    BlackCalculator black(payoff, forward, stdDev, rDiscount);
    money_account += black.value();
    Real delta = black.delta(stock);
    Real stockAmount = delta;
    money_account -= stockAmount*stock;
    
    // hedging //
    for (Size step = 0; step<n-1;step++)
    {
        t += dt; // time passes
        Real cashdiff = delta-stockAmount;
        money_account *= std::exp(r_*dt);
        stock = path[step+1];
        rDiscount = std::exp(-r_*(maturity_-t));
        qDiscount = std::exp(-q_*(maturity_-t));
        forward = stock*qDiscount/rDiscount;
        stdDev=std::sqrt(sigma_*sigma_*(maturity_-t));
        BlackCalculator black(payoff, forward, stdDev, rDiscount);
        delta = black.delta(stock);
        
        if (cashdiff >=0)
        {
            money_account = money_account -(delta-stockAmount)*stock;
        }
        else if (cashdiff < 0)
        {
            money_account = money_account -(delta-stockAmount)*stock * (1-transactionCost_);
        }
         // rehedging
        stockAmount = delta;
    }
        // at maturity //
        money_account *= std::exp(r_*dt);
        stock = path[n];
        Real optionPayoff = PlainVanillaPayoff(type_, strike_)(stock);
        money_account -= optionPayoff;
        money_account += stockAmount *stock;
        return money_account;
        
}
    
    
 

void ReplicationError::optimalHedging(Size maxDt)
{
	vector<pair<Real, Real>> objectInput;
	vector<Real> objectFunc;
	Size startDt = 1;

	for (Size dt = startDt; dt <= maxDt; dt++)
	{
		compute(dt, nSamples);
		objectInput.push_back(make_pair(PLMean, PLStddev));
		objectFunc.push_back(objectInput[dt-1].first + objectInput[dt-1].second);
	}

    Size min_index;
    min_index = min_element(objectFunc.begin(), objectFunc.end()) - objectFunc.begin()- startDt + 1;
    optimalTimes = min_index;
	std::cout << "Optimal Hedge Numbers / Counts : " << optimalTimes << std::endl;
}


