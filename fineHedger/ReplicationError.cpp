#include "ReplicationError.hpp"
#include <stdio.h>
#include <vector>
#include <cmath>
#include <random>
#include <cstdlib>

using namespace QuantLib;
using namespace std;


void ReplicationError::compute(Size nTimeSteps, Size nSamples)
{
    QL_REQUIRE(nTimeSteps > 0, "the number of steps must be > 0");
    initialHedge();
    // if interval set,
    // Time tau = maturity_ / nTimeSteps;

    Calendar calendar = SouthKorea();
    Date today = Date::todaysDate();
    DayCounter dayCount = Actual365Fixed();
    Handle<Quote> stateVariable(ext::shared_ptr<Quote>(new SimpleQuote(s0_)));
    Handle<YieldTermStructure> riskFreeRate(ext::shared_ptr<YieldTermStructure>(new FlatForward(today, u_, dayCount)));
    Handle<YieldTermStructure> dividendYield(ext::shared_ptr<YieldTermStructure>(new FlatForward(today, 0.0, dayCount)));
    Handle<BlackVolTermStructure> volatility(ext::shared_ptr<BlackVolTermStructure>(new BlackConstantVol(today, calendar, sigma_, dayCount)));
    //Handle<BlackVolTermStructure> volatility(ext::shared_ptr<BlackVolTermStructure>(new BlackConstantVol(today, calendar, simulSigma_, dayCount)));
    ext::shared_ptr<StochasticProcess1D> diffusion(new BlackScholesMertonProcess(stateVariable, dividendYield, riskFreeRate, volatility));

    //
    PseudoRandom::rsg_type rsg = PseudoRandom::make_sequence_generator(nTimeSteps, 1);
    bool brownianBridge = false;
    typedef SingleVariate<PseudoRandom>::path_generator_type generator_type;

    //ext::shared_ptr<generator_type> myPathGenerator(new generator_type(diffusion, maturity_, nTimeSteps, rsg, brownianBridge));
    ext::shared_ptr<generator_type> myPathGenerator(new generator_type(diffusion, maturity_, nTimeSteps, rsg, brownianBridge));

    ext::shared_ptr<PathPricer<Path> > myPathPricer(new ReplicationPathPricer(payoff_.optionType(), payoff_.strike(),
        r_, q_, maturity_, HedgeVol_, transactionCost_));
    Statistics statisticsAccumulator;
    MonteCarloModel<SingleVariate, PseudoRandom> MCSimulation(myPathGenerator, myPathPricer, statisticsAccumulator, false);
    MCSimulation.addSamples(nSamples);

    PLMean = MCSimulation.sampleAccumulator().mean();
    PLStddev = MCSimulation.sampleAccumulator().standardDeviation();
    PLSkew = MCSimulation.sampleAccumulator().skewness();
    PLKurt = MCSimulation.sampleAccumulator().kurtosis();

    // Derman and Kamil //
    theorStD = std::sqrt(M_PI / 4 / nTimeSteps)*vega_*sigma_;

    // // //
    std::vector<pair<Real, Real> > temp1;
    temp1 = MCSimulation.sampleAccumulator().data();

    for (Size i = 0; i < temp1.size(); i++)
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
    QL_REQUIRE(n > 0, "the path cannot be empty");

    Time dt = maturity_ / n;
    //Rate stockDividendYield = q_;

    // simulation
    Time t = 0;
    Real stock = path.front();
    Real money_account = 0.0;

    DiscountFactor rDiscount = std::exp(-r_ * maturity_);
    DiscountFactor qDiscount = std::exp(-q_ * maturity_);
    Real forward = stock * qDiscount / rDiscount;
    Real stdDev = std::sqrt(sigma_*sigma_*maturity_);
    //Real stdDev = std::sqrt(sigma_*sigma_*maturity_); // Input Vol
    ext::shared_ptr<StrikedTypePayoff> payoff(new PlainVanillaPayoff(type_, strike_));


    BlackCalculator black(payoff, forward, stdDev, rDiscount);
    money_account += black.value();
    Real delta = black.delta(stock);
    Real stockAmount = delta;
    money_account -= stockAmount * stock;

    // hedging //
    for (Size step = 0; step < n - 1; step++)
    {
        t += dt; // time passes
        Real cashdiff = delta - stockAmount;
        money_account *= std::exp(r_*dt);
        stock = path[step + 1];
        rDiscount = std::exp(-r_ * (maturity_ - t));
        qDiscount = std::exp(-q_ * (maturity_ - t));
        forward = stock * qDiscount / rDiscount;
        stdDev = std::sqrt(sigma_*sigma_*(maturity_ - t));
        BlackCalculator black(payoff, forward, stdDev, rDiscount);
        delta = black.delta(stock);

        if (cashdiff >= 0)
        {
            money_account = money_account - (delta - stockAmount)*stock;
        }
        else if (cashdiff < 0)
        {
            money_account = money_account - (delta - stockAmount)*stock * (1 - transactionCost_);
        }
        // rehedging
        stockAmount = delta;
    }
    // at maturity //
    money_account *= std::exp(r_*dt);
    stock = path[n];
    Real optionPayoff = PlainVanillaPayoff(type_, strike_)(stock);
    money_account -= optionPayoff;
    money_account += stockAmount * stock;
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
        objectFunc.push_back(objectInput[dt - 1].first + objectInput[dt - 1].second);
    }

    Size min_index;
    min_index = min_element(objectFunc.begin(), objectFunc.end()) - objectFunc.begin() - startDt + 1;
    optimalTimes = min_index;
    std::cout << "Optimal Hedge Numbers / Counts : " << optimalTimes << std::endl;
}



double gaussian_box_muller() {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;

  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}


void ReplicationError::computePnL(Size nTimeSteps, Size nSamples)
{
    //std::random_device rd;
    std::vector<vector<Real>> sTHedgeVol(nTimeSteps+1, vector<Real>(nSamples+1, 0.000000));
    std::vector<vector<Real>> sTdrift(nTimeSteps + 1, vector<Real>(nSamples+1, 0.000000));
    std::vector<vector<Real>> sK(nTimeSteps + 1, vector<Real>(nSamples+1, 0.000000));
    std::vector<vector<Real>> tempDrift(nTimeSteps + 1, vector<Real>(nSamples+1, 0.000000));
    std::vector<vector<Real>> DeltaHedgeVol(nTimeSteps+1, vector<Real>(nSamples+1, 0.000000));
    std::vector<vector<Real>> DeltaInputVol(nTimeSteps+1, vector<Real>(nSamples+1, 0.000000));
    std::vector<vector<Real>> GammaInputVol(nTimeSteps+1, vector<Real>(nSamples+1, 0.000000));
    std::vector<vector<Real>> pnlProcess(nTimeSteps+1, vector<Real>(nSamples+1, 0.000000));
    std::vector<vector<Real>> pnlResult(nTimeSteps+1, vector<Real>(nSamples+1, 0.000000));
    
    Real HedgePrice = BSMPrice(HedgeVol_, s0_, maturity_);
    Real InputPrice = BSMPrice(sigma_, s0_, maturity_);
    //cout << HedgePrice << endl;
    //cout << InputPrice << endl;
    double upperBound = HedgePrice - InputPrice;
    double constInputDelta = BSMDelta(sigma_, s0_, maturity_);
    //cout << constInputDelta << endl;
    double constHedgeDelta=BSMDelta(HedgeVol_, s0_, maturity_);
    //cout << constHedgeDelta << endl;
    double constInputGamma = BSMGamma(sigma_, s0_, maturity_);
    //cout << constInputGamma << endl;
    double driftTerm = (0.5 * constInputGamma * s0_ * s0_ * (HedgeVol_ * HedgeVol_ - sigma_ * sigma_));// +(constInputD = 0.0elta - constHedgeDelta)*(u_ - r_ + q_) * s0_);
    
    std::vector<vector<double>> lowerBounds(nTimeSteps+1, vector<double>(nSamples+1, 0.000000));
    Decimal dt;
    dt = 1.000000 / nTimeSteps;
    
    std::mt19937_64 gen;
    std::normal_distribution<double> engine(0.0, 1.0);
    gen.seed(1);
    double rnStock1;
    rnStock1 = engine(gen);
    
    for (double i = 1; i <= nTimeSteps; i++)
    {
        //cout << constInputDelta << endl;
        for (double j = 1; j <= nSamples; j++)
        {
            double rnStock;
            rnStock = gaussian_box_muller();
            //cout << constInputDelta << endl;
            Real dtCum;
            Real Temp;
            Temp = (i) * dt * maturity_;
            dtCum = Temp;
            Real timeToMaturity = maturity_ - dtCum;
            
            sTHedgeVol[0][j] = s0_;
            sTHedgeVol[i][j] = (sTHedgeVol[i - 1][j] * std::exp((u_ -q_ - HedgeVol_ * HedgeVol_ * 0.5)*dt + HedgeVol_ * std::sqrt(dt)*rnStock));
            //sTHedgeVol[i][j] = sTHedgeVol[i - 1][j]+ sTHedgeVol[i - 1][j]* (u_)*dt + sTHedgeVol[i - 1][j]* HedgeVol_ * std::sqrt(dt)*rnStock*(dt*dt-1);
            sTdrift[0][j] = s0_;
            sTdrift[i][j] = (sTdrift[i - 1][j] * std::exp((u_ - r_ + q_ - HedgeVol_ * HedgeVol_ * 0.5)*dt + HedgeVol_ * std::sqrt(dt)*rnStock));

            DeltaHedgeVol[i][j] = BSMDelta(HedgeVol_, sTHedgeVol[i][j], timeToMaturity);
            DeltaInputVol[i][j] = BSMDelta(sigma_, sTHedgeVol[i][j], timeToMaturity);
            GammaInputVol[i][j] = BSMGamma(sigma_, sTHedgeVol[i][j], timeToMaturity);
            
            //double driftTerm2 = 0.5* GammaInputVol[i-1][j] * sTHedgeVol[i][j] *sTHedgeVol[i][j]
            //* (HedgeVol_ * HedgeVol_ - sigma_ * sigma_);
            double driftTerm2 = 0.5* GammaInputVol[i-1][j] * sTHedgeVol[i][j] *sTHedgeVol[i][j] * (HedgeVol_ * HedgeVol_ - sigma_ * sigma_);
            //double  (constInputDelta - constHedgeDelta) *  (u_ - r_ + q_);
            //double sigmaTerm = (constInputDelta - constHedgeDelta) * HedgeVol_;
            
            pnlProcess[i][j] = pnlProcess[i - 1][j] + dt * driftTerm2 + (DeltaInputVol[i-1][j] - DeltaHedgeVol[i-1][j]) * (sTdrift[i][j]- sTdrift[i-1][j]);
            //(sTdrift[i][j] - sTdrift[i - 1][j]);//sTdrift[i-1][j-1] * (exp((u_ - r_ + q_ - HedgeVol_ * HedgeVol_ * 0.5)*dt + HedgeVol_ * sqrt(dt)*rnStock)-1);
            //- sTdrift[i - 1][j]);//sTdrift[i-1][j-1] * (exp((u_ - r_ + q_ - HedgeVol_ * HedgeVol_ * 0.5)*dt + HedgeVol_ * sqrt(dt)*rnStock)-1);

            pnlResult[i][j] = pnlProcess[i][j] * std::exp(-r_ *(timeToMaturity));
            
            lowerBounds[i][j] = upperBound - strike_ * 2 * std::exp(-r_* (timeToMaturity)) * (std_cum_norm(0.5 * (HedgeVol_ - sigma_) * std::sqrt(timeToMaturity)) - 0.5);

            std::cout << "Time at : " << i << " , Path at : " << j << " , Cumm PnL : " << pnlResult[i][j] <<
                " , Upper Bound : " << upperBound << " , Lower Bound : " << lowerBounds[i][j] <<
                std::endl;
        }
        
    }
}
