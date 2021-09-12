/*

void ReplicationError::computePnL(Size nSim)
{
    DiscountFactor rDiscount;// = std::exp(-r_*maturity_);
    DiscountFactor qDiscount;// = std::exp(-q*maturity_);;
    
    Real optionPayoff;// = PlainVanillaPayoff(type_, strike_)(stock);
    Real forward;//= s0_*qDiscount/rDiscount;
    Real stdDev;// = std::sqrt(sigma_*sigma_*maturity_);
    
    ext::shared_ptr<StrikedTypePayoff> payoff(new PlainVanillaPayoff(type_, strike_));

    vega_ = black.vega(maturity_);
    OptionValue_ = black.value();
    
    std::mt19937_64 gen;
    std::normal_distribution<double> engine(0.0, 1.0);
    gen.seed(1);
    
    vector<vector<Real>> PNLPaths;
    vector<vector<Real>> S;
    Real UpperBound;
    vector<Real> LowerBound;
    
    PNLPaths[0] = vector<Real>(0);
    S[0] = vector<Real>(s0_);
    
    Real gammaIV;
    Real DeltaIV;
    Real DeltaHG;
    
    for(Time i = 1; i < nTimeSteps; i++)
    {
        double df = std::exp(-r_ * i / nTimeSteps * maturity_);
        for (Size j = 0; j< nSamples;j++)
        {
            double e = engine(gen);
            S[i][j] = S[0][0] * std::exp((r_-q_-0.5*simulSigma_*simulSigma_) * i / nTimeSteps * maturity_)
            + simulSigma_ * std::sqrt(i / nTimeSteps * maturity_) * e;
            
            
        
            BlackCalculator black(payoff, forward, stdDev, rDiscount);
            
            
            PNLPaths[i][j] = PNLPaths[i-1][j] + i / nTimeSteps * (0.5* gammaIV * S[i][j] * S[i][j] * (HedgeVol_ * HedgeVol_ - simulSigma_ * simulSigma_) + (DeltaIV - DeltaHG) * S[i][j] * (u_ - r_ + q_)) +
            (DeltaIV - DeltaHG) * HedgeVol_ * S[i][j] * e;
        }
    }
}

*/
