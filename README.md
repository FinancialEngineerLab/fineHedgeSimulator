# fineHedger

You need to read the Result.png

This Program is used for discrete dynamic hedging simulation to analyize impact of hedging numbers as steps, transcation cost, mismatch of volatility and mismatch of drift term.

2021-07-30 : hedging numbers when paying transaction cost in delta hedging

2021-08-02 : mismatch volatility and drift term, respectively
 > What if there is mismatched of volatility (i.e Real vs Implied) or drift (i.e Risk Free Rate vs Mean)

2021-08-10 : Optimal Hedging Trading Numbers with Objective Functions as Optimization
 > you can set up the Utility function for optimization in hedging performance,
 > then this program give the Optimal Hedging Frequency satisfied with the utility function as well as minizing transaction cost

2021-08-31 : Hedging Real vol vs Pricing Implied vol with Upper and Lower Bound Level
 > When Realized Volatility is larger than Implied Volatlity, then long call option and delta hedging,
 > the Profit and Loss Process (SDE) can be described as
 > convergence to Upper bound with difference of two options (+) following to decay of time maturity
 > and Lower bound can be calculated in formula of Vol difference and time decay
 > and the Profit and Loss Process (SDE) decomposed into two major term
 > first, Gamma (Gamma + Gamma from theta)
 > Second, Stock Path in drift term manipulated
2021-08-? : Gamma Trading


main refernece:

Emanuel Derman's
Ch4. The P&L of hedged option strategy in the BSM world, Volatility Smile (2016)

Ch5. The Effect of Discrete Hedging on P&L, Volatility Smile (2016)

Ch6. The Effect of Transaction cost on P&L, Volatility Smile (2016)

Derman & Kamal's
(Goldman Sachs Equity Derivatives Research) Research Note: "When You Cannot Hedge Continuously: The Corrections to Black-Scholes"
    http://www.ederman.com/emanuelderman/GSQSpapers/when_you_cannot_hedge.pdf

needs: 
boost library
quantlib library (1.21 ~)

environment: xcode
