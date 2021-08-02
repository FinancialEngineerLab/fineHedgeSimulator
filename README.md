# fineHedger

This Program is used for discrete dynamic hedging simulation to analyize impact of hedging numbers as steps, transcation cost, mismatch of volatility and mismatch of drift term.

2021-07-30 : hedging numbers transaction cost

2021-08-?? : mismatch volatility and drift term, respectively

main refernece:

Emanuel Derman's
Ch5. The Effect of Discrete Hedging on P&L, Volatility Smile (2016)

Ch6. The Effect of Transaction cost on P&L, Volatility Smile (2016)


Derman & Kamal's
(Goldman Sachs Equity Derivatives Research) Research Note: "When You Cannot Hedge Continuously: The Corrections to Black-Scholes"
    http://www.ederman.com/emanuelderman/GSQSpapers/when_you_cannot_hedge.pdf


Sample Result


Option Value: 2.3419

Samples |   Trades |     Mean |   Stddev |     formaula | Skewness | Kurtosis

------------------------------------------------------------------------------
 

Hedge : 21hedging number

   50000 |       21 |   -0.044 |    0.427 |        0.445 |   -0.299 |    1.439

 

 

Hedge : 84hedging number

   50000 |       84 |   -0.044 |    0.220 |        0.222 |   -0.204 |    1.684

 

 

*********** Mismatch of Drift *************

Hedge : 84hedging number

   50000 |       84 |   -0.044 |    0.220 |        0.222 |   -0.201 |    1.678

 

 

*********** Mismatch of Volatility *************

Hedge : 84hedging number

   50000 |       84 |   -2.350 |    1.243 |        0.222 |   -0.655 |    0.052

 

 

******* Mismatch of Volatility & Drift Term *******

Hedge : 84hedging number

   50000 |       84 |   -2.349 |    1.244 |        0.222 |   -0.654 |    0.040



needs: 
boost library
quantlib library (1.21 ~)

environment: xcode
