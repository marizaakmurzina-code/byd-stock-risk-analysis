# BYD Stock Risk Analysis

**Financial Risk Analysis using R | Academic Project**

## Project Overview

This project analyzes the historical stock returns and financial risks of BYD Company Limited using statistical methods and probability distributions.

**Dataset:** Yahoo Finance, 2022–2026  
**Observations:** 969 daily stock prices  
**Tools:** R, RStudio, quantmod, fitdistrplus, sn, moments

## Methodology

- Data collection and preprocessing
- Exploratory Data Analysis (EDA)
- Daily logarithmic returns calculation
- Descriptive statistics and data visualization
- Probability distribution fitting
- Model comparison using AIC and goodness-of-fit tests
- Estimation of extreme loss probabilities

## Key Findings

- Daily volatility was approximately 1.99%.
- The largest observed daily price decline was approximately 10%.
- Skewed Student's t provided the best fit among the compared return distributions according to AIC.
- The Normal distribution did not adequately capture the heavy tails of stock returns.

## Conclusion

The analysis demonstrates the importance of selecting appropriate statistical distributions for financial risk assessment. Models accounting for heavy tails provide a more suitable description of extreme stock returns than a standard Normal distribution.
