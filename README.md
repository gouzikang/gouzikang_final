# A New Regularization Method for High-dimensional Portfolio Selection
RUC course project: Advanced Statistical Analysis in Spring 2023.

# File description

Folder **Code** contains the main code for this project, where **data processing.ipynb** helps us to collate raw data. **admm.R**, **qpsolve.R** contain two models to tackle portfolio allocation problem via simulation data;  **admm_empirical.R** and **qpsolve_empirical.R** use the S&P 500 index to tackle portfolio allocation problem.

**SP 500 Stock Prices 2014-2017.csv** is the main data I use in my project. The data can be obtained via  https://www.kaggle.com/datasets/gauravmehta13/sp-500-stock-prices.


File **prensentation.pdf** is the slides of my presentation. 
File **report.pdf** is the final report.

# Conclusions

We conduct exploratory data analysis on stocks in S&P 500 index that have complete records from 2014 to 2016 and propose a new model to tackle the portfolio allocation problem. Our method has better performance in both simulation and empirical analysis, it has a better sharpe ratio, and selects fewer assets in the portfolio.

# Acknowledgement

The author thanks Prof. Ma Wei for his insightful comments and valuable suggestions that improve my final report. The author also thanks the provider of the dataset that allows my analysis to proceed.