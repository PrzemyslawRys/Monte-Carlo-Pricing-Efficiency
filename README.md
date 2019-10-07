# Improving efficiency of Monte Carlo option pricing

The project contains the implementation of different features improving the efficiency of Monte Carlo option pricing in the Black-Scholes framework. The results of the tests are included in the RMarkdown report. Some of the considered features were based on efficient techniques of R programming, such as vectorization and others on statistical properties of Monte Carlo methods, mainly variance reduction methods, such as antithetic variables, control variables techniques and low-discrepancy series. Those methods make the procedure of derivative pricing more complex and time-consuming, which means that the same number of simulations last longer. Nevertheless, advantages in variance reduction were strong enough to recompensate that effect. Considered features need much fewer simulations to be accurate, thus obtaining results of the same quality requires, in fact, less time, than in classic approach.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Please check the prerequisites below and make sure to install all packages and libraries before running the code.

- *fOptions* package for *runif.halton* function, generating Halton sequence, Quasi Monte Carlo method is based on.
- *rbenchmark* package for measuring running time of R code,

### Input data structure

The input data structure is not required. Methods in this project are suited for derivatives pricing problem for predefined parameters. If you would like to use them for real financial assets, you should calibrate Black-Scholes model to market data. Calibrating methods are not part of this project.

### How to reproduce results?

In order to reproduce results you please start with running script *script-testEfficiency.R* from main folder, and then knit report from *report.Rmd* to produce .pdf report.

## Authors

* **Przemysław Ryś** - [see Github](https://github.com/PrzemyslawRys), [see Linkedin](https://www.linkedin.com/in/przemyslawrys/)

Codes were prepared for the Advanced R Programming course on the Faculty of Economic Sciences, University of Warsaw coordinated by dr. Piotr Wójcik.