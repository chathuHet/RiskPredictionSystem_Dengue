# RiskPredictionSystem_Dengue
This folder contains the code and data associated to the paper following.

**Title:** Risk prediction system for dengue transmission based on high resolution weather data, PLoS ONE, accepted for publication. 2018. 

**Authors:** C. Hettiarachchige, S. von Cavallar,  T.M. Lynar, R.I. Hickson, M. Gambhir

## Abstract
**Background**
Dengue is the fastest spreading vector-borne viral disease, resulting in an estimated 390
million infections annually. Precise prediction of many attributes related to dengue is
still a challenge due to the complex dynamics of the disease. Important attributes to
predict include: the risk of and risk factors for an infection; infection severity; and the
timing and magnitude of outbreaks. In this work, we build a model for predicting the
risk of dengue transmission using high-resolution weather data. The level of dengue
transmission risk depends on the vector density, hence we predict risk via vector
prediction.

**Methods and Findings**
We make use of surveillance data on Aedes aegypti larvae collected by the Taiwan
Centers for Disease Control as part of the national routine entomological surveillance of
dengue, and weather data simulated using the IBM's Containerized Forecasting
Workflow, a high spatial- and temporal-resolution forecasting system. We propose a two
stage risk prediction system for assessing dengue transmission via Aedes aegypti
mosquitoes. In stage one, we perform a logistic regression to determine whether larvae
are present or absent at the locations of interest using weather attributes as the
explanatory variables. The results are then aggregated to an administrative division,
with presence in the division determined by a threshold percentage of larvae positive
locations resulting from a bootstrap approach. In stage two, larvae counts are estimated
for the predicted larvae positive divisions from stage one, using a zero-inflated negative
binomial model. This model identifies the larvae positive locations with 71% accuracy
and predicts the larvae numbers producing a coverage probability of 98% over 95%
nominal prediction intervals. This two-stage model improves the overall accuracy of
identifying larvae positive locations by 29%, and the mean squared error of predicted
larvae numbers by 9.6%, against a single-stage approach which uses a zero-inflated
binomial regression approach.

**Conclusions**
We demonstrate a risk prediction system using high resolution weather data can provide
valuable insight to the distribution of risk over a geographical region. The work also
shows that a two-stage approach is beneficial in predicting risk in non-homogeneous
regions, where the risk is localised.

## File organisation:

* The dataset used for analysis is in the file mosdf.rds.

* The R markdown file analysis.rmd creates a report (analysis.html) which contains all the output in the paper.

* We also present this code paper_code.R file.

* All the functions required for running the analysis is saved in functions.R. We call this file in both the analysis.rmd and paper_code.r by ```source(here::here("functions.R"))```.

* The shapefile for Taiwan used for plotting Fig 3 is saved as GADM_2.8_TWN_adm2.rds.

* We have saved the dataframe resulting from the bootstarp approach in pos_percentage.rds to save the time in running the code.
If one wants to run the bootstrap approach, it can be done by uncommenting the relevant lines (lines 141-142 in analysis.rmd or lines 105-106 in paper_code.r).


