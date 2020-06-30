# US_NPI_Re

This repo contains data and code to replicate the findings presented in "Quantifying the impact of US state non-pharmaceutical interventions on COVID-19 transmission."

# csv files

csv files contain county level data necessary to run "to_report.R"

# R code files 

"Re_cases.R" provides code necessary to estimate Re from case or death data without stohcastic back forecasting.

"Re_rev.R" provides code necessary to estimate R0 from the exponential growth phase of the epidemic

"Re_lag.R" provides code necessary to estimate Re using stochastic back-forecasting to reconstruct incidence data

"to_report.R" contains the code for the comparison of Re and R0 estimates with county level data

# state actions

More detailed information on state actions is available in a google spreadsheet, though it will not updated after June, 2020:
https://docs.google.com/spreadsheets/d/17-cobYm1d6crJ3zPaOXLysX3kBZTez_dGP2tYI0p4aw/edit?usp=sharing
