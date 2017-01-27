This project aims to create and refine a statistical graphic and Shiny app for the purposes of making individual PK/PD predictions using oppportunistic/scavenged samples. The graphic includes statistical summaries, including credible bands and intervals. The methods used to generate these summaries will be evaluated and revised to achieve 1) approximately correct operating characteristics and 2) minimize computing latency. The products of this work will be used to facilitate an therapeutic drug monitoring (TDM) hospital intervention.

### Dependencies (R packages)
- 'rhandsontable'
- 'shinydashboard'

### Files:

- 'README.md': this README
- 'model.R': two-compartment PK model for multiple infusion dosing
- 'Bayes.R': Bayesian framework and plotting prior/posterior summaries
- 'server.R': Shiny application server
- 'ui.R': Shiny application UI


### Resources:

http://jrowen.github.io/rhandsontable/

### Comments:
- Graphics produced by application look much better when the R package 'Cairo' is available. On some machines (notable Linux), this may require installation of external libraries (e.g., libcairo2-dev, libxt-dev)
