# Seasonal movements and distributions of radiotagged burbot in the mainstem Tanana River, 2018-2021
An analysis of the movements of instrumented burbot in the Tanana River and tributaries.  THIS IS WHERE WE SHOULD GIVE A BIT MORE BACKGROUND

## Operational Plan

The current operational can be accessed here:
http://www.adfg.alaska.gov/FedAidPDFs/ROP.SF.3F.2018.10.pdf

## Folder structure 

### data/

* `TananaBurbotTelem.csv` is the full telemetry dataset, saved as a .csv file

* `tyb.Rdata` is an R workspace giving the river shapefile expressed as a rivernetwork object (for use within the 'riverdist' package)

### R/

* `tyb_analysis` 
  * Reads & processes data
  * Runs analyses
    * Length composition
      * means, proportions (binned)
      * also by river section - evidence of differences found
    * Identifying areas of relatively high density
      * kernel density (network)
      * kernel density (linear)
    * Proportions in each river section
      * per survey
      * per survey, for each tagged stock
    * Movement by individual - minimum homerange and minimum cumulative distance
      * also relationships between this and other individual-level variables
    * Distances between sequential surveys
    * Overlap in spatial extent 
    * Estimation of survival parameters between survey events
     * Bayesian state-space model in which un-observed survival is modeled as unobserved parameter


