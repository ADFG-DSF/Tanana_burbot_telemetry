# Seasonal movements and distributions of radiotagged burbot in the mainstem Tanana River, 2018-2021
An analysis of the movements of instrumented burbot in the Tanana River and tributaries.  309 burbot were radiotagged within the three sections of the Tanana River (Lower, Middle, and Upper, as defined by the Statewide Harvest Survey) and spatial patterns and movements were monitored over seasonal telemetry flights between January 2019 and January 2021.

The objectives of the study were to 
1)	estimate the proportion of fish located in each river section during each survey,
2)	describe the spatial distributions and migrations of fish during surveys and seasons, and
3)	identify probable spawning areas during late January.


## Operational Plan

The current operational can be accessed here:
http://www.adfg.alaska.gov/FedAidPDFs/ROP.SF.3F.2018.10.pdf

## Folder structure 

### data/

* Three iterations of the telemetry dataset are preserved here, reflecting subsequent changes.      
  * **`TananaBurbotTelem3.csv` is the most recent**
  * `TananaBurbotTelem2.csv`
  * `TananaBurbotTelem.csv`

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
     
### tables/

The source data from several of the report tables are exported as .csv files here

### figures/

Several report figures are exported as .jpg here, as well as alternate versions.


