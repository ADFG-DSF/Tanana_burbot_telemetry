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


