Readme for the 'biodiversity_ES_ensembles' project. 

Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2024-06-05
Last update: 2025-05-08

Raw data requirements
---------------------
-- Biodiversity (species richness) data
    |- RData files downloaded Cai, L. et al. Global models and predictions of plant diversity based on advanced machine learning techniques. New Phytologist 237, 1432–1445 (2023). These were downloaded at the 7,774 km2 resolution for all of the files that was not the ensemble.
-- Ecosystem services data
	|- Obtained from Willcock, S. et al. Model ensembles of ecosystem services fill global certainty and capacity gaps. Science Advances 9, eadf5492 (2023). The native output from all the models was obtained from the authors of that paper, as was the 1 km2 normalised data. 
---------------------

(final) File architecture
---------------------
-- data
    |- raw
    |- checks
    |- intermediate_outputs
		|- biodiversity
			|- normalised
				|- res[resolution]
					|- continents
						|- normalised plant species richness data at the continental level
					|- countries
						|- normalised plant species richness data at the country level
		|- continents
			|- outline polygons
			|- template rasters
		|- countries
			|- outline polygons
			|- template rasters
		|- ecosystem_services
			|- normalised
				|- res[resolution]
					|- continents
						|- normalised ecosystem service data at the continental level
					|- countries
						|- normalised ecosystem service data at the country level			
-- results
	|- figures
	|- tables
-- code
	|- r_scripts
		|- all scripts for the analysis and visualisation
---------------------

File origins
---------------------
-- create_continents.R
    |- data/intermediate_outputs/continents/[continent]_final.gpkg
-- create_countries.R    
    |- data/intermediate_outputs/countries/[country]_final.gpkg
-- Cai_data_download.R
	|- data/raw/vascular_plants/Cai_2023/sr_[model type]_Prediction_[resolution]_wgs.gpkg
-- template_rasters.R
	|- data/[continent]_biodiversity_[resolution]_weighting.csv
	|- data/intermediate_outputs/continents/[continent][resolution]template.tif
-- data_to_continent_extent.R
	|- data/intermediate_outputs/biodiversity/normalised/res100/continents/sr_[model type]_norm_[continent].tif
	|- data/intermediate_outputs/ecosystem_services/normalised/res100/continents/[ecosystem service]/[model type]_norm_[continent].tif
-- correlation_test.R
	|- results/tables/correlation_summary_BD.csv
-- hotspot_analysis.R
	|- results/tables/[ecosystem service]100[continent]_sr_rw[row number].csv
	|- results/tables/[ecosystem service]100[continent]_sr_END.csv
-- uncert_graphs.R
	|- results/figures/combined/combinedES[date].pdf
	|- results/figures/single_item/singleES[date].pdf
	|- results/figures/single_item/singleES[ecosystem service]_[continent].png
-- map_creation_final.R
	|- results/figures/maps/spatial_agreement[date].pdf
	|- results/figures/maps/spatAgree_maps_graphs[date].pdf
-- paper_graphs.R
	|- results/tables/final_results.csv
	|- results/tables/summExitResults.csv
	|- results/tables/ExitResultsShort.csv
	|- results/tables/SI_table_[ecosystem service][date].csv
	|- results/figures/maps/Maps_[ecosystem service]_SR_[continent].pdf
	|- results/figures/maps/[ecosystem service]_[percentile cut-off]_[continent][date].pdf
	|- results/figures/singlePercGraphs_[percentile cut-off][date].pdf
	|- results/figures/SI_maps_SpatAgree.pdf
	|- results/figures/Spatial_Agreement_figures_[date].pdf
	|- results/figures/SI_graphs[ecosystem service][date].pdf
---------------------