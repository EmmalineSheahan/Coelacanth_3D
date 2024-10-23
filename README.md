# Coelacanth_3D
A code repository for creating 3 dimensional and 2 dimensional ecological niche models of the coelacanth, Latimeria chalumnae

Important scripts:
1. coelacanth_full_occ_cleaning_and_thinning.R: cleaning and thinning of 2D and 3D occurrence records
2. formatting_env_data.R: formatting environmental data retrieved from WOA and the Copernicus Marine Service in both 3D and 2D
3. 3D_accessible_area.R: creation of accessible area polygons on a per depth slice basis for 3D modelling
4. 2D_accessible_area.R: creation of singular 2D accessible area polygon for 2D ENM
5. all_functions_3D_sdm.R: script containing all functions necessary for 3D ecological niche modelling. this script is sourced in 3D_sdm_coelacanth.R
6. 2D_sdm_coelacanth.R: the script where the 2D models for chalumnae is produced
7. 3D_sdm_coelacanth.R: the script where the 3D models for chalumnae is produced
8. l_menadoensis_occ_cleaning.R: where the 2D and 3D occurrence points for menadoensis are cleaned
9. projecting_to_indopacific.R: where the 2D and 3D models are projected onto the indopacific
10. results_plotting.R: where the primary figures are made
11. chloropleth_maps.R: where the figures containing the bivariate chloropleth maps are made
12. depth_slice_plots.R: where the depth slice plots, figures 10 and 11, are made

Supplemental scripts:
1. mess_analysis.R: where the MESS analyses in 2D and 3D are conducted
2. recipricol_model.R: where the 3D ENM pipeline is run for menadoensis
3. reciprocal_model_projection.R: where the resulting menadoensis model is projected back onto the Mozambique channel

One will want to be certain to rewrite code as needed to fit one's working directory, as well as ensure one has all of the necessary folders created. 
