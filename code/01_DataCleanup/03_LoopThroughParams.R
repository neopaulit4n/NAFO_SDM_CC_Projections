
# Looping through periods/SSPs/VME groups

# Load data ----
source("code/01_LoadData.R")


# Create table of each combination of iterations to save variable selection results ----
vme_all <- unique(resp_df$VME_Group)
period_all <- c("P1", "P2", "P3", "P4")
ssp_all <- unique(cmip_df_period_ssp$ssp)

var_select_df <- expand_grid(vmeoi = vme_all,
                       poi = period_all,
                       sspoi = ssp_all) %>%
  arrange(vmeoi, poi, sspoi) %>%
  # Add columns for CMIP variables
  cbind(matrix(NA, nrow = nrow(.), ncol = length(cmip_vars))) %>%
  set_names(c("vmeoi", "poi", "sspoi", cmip_vars))


# Initialise fold metrics dataframe to save results from each fold of each iteration
fold_metrics_summary_df <- data.frame(VME_Group = character(),
                                     Period = character(),
                                     SSP = character(),
                                     Subsampled = logical(),
                                     SelVars = logical(),
                                     metric = character(),
                                     mean_value = numeric(),
                                     sd_value = numeric())

# Save raster prediction outputs for each iteration in a list
raster_output_list <- list()


# Loop through each combination of VME group, period, SSP, and subsampling option ----
loop_seed <- 411
for (vmeoi in "black_corals") {  # vme_all
  for (poi in period_all) {
    for (sspoi in ssp_all) {
      # Call modelling code for this combination of parameters
      cat("Running modelling for VME group:", vmeoi, 
          "Period:", poi, 
          "SSP:", sspoi,
          "\n")
      loop_seed <- loop_seed + 1
      source("code/03_Modelling.R")
    }
  }
}


# List the output raster names for each VME ----

sdm2024_raster_output_df <- data.frame(VME_group = vme_all,
                                      MaxClass = c("Black Coral/BlackCorals_raster_output_sens_spe/BlackCorals_raster_output_maxclass.tif",
                                                   "Boltenia/Boltenia_Model_Outputs/Boltenia_SensSpec/Boltenia_SensSpecSENSSPEC_raster_output_maxclass.tif",
                                                   "Bryozoans/Bryozoans_Model_Outputs/Bryozoans_SensSpec/Bryozoans_SensSpecSENSSPEC_raster_output_maxclass.tif",
                                                   "Large Gorgonians/LG_Model_Outputs/LG_SensSpec/lg_SensSpecSENSSPEC_raster_output_maxclass.tif",
                                                   "Sea Pens/Seapens_raster_output_maxclass.tif",
                                                   "Small Gorgonians/SG_Model_Outputs/SG_FG_SensSpec/sg_fg_SensSpecSENSSPEC_raster_output_maxclass.tif",
                                                   "Sponges/SpongesFG_v2/SpongesFG_v2_raster_output_maxclass.tif",
                                                   "Sponges/Sponges_raster_output_maxclass.tif"),
                                      MaxClassF = c("Black Coral/BlackCorals_raster_output_sens_spe/BlackCorals_raster_output_maxclassf.tif",
                                                   "Boltenia/Boltenia_Model_Outputs/Boltenia_SensSpec/Boltenia_SensSpecSENSSPEC_raster_output_maxclassf.tif",
                                                   "Bryozoans/Bryozoans_Model_Outputs/Bryozoans_SensSpec/Bryozoans_SensSpecSENSSPEC_raster_output_maxclassf.tif",
                                                   "Large Gorgonians/LG_Model_Outputs/LG_SensSpec/lg_SensSpecSENSSPEC_raster_output_maxclassf.tif",
                                                   "Sea Pens/Seapens_raster_output_maxclassf.tif",
                                                   "Small Gorgonians/SG_Model_Outputs/SG_FG_SensSpec/sg_fg_SensSpecSENSSPEC_raster_output_maxclassf.tif",
                                                   "Sponges/SpongesFG_v2/SpongesFG_v2_raster_output_maxclassf.tif",
                                                   "Sponges/Sponges_raster_output_maxclassf.tif"),
                                      MaxClassAvgProb = c("Black Coral/BlackCorals_raster_output_sens_spe/BlackCorals_raster_output_maxclassaveprob.tif",
                                                   "Boltenia/Boltenia_Model_Outputs/Boltenia_SensSpec/Boltenia_SensSpecSENSSPEC_raster_output_maxclassaveprob.tif",
                                                   "Bryozoans/Bryozoans_Model_Outputs/Bryozoans_SensSpec/Bryozoans_SensSpecSENSSPEC_raster_output_maxclassaveprob.tif",
                                                   "Large Gorgonians/LG_Model_Outputs/LG_SensSpec/lg_SensSpecSENSSPEC_raster_output_maxclassaveprob.tif",
                                                   "Sea Pens/Seapens_raster_output_maxclassaveprob.tif",
                                                   "Small Gorgonians/SG_Model_Outputs/SG_FG_SensSpec/sg_fg_SensSpecSENSSPEC_raster_output_maxclassaveprob.tif",
                                                   "Sponges/SpongesFG_v2/SpongesFG_v2_raster_output_maxclassaveprob.tif",
                                                   "Sponges/Sponges_raster_output_maxclassaveprob.tif"),
                                      CombConf = c("Black Coral/BlackCorals_raster_output_sens_spe/BlackCorals_raster_output_combconf.tif",
                                                   "Boltenia/Boltenia_Model_Outputs/Boltenia_SensSpec/Boltenia_SensSpecSENSSPEC_raster_output_combconf.tif",
                                                   "Bryozoans/Bryozoans_Model_Outputs/Bryozoans_SensSpec/Bryozoans_SensSpecSENSSPEC_raster_output_combconf.tif",
                                                   "Large Gorgonians/LG_Model_Outputs/LG_SensSpec/lg_SensSpecSENSSPEC_raster_output_combconf.tif",
                                                   "Sea Pens/Seapens_raster_output_combconf.tif",
                                                   "Small Gorgonians/SG_Model_Outputs/SG_FG_SensSpec/sg_fg_SensSpecSENSSPEC_raster_output_combconf.tif",
                                                   "Sponges/SpongesFG_v2/SpongesFG_v2_raster_output_combconf.tif",
                                                   "Sponges/Sponges_raster_output_combconf.tif"),
                                      CVSum = c("Black Coral/BlackCorals_raster_output_sens_spe/BlackCoralsVME_raster_output_cvsum.tif",
                                                   "Boltenia/Boltenia_Model_Outputs/Boltenia_SensSpec/Boltenia_SensSpecSENSSPECVME_raster_output_cvsum.tif",
                                                   "Bryozoans/Bryozoans_Model_Outputs/Bryozoans_SensSpec/Bryozoans_SensSpecSENSSPECVME_raster_output_cvsum.tif",
                                                   "Large Gorgonians/LG_Model_Outputs/LG_SensSpec/lg_SensSpecSENSSPECVME_raster_output_cvsum.tif",
                                                   "Sea Pens/SeapensVME_raster_output_cvsum.tif",
                                                   "Small Gorgonians/SG_Model_Outputs/SG_FG_SensSpec/sg_fg_SensSpecSENSSPECVME_raster_output_cvsum.tif",
                                                   "Sponges/SpongesFG_v2/SpongesFG_v2VME_raster_output_cvsum.tif",
                                                   "Sponges/SpongesVME_raster_output_cvsum.tif"))



