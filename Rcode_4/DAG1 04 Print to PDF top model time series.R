# now there is one more step, where I need to go back to the raw data and plot the time series that went into each model. 
# So for example, for DAG1C_long_topmodel, there were 7 nodes total. I want to plot those 7 time series. 
# The names that went into the model can be found in the guild.dfas1 dataframe. There should be an X at the beginning of each one, because it is a column name. I think we already had the X in there. 
# Can you use the same logic you used to find the right names and make a page for all of the time series that went into the final model?
  
source("LisaXP/functions/generate_model_time_series_plot_fxn")


# Define output PDF destination path
library(Cairo)
CairoPDF("LisaXP/outputs_4/DAG1_Raw_TimeSeries_Plots.pdf", width = 16, height = 12)

model_bases <- c("DAG1A_short", "DAG1A_long", "DAG1B_long", "DAG1B_short", "DAG1C_long", "DAG1C_short")

for (m_base in model_bases) {
  
  # ----------------------------------------------------
  # PART A: Standard Topmodel Time Series
  # ----------------------------------------------------
  model_name <- paste0(m_base, "_topmodel")
  
  if (exists(model_name)) {
    message(paste("Generating Time Series PDF page for:", model_name))
    fit_obj <- get(model_name)
    
    plot_ts_full <- generate_model_time_series_plot_fxn(model_name, fit_obj, var_lookup_NCC_AK, my_pretty_names)
    print(plot_ts_full)
  }
  
  # ----------------------------------------------------
  # PART B: Reduced Topmodel Time Series
  # ----------------------------------------------------
  reduced_name <- paste0(m_base, "_topmodel_reduced")
  
  if (exists(reduced_name)) {
    message(paste("Generating Time Series PDF page for:", reduced_name))
    fit_obj_red <- get(reduced_name)
    
    plot_ts_red <- generate_model_time_series_plot_fxn(reduced_name, fit_obj_red, var_lookup_NCC_AK, my_pretty_names)
    print(plot_ts_red)
  }
}

dev.off()
message("Truncation-aligned Time Series PDF processing complete!")