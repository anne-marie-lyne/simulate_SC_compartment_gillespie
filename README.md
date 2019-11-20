# simulate_SC_compartment_gillespie
Use Gillespie algorithm to simulate a stem cell compartment in which cells can self renew and differentiate whilst accumulating mutations

Main file is gillespie_simulate_SC_compartment.R in which simulation parameters can be set (mutation rate, number of cells, selective advantage etc)

compute_index_function.R is called by gillespie_simulate_SC_compartment.R and computes the linearity index at the final time point for each independent replicate

make_muller_plots.R can be used to plot figures showing mutation accumulation for specified replicates
