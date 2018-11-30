load('SIM_OUTPUT/SIM_Wr0_nicheBreadth1_20181129_130519_989818.rda')

pdf('SIM_OUTPUT/Fig_coenoclines.pdf')
plot_coenoclines(sim_result = sim.result)
dev.off()
