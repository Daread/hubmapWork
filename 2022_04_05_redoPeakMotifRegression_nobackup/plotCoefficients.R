



library(ggplot2)





csvPath = "./plots/Vascular_Endothelium_fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult/Vascular_Endothelium_SexM_Table.csv"


myCSV = read.csv(csvPath)



png('./plots/VascEndothPlot.png', res=200, height =1200, width=1400)

myPlot = ggplot(myCSV, aes(x=gene, y=coefficientValue)) + 
	geom_col() + ylab("Coefficient (Positive -> Higher in Men)") + 
	xlab("TF Motif")

plot(myPlot)
dev.off()










