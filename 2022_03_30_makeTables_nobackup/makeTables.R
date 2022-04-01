
# install.packages("gt")

library(gt)
library(tidyverse)

library(webshot)
# webshot::install_phantomjs()

# Read in a csv giving sample data
sampleTable = "./inputFiles/Heart_Sample_Data.csv"

# Output as an aesthetically nice table
sampleDF = read.csv(sampleTable)

myTable = sampleDF %>% gt() %>%
		   tab_header(
		   	title="Samples Processed for sci RNA- and ATAC-Seq"
		   	)

tableFile = "./outputs/Heart_Sample_Data_Table.png"


gtsave(myTable, tableFile)

# png()
# print(myTable)
# dev.off()




