##### GEP #####
library(openxlsx)
ttest_path <- "/Users/benho/Desktop/GoogleDrive_BS/benho/Projects/NanostringProject/Projects/210811_GEP/GEP_Ttests_consensus_grouping_2017-07-24_For_Annie_v2.xlsx"

SHH <- read.xls(ttest_path, sheet = "1vOthers")
TYR <- read.xls(ttest_path, sheet = "2vOthers")
MYC <- read.xls(ttest_path, sheet = "3vOthers")


#selection criteria
GEP_data <- SHH
pval_thres <- 0.05
fc_thres <- 3
type <- "up"


n_ori <- nrow(GEP_data)
GEP_data <- GEP_data[GEP_data$adj.P.Val < pval_thres,] # filter by pval first

if(type == "up"){
  GEP_data <- GEP_data[GEP_data$FoldChange > fc_thres,]
} else if (type == "down"){
  GEP_data <- GEP_data[GEP_data$FoldChange < -fc_thres,]
} else if (type == "both"){
  GEP_data_up <- GEP_data[GEP_data$FoldChange > fc_thres,]
  GEP_data_down <- GEP_data[GEP_data$FoldChange < -fc_thres,]
  GEO_data <- rbind(GEP_data_up,GEP_data_down)
}

n_new <- nrow(GEP_data)

print(paste0("[MSG] Ori: ",n_ori, " New: ", n_new))
