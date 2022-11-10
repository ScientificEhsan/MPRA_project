# This R script helps download UKBB dataset by tidying the UKBB Cancer data xlsx
# Prepared by Shafqat Ehsan

library(tidyverse)
library(readxl)
library(dplyr)
UKBB_Cancers <- read_excel("MPRA/UKBB Cancers.xlsx") #reads the UKBB file
head(UKBB_Cancers)
colnames(UKBB_Cancers)

# str(UKBB_Cancers)

aws_link_column = UKBB_Cancers[,c(2,39,41)] #filters everything except the columns Phenocode, filename and aws_link
head(aws_link_column)
colnames(aws_link_column)
view(aws_link_column)
df=UKBB_Cancers[(1:29),c(2,39,41)] #filers for C02-C076
view(df)
write.csv(df, file="C76_AWS_Dowloadlink.csv") #creates new file for the filtered data for C02-C076


for (i in 1:ncol(df)){
                      download.file(url = df$column_1[i], 
                      destfile = paste0('files', df$column_ID_number[i], 
                                  '.tar.gz'), method = 'curl')}
