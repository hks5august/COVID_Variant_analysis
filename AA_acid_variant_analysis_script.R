library(ggplot2)
library(reshape)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(stringr)

library(gt)




seq <- 'ORF1a_protein'

#define start position minus 1
start_pos_1 = 265
Protein <- 'ORF1a_Protein_AA_report.txt'
add = 0.22


report_name <- toString(sub(".txt", ".pdf", Protein))


#Define Nucleotide variant positions
pos <- factor(c(353,564,625,683,845,896,913,1059,1284,1482,1944,2110,2992,3037,3261,3267,3483,3485,4069,4105,4150,4212,4486,5388,5852,5895,5986,6166,6954,7042,7332,7834,7984,8083,8123,8127,8179,8917,9043,9120,9389,9929,10030,10039,10153,10195,10274,10319,10332,10369,10793,10802,10901,11083,11130,11198,11296,11598,11747,11782,12357,12484,12525,12708,12747,13130,13261))

mydata1 <- read.table('AA_input.txt', sep="\t", header=TRUE, check.names = F)

#Extract only mutated positions
mydata2 <- mydata1[mydata1$POS %in% pos,]
head(mydata2)
new_dfT <- as.data.frame(mydata2)

msa1sel <- new_dfT
dim(msa1sel)
head(msa1sel[1])

#select only data where any sample has a variant to create a report
#new_df <- msa1sel[!(msa1sel[2] == msa1sel[3] & msa1sel[2] == msa1sel[4] & msa1sel[2] == msa1sel[5] & msa1sel[2] == msa1sel[6] & msa1sel[2] == msa1sel[7] & msa1sel[2] == msa1sel[8] & msa1sel[2] == msa1sel[9] & msa1sel[2] == msa1sel[10] & msa1sel[2] == msa1sel[11]),] 
new_df <- msa1sel[!(msa1sel[2] == msa1sel[2] & msa1sel[2] == msa1sel[3] & msa1sel[2] == msa1sel[4] & msa1sel[2] == msa1sel[5] & msa1sel[2] == msa1sel[6] & msa1sel[2] == msa1sel[7] & msa1sel[2] == msa1sel[8] & msa1sel[2] == msa1sel[9] & msa1sel[2] == msa1sel[10] & msa1sel[2] == msa1sel[11]& msa1sel[2] == msa1sel[12]& msa1sel[2] == msa1sel[13]& msa1sel[2] == msa1sel[14]& msa1sel[2] == msa1sel[15]& msa1sel[2] == msa1sel[16]& msa1sel[2] == msa1sel[17]& msa1sel[2] == msa1sel[18]& msa1sel[2] == msa1sel[19]& msa1sel[2] == msa1sel[20]& msa1sel[2] == msa1sel[21]& msa1sel[2] == msa1sel[22]& msa1sel[2] == msa1sel[23]& msa1sel[2] == msa1sel[24]),] 


#create column no. and row no.
x <- ncol(new_df)
x
y<- nrow(new_df)
y
#create an empty data frame to store the new data
df_new1 <- data.frame(matrix(0, ncol = x+1, nrow = y))
colnames(df_new1) <- colnames(new_df)

#define i value
i=3

#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 3:x) {
  
  df_new1[,i] <- ifelse(new_df[,2]== new_df[,i],"-", paste0(new_df[,2],">",new_df[,i]))
}


#Remove reference genome column and last column (NA column) from the dataframe 
df_new2 <- select(df_new1,-1, -2, -ncol(df_new1))
head(df_new2)

NTD_pos <- rownames(new_df[1])
#new_pos <- as.data.frame(t(row.names(new_df)))
NTD_pos
#Compute AA position
aa_pos1 <- (as.numeric(NTD_pos) - start_pos_1)/3 +add
AA_pos <- as.data.frame(round(aa_pos1,0))
positions <- cbind(NTD_pos, AA_pos)
colnames(positions) <- c( "Nucleotide position", "Amino Acid position")

#Add position column
final_df <- cbind(positions, df_new2)
final_df
#Check whether output is correct
head(final_df)

#Write into file
write.table(final_df, file= paste(seq, "_AA_variant_info.txt"), row.names = F, sep ="\t")

#Visualize
final_df %>% gt()


#summary table

main.title <- paste("Summary of AA mutations in", seq)
mutation_summary <- ggtexttable(final_df, rows = NULL, theme = ttheme("light", padding = unit(c(3, 3), "mm"), tbody.style = tbody_style(fill = "white", size = 9), colnames.style = colnames_style(fill = "white", size = 10, color = "Black", rot=90)))
mutation_summary1 <- mutation_summary %>% tab_add_title(text = main.title, face = "bold",  padding = unit(5, "line"))

mutation_summary1

#Save into a pdf file

grid <- grid.arrange(mutation_summary1, nrow=1)
ggsave(grid, file=report_name, height = 15 , width = 12, device = pdf,  dpi = 300) ## size of pdf file can be changed depending upon the no. of variants







