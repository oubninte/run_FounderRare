

##**********************************************************************************
#*****************Distribution of Smsg ***************************
#********************************************************************************

rm(list = ls());

library(dplyr)
library(FounderRare)
library(GenomicRanges)
library(ggplot2)
library(hrbrthemes)
library(cowplot)

#Chemin.SmsgH0="/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/Real_Vs_Sims/New_RealData2022/Smsg_H0/Smsg/"
Chemin.SmsgH0="/lustre09/project/6033529/genealogy_sims/results/Samir/Real_Vs_Sims/Smsg_H0_B/Smsg/"
Chemin.SmsgRealData="/lustre09/project/6033529/genealogy_sims/results/Samir/Real_Vs_Sims/New_RealData2022/Smsg_RealData/Smsg/"
Chemin.results="/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/Real_Vs_Sims/New_RealData2022/"
#******************************* Import data **********************************

# Build dataset with different distributions 

Smsg.global=read.csv2(paste0(Chemin.SmsgH0,"Smsg_chr10rep100.csv"), header=T)
Smsg.global=Smsg.global[0,]

for (chr in 1:22) {
  for (rep in 1:100) {
    SmsgH0.chrom =read.csv(paste0(Chemin.SmsgH0,"Smsg_chr",chr,"rep",rep,".csv"), header=T)
    SmsgH0.chrom$Position=(SmsgH0.chrom$End+SmsgH0.chrom$End)/2
    SmsgH0.chrom$CAT=paste0("H0.rep", rep)
    Smsg.global=rbind(Smsg.global, SmsgH0.chrom)
  }
  SmsgRealD.chrom =read.csv(paste0(Chemin.SmsgRealData,"Smsg_chr",chr,".csv"), header=T)
  SmsgRealD.chrom$Position=(SmsgRealD.chrom$Start+SmsgRealD.chrom$End)/2
  SmsgRealD.chrom$CAT="RealData"
  Smsg.global=rbind(Smsg.global, SmsgRealD.chrom)
}

Smsg.global=Smsg.global[Smsg.global$Smsg>0,]
Smsg.global=Smsg.global %>%
  group_by(chr, SG, CAT) %>%
  slice_max(Smsg, n = 1, with_ties = FALSE) 

# Calculate the summary statistics for 'length' within each category
summary_data.global <- Smsg.global %>%
  mutate(case=if_else(CAT=="RealData", "Real.Data", "Sims.Data") )%>%
  group_by(case) %>%
  summarise(Num = length(Smsg ),
            Min = min(Smsg ),
            Q1 = quantile(Smsg , 0.25),
            Median = median(Smsg ),
            Mean = mean(Smsg ),
            Q3 = quantile(Smsg , 0.75),
            Max = max(Smsg ))


summary_data.global.byRep <- Smsg.global %>%
  mutate(case=if_else(CAT=="RealData", "Real.Data", CAT) )%>%
  group_by(case) %>%
  summarise(Num = length(Smsg ),
            Min = min(Smsg ),
            Q1 = quantile(Smsg , 0.25),
            Median = median(Smsg ),
            Mean = mean(Smsg ),
            Q3 = quantile(Smsg , 0.75),
            Max = max(Smsg ))
# Print the summary data all segment
print(summary_data.global)
print(summary_data.global.byRep)

#******************************* Visualisation **********************************
Randomrep=c(paste0("H0.rep", sample(1:100, 10)), "RealData")
colorscat=generate_random_colors((length(Randomrep)))

p.Smsg.AllChr.box <- Smsg.global[Smsg.global$CAT %in% Randomrep, ] %>%
  ggplot( aes(x=CAT, y=Smsg, fill=CAT)) +
  geom_violin()+
  geom_boxplot( color=generate_random_colors(1), alpha=0.6, position = 'identity', show.legend = F, width=.15) +
  theme_ipsum() +
  #labs(title = "Smsg Distribution Across All Chromosomes: A Comparative of Ten Random Replicates and Real Data")+
  theme(plot.title = element_text(hjust = 0.5,face = "plain"),plot.margin = margin(1.5,1.5,1.5,1.5),axis.title.x = element_text(face="bold"), # Bold the x-axis label
        axis.title.y = element_text(face="bold"))+
  labs(x = "Data")

p.Smsg.AllChr.hist=Smsg.global[Smsg.global$CAT %in% Randomrep, ] %>%
  ggplot( aes(x=Smsg, fill=CAT)) +
  geom_histogram( color=generate_random_colors(1), alpha=0.6, position = 'identity', bins = 18) +
  scale_fill_manual(values=colorscat) +
  #geom_histogram(binwidth = 30)+
  theme_ipsum() +
  #labs(title =paste0("Smsg Distribution Across All Chromosomes: A Comparative of Ten Random Replicates and Real Data") )+
  theme(plot.title = element_text(hjust = 0.5,face = "plain"),plot.margin = margin(1.5,1.5,1.5,1.5),axis.title.x = element_text(face="bold"), # Bold the x-axis label
        axis.title.y = element_text(face="bold"))+
  labs(x = "Smsg") 

p.Smsg.AllChr.vio <- Smsg.global[Smsg.global$CAT %in% Randomrep, ] %>%
  ggplot( aes(x=CAT, y=Smsg, fill=CAT)) +
  geom_violin()+
  geom_boxplot( color=generate_random_colors(1), alpha=0.6, position = 'identity', show.legend = F, width=.15) +
  theme_ipsum() +
  #labs(title = "Smsg Distribution Across All Chromosomes: A Comparative of Ten Random Replicates and Real Data")+
  theme(plot.title = element_text(hjust = 0.5,face = "plain"), plot.margin = margin(1.5,1.5,1.5,1.5),axis.title.x = element_text(face="bold"), # Bold the x-axis label
        axis.title.y = element_text(face="bold"))+
  labs(x = "Data")

plot = plot_grid(p.Smsg.AllChr.box, p.Smsg.AllChr.hist , ncol = 1 ) #p.Smsg.AllChr.vio

# Add a single title to the combined plot
title <- ggdraw() + draw_label("Smsg Distribution Across All Chromosomes: A Comparative of Ten Random Replicates and Real Data", 
                               fontface='bold', x = 0.5, hjust = 0.5, vjust = 1.5) # Center the title and make it bold
combined_plot<- plot_grid(title, plot, ncol=1, rel_heights=c(0.05, 0.95)) # Combine title and plot
combined_plot

ggsave(filename = paste0(Chemin.results, "Smsg_Distribution.png"), plot = combined_plot,   width = 16, height = 20, dpi = 200)



# manhattan plot
library(qqman)
Smsg.global$SNP=paste0(Smsg.global$SG, Smsg.global$cd)

maxH0=max(Smsg.global$Smsg[Smsg.global$CAT!="RealData"])
maxReal=max(Smsg.global$Smsg[Smsg.global$CAT=="RealData"])
Seuil=c(0.05)

m_by_rep <- Smsg.global[Smsg.global$CAT!="RealData",] %>%
  group_by(CAT) %>%
  summarise(m_max = max(Smsg, na.rm = TRUE))

SeuilQ=quantile(m_by_rep$m_max, 1-Seuil)

png(paste0(Chemin.results, "Manhattan_plot_of_Smsg_AllAff.png"), width = 2*800, height = 2*400)  # Specify the file name and dimensions

# Get the maximum x value
max_x <- max(Smsg.global[Smsg.global$CAT=="RealData",]$Position)
# Add some space for the labels
xlim <- c(0, 11*max_x * 1.1)

# Create the plot with the adjusted x limits
manhattan(Smsg.global[Smsg.global$CAT=="RealData",], chr="chr", bp="Position", snp="SNP", p="Smsg" , ylim = c(0,max(maxH0, maxReal)), logp=FALSE, 
          col = c("#69b3a2", "#404080"), xlab="", ylab = "", cex = 0.8, xlim=xlim, suggestiveline=0,   genomewideline=0, annotateTop = TRUE)

abline(h =SeuilQ, col = "red", lty = 2)

# Add labels for each line
text(x = 10*max_x , y = 1.01*SeuilQ[1], labels = paste(Seuil[1]*100, "% threshold"), pos = 4)

# Modify the title and axis labels
title("Visualizing Smsg Across Chromosomes for Real Data", font.main=4, cex.main=1) # Reduce the size of the title and make it bold
mtext("Chromosome", side=1, line=2.5, font=2) # Bold the x-axis label
mtext("Smsg", side=2, line=2.5, font=2) # Bold the y-axis label

dev.off() 



png(paste0(Chemin.results, "QQplot_Smsg_AllChr.png"), width = 1.3*800, height = 1.3*600)  # Specify the file name and dimensions
qqplot(Smsg.global$Smsg[Smsg.global$CAT!="RealData"], Smsg.global$Smsg[Smsg.global$CAT == "RealData"],
       ylab = "",       xlab = "",        pch = 20,       col = "black")
abline(a = 0, b = 1, col = "red")

title("QQPlot Comparing Smsg Across All Chromosomes: Real Data vs. Simulated Data", font.main=4, cex.main=0.85) # Reduce the size of the title and make it bold
mtext("Simulated Data", side=1, line=2.5, font=2) # Bold the x-axis label
mtext("Real Data", side=2, line=2.5, font=2) # Bold the y-axis label

dev.off() 


#QQplot with the size of point reflecting the count of 
#the size of the points in the plot corresponds to the number of observations at each location

# Sort the data
realdata <- sort(Smsg.global$Smsg[Smsg.global$CAT=="RealData"])

# Calculate the quantiles of the data
realdata_quantiles <- quantile(realdata, probs = seq(0, 1, length.out = length(realdata)))

# Generate the theoretical quantiles
Simadata <- sort(Smsg.global$Smsg[Smsg.global$CAT!="RealData"])

# Calculate the quantiles of the data
Simadata_quantiles <- quantile(Simadata, probs = seq(0, 1, length.out = length(realdata)))

library(ggplot2)

# Assuming Simadata_quantiles and realdata_quantiles are your data
df <- data.frame(Simadata_quantiles, realdata_quantiles)

ggplot(df, aes(x = Simadata_quantiles, y = realdata_quantiles)) +
  geom_count() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Add line y=x
  guides(size = "none") +  # Hide the legend
  theme_minimal() +
  labs(x = "Simulated Data", y = "Real Data", title = "QQPlot Comparing Smsg Across All Chromosomes: Real Data vs. Simulated Data")

#signal ...
realdf=Smsg.global[Smsg.global$CAT=="RealData",]
tmp=realdf[realdf$Smsg==max(realdf$Smsg),]








