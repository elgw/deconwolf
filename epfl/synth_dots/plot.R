library(ggplot2)
library(plyr)

files = c("tab1.csv", "tab2.csv", "tab3.csv", "tab4.csv", "tab5.csv")

tab = data.frame()
for(f in files)
{
  tab1 = read.csv(f, header=T, sep=',');
  tab  <-  rbind(tab, tab1)
}


tabg = read.csv("tab6.csv", header=T, sep=',');
# tabg$DWnDotsIsLMAX by definition this is 100%

summary(tab)
dim(tab)


names(tab)
tab$file = factor(tab$file)
tab$file

fs  <- tab$file
fs = mapvalues(fs, from = c("50_sdots.tif", "dl2_50_sdots.tif", "dl2_100_sdots.tif", "dl2_200_sdots.tif", "sdots.tif"), 
               to=c("DW@50", "Dela2@50", "DeLa2@100", "DeLa2@200", "Synth"))
tab$Dset = fs

P = ggplot(tab, aes(x=N, y=DWnDotsIsLMAX, color=Dset))+geom_line()
P = P + xlab("dots") + ylab("percentage of detected dots")
P = P + labs(colour = "Type")
ggsave('plots/perc_detected.png')
ggsave('plots/perc_detected.pdf')
P

Vnuc = 374
Vsim = 256*256*40*0.13*0.13*0.13
tab$PerNuc = Vnuc/Vsim*tab$N

P = ggplot(tab, aes(x=PerNuc, y=DWnDotsIsLMAX, color=Dset))+geom_line()
P = P + xlab("dots per nuclei") + ylab("percentage of detected dots")
P = P + labs(colour = "Type")
ggsave('plots/perc_detected_nuc.png')
ggsave('plots/perc_detected_nuc.pdf')
P

