library(ggplot2)
library(plyr)
require('purrr')
#install.packages('sjmisc')
#install.packages('glue')
#install.packages('plyr')
library('sjmisc')

# Define files to load
files = c("newtab1.csv", "newtab2.csv", "newtab3.csv", "newtab4.csv", "newtab5.csv", "newtab6.csv", "newtab7.csv", "newtab8.csv", "newtab9.csv", "newtab10.csv")

# Load the files
tab = data.frame()
for(f in files)
{
  tab1 = read.csv(f, header=T, sep=',');
  tab  <-  rbind(tab, tab1)
}

summary(tab)
dim(tab)


names(tab)
tab$file = factor(tab$file)
tab$file

fs  <- tab$file
fs = mapvalues(fs, from = 
c("dl2_50_sdots.tif", "dl2_100_sdots.tif", "dl2_200_sdots.tif", "sdots.tif", "25_sdots.tif", "50_sdots.tif", "100_sdots.tif", "125_sdots.tif", "150_sdots.tif"), 
  to=c("DeLa2@50", "DeLa2@100", "DeLa2@200", "Synth", "DW@25", "DW@50", "DW@100", "DW@125", "DW@150"))
tab$Dset = fs

getMethod <- function(mstring)
{
#  browser();
  mstr = as.character(mstring);
  if(str_contains(mstr, 'Ground'))
  {
    return('GT');
  }
  if(str_contains(mstr, 'DW'))
  {
    return('DW');
  }
  if(str_contains(mstr, 'DeLa'))
  {
    return('DeLa2');
  }
  if(str_contains(mstr, 'Synth'))
  {
  return('Synth');
  }
  browser();
  return('Unknown');
}


tab$Class = as.factor(map_chr(tab$Dset, getMethod));




P = ggplot(tab, aes(x=N, y=DWnDotsIsLMAX, color=Dset))+geom_line(aes(linetype=Class))
P = P + xlab("dots") + ylab("percentage of detected dots")
P = P + labs(colour = "Type")
ggsave('plots/perc_detected.png')
ggsave('plots/perc_detected.pdf')
P

Vnuc = 374
Vsim = 256*256*40*0.13*0.13*0.13
tab$PerNuc = Vnuc/Vsim*tab$N

P = ggplot(tab, aes(x=PerNuc, y=DWnDotsIsLMAX, color=Dset))+geom_line(aes(linetype=Class))
P = P + xlab("dots per nuclei") + ylab("percentage of detected dots")
P = P + labs(colour = "Type")
ggsave('plots/perc_detected_nuc.png')
ggsave('plots/perc_detected_nuc.pdf')
P

P = ggplot(tab, aes(x=PerNuc, y=DWnDotsIsLMAX, color=Dset))+geom_line(aes(linetype=Class))
P = P + xlab("dots per nuclei") + ylab("percentage of detected dots")
P = P + labs(colour = "Type")+xlim(0, 550)
ggsave('plots/perc_detected_nuc_zoom.png')
ggsave('plots/perc_detected_nuc_zoom.pdf')
P


P = ggplot(tab, aes(x=PerNuc, y=mse, color=Dset))+geom_line(aes(linetype=Class))
P = P + xlab("dots per nuclei") + ylab("MSE")
P = P + labs(colour = "Type")
P
ggsave('plots/mse_detected_nuc.png')
ggsave('plots/mse_detected_nuc.pdf')
P

