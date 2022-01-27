data<-read.csv("210416_motif_time_25Hz_v5.csv")
dim(data)
#[1] 1197824      41

summary(data[,1:11])
#mouse: 2-126, 32 mice
#trial: 32 trial
#genotype:2
#sex:2
#community:A~H
#motif: 0~29, 30 motifs

# Mouse            Trial          Genotype           Sex     
# Min.   :  2.00   Min.   : 1.00   Min.   :0.0000   Min.   :0.0  
# 1st Qu.: 31.75   1st Qu.: 8.75   1st Qu.:0.0000   1st Qu.:0.0  
# Median : 70.00   Median :16.50   Median :1.0000   Median :0.5  
# Mean   : 70.09   Mean   :16.50   Mean   :0.5625   Mean   :0.5  
# 3rd Qu.:109.25   3rd Qu.:24.25   3rd Qu.:1.0000   3rd Qu.:1.0  
# Max.   :126.00   Max.   :32.00   Max.   :1.0000   Max.   :1.0  
# 
# Seconds         Binned150sec     Binned300sec     Community     
# Min.   :   0.04   Min.   : 1.000   Min.   :1.000   A      :283119  
# 1st Qu.: 374.35   1st Qu.: 3.000   1st Qu.:2.000   B      :232705  
# Median : 748.66   Median : 5.000   Median :3.000   E      :172681  
# Mean   : 748.66   Mean   : 5.499   Mean   :2.999   F      :167272  
# 3rd Qu.:1122.97   3rd Qu.: 8.000   3rd Qu.:4.000   D      :128700  
# Max.   :1497.28   Max.   :10.000   Max.   :5.000   C      :123196  
# (Other): 90151  
# Motif         Speedbelly          SpeedNose       
# Min.   : 0.00   Min.   :0.0000000   Min.   :0.000000  
# 1st Qu.: 8.00   1st Qu.:0.0007889   1st Qu.:0.002275  
# Median :15.00   Median :0.0026887   Median :0.007650  
# Mean   :14.48   Mean   :0.0062263   Mean   :0.010976  
# 3rd Qu.:22.00   3rd Qu.:0.0074159   3rd Qu.:0.015924  
# Max.   :29.00   Max.   :0.2878697   Max.   :0.537502  
sum( data[2:nrow(data),5]-data[1:(nrow(data)-1),5] <0)
#[1] 958
length(unique(data[,1]))
#[1] 32 #time not ordered?

colnames(data)
#[1] "Mouse"        "Trial"        "Genotype"     "Sex"          "Seconds"     
#[6] "Binned150sec" "Binned300sec" "Community"    "Motif"        "Speedbelly"  
#[11] "SpeedNose" 

data=data[order(data[,1],data[,5]),]

mouse_type=table(data[,1])

for(m in sort(names(mouse_type))){
  
  data_sub=data[data[,1]==m,]
  
  index_end=c(which(data_sub[2:nrow(data_sub),9]!=data_sub[1:(nrow(data_sub)-1),9]),nrow(data_sub))
  index_end=sort(unique(index_end))
  
  index_start=c(1,index_end[-length(index_end)]+1)
  index_start=sort(unique(index_start))
  
  data_sub_end=data_sub[index_end,]
  data_sub_start=data_sub[index_start,]
  which(data_sub_end[,5]!=data_sub_start[,5])[1]
  head(index_start,30)
  head(index_end,30)
  cbind(head(data_sub[index_start,c(1,5,9)],30),head(data_sub[index_end,c(1,5,9)],30))
  cbind(tail(data_sub[index_start,c(1,5,9)],30),tail(data_sub[index_end,c(1,5,9)],30))
  
  lingering_time=data_sub[index_end,5]-data_sub[index_start,5]
  lingered_motif=data_sub[index_start,9]
  summary(lingering_time)
  
  #Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  #-1496.680     0.080     0.160     0.271     0.320    20.760 
  sum(lingering_time<0)
  
  
  #per motif per mouse
  
  unique_lingered_motif=sort(unique(lingered_motif))
 
  pdf(paste0("mouse",m,"_perMotif_perMouse_lingering_time_hist.pdf"))
  for(i in 1:length(unique_lingered_motif)) {
    print(hist(lingering_time[lingered_motif==unique_lingered_motif[i]],breaks=100,main=unique_lingered_motif[i]))
  
  }
  dev.off()
  
}



cbind((data_sub[(index_start[27]-5):(index_start[27]+5),c(1,5,9)]),(data_sub[(index_start[27]-5):(index_start[27]+5),c(1,5,9)]))





genotype=table(data[,3])

for(m in sort(names(genotype))){
  
  data_sub=data[data[,3]==m,]
  
  index_end=c(which(data_sub[2:nrow(data_sub),9]!=data_sub[1:(nrow(data_sub)-1),9]),nrow(data_sub))
  index_end=sort(unique(index_end))
  
  index_start=c(1,index_end[-length(index_end)]+1)
  index_start=sort(unique(index_start))
  
  data_sub_end=data_sub[index_end,]
  data_sub_start=data_sub[index_start,]
  which(data_sub_end[,5]!=data_sub_start[,5])[1]
  head(index_start,30)
  head(index_end,30)
  cbind(head(data_sub[index_start,c(1,5,9)],30),head(data_sub[index_end,c(1,5,9)],30))
  cbind(tail(data_sub[index_start,c(1,5,9)],30),tail(data_sub[index_end,c(1,5,9)],30))
  
  lingering_time=data_sub[index_end,5]-data_sub[index_start,5]
  lingered_motif=data_sub[index_start,9]
  
  lingered_motif=lingered_motif[-which(lingering_time<0)]
  lingering_time=lingering_time[-which(lingering_time<0)]
  
  summary(lingering_time)
  
  #Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  #-1496.680     0.080     0.160     0.271     0.320    20.760 
  sum(lingering_time<0)
  
  
  #per motif per mouse
  
  unique_lingered_motif=sort(unique(lingered_motif))
  
  pdf(paste0("genotype",m,"_perMotif_perMouse_lingering_time_hist.pdf"))
  for(i in 1:length(unique_lingered_motif)) {
    print(hist(lingering_time[lingered_motif==unique_lingered_motif[i]],breaks=100,main=unique_lingered_motif[i]))
    
  }
  dev.off()
  
}




#per motif per genotype
