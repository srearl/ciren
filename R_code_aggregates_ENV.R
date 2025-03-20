#Merger of Hannah's code and Ro's code to handle aggregates data. Modified AEM 3/12/25. Have not yet populated teh env data....
#load needed packages
packages <- c("tidyr","dplyr","readxl", "readr", "stringr","devtools", "data.table","ggplot2","LaCroixColoR","readxl", 
              "ggpubr", "gridExtra", "RColorBrewer","colorspace","OneR", "lookup",
              "FactoMineR", "factoextra", "gginnards", "cellWise", "corrplot", 
              "vegan", "morphr", "purrr", "imager", "ggrepel", "cowplot", "Nmisc", "grid")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
#load ecotaxa import file from location on computer  (.xls version)
ecotaxa_aggregates <- read_delim("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Aggregates_March_25_raw.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

#pull lines of date that we want
sub1<-ecotaxa_aggregates[,c("object_id","object_depth_min","object_depth_max","object_annotation_hierarchy",
                        "object_area","object_major",	"object_minor","object_esd","sample_tot_vol","acq_sub_part", 
                        "object_mean", "object_mode", "object_min", "object_max", "object_angle", "object_circ.", 
                        "object_kurt", "object_fractal", "object_nb1", "object_symetrievc", "object_fcons", 
                        "object_thickr", "object_elongation", "object_meanpos", "object_sr", 
                        "object_feretareaexc", "object_perimferet", "object_cdexc")]

#label lines of data that we want
names(sub1)=c("object_id","Min_depth","Max_depth","Taxa","area","major","minor",
              "ESD","Tow_Vol","Sub_part","mean", "mode", "min", "max", "angle", "circ.", 
              "kurt", "fractal", "nb1", "symetrievc", "fcons", 
              "thickr", "elongation", "meanpos", "sr", 
              "feretareaexc", "perimferet", "cdexc")

#The "object_ID" (now called Label) typically varies among Ecotaxa users. Generally it contains information about the contains cruise, tow, net, and size fraction
#This information is typically (hopefully!) separated by some special character. The following code infers these splits and labels the information associated with each particle
#You will need to change the label order if your data is included differently.
#sub1$num<-gsub(".*_","",sub1$object_id)
#sub1<-separate(sub1,object_id,
#               c("cruise","moc","net","fraction"), extra="drop", remove=FALSE)

# Reads the Hydrography file
env_data <- read_csv("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Aggregates_MOCNESS_net_hydrography.csv") 
head(env_data)

## Data Handling ####

# Note: Modify the column names and keys in steps below according to the format of your data
# Process data:  
sub1 <- sub1 %>%
  # Split 'object_id' into multiple columns
  separate(object_id, c("cruise", "moc", "net", "fraction"), extra = "drop", remove = FALSE) %>% 
  # Calculate density and abundance.m2
  mutate(density = Sub_part / Tow_Vol,
         abundance.m2 = (Sub_part / Tow_Vol) * 
           (Max_depth - Min_depth))

# Rename columns in environmental data (for simplicity) 
env_data1 <- env_data %>%
  rename(
    #    Temp = `temp`,
    #    Sal = `sal`,
    #    Fluor = `avg fluor`,
    #    O2 = `avg o2`,
    Depth_max = `max_depth`,
    Depth_min = `min_depth`,
    D.N = `D/N`
  ) %>%
  # Add 'n' prefix to Net column to match the format in data_p
  mutate(Cruise = tolower(Cruise),
         net = paste0("n", net), 
         tow = paste0("m", tow)) 

# Merge environmental data with processed data
sub1_p <- sub1 %>%
  left_join(env_data1, by = c("cruise" = "Cruise", "moc" = "tow", "net" = "net")) %>%
  # Remove 'object_' prefix from column names for simplicity
  rename_with(~ gsub("^object_", "", .x)) 

#remove non-living
a<-"not-living"
sub1<-filter(sub1_p,!grepl(a, Taxa))

sub1$area_mm2 <- sub1$"area" * ifelse(sub1$DPI == 4800, 0.00002809, 
                                    ifelse(sub1$DPI == 2400, 0.000112, NA))
sub1$major_mm <- sub1$"major" * ifelse(sub1$DPI == 4800, .0053, 
                                      ifelse(sub1$DPI == 2400, .010583, NA))
sub1$minor_mm <- sub1$"minor" * ifelse(sub1$DPI == 4800, .0053, 
                                        ifelse(sub1$DPI == 2400, .010583, NA))
#sub1$feret_mm<-as.numeric(sub1$"feret"*0.010583)
sub1$esd_mm <- sub1$"ESD" * ifelse(sub1$DPI == 4800, .0053, 
                                       ifelse(sub1$DPI == 2400, .010583, NA))
sub1$vol<-(4/3)*pi*((sub1$minor_mm*0.5)^2)*(sub1$major_mm/2)
sub1$hdif<-(sub1$Max_depth-sub1$Min_depth)
sub1$split<-(1/(sub1$Sub_part))
sub1$cruise_moc_net<-as.factor(paste(sub1$cruise,sub1$moc,sub1$net, sep="_"))



#The split information can be added separately and indexed using these lines
#scan_split_index <- read_excel("C:/Users/amy.maas/BIOS Dropbox/Amy Maas/ZoopGroup_LAB/Projects/EXPORTS/Cruise 2018/Zooscan/R analysis 2024/scan_split_index.xlsx")
#sub1$tow_net_bin<-as.factor(paste(sub1$tow_net,sub1$fraction, sep="_"))
#sub1$sub_part_correct<-as.factor(as.character(with(lookup, scan_split_index$Acq_sub[match(sub1$tow_net_bin,scan_split_index$Fraction)])))
#sub1$split_correct<-(1/as.numeric((sub1$sub_part_correct)))


sub1$DW_mg<-(0.055*(sub1$vol))
#These are teh dw to volume conversions from Maas et al. 2021 from BATS.
#If you are in a different ocean basin or want different volume or mass estimates you need to change
sub1$DW_mg<-with(sub1, ifelse(Taxa %like% "Calanoida", 0.055*vol,
                           ifelse(Taxa %like% "Chaetognatha", 0.013*vol,
                                  ifelse(Taxa %like% "Ostracoda", 0.052*vol,
                                         ifelse(Taxa %like% "Thecosomata",0.1913*vol,
                                                ifelse(Taxa %like% "Amphipoda",0.034*vol,
                                                       ifelse(Taxa %like% "Euphausiacea",0.027*vol,
                                                              ifelse(Taxa %like% "Poecilostomatoida", 0.074*vol,
                                                                     ifelse(Taxa %like% "Foraminifera",0.142*vol,
                                                                            ifelse(Taxa %like% "Decapoda",0.034*vol,
                                                                                   0.055*vol))))))))))
#This makes everything a copepod breathing at depth where the temperature is 15*C.
sub1$O2_umol<-(exp(-0.339+(0.801*log(sub1$DW)))+0.069*(15))/22.4 
#This makes everything a copepod breathing at local depth/temp. You need to have loaded the temperature into the info above to use this.
#sub2$O2_umol<-(exp(-0.339+(0.801*log(sub2$DW)))+0.069*(sub2$temp))/22.4 
#This takes into account taxonomy breathing at local depth/temp You need to have loaded the temperature into the info above to use this.
#sub1$O2_umol<-with(sub1, ifelse(Taxa %in% "Copepoda", ((exp(-0.399+(0.801*log(DW)))+0.069*(temp))/22.4),
#                                ifelse(Taxa %in% "Chaetognatha", ((exp(-0.173+(0.805*log(DW)))+0.068*(temp))/22.4),
#                                       ifelse(Taxa %in% "Amphipoda", ((exp(0.407+(0.743*log(DW)))+0.037*(temp))/22.4),
#                                              ifelse(Taxa %in% "Euphausiacea",((exp(0.392+(0.753*log(DW)))+0.046*(temp))/22.4),
#                                                     ifelse(Taxa %in% "Mollusca",((exp(-0.56+(0.82*log(DW)))+0.046*(temp))/22.4),
#                                                            ((exp(-0.399+(0.801*log(DW)))+0.069*(temp))/22.4)))))))

sub1$CO2_umol<-(sub1$O2_umol)*0.87 #Converts between O2 and CO2 using a general RQ 
#sub1$CO2<-with(sub1, ifelse(Taxa %in% "Copepoda", O2_umol*0.87,
#                            ifelse(Taxa %in% "Chaetognatha", O2_umol*1.35,
#                                   ifelse(Taxa %in% "Amphipoda", O2_umol*1.35,
#                                          ifelse(Taxa %in% "Euphausiacea",O2_umol*1.35,
#                                                 ifelse(Taxa %in% "Mollusca",O2_umol*0.94,
#                                                        O2_umol*0.87))))))

##We need to build these functions Euphausiacea, Pleuromamma, Thecosomata. 
#I told it to estimate no poop from everything else so we can get a fast summary
sub1$FP<-with(sub1, ifelse(Taxa %in% "Euphausiacea", 1*(#EXPERIMENTAL RESULTS),
                           ifelse(Taxa %in% "Pleuromamma", O2_umol*1.35,
                                   ifelse(Taxa %in% "Thecosomata", O2_umol*1.35,
                                          0)))))

summary(sub1)

write.csv(sub1,file=paste
          ("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Aggregates_March25_live.csv", sep=""),row.names=F)


#To include all of the selected group please run: 
#Pleuromamma<-M_filt%>%filter(str_detect(as.character(cruise),x))
#To exclude all of the selected group please run: 
M_filt<-sub1%>%filter(!grepl("ae2306", sub1$cruise))


B<-as.character(c(seq(0.25, 27.75, by=0.25)))
M_filt$bin<-cut(M_filt$esd_mm, breaks=c(seq(0.25, 28, by=0.25)), labels=B)
M_filt$bin<-factor(M_filt$bin)
M_filt$Nbin<-as.numeric(as.character(M_filt$bin))
M_filt$net<-as.factor(M_filt$net)

summary(M_filt)
write.csv(M_filt,file=paste
          ("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Aggregates_March25_live_bin.csv", sep=""),row.names=F)

summary<-M_filt%>%group_by(cruise,moc,net,fraction,bin)%>%
  summarize(count=n(),depth=(mean(Min_depth)+mean(Max_depth))/2,max_depth=(Max_depth),hdif=median(hdif),split=median(split),
            freq=count/split,tv=mean(Tow_Vol), Density_m3=freq/tv, bin2=mean(Nbin),
            NBV_m3=(sum(vol)/split/tv),BM_m3=(sum(DW)/split/tv), oxy_m3=(sum(O2_umol))/tv/split,
            CO2_m3=(sum(CO2)/tv/split),Abundance_m2=(freq/tv*hdif),NBV_m2=NBV_m3*hdif,
            BM_m2=BM_m3*hdif, oxy_m2=oxy_m3*hdif,CO2_m2=CO2_m3*hdif) 

write.csv(summary,file=paste
          ("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Aggregates_March25_bins.csv", sep=""),row.names=F)


##Analysis for Yuuki Fecal Pellets

# FILTER FOR CLASSIFICATION LEVEL; CAN INSERT ANY VARIATION OF CLASSIFICATION LEVELS
# this is where you can run species specific information or just live/dead.
# Note that this filter is case sensitive and inclusive (living will also pick up non-living)
# For this analysis we are doing Euphausiacea, Pleuromamma, Thecosomata. For Thecosomata we should put a size threshold....
#x<-"Pleuromamma"
#To include all of the selected group please run: 
#Pleuromamma<-M_filt%>%filter(str_detect(as.character(Taxa),x))
#To exclude all of the selected group please run: 
#Taxa<-filter(Tow_for_R,!grepl(x, Taxa))

## Create a folder name for this run (e.g. Euphausiacea, Pleuromamma, Thecosomata)
#descr<-("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Pleuromamma")
#dir.create(descr)
#write.csv(summary,file=paste
          ("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Pleuromamma/Pleuromamma_bins.csv", sep=""),row.names=F)

###PLOTTING

M_filt <- read_delim("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Aggregates_Aug24_bins.csv")
#I don't know how to do all the plots at once - I only know how to filter for each cruise and run independently
x<-"m22"
sub3<-M_filt%>%filter(str_detect(as.character(moc),x))

#save the specific cruise
write.csv(sub3,file=paste
          ("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Aggregates_m22_bins.csv", sep=""),row.names=F)

### BIOMASS SuMMARY
ggplot(data=(sub3), aes(x=net, y=BM_m3, fill="#D72000", alpha=fraction))+geom_bar(stat="identity")+
  scale_fill_manual(values=c("#D72000"))+
  labs(x="Net", y=expression("Biomass (Dry Weight)" ~(mg~m^-3)), title=paste(x,"Total Biomass by Net"))+
  theme(plot.title=element_text(face="bold",hjust=0.5))+
  scale_alpha_discrete(range=c(0.5,1))+
  guides(fill=FALSE)+
  ylim(0,50000)
ggsave(filename=paste(x,"_DW3_bynet.png", sep="_"),path="C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/", width=6, height=6, units="in", dpi=300)

### OXYGEN USE SUMMARY (You will want to change this to FP mg m-2)
ggplot(data=(summary), aes(x=net, y=oxy_m2, fill="#FFAD0A", alpha=fraction))+geom_bar(stat="identity")+
  scale_fill_manual(values="#FFAD0A")+
  labs(x="Net", y=expression(mu*mol~O[2]*m^-2*h^-1), title=paste(H,"Apparent Oxygen Utilization by Net"))+
  theme(plot.title=element_text(face="bold",hjust=0.5, size=11))+
  scale_alpha_discrete(range=c(0.5,1))+
  guides(fill=FALSE)+
  ylim(0,8000)
ggsave(filename=paste(Tow,"AOU_bynet.png", sep="_"),path="C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/", width=6, height=6, units="in", dpi=300)


## DAY-NIGHT PLOTS #############################
###############################################


#SET YOUR PATH HERE, PICKING THE DAY OF A PAIR
day_bin<-read.csv("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Aggregates_m22_bins.csv")
#SET YOUR PATH HERE, PICKING THE NIGHT OF A PAIR
night_bin<-read.csv("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/Aggregates_m23_bins.csv")

I="June 2023" #SET THE NAME OF THE PAIR

##Create a new folder for the pair
dir.create(paste("C:/Users/amy.maas/Desktop/MOC_Ecotaxa_Analysis/Aggregates_July24/PairedNets/",I, sep=""))

#FILTER BIN SIZES FOR BOTH DAY AND NIGHT TO SELECTED RANGE
#day_bin_sub<-day_bin%>%filter(bin2 >=1&bin2<=20)
day<-day_bin%>%group_by(net)%>%summarize(Tot_BV_m3=sum(NBV_m3),Tot_BM_m3=sum(BM_m3), Tot_Ox_m3=sum(oxy_m3),
                                             Tot_CO2_m3=sum(CO2_m3),Tot_BV_m2=sum(NBV_m2),
                                             Tot_BM_m2=sum(BM_m2), Tot_Ox_m2=sum(oxy_m2),
                                             Tot_CO2_m2=sum(CO2_m2),med_depth=mean(depth))
#night_bin_sub<-night_bin%>%filter(bin2 >=1&bin2<=20)
night<-night_bin_sub%>%group_by(net)%>%summarize(Tot_BV_m3=sum(NBV_m3),Tot_BM_m3=sum(BM_m3), Tot_Ox_m3=sum(oxy_m3),
                                                 Tot_CO2_m3=sum(CO2_m3),Tot_BV_m2=sum(NBV_m2),
                                                 Tot_BM_m2=sum(BM_m2), Tot_Ox_m2=sum(oxy_m2),
                                                 Tot_CO2_m2=sum(CO2_m2),med_depth=mean(depth))


# If you are missing a net (like in M11) use this text to block it out
# day<-filter(day,as.factor(net)!='n8')
# day$net<-as.factor(day$net)
# night<-filter(night,as.factor(net)!='n8')
# night$net<-as.factor(night$net)

#This is the standard code, but use it always
DayNight<-data.frame("BV_Mig"=c(abs(day$Tot_BV_m2-night$Tot_BV_m2)), "BV_Res"=do.call(pmin,(as.data.frame(cbind(day$Tot_BV_m2,night$Tot_BV_m2)))),
                     "BM_Mig"=c(abs(day$Tot_BM_m2-night$Tot_BM_m2)), "BM_Res"=do.call(pmin,(as.data.frame(cbind(day$Tot_BM_m2,night$Tot_BM_m2)))),
                     "Ox_Mig"=c(abs(day$Tot_Ox_m2-night$Tot_Ox_m2)), "Ox_Res"=do.call(pmin,(as.data.frame(cbind(day$Tot_Ox_m2,night$Tot_Ox_m2)))),
                     "CO2_Mig"=c(abs(day$Tot_CO2_m2-night$Tot_CO2_m2)),"CO2_Mig"=do.call(pmin,(as.data.frame(cbind(day$Tot_CO2_m2,night$Tot_CO2_m2)))),
                     "BM_DVM"=c((day$Tot_BM_m2-night$Tot_BM_m2)),"BM_day"=c(day$Tot_BM_m2),"BM_night"=c(night$Tot_BM_m2),
                     "BV_DVM"=c((day$Tot_BV_m2-night$Tot_BV_m2)),
                     "Med_Depth"=as.factor(apply(as.data.frame(cbind(day$med_depth,night$med_depth)),1,FUN=mean)),
                     "Net"=c(day$net))
# This is the standard code, skip below if you have a missing net
DN2<-data.frame("Net"=as.factor(c(1:8,1:8)),"M_R"=as.factor(c(rep.int("M",8), rep.int("R",8))),"BM"=c(DayNight$BM_Mig,DayNight$BM_Res), "BV"=c(DayNight$BV_Mig,DayNight$BV_Res),
                "Ox"=c(DayNight$Ox_Mig,DayNight$Ox_Res), "CO2"=c(DayNight$CO2_Mig,DayNight$CO2_Res),
                "Med_Depth"=as.factor(rep.int(DayNight$Med_Depth, 2)))

# Use this code if you are skipping a net
#DN2<-data.frame("Net"=as.factor(c(1:7,1:7)),"M_R"=as.factor(c(rep.int("M",7), rep.int("R",7))),"BM"=c(DayNight$BM_Mig,DayNight$BM_Res), "BV"=c(DayNight$BV_Mig,DayNight$BV_Res),
#               "Ox"=c(DayNight$Ox_Mig,DayNight$Ox_Res),"CO2"=c(DayNight$CO2_Mig,DayNight$CO2_Res),
#              "Med_Depth"=as.factor(rep.int(DayNight$Med_Depth, 2)))

write.csv(DN2,file=paste("Output/PairedNets/",I,"/",I,"_DayNightNETS.csv",sep=""),row.names=FALSE)

d_labs=c("900","700","550","400","300","200","50","0")
## 3/4    c("900","750","600","450","300","225","50","0")
## 8/9    c("850","700","550","400","275","200","50","0")
## 10/11  c("850","700","550","400","250","175","50","0")
## 12/13  c("900","700","550","400","300","200","50","0")


plot1<-ggplot(data=DN2, aes(x=Med_Depth, y=BM,fill="grey7", alpha=M_R))+
  scale_fill_manual(values=("grey7"))+
  scale_alpha_discrete(range=c(0.5,1), labels=c("Migratory","Resident"), drop=FALSE)+
  geom_col(position="stack", na.rm=FALSE)+coord_flip()+
  labs(x="Minimum Net Depth (m)", y=expression("mg Biomass"~m^-2),
       title="Dry Weight Biomass", alpha="")+
  theme(plot.title=element_text(face="bold",hjust=0.5, size=12))+
  guides(fill=FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=12),axis.text.y= element_text(size=12, vjust=-1),
        axis.ticks.y=(element_blank()),
        axis.title=element_text(size=14),strip.text=element_text(size=14, face="bold"),
        legend.text=element_text(size=11),
        legend.position='bottom')+
  scale_x_discrete(limits=c("1000",rev(levels(DN2$Med_Depth))),labels=rev(c(rev(d_labs),"1000")), drop=FALSE)+
  scale_y_continuous(limits=c(0,500))
plot2<-ggplot(data=DN2, aes(x=Med_Depth, y=Ox, fill="grey7", alpha=M_R))+
  scale_fill_manual(values=("grey7"))+scale_alpha_discrete(range=c(0.5,1), labels=c("Migratory","Resident"))+
  geom_col(position="stack")+coord_flip()+
  labs(x="Minimumn Net Depth (m)", y=expression(mu*mol~O[2]~m^-2*h^-1),
       title="Apparent Oxygen Usage", alpha="")+
  guides(fill=FALSE, alpha=FALSE)+
  theme(plot.title=element_text(face="bold",hjust=0.5, size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=12),axis.text.y= element_text(size=12, vjust=-1),
        axis.ticks.y=(element_blank()),
        axis.title=element_text(size=14),strip.text=element_text(size=14, face="bold"),
        legend.text=element_text(size=11),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(c(0,4000))+
  scale_x_discrete(limits=c("1000",rev(levels(DN2$Med_Depth))),labels=rev(c(rev(d_labs),"1000")), drop=FALSE)

plot3<-ggplot(data=DayNight, aes(x=Med_Depth, y=BM_DVM, fill="grey7"))+
  geom_col()+coord_flip()+
  labs(x="Minimum Net Depth (m)", y=expression("Biomass"~mg~m^-2),
       title="Biomass (Dry Weight)")+
  scale_fill_manual(values="grey7")+
  theme(plot.title=element_text(face="bold",hjust=0.5, size=12))+
  scale_x_discrete(limits=c("1000",rev(levels(DN2$Med_Depth))),labels=rev(c(rev(d_labs),"1000")), drop=FALSE)+
  ylim(c(-300,100))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=12),axis.text.y= element_text(size=12, vjust=-1),
        axis.ticks.y=(element_blank()),
        axis.title=element_text(size=14),strip.text=element_text(size=14, face="bold"),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill=FALSE)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(plot1)

p1<-grid.arrange(arrangeGrob(plot1+theme(legend.position="none"),plot2,plot3, ncol=3, top=paste(I,"Migrators and Resident Zooplankton")),
                 mylegend, nrow=2,heights=c(10, 1))

ggsave(p1,filename=paste(I,"DNPlot.png", sep="_"),path=paste("Output/PairedNets/",I, sep=""), width=8, height=7, units="in", dpi=300)


### MIGRATION HEATMAPS #########################################################
################################################################################

#SET YOUR PATH HERE, PICKING THE DAY OF A PAIR
day2<-read.csv("Output/NetFilter/M13_delT_allTaxa/M13_bins.csv")
#SET YOUR PATH HERE, PICKING THE NIGHT OF A PAIR
night2<-read.csv("Output/NetFilter/M12_delT_allTaxa/M12_bins.csv")

# M3 July 2016 (Night) 
# M4 July 2016 (Day)       
# M8 July 2017 (Day)
# M9 July 2017 (Night)
# M10 July 2018 (Day)
# M11 July 2018 (Night)
# M12 Oct 2018 (Night)
# M13 Oct 2018 (Day)

J<-"Oct 2018" #SET THE NAME OF THE PAIR
#SET THE DEPTH INTERVALS FOR THE PAIR
d_labs=c("900","700","550","400","300","200","50","0")
## 3/4    c("900","750","600","450","300","225","50","0")
## 8/9    c("850","700","550","400","275","200","50","0")
## 10/11  c("850","700","550","400","250","175","50","0")
## 12/13  c("900","700","550","400","300","200","50","0")


day2$bin<-as.factor(day2$bin)
night2$bin<-as.factor(night2$bin)
dn.hm<-merge(day2,night2, by.x=c("net","bin","fraction"),by.y=c("net","bin","fraction"),all.x=TRUE,all.y=TRUE)
head(dn.hm)
summary(dn.hm)
dn.hm[is.na(dn.hm)]<-0
dn.hm2<-dn.hm%>%group_by(net,bin)%>%summarize(BVm3_day=sum(NBV_m3.x), BVm3_night=sum(NBV_m3.y),
                                              BVm2_day=sum(NBV_m2.x), BVm2_night=sum(NBV_m2.y),
                                              BMm2_day=sum(BM_m2.x), BMm2_night=sum(BM_m2.y),
                                              Oxm2_day=sum(oxy_m2.x),Oxm2_night=sum(oxy_m2.y),
                                              CO2m2_day=sum(CO2_m2.x),CO2m2_night=sum(CO2_m2.y))

#write.csv(dn.hm2,file=paste("Output/PairedNets/",I,"/",J,"_DayNight_dnhm2.csv",sep=""),row.names=FALSE)

dn.hm3<-dn.hm2%>%mutate(BV_m3=(BVm3_day-BVm3_night), BV_m2=(BVm2_day-BVm2_night),BM_m2=(BMm2_day-BMm2_night),
                        Ox=(Oxm2_day-Oxm2_night),CO2=(CO2m2_day-CO2m2_night)) %>%select(,c(1:2,13:17))  

#write.csv(dn.hm3,file=paste("Output/PairedNets/",I,"/",J,"_DayNight_dnhm3.csv",sep=""),row.names=FALSE)
dn.hm3<-as.data.frame(dn.hm3)
dn.hm3$binnum<-as.numeric(as.character(dn.hm3$bin))
dn.hm4<-filter(dn.hm3,binnum>0.01&binnum<100)

#Use these for July 2018 to drop the data from net 8
#dn.hm4<-filter(dn.hm4, net!='n8')
#dn.hm4$net<-as.factor(dn.hm4$net)
#levels(dn.hm4$net)<-c(levels(dn.hm4$net),'n8')


write.csv(dn.hm4,file=paste("Output/PairedNets/",I,"/",J,"_DayNightBINS.csv",sep=""),row.names=FALSE)


BV.heatmap<-ggplot(data=dn.hm4, mapping=aes(x=bin, y=net, fill=BV_m2, color=""))+
  geom_tile()+
  xlab(label="Size Class (mm^3)")+
  ylab(label="Minimum Net Depth (m)")+
  scale_fill_continuous_divergingx(na.value="gray45",limits=c(-1000,1000),palette = 'RdBu',
                                   rev=TRUE, mid =0, l3 = 0, p1=0.4, p3 = .4, p4 = .5)+
  scale_x_discrete(breaks=c("0.01", "0.1", "1","10", "100"), 
                   labels=c("0.01", "0.1", "1","10", "100"), drop=TRUE)+
  coord_cartesian(clip='off')+
  scale_y_discrete(labels=d_labs, drop=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.title= element_text(size=14),
        axis.text.y=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust = 0.5,face="bold",size=16))+
  ggtitle(paste(J,'(Day-Night) Plankton Biovolume Shift'))+
  labs(fill=expression("Biovolume"~(mm^3/m^2)),
       x=expression("Size Class"~mm^3),
       color=expression("<-1000"~mm^3/m^2))+
  scale_colour_manual(values=NA) +   
  guides(fill=guide_colorbar(order=1))+
  guides(color=guide_legend(order=2, override.aes=list(fill="gray45")))
BV.heatmap
ggsave(BV.heatmap,filename=paste(J,"BV_Shift.png", sep="_"),
       path=paste("Output/PairedNets/",I, sep=""), width=8, height=5, units="in", dpi=300)

CO2.heatmap<-ggplot(data=dn.hm4, mapping=aes(x=bin, y=net, fill=CO2, color=""))+
  geom_tile()+
  xlab(label="Size Class (mm^3)")+
  ylab(label="Minimum Net Depth (m)")+
  scale_fill_continuous_divergingx(na.value="gray45",limits=c(-80,80),palette = 'RdBu',
                                   rev=TRUE, mid =0, l3 = 0, p1=0.4, p3 = .4, p4 = .5)+
  scale_x_discrete(breaks=c("0.01", "0.1", "1","10", "100"), 
                   labels=c("0.01", "0.1", "1","10", "100"), drop=TRUE)+
  coord_cartesian(clip='off')+
  scale_y_discrete(labels=d_labs, drop=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.title= element_text(size=14),
        axis.text.y=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust = 0.5,face="bold",size=16))+
  ggtitle(paste(J,'(Day-Night) Plankton Respiratory Shift'))+
  labs(fill=expression("CO2"~(mu*mol~m^-2*h^-1)),
       x=expression("Size Class"~mm^3),
       color=expression("<-80"~mu*mol~m^-2*h^-1))+
  scale_colour_manual(values=NA) +   
  guides(fill=guide_colorbar(order=1))+
  guides(color=guide_legend(order=2, override.aes=list(fill="gray45")))
CO2.heatmap
ggsave(CO2.heatmap,filename=paste(J,"CO2_Shift.png", sep="_"),
       path=paste("Output/PairedNets/",I, sep=""), width=8, height=5, units="in", dpi=300)



