# Packages and functions ----------------------------------------------------------------
library(dplyr)
library(tibble)
library(DESeq2)
library(vegan)
library(conflicted)
conflict_prefer_all("dplyr")
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(FUNGuildR)
library(ggh4x)

my_shannon <- function(x) {
  relabund <- x[x>0]/sum(x)
  -sum(relabund*log(relabund))
}
rarefy_and_average <- function(dataframe, target_reads, num_iterations) {
  rarefied_dataframes <- list()
  for (i in 1:num_iterations) {
    rarefied_df <- rrarefy(dataframe, sample = target_reads)
    rarefied_dataframes[[i]] <- rarefied_df
  }
  
  # Calculate average rarefied counts
  average_rarefied_df <- Reduce(`+`, rarefied_dataframes) / num_iterations
  return(average_rarefied_df)
}



# Soil Biodiversity Profiles ----------------------------------------------

## Fungi ------------------------------------------------------------------

### Dataframe setup ---------------------------------------------------------

manag_ITS_data_soil<-read.csv("path/to/fungal_dataset.csv")%>%
  select(-grep("Run1", names(.), value = TRUE))
 
plotsite_metadata_soil<-data.frame(cbind(c(names(manag_ITS_data_soil)[c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2"),
                                                                        which(names(manag_ITS_data_soil)=="Sample_F01_01_gITS7.ITS4__Run3"):which(names(manag_ITS_data_soil)=="Sample_F03_05_gITS7.ITS4__Run3"))]),
                                         c(rep("034", 6),
                                           rep("050", 6),
                                           rep("052", 6),
                                           rep("081", 6),
                                           rep("083", 6),
                                           rep("095", 6),
                                           rep("132", 6),
                                           rep("F01", 6),
                                           rep("F02", 6),
                                           rep("F03", 5)),
                                         c(rep("Reduced", 18),
                                           rep("Intermediate", 18),
                                           rep("Enhanced", 6),
                                           rep("Forest", 17)))) %>%
  rename("Sample" = "X1",
         "Plot" = "X2",
         "Understory" = "X3") %>%
  mutate(Plot=factor(Plot,
                     levels=c("034", "050", "052",
                              "081", "083", "095",
                              "132", "F01", "F02", "F03")),
         Understory=factor(Understory,
                           levels=c("Reduced",
                                    "Intermediate",
                                    "Enhanced",
                                    "Forest"))) %>%
  mutate(Biomass=c(read.csv("path/to/Fieldwork metadata.csv")$Biomass.mean))%>%
  mutate(NMDS_centroid=NA)%>%
  mutate(NMDS_centroid=ifelse(Understory=="Reduced" | Understory =="Intermediate",
                              "Y",
                              "N"))

for (e in c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2"))){
  manag_ITS_data_soil[,e]<-c(manag_ITS_data_soil[,e]-manag_ITS_data_soil$Sample_NegITS_gITS7.ITS4__Run2)
}
for (e in c(which(names(manag_ITS_data_soil)=="Sample_F01_01_gITS7.ITS4__Run3"):which(names(manag_ITS_data_soil)=="Sample_F03_05_gITS7.ITS4__Run3"))){
  manag_ITS_data_soil[,e]<-c(manag_ITS_data_soil[,e]-manag_ITS_data_soil$Sample_ForNeg_gITS7.ITS4__Run3)
}

for (e in c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_F03_05_gITS7.ITS4__Run3"))){
  manag_ITS_data_soil[which(manag_ITS_data_soil[,e]<0),e]<-0
}
for (e in c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_F03_05_gITS7.ITS4__Run3"))){
  manag_ITS_data_soil[which(manag_ITS_data_soil[,e]==1),e]<-0
}

manag_ITS_data_soil<-manag_ITS_data_soil[,
                                         -c(which(names(manag_ITS_data_soil)=="Sample_NegITS_gITS7.ITS4__Run2"),
                                            which(names(manag_ITS_data_soil)=="Sample_ForNeg_gITS7.ITS4__Run3"))]

non_occur_ESV<-c()
for (e in c(1:nrow(manag_ITS_data_soil))){
  non_occur_ESV<-c(non_occur_ESV,
                   sum(manag_ITS_data_soil[e,
                                           c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_F03_05_gITS7.ITS4__Run3"))]))
}
manag_ITS_data_soil<-manag_ITS_data_soil[-c(which(non_occur_ESV==0)),]
manag_ITS_data_soil<-subset(manag_ITS_data_soil,
                            manag_ITS_data_soil$kingdom=="Fungi")


### Rarefaction and normalisation -------------------------------------------
manag_ITS_soil_rarefaction<-manag_ITS_data_soil[,
                                                c(1:which(names(manag_ITS_data_soil)=="Sample_F03_05_gITS7.ITS4__Run3"))] %>%
  remove_rownames(.) %>%
  column_to_rownames(.,
                     var="ESV") %>%
  t(.) %>%
  as.data.frame(.) %>%
  filter(!row_number() %in% c(1:3))%>%
  mutate(across(everything(), as.numeric))

manag_ITS_soil_normal<-
  DESeqDataSetFromMatrix(countData = manag_ITS_soil_rarefaction %>%
                           t() %>%
                           as.data.frame(),
                         colData = plotsite_metadata_soil %>%
                           column_to_rownames("Sample"),
                         design = ~Understory) %>%
  DESeq2::estimateSizeFactors(.,type="poscounts")%>%
  #DESeq(.) %>% 
  DESeq2::counts(.,
                 normalized=TRUE)%>%
  as.data.frame() %>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_ITS_data_soil%>%
              select(-c("Sample_034_01_gITS7.ITS4__Run2":"Sample_F03_05_gITS7.ITS4__Run3")),
            by=c("ESV"))

set.seed(131)
manag_ITS_soil_rarefied<-manag_ITS_soil_rarefaction %>%
  rarefy_and_average(.,
                     min(sapply(manag_ITS_data_soil[c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_F03_05_gITS7.ITS4__Run3"))], 
                                sum)),
                     500) %>% 
  t()%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_ITS_data_soil%>%
              select(-c("Sample_034_01_gITS7.ITS4__Run2":"Sample_F03_05_gITS7.ITS4__Run3")),
            by=c("ESV"))


### Making taxonomic assignment user-friendly -------------------------------
manag_ITS_soil_normal$phylum[which(manag_ITS_soil_normal$phylum=="unidentified")]<-"unclassified_Fungi"
manag_ITS_soil_normal$class[which(manag_ITS_soil_normal$phylum=="unclassified_Fungi")]<-manag_ITS_soil_normal$phylum[which(manag_ITS_soil_normal$phylum=="unclassified_Fungi")] 
manag_ITS_soil_normal$class[c(grep("unidentified", manag_ITS_soil_normal$class))]<-paste("unclassified",
                                                                                         manag_ITS_soil_normal$phylum[c(grep("unidentified", manag_ITS_soil_normal$class))],
                                                                                         sep="_")
manag_ITS_soil_normal$class[c(grep("Incertae_sedis", manag_ITS_soil_normal$class))]<-paste("unclassified",
                                                                                           manag_ITS_soil_normal$phylum[c(grep("Incertae_sedis", manag_ITS_soil_normal$class))],
                                                                                           sep="_")
manag_ITS_soil_normal$order[which(is.na(manag_ITS_soil_normal$order))]<-manag_ITS_soil_normal$class[which(is.na(manag_ITS_soil_normal$order))]
manag_ITS_soil_normal$order[c(grep("unidentified", manag_ITS_soil_normal$order))]<-paste("unclassified",
                                                                                         manag_ITS_soil_normal$class[c(grep("unidentified", manag_ITS_soil_normal$order))],
                                                                                         sep="_")
manag_ITS_soil_normal$order[c(grep("Incertae_sedis", manag_ITS_soil_normal$order))]<-paste("unclassified",
                                                                                           manag_ITS_soil_normal$class[c(grep("Incertae_sedis", manag_ITS_soil_normal$order))],
                                                                                           sep="_")
manag_ITS_soil_normal$order[c(grep("unclassified_unclassified", manag_ITS_soil_normal$order))]<-manag_ITS_soil_normal$class[c(grep("unclassified_unclassified", manag_ITS_soil_normal$order))]
manag_ITS_soil_normal$family[which(is.na(manag_ITS_soil_normal$family))]<-manag_ITS_soil_normal$order[which(is.na(manag_ITS_soil_normal$family))]
manag_ITS_soil_normal$family[c(grep("unidentified", manag_ITS_soil_normal$family))]<-paste("unclassified",
                                                                                           manag_ITS_soil_normal$order[c(grep("unidentified", manag_ITS_soil_normal$family))],
                                                                                           sep="_")
manag_ITS_soil_normal$family[c(grep("Incertae_sedis", manag_ITS_soil_normal$family))]<-paste("unclassified",
                                                                                             manag_ITS_soil_normal$order[c(grep("Incertae_sedis", manag_ITS_soil_normal$family))],
                                                                                             sep="_")
manag_ITS_soil_normal$family[c(grep("unclassified_unclassified", manag_ITS_soil_normal$family))]<-manag_ITS_soil_normal$order[c(grep("unclassified_unclassified", manag_ITS_soil_normal$family))]
manag_ITS_soil_normal$genus[which(is.na(manag_ITS_soil_normal$genus))]<-manag_ITS_soil_normal$family[which(is.na(manag_ITS_soil_normal$genus))]
manag_ITS_soil_normal$genus[c(grep("unidentified", manag_ITS_soil_normal$genus))]<-paste("unclassified",
                                                                                         manag_ITS_soil_normal$family[c(grep("unidentified", manag_ITS_soil_normal$genus))],
                                                                                         sep="_")
manag_ITS_soil_normal$genus[c(grep("Incertae_sedis", manag_ITS_soil_normal$genus))]<-paste("unclassified",
                                                                                           manag_ITS_soil_normal$family[c(grep("Incertae_sedis", manag_ITS_soil_normal$genus))],
                                                                                           sep="_")
manag_ITS_soil_normal$genus[c(grep("unclassified_unclassified", manag_ITS_soil_normal$genus))]<-manag_ITS_soil_normal$family[c(grep("unclassified_unclassified", manag_ITS_soil_normal$genus))]



manag_ITS_data_soil$phylum[which(manag_ITS_data_soil$phylum=="unidentified")]<-"unclassified_Fungi"
manag_ITS_data_soil$class[which(manag_ITS_data_soil$phylum=="unclassified_Fungi")]<-manag_ITS_data_soil$phylum[which(manag_ITS_data_soil$phylum=="unclassified_Fungi")] 
manag_ITS_data_soil$class[c(grep("unidentified", manag_ITS_data_soil$class))]<-paste("unclassified",
                                                                                     manag_ITS_data_soil$phylum[c(grep("unidentified", manag_ITS_data_soil$class))],
                                                                                     sep="_")
manag_ITS_data_soil$class[c(grep("Incertae_sedis", manag_ITS_data_soil$class))]<-paste("unclassified",
                                                                                       manag_ITS_data_soil$phylum[c(grep("Incertae_sedis", manag_ITS_data_soil$class))],
                                                                                       sep="_")
manag_ITS_data_soil$order[which(is.na(manag_ITS_data_soil$order))]<-manag_ITS_data_soil$class[which(is.na(manag_ITS_data_soil$order))]
manag_ITS_data_soil$order[c(grep("unidentified", manag_ITS_data_soil$order))]<-paste("unclassified",
                                                                                     manag_ITS_data_soil$class[c(grep("unidentified", manag_ITS_data_soil$order))],
                                                                                     sep="_")
manag_ITS_data_soil$order[c(grep("Incertae_sedis", manag_ITS_data_soil$order))]<-paste("unclassified",
                                                                                       manag_ITS_data_soil$class[c(grep("Incertae_sedis", manag_ITS_data_soil$order))],
                                                                                       sep="_")
manag_ITS_data_soil$order[c(grep("unclassified_unclassified", manag_ITS_data_soil$order))]<-manag_ITS_data_soil$class[c(grep("unclassified_unclassified", manag_ITS_data_soil$order))]
manag_ITS_data_soil$family[which(is.na(manag_ITS_data_soil$family))]<-manag_ITS_data_soil$order[which(is.na(manag_ITS_data_soil$family))]
manag_ITS_data_soil$family[c(grep("unidentified", manag_ITS_data_soil$family))]<-paste("unclassified",
                                                                                       manag_ITS_data_soil$order[c(grep("unidentified", manag_ITS_data_soil$family))],
                                                                                       sep="_")
manag_ITS_data_soil$family[c(grep("Incertae_sedis", manag_ITS_data_soil$family))]<-paste("unclassified",
                                                                                         manag_ITS_data_soil$order[c(grep("Incertae_sedis", manag_ITS_data_soil$family))],
                                                                                         sep="_")
manag_ITS_data_soil$family[c(grep("unclassified_unclassified", manag_ITS_data_soil$family))]<-manag_ITS_data_soil$order[c(grep("unclassified_unclassified", manag_ITS_data_soil$family))]
manag_ITS_data_soil$genus[which(is.na(manag_ITS_data_soil$genus))]<-manag_ITS_data_soil$family[which(is.na(manag_ITS_data_soil$genus))]
manag_ITS_data_soil$genus[c(grep("unidentified", manag_ITS_data_soil$genus))]<-paste("unclassified",
                                                                                     manag_ITS_data_soil$family[c(grep("unidentified", manag_ITS_data_soil$genus))],
                                                                                     sep="_")
manag_ITS_data_soil$genus[c(grep("Incertae_sedis", manag_ITS_data_soil$genus))]<-paste("unclassified",
                                                                                       manag_ITS_data_soil$family[c(grep("Incertae_sedis", manag_ITS_data_soil$genus))],
                                                                                       sep="_")
manag_ITS_data_soil$genus[c(grep("unclassified_unclassified", manag_ITS_data_soil$genus))]<-manag_ITS_data_soil$family[c(grep("unclassified_unclassified", manag_ITS_data_soil$genus))]


### Ordination of community composition -------------------------------------
manag_ITS_soil_NMDS<-manag_ITS_soil_normal%>%
  select("ESV", "Sample_034_01_gITS7.ITS4__Run2":"Sample_F03_05_gITS7.ITS4__Run3") %>%
  column_to_rownames("ESV") %>%
  t() %>%
  as.data.frame()

set.seed(19990825)
manag_ITS_soil_NMDS_plot<-as.matrix(manag_ITS_soil_NMDS) %>%
  vegdist(., method = "bray")%>%
  metaMDS(., trymax = 500) %>%
  scores(.) %>%
  as.data.frame() %>%
  mutate("Understory" = plotsite_metadata_soil$Understory,
         "Plot" = plotsite_metadata_soil$Plot,
         .before="NMDS1") %>%
  mutate(Understory = factor(Understory,
                             levels=c("Reduced",
                                      "Intermediate",
                                      "Enhanced",
                                      "Forest")),
         Plot = factor(Plot,
                       levels=c("034", "050", "052",
                                "081", "083", "095",
                                "132", "F01", "F02", "F03")))

polygon_data_ITS_soil<-manag_ITS_soil_NMDS_plot[manag_ITS_soil_NMDS_plot %>%
                                                  group_by(Understory) %>%
                                                  reframe(chull(NMDS1,NMDS2)) %>%
                                                  rename("Rows" = "chull(NMDS1, NMDS2)") %>%
                                                  mutate(Rows=ifelse(Understory=="Intermediate",
                                                                     Rows+18,
                                                                     Rows),
                                                         Rows=ifelse(Understory=="Enhanced",
                                                                     Rows+36,
                                                                     Rows),
                                                         Rows=ifelse(Understory=="Forest",
                                                                     Rows+42,
                                                                     Rows)) %>%
                                                  .$Rows,]

set.seed(324)
manag_ITS_soil_NMDS_envfit<-envfit(manag_ITS_soil_NMDS_plot[,c("NMDS1",
                                                               "NMDS2")], 
                                   plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                   permutations=10000, 
                                   na.rm=T) %>%
  scores(., "vectors") %>%
  as.data.frame(.) %>%
  mutate("PVal"=c(envfit(manag_ITS_soil_NMDS_plot[,c("NMDS1",
                                                     "NMDS2")], 
                         plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                         permutations=10000, 
                         na.rm=T)$vectors$pval),
         "RSq"=c(envfit(manag_ITS_soil_NMDS_plot[,c("NMDS1",
                                                    "NMDS2")], 
                        plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                        permutations=10000, 
                        na.rm=T)$vectors$r))%>%
  filter(PVal<0.05)%>%
  mutate(NMDS1start=envfit(manag_ITS_soil_NMDS_plot[,c("NMDS1",
                                                       "NMDS2")], 
                           plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                           permutations=10000, 
                           na.rm=T) %>%
           scores(., "factors") %>%
           as.data.frame()%>%
           select(NMDS1)%>%
           as.numeric(.),
         NMDS2start=envfit(manag_ITS_soil_NMDS_plot[,c("NMDS1",
                                                       "NMDS2")], 
                           plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                           permutations=10000, 
                           na.rm=T) %>%
           scores(., "factors") %>%
           as.data.frame()%>%
           select(NMDS2)%>%
           as.numeric(.),
         NMDS1end=NMDS1+(envfit(manag_ITS_soil_NMDS_plot[,c("NMDS1",
                                                            "NMDS2")], 
                                plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                permutations=10000, 
                                na.rm=T) %>%
                           scores(., "factors") %>%
                           as.data.frame()%>%
                           select(NMDS1)%>%
                           as.numeric(.)),
         NMDS2end=NMDS2+(envfit(manag_ITS_soil_NMDS_plot[,c("NMDS1",
                                                            "NMDS2")], 
                                plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                permutations=10000, 
                                na.rm=T) %>%
                           scores(., "factors") %>%
                           as.data.frame()%>%
                           select(NMDS2)%>%
                           as.numeric(.)))

ggplot(data=manag_ITS_soil_NMDS_plot)+
  geom_point(size=5, shape=21,colour="black",
             aes(x=NMDS1, y=NMDS2, fill=Understory))+
  geom_segment(data = manag_ITS_soil_NMDS_envfit, 
               aes(x = NMDS1start, y = NMDS2start, 
                   xend = NMDS1end, yend = NMDS2end), 
               arrow = arrow(length = unit(0.1, "inches")),
               col="grey30", linewidth=1, alpha=1,) +
  geom_text_repel(data = manag_ITS_soil_NMDS_envfit, 
                  aes(label = rownames(manag_ITS_soil_NMDS_envfit),
                      x=NMDS1end, y=NMDS2end), 
                  box.padding = 1,       
                  point.padding = 0.5,     
                  segment.color = 'grey50',
                  max.overlaps = Inf,
                  size=5,
                  vjust=-2,
                  hjust=-0.5) + 
  geom_polygon(data = polygon_data_ITS_soil, 
               aes(x = NMDS1, y = NMDS2,
                   fill=Understory),
               alpha=0.25,
               show.legend = F)+
  theme_bw()+
  theme(axis.text = element_text(size=),
        axis.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_colour_manual(name = "Treatment", values= c("Orange", "#DC267F",
                                                    "steelblue1", "#2E8B57"))+
  scale_fill_manual(name = "Treatment", values= c("Orange", "#DC267F",
                                                  "steelblue1", "#2E8B57"))



### Soil mycorrhizal ordination ---------------------------------------------
Taxonomy<-c()
for (e in c(1:nrow(manag_ITS_data_soil))){
  Taxonomy<-c(Taxonomy,
              paste(manag_ITS_data_soil$kingdom[e],
                    manag_ITS_data_soil$phylum[e],
                    manag_ITS_data_soil$class[e],
                    manag_ITS_data_soil$order[e],
                    manag_ITS_data_soil$family[e],
                    manag_ITS_data_soil$genus[e],
                    manag_ITS_data_soil$species[e],
                    sep=";"))
}
manag_ITS_soil_funguild<-cbind(manag_ITS_data_soil[,
                                                   c(which(names(manag_ITS_data_soil)=="ESV"),
                                                     which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_F03_05_gITS7.ITS4__Run3"))],
                               Taxonomy) %>%
  funguild_assign(.) %>%
  mutate(trophicMode = if_else(is.na(trophicMode) |
                                 confidenceRanking=="Possible" |
                                 trophicMode=="Pathotroph-Saprotroph" |
                                 trophicMode=="Pathotroph-Saprotroph-Symbiotroph" |
                                 trophicMode=="Saprotroph-Symbiotroph" |
                                 trophicMode=="Pathotroph-Symbiotroph" |
                                 trophicMode==" Pathotroph-Pathotroph-Saprotroph" |
                                 trophicMode==" Saprotroph-Pathotroph-Saprotroph",
                               "Unknown",
                               trophicMode),
         trophicMode = if_else(trophicMode==" Saprotroph",
                               "Saprotroph",
                               trophicMode)) %>%
  mutate(guild = if_else(is.na(guild) |
                           trophicMode=="Unknown",
                         "Unknown",
                         guild)) %>%
  mutate(guild = if_else(guild == "Lichenized-Wood Saprotroph",
                         "Lichenized", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Plant Pathogen|",
                         "Plant Pathogen", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Wood Saprotroph|",
                         "Wood Saprotroph", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Epiphyte|",
                         "Epiphyte", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Plant Saprotroph|",
                         "Plant Saprotroph", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Fungal Parasite|",
                         "Fungal Parasite", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Lichen Parasite|",
                         "Lichen Parasite", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Ectomycorrhizal|",
                         "Ectomycorrhizal", 
                         guild)) %>%
  mutate(guild = if_else(guild == "NULL" |
                           guild == "Plant Pathogen-Undefined Parasite-Undefined Saprotroph" |
                           guild == "|Animal Parasite|-Animal Pathogen-Arthropod Parasite",
                         "Undefined Pathotroph", 
                         guild)) %>%
  mutate(guild = if_else(guild == "Dung Saprotroph-Undefined Saprotroph-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Plant Saprotroph-Wood Saprotroph" |
                           guild == "Undefined Saprotroph-Wood Saprotroph" |
                           guild == "Leaf Saprotroph-Soil Saprotroph" |
                           guild == "Litter Saprotroph-Soil Saprotroph-Wood Saprotroph" |
                           guild == "Soil Saprotroph-Undefined Saprotroph" |
                           guild == "Dung Saprotroph-Undefined Saprotroph" |
                           guild == "Plant Saprotroph-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Plant Saprotroph" |
                           guild == "Dung Saprotroph-Soil Saprotroph" |
                           guild == "Fungal Parasite-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Plant Saprotroph-Soil Saprotroph" |
                           guild == "|Plant Saprotroph|-Undefined Saprotroph" |
                           guild == "Undefined Saprotroph-|Wood Saprotroph|" |
                           guild == "Plant Saprotroph-|Wood Saprotroph|" |
                           guild == "Plant Saprotroph-Undefined Saprotroph-|Wood Saprotroph|" |
                           guild == "Dung Saprotroph-Plant Saprotroph-|Wood Saprotroph|" |
                           guild == "Plant Saprotroph-|Undefined Saprotroph|" |
                           guild == "|Plant Saprotroph|-Wood Saprotroph" |
                           guild == "|Undefined Saprotroph|" |
                           guild == "Dung Saprotroph-Fungal Parasite-Undefined Saprotroph",
                         "Undefined Saprotroph",
                         guild))%>%
  left_join(.,
            manag_ITS_data_soil %>%
              select(ESV, genus),
            by="ESV")%>%
  mutate(genus=ifelse(is.na(genus),
                      "Unknown",
                      genus),
         guild=ifelse(genus=="Fusarium" | genus=="Ganoderma" | genus=="Corticium" |
                        genus=="Ceratobasidium" | genus=="Lasiodiplodia" | genus=="Phyllosticta" | 
                        genus=="Mycosphaerella" | genus=="Melanconium" | genus=="Cercospora" | 
                        genus=="Alternaria" | genus=="Graphiola" | genus=="Rhizoctonia" | genus=="Graphium" |
                        genus=="Bipolaris" | genus=="Chalara" | genus=="Ceratocystis",
                      "Plant Pathogen",
                      guild),
         trophicMode=ifelse(genus=="Fusarium" | genus=="Ganoderma" | genus=="Corticium" |
                              genus=="Ceratobasidium" | genus=="Lasiodiplodia" | genus=="Phyllosticta" | 
                              genus=="Mycosphaerella" | genus=="Melanconium" | genus=="Cercospora" | 
                              genus=="Alternaria" | genus=="Graphiola" | genus=="Rhizoctonia" | genus=="Graphium" |
                              genus=="Bipolaris" | genus=="Chalara" | genus=="Ceratocystis",
                            "Pathotroph",
                            trophicMode))

manag_ITS_soil_funguild_normal<-manag_ITS_soil_funguild %>%
  filter(trophicMode!="Unknown") %>%
  select(ESV,
         Sample_034_01_gITS7.ITS4__Run2:Sample_F03_05_gITS7.ITS4__Run3) %>%
  column_to_rownames("ESV")%>%
  DESeqDataSetFromMatrix(countData = .,
                         colData = plotsite_metadata_soil %>%
                           column_to_rownames("Sample"),
                         design = ~Understory) %>%
  DESeq2::estimateSizeFactors(.,type="poscounts")%>%
  #DESeq(.) %>% 
  DESeq2::counts(.,
                 normalized=TRUE)%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_ITS_soil_funguild %>%
              select(-c(Sample_034_01_gITS7.ITS4__Run2:Sample_F03_05_gITS7.ITS4__Run3)),
            by="ESV")

manag_ITS_soil_AMF<-manag_ITS_soil_funguild_normal %>%
  filter(guild=="Arbuscular Mycorrhizal") %>%
  select(ESV,
         Sample_034_01_gITS7.ITS4__Run2:Sample_F03_05_gITS7.ITS4__Run3) %>%
  left_join(.,
            manag_ITS_soil_funguild[,c(which(names(manag_ITS_soil_funguild)=="ESV"),
                                       which(names(manag_ITS_soil_funguild)=="guild"))],
            by="ESV") %>%
  left_join(.,
            manag_ITS_soil_normal[,c(which(names(manag_ITS_soil_normal)=="ESV"),
                                     which(names(manag_ITS_soil_normal)=="class"),
                                     which(names(manag_ITS_soil_normal)=="order"),
                                     which(names(manag_ITS_soil_normal)=="family"),
                                     which(names(manag_ITS_soil_normal)=="genus"))],
            by="ESV") %>%
  select(-guild) %>%
  group_by(family) %>%
  summarise(across(matches("^Sample_"), c(abundance=sum))) %>%
  pivot_longer(cols = c(2:ncol(.)),
               names_to = "Sample",
               values_to = "Abundance")%>%
  mutate(Sample=substr(Sample, 1, nchar(Sample) - 10)) %>%
  left_join(.,
            manag_ITS_soil_funguild %>%
              select(-c("Sample_034_01_gITS7.ITS4__Run2":"Sample_F03_05_gITS7.ITS4__Run3")) %>%
              left_join(.,
                        manag_ITS_soil_rarefied %>%
                          select(ESV,
                                 Sample_034_01_gITS7.ITS4__Run2:Sample_F03_05_gITS7.ITS4__Run3),
                        by="ESV") %>%
              relocate(ESV, Sample_034_01_gITS7.ITS4__Run2:Sample_F03_05_gITS7.ITS4__Run3,
                       Taxonomy:citationSource)%>%
              filter(guild=="Arbuscular Mycorrhizal") %>%
              select(ESV,
                     Sample_034_01_gITS7.ITS4__Run2:Sample_F03_05_gITS7.ITS4__Run3) %>%
              left_join(.,
                        manag_ITS_soil_funguild[,c(which(names(manag_ITS_soil_funguild)=="ESV"),
                                                   which(names(manag_ITS_soil_funguild)=="guild"))],
                        by="ESV") %>%
              left_join(.,
                        manag_ITS_soil_normal[,c(which(names(manag_ITS_soil_normal)=="ESV"),
                                                 which(names(manag_ITS_soil_normal)=="class"),
                                                 which(names(manag_ITS_soil_normal)=="order"),
                                                 which(names(manag_ITS_soil_normal)=="family"),
                                                 which(names(manag_ITS_soil_normal)=="genus"))],
                        by="ESV") %>%
              select(-guild) %>%
              group_by(family) %>%
              summarise(across(matches("^Sample"),
                               c(richness=~my_shannon(.)))) %>%
              pivot_longer(cols = c(2:ncol(.)),
                           names_to = "Sample",
                           values_to = "Richness")%>%
              mutate(Sample=substr(Sample, 1, nchar(Sample) - 9)),
            by=c("family", "Sample"))%>%
  mutate(Understory=c(rep(plotsite_metadata_soil$Understory,
                          length(unique(family)))),
         .before = "Abundance") %>%
  pivot_longer(cols = c(Abundance:Richness),
               names_to = "Metric",
               values_to = "Value")

set.seed(43123123)
manag_ITS_soil_AMF_NMDS_plot<-as.matrix(manag_ITS_soil_funguild_normal %>%
                                          filter(guild=="Arbuscular Mycorrhizal") %>%
                                          column_to_rownames("ESV") %>%
                                          select(Sample_034_01_gITS7.ITS4__Run2:Sample_F03_05_gITS7.ITS4__Run3) %>%
                                          t() %>%
                                          as.data.frame()) %>%
  vegdist(., method = "bray") %>%
  metaMDS(., trymax = 500) %>%
  scores(.) %>%
  as.data.frame() %>%
  mutate("Understory" = plotsite_metadata_soil$Understory,
         .before="NMDS1") %>%
  mutate(Understory = factor(Understory,
                             levels=c("Reduced",
                                      "Intermediate",
                                      "Enhanced",
                                      "Forest")))
polygon_data_ITS_soil_AMF<-manag_ITS_soil_AMF_NMDS_plot[manag_ITS_soil_AMF_NMDS_plot %>%
                                                          group_by(Understory) %>%
                                                          reframe(chull(NMDS1,NMDS2)) %>%
                                                          rename("Rows" = "chull(NMDS1, NMDS2)") %>%
                                                          mutate(Rows=ifelse(Understory=="Intermediate",
                                                                             Rows+18,
                                                                             Rows),
                                                                 Rows=ifelse(Understory=="Enhanced",
                                                                             Rows+36,
                                                                             Rows),
                                                                 Rows=ifelse(Understory=="Forest",
                                                                             Rows+42,
                                                                             Rows)) %>%
                                                          .$Rows,]


set.seed(42342)
manag_ITS_soil_AMF_NMDS_envfit<-envfit(manag_ITS_soil_AMF_NMDS_plot[,c("NMDS1",
                                                                       "NMDS2")], 
                                       manag_ITS_soil_AMF %>%
                                         filter(Metric=="Abundance") %>%
                                         pivot_wider(names_from = family,
                                                     values_from = Value) %>%
                                         select(Acaulosporaceae:unclassified_Paraglomerales), 
                                       permutations=10000, 
                                       na.rm=T) %>%
  scores(., "vectors") %>%
  as.data.frame(.) %>%
  mutate("PVal"=c(envfit(manag_ITS_soil_AMF_NMDS_plot[,c("NMDS1",
                                                         "NMDS2")], 
                         manag_ITS_soil_AMF %>%
                           filter(Metric=="Abundance") %>%
                           pivot_wider(names_from = family,
                                       values_from = Value) %>%
                           select(Acaulosporaceae:unclassified_Paraglomerales), 
                         permutations=10000, 
                         na.rm=T)$vectors$pval),
         "RSq"=c(envfit(manag_ITS_soil_AMF_NMDS_plot[,c("NMDS1",
                                                        "NMDS2")], 
                        manag_ITS_soil_AMF %>%
                          filter(Metric=="Abundance") %>%
                          pivot_wider(names_from = family,
                                      values_from = Value) %>%
                          select(Acaulosporaceae:unclassified_Paraglomerales), 
                        permutations=10000, 
                        na.rm=T)$vectors$r)) %>%
  filter(PVal<0.05) %>%
  rownames_to_column("Family")%>%
  mutate(Family=ifelse(Family=="unclassified_Diversisporales",
                       "Order Diversisporales",
                       Family),
         Family=ifelse(Family=="unclassified_Glomerales",
                       "Order Glomerales",
                       Family),
         Family=ifelse(Family=="unclassified_Glomeromycetes",
                       "Class Glomeromycetes",
                       Family),
         Family=ifelse(Family=="unclassified_Glomeromycota",
                       "Phylum Glomeromycota",
                       Family),
         Family=ifelse(Family=="unclassified_Paraglomerales",
                       "Order Paraglomerales",
                       Family))%>%
  column_to_rownames("Family")

ggplot(data=manag_ITS_soil_AMF_NMDS_plot)+
  geom_point(size=5, aes(x=NMDS1, y=NMDS2, colour=Understory))+
  #stat_ellipse(level=0.5,
  #             linewidth=1.5,
  #             linetype=4,
  #             aes(x=NMDS1, y=NMDS2, col=Understory))+
  geom_polygon(data = polygon_data_ITS_soil_AMF, 
               aes(x = NMDS1, y = NMDS2,
                   fill=Understory),
               alpha=0.25,
               show.legend = F)+
  geom_segment(data = manag_ITS_soil_AMF_NMDS_envfit, 
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.1, "inches")),
               col="grey30", linewidth=1, alpha=1,) +
  geom_text_repel(data = manag_ITS_soil_AMF_NMDS_envfit, 
                  aes(label = gsub(" ", "\n", rownames(manag_ITS_soil_AMF_NMDS_envfit)),
                      x = NMDS1, y = NMDS2), 
                  box.padding = 1.5,       
                  point.padding = 0.5,     
                  segment.color = 'grey50',
                  max.overlaps = Inf,
                  size=7) +    
  theme_bw()+
  theme(axis.text = element_text(size=),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_colour_manual(name = "Treatment", values= c("Orange", "#DC267F",
                                                    "steelblue1",
                                                    "#2E8B57"))+
  scale_fill_manual(name = "Treatment", values= c("Orange", "#DC267F",
                                                  "steelblue1",
                                                  "#2E8B57"))


## Bacteria ----------------------------------------------------------------

### Dataframe setup ---------------------------------------------------------

manag_16s_data_soil<-read.csv("Available data/bacterial_dataset.csv")%>%
  select(-grep("Run1", names(.), value = TRUE))

plotsite_metadata_soil<-data.frame(cbind(c(names(manag_16s_data_soil)[c(which(names(manag_16s_data_soil)=="Sample_034_01_515F.806R__Run2"):which(names(manag_16s_data_soil)=="Sample_132_06_515F.806R__Run2"),
                                                                        which(names(manag_16s_data_soil)=="Sample_F01_01_515F.806R__Run3"):which(names(manag_16s_data_soil)=="Sample_F03_05_515F.806R__Run3"))]),
                                         c(rep("034", 6),
                                           rep("050", 6),
                                           rep("052", 6),
                                           rep("081", 6),
                                           rep("083", 6),
                                           rep("095", 6),
                                           rep("132", 6),
                                           rep("F01", 6),
                                           rep("F02", 6),
                                           rep("F03", 5)),
                                         c(rep("Reduced", 18),
                                           rep("Intermediate", 18),
                                           rep("Enhanced", 6),
                                           rep("Forest", 17)))) %>%
  rename("Sample" = "X1",
         "Plot" = "X2",
         "Understory" = "X3") %>%
  mutate(Plot=factor(Plot,
                     levels=c("034", "050", "052",
                              "081", "083", "095",
                              "132", "F01", "F02", "F03")),
         Understory=factor(Understory,
                           levels=c("Reduced",
                                    "Intermediate",
                                    "Enhanced",
                                    "Forest"))) %>%
  mutate(Biomass=c(read.csv("Available data/Fieldwork metadata.csv")$Biomass.mean))%>%
  mutate(NMDS_centroid=NA)%>%
  mutate(NMDS_centroid=ifelse(Understory=="Reduced" | Understory =="Intermediate",
                              "Y",
                              "N"))

for (e in c(which(names(manag_16s_data_soil)=="Sample_034_01_515F.806R__Run2"):which(names(manag_16s_data_soil)=="Sample_132_06_515F.806R__Run2"))){
  manag_16s_data_soil[,e]<-c(manag_16s_data_soil[,e]-manag_16s_data_soil$Sample_Neg16s_515F.806R__Run2)
}
for (e in c(which(names(manag_16s_data_soil)=="Sample_F01_01_515F.806R__Run3"):which(names(manag_16s_data_soil)=="Sample_F03_05_515F.806R__Run3"))){
  manag_16s_data_soil[,e]<-c(manag_16s_data_soil[,e]-manag_16s_data_soil$Sample_ForNeg_515F.806R__Run3)
}
for (e in c(which(names(manag_16s_data_soil)=="Sample_034_01_515F.806R__Run2"):which(names(manag_16s_data_soil)=="Sample_F03_05_515F.806R__Run3"))){
  manag_16s_data_soil[which(manag_16s_data_soil[,e]<0),e]<-0
}
for (e in c(which(names(manag_16s_data_soil)=="Sample_034_01_515F.806R__Run2"):which(names(manag_16s_data_soil)=="Sample_F03_05_515F.806R__Run3"))){
  manag_16s_data_soil[which(manag_16s_data_soil[,e]==1),e]<-0
}
manag_16s_data_soil<-manag_16s_data_soil[,
                                         -c(which(names(manag_16s_data_soil)=="Sample_Neg16s_515F.806R__Run2"),
                                            which(names(manag_16s_data_soil)=="Sample_ForNeg_515F.806R__Run3"))]

non_occur_ESV<-c()
for (e in c(1:nrow(manag_16s_data_soil))){
  non_occur_ESV<-c(non_occur_ESV,
                   sum(manag_16s_data_soil[e,
                                           c(which(names(manag_16s_data_soil)=="Sample_034_01_515F.806R__Run2"):which(names(manag_16s_data_soil)=="Sample_F03_05_515F.806R__Run3"))]))
}
manag_16s_data_soil<-manag_16s_data_soil[-c(which(non_occur_ESV==0)),]

manag_16s_data_soil<-subset(manag_16s_data_soil,
                            manag_16s_data_soil$domain=="Bacteria")

manag_16s_data_soil<-manag_16s_data_soil[-c(which(manag_16s_data_soil$order=="Chloroplast"),
                                            which(manag_16s_data_soil$family=="Mitochondria")),]


### Rarefaction and normalisation -------------------------------------------
manag_16s_soil_rarefaction<-manag_16s_data_soil[,
                                                c(1:which(names(manag_16s_data_soil)=="Sample_F03_05_515F.806R__Run3"))] %>%
  remove_rownames(.) %>%
  column_to_rownames(.,
                     var="ESV") %>%
  t(.) %>%
  as.data.frame(.) %>%
  filter(!row_number() %in% c(1:2))%>%
  mutate(across(everything(), as.numeric))

manag_16s_soil_normal<-
  DESeqDataSetFromMatrix(countData = manag_16s_soil_rarefaction%>%
                           t() %>%
                           as.data.frame(),
                         colData = plotsite_metadata_soil %>%
                           column_to_rownames("Sample"),
                         design = ~Understory)%>%
  DESeq2::estimateSizeFactors(.,type="poscounts")%>%
  #DESeq(.) %>% 
  DESeq2::counts(.,
                 normalized=TRUE)%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_16s_data_soil%>%
              select(-c("Sample_034_01_515F.806R__Run2":"Sample_F03_05_515F.806R__Run3")),
            by=c("ESV"))

set.seed(2432)
manag_16s_soil_rarefied<-manag_16s_soil_rarefaction %>%
  rarefy_and_average(.,
                     min(sapply(manag_16s_data_soil[c(which(names(manag_16s_data_soil)=="Sample_034_01_515F.806R__Run2"):which(names(manag_16s_data_soil)=="Sample_F03_05_515F.806R__Run3"))], 
                                sum)),
                     500) %>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_16s_data_soil%>%
              select(-c("Sample_034_01_515F.806R__Run2":"Sample_F03_05_515F.806R__Run3")),
            by=c("ESV"))

### Making taxonomic assignment user-friendly -------------------------------

manag_16s_soil_normal$class[which(is.na(manag_16s_soil_normal$class))]<-paste("unclassified",
                                                                              manag_16s_soil_normal$phylum[which(is.na(manag_16s_soil_normal$class))],
                                                                              sep="_")
manag_16s_soil_normal$class[c(grep("unidentified", manag_16s_soil_normal$class))]<-paste("unclassified",
                                                                                         manag_16s_soil_normal$phylum[c(grep("unidentified", manag_16s_soil_normal$class))],
                                                                                         sep="_")
manag_16s_soil_normal$class[c(grep("Incertae Sedis", manag_16s_soil_normal$class))]<-paste("unclassified",
                                                                                           manag_16s_soil_normal$phylum[c(grep("Incertae Sedis", manag_16s_soil_normal$class))],
                                                                                           sep="_")
manag_16s_soil_normal$class[c(grep("unclassified_unclassified", manag_16s_soil_normal$class))]<-manag_16s_soil_normal$phylum[c(grep("unclassified_unclassified", manag_16s_soil_normal$class))]

manag_16s_soil_normal$order[which(is.na(manag_16s_soil_normal$order))]<-paste("unclassified",
                                                                              manag_16s_soil_normal$class[which(is.na(manag_16s_soil_normal$order))],
                                                                              sep="_")
manag_16s_soil_normal$order[c(grep("unidentified", manag_16s_soil_normal$order))]<-paste("unclassified",
                                                                                         manag_16s_soil_normal$class[c(grep("unidentified", manag_16s_soil_normal$order))],
                                                                                         sep="_")
manag_16s_soil_normal$order[c(grep("Incertae Sedis", manag_16s_soil_normal$order))]<-paste("unclassified",
                                                                                           manag_16s_soil_normal$class[c(grep("Incertae Sedis", manag_16s_soil_normal$order))],
                                                                                           sep="_")
manag_16s_soil_normal$order[c(grep("unclassified_unclassified", manag_16s_soil_normal$order))]<-manag_16s_soil_normal$class[c(grep("unclassified_unclassified", manag_16s_soil_normal$order))]

manag_16s_soil_normal$family[which(is.na(manag_16s_soil_normal$family))]<-paste("unclassified",
                                                                                manag_16s_soil_normal$order[which(is.na(manag_16s_soil_normal$family))],
                                                                                sep="_")
manag_16s_soil_normal$family[c(grep("unidentified", manag_16s_soil_normal$family))]<-paste("unclassified",
                                                                                           manag_16s_soil_normal$order[c(grep("unidentified", manag_16s_soil_normal$family))],
                                                                                           sep="_")
manag_16s_soil_normal$family[c(grep("Unknown", manag_16s_soil_normal$family))]<-paste("unclassified",
                                                                                      manag_16s_soil_normal$order[c(grep("Unknown", manag_16s_soil_normal$family))],
                                                                                      sep="_")
manag_16s_soil_normal$family[c(grep("Incertae Sedis", manag_16s_soil_normal$family))]<-paste("unclassified",
                                                                                             manag_16s_soil_normal$order[c(grep("Incertae Sedis", manag_16s_soil_normal$family))],
                                                                                             sep="_")

manag_16s_soil_normal$family[c(grep("unclassified_unclassified", manag_16s_soil_normal$family))]<-manag_16s_soil_normal$order[c(grep("unclassified_unclassified", manag_16s_soil_normal$family))]
manag_16s_soil_normal$genus[which(is.na(manag_16s_soil_normal$genus))]<-paste("unclassified",
                                                                              manag_16s_soil_normal$family[which(is.na(manag_16s_soil_normal$genus))],
                                                                              sep="_")
manag_16s_soil_normal$genus[c(grep("unidentified", manag_16s_soil_normal$genus))]<-paste("unclassified",
                                                                                         manag_16s_soil_normal$family[c(grep("unidentified", manag_16s_soil_normal$genus))],
                                                                                         sep="_")
manag_16s_soil_normal$genus[c(grep("Unknown", manag_16s_soil_normal$genus))]<-paste("unclassified",
                                                                                    manag_16s_soil_normal$family[c(grep("Unknown", manag_16s_soil_normal$genus))],
                                                                                    sep="_")
manag_16s_soil_normal$genus[c(grep("Incertae Sedis", manag_16s_soil_normal$genus))]<-paste("unclassified",
                                                                                           manag_16s_soil_normal$family[c(grep("Incertae Sedis", manag_16s_soil_normal$genus))],
                                                                                           sep="_")

manag_16s_soil_normal$genus[c(grep("unclassified_unclassified", manag_16s_soil_normal$genus))]<-manag_16s_soil_normal$family[c(grep("unclassified_unclassified", manag_16s_soil_normal$genus))]

manag_16s_data_soil$class[which(is.na(manag_16s_data_soil$class))]<-paste("unclassified",
                                                                          manag_16s_data_soil$phylum[which(is.na(manag_16s_data_soil$class))],
                                                                          sep="_")
manag_16s_data_soil$class[c(grep("unidentified", manag_16s_data_soil$class))]<-paste("unclassified",
                                                                                     manag_16s_data_soil$phylum[c(grep("unidentified", manag_16s_data_soil$class))],
                                                                                     sep="_")
manag_16s_data_soil$class[c(grep("Incertae Sedis", manag_16s_data_soil$class))]<-paste("unclassified",
                                                                                       manag_16s_data_soil$phylum[c(grep("Incertae Sedis", manag_16s_data_soil$class))],
                                                                                       sep="_")
manag_16s_data_soil$class[c(grep("unclassified_unclassified", manag_16s_data_soil$class))]<-manag_16s_data_soil$phylum[c(grep("unclassified_unclassified", manag_16s_data_soil$class))]

manag_16s_data_soil$order[which(is.na(manag_16s_data_soil$order))]<-paste("unclassified",
                                                                          manag_16s_data_soil$class[which(is.na(manag_16s_data_soil$order))],
                                                                          sep="_")
manag_16s_data_soil$order[c(grep("unidentified", manag_16s_data_soil$order))]<-paste("unclassified",
                                                                                     manag_16s_data_soil$class[c(grep("unidentified", manag_16s_data_soil$order))],
                                                                                     sep="_")
manag_16s_data_soil$order[c(grep("Incertae Sedis", manag_16s_data_soil$order))]<-paste("unclassified",
                                                                                       manag_16s_data_soil$class[c(grep("Incertae Sedis", manag_16s_data_soil$order))],
                                                                                       sep="_")
manag_16s_data_soil$order[c(grep("unclassified_unclassified", manag_16s_data_soil$order))]<-manag_16s_data_soil$class[c(grep("unclassified_unclassified", manag_16s_data_soil$order))]

manag_16s_data_soil$family[which(is.na(manag_16s_data_soil$family))]<-paste("unclassified",
                                                                            manag_16s_data_soil$order[which(is.na(manag_16s_data_soil$family))],
                                                                            sep="_")
manag_16s_data_soil$family[c(grep("unidentified", manag_16s_data_soil$family))]<-paste("unclassified",
                                                                                       manag_16s_data_soil$order[c(grep("unidentified", manag_16s_data_soil$family))],
                                                                                       sep="_")
manag_16s_data_soil$family[c(grep("Unknown", manag_16s_data_soil$family))]<-paste("unclassified",
                                                                                  manag_16s_data_soil$order[c(grep("Unknown", manag_16s_data_soil$family))],
                                                                                  sep="_")
manag_16s_data_soil$family[c(grep("Incertae Sedis", manag_16s_data_soil$family))]<-paste("unclassified",
                                                                                         manag_16s_data_soil$order[c(grep("Incertae Sedis", manag_16s_data_soil$family))],
                                                                                         sep="_")

manag_16s_data_soil$family[c(grep("unclassified_unclassified", manag_16s_data_soil$family))]<-manag_16s_data_soil$order[c(grep("unclassified_unclassified", manag_16s_data_soil$family))]
manag_16s_data_soil$genus[which(is.na(manag_16s_data_soil$genus))]<-paste("unclassified",
                                                                          manag_16s_data_soil$family[which(is.na(manag_16s_data_soil$genus))],
                                                                          sep="_")
manag_16s_data_soil$genus[c(grep("unidentified", manag_16s_data_soil$genus))]<-paste("unclassified",
                                                                                     manag_16s_data_soil$family[c(grep("unidentified", manag_16s_data_soil$genus))],
                                                                                     sep="_")
manag_16s_data_soil$genus[c(grep("Unknown", manag_16s_data_soil$genus))]<-paste("unclassified",
                                                                                manag_16s_data_soil$family[c(grep("Unknown", manag_16s_data_soil$genus))],
                                                                                sep="_")
manag_16s_data_soil$genus[c(grep("Incertae Sedis", manag_16s_data_soil$genus))]<-paste("unclassified",
                                                                                       manag_16s_data_soil$family[c(grep("Incertae Sedis", manag_16s_data_soil$genus))],
                                                                                       sep="_")

manag_16s_data_soil$genus[c(grep("unclassified_unclassified", manag_16s_data_soil$genus))]<-manag_16s_data_soil$family[c(grep("unclassified_unclassified", manag_16s_data_soil$genus))]

### Ordination of community composition -------------------------------------
manag_16s_soil_NMDS<-manag_16s_soil_normal %>%
  select("ESV", "Sample_034_01_515F.806R__Run2":"Sample_F03_05_515F.806R__Run3") %>%
  column_to_rownames("ESV") %>%
  t()%>%
  as.data.frame()

set.seed(43843)
manag_16s_soil_NMDS_plot<-as.matrix(manag_16s_soil_NMDS) %>%
  vegdist(., method = "bray") %>%
  metaMDS(.) %>%
  scores(.) %>%
  as.data.frame() %>%
  mutate("Understory" = plotsite_metadata_soil$Understory,
         "Plot" = plotsite_metadata_soil$Plot,
         .before="NMDS1") %>%
  mutate(Understory = factor(Understory,
                             levels=c("Reduced",
                                      "Intermediate",
                                      "Enhanced",
                                      "Forest")),
         Plot = factor(Plot,
                       levels=c("034", "050", "052",
                                "081", "083", "095",
                                "132", "F01", "F02",
                                "F03")))

polygon_data_16s_soil<-manag_16s_soil_NMDS_plot[manag_16s_soil_NMDS_plot %>%
                                                  group_by(Understory) %>%
                                                  reframe(chull(NMDS1,NMDS2)) %>%
                                                  rename("Rows" = "chull(NMDS1, NMDS2)") %>%
                                                  mutate(Rows=ifelse(Understory=="Intermediate",
                                                                     Rows+18,
                                                                     Rows),
                                                         Rows=ifelse(Understory=="Enhanced",
                                                                     Rows+36,
                                                                     Rows),
                                                         Rows=ifelse(Understory=="Forest",
                                                                     Rows+42,
                                                                     Rows)) %>%
                                                  .$Rows,]

set.seed(42342)
manag_16s_soil_NMDS_envfit<-envfit(manag_16s_soil_NMDS_plot[,c("NMDS1",
                                                               "NMDS2")], 
                                   plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                   permutations=10000, 
                                   na.rm=T) %>%
  scores(., "vectors") %>%
  as.data.frame(.) %>%
  mutate("PVal"=c(envfit(manag_16s_soil_NMDS_plot[,c("NMDS1",
                                                     "NMDS2")], 
                         plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                         permutations=10000, 
                         na.rm=T)$vectors$pval),
         "RSq"=c(envfit(manag_16s_soil_NMDS_plot[,c("NMDS1",
                                                    "NMDS2")], 
                        plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                        permutations=10000, 
                        na.rm=T)$vectors$r))%>%
  filter(PVal<0.05)%>%
  mutate(NMDS1start=envfit(manag_16s_soil_NMDS_plot[,c("NMDS1",
                                                       "NMDS2")], 
                           plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                           permutations=10000, 
                           na.rm=T) %>%
           scores(., "factors") %>%
           as.data.frame()%>%
           select(NMDS1)%>%
           as.numeric(.),
         NMDS2start=envfit(manag_16s_soil_NMDS_plot[,c("NMDS1",
                                                       "NMDS2")], 
                           plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                           permutations=10000, 
                           na.rm=T) %>%
           scores(., "factors") %>%
           as.data.frame()%>%
           select(NMDS2)%>%
           as.numeric(.),
         NMDS1end=NMDS1+(envfit(manag_16s_soil_NMDS_plot[,c("NMDS1",
                                                            "NMDS2")], 
                                plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                permutations=10000, 
                                na.rm=T) %>%
                           scores(., "factors") %>%
                           as.data.frame()%>%
                           select(NMDS1)%>%
                           as.numeric(.)),
         NMDS2end=NMDS2+(envfit(manag_16s_soil_NMDS_plot[,c("NMDS1",
                                                            "NMDS2")], 
                                plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                permutations=10000, 
                                na.rm=T) %>%
                           scores(., "factors") %>%
                           as.data.frame()%>%
                           select(NMDS2)%>%
                           as.numeric(.)))

ggplot(data=manag_16s_soil_NMDS_plot)+
  geom_point(size=5, shape=21,colour="black",
             aes(x=NMDS1, y=NMDS2, fill=Understory))+
  geom_segment(data = manag_16s_soil_NMDS_envfit, 
               aes(x = NMDS1start, y = NMDS2start, 
                   xend = NMDS1end, yend = NMDS2end), 
               arrow = arrow(length = unit(0.1, "inches")),
               col="grey30", linewidth=1, alpha=1,) +
  geom_text_repel(data = manag_16s_soil_NMDS_envfit, 
                  aes(label = rownames(manag_16s_soil_NMDS_envfit),
                      x=NMDS1end, y=NMDS2end), 
                  box.padding = 1,       
                  point.padding = 0.5,     
                  segment.color = 'grey50',
                  max.overlaps = Inf,
                  size=5,
                  vjust=-2) + 
  geom_polygon(data = polygon_data_16s_soil, 
               aes(x = NMDS1, y = NMDS2,
                   fill=Understory),
               alpha=0.25,
               show.legend = F)+
  theme_bw()+
  theme(axis.text = element_text(size=),
        axis.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_colour_manual(name = "Treatment", values= c("Orange", "#DC267F",
                                                    "steelblue1", "#2E8B57"))+
  scale_fill_manual(name = "Treatment", values= c("Orange", "#DC267F",
                                                  "steelblue1", "#2E8B57"))

## Animalia ----------------------------------------------------------------

### Dataframe setup ---------------------------------------------------------
manag_18s_data_soil<-read.csv("path/to/animal_dataset.csv")

plotsite_metadata_soil<-data.frame(cbind(c(names(manag_18s_data_soil)[c(which(names(manag_18s_data_soil)=="Sample_034_01_M620F.M1041F.RC__Run2"):which(names(manag_18s_data_soil)=="Sample_132_06_M620F.M1041F.RC__Run2"),
                                                                        which(names(manag_18s_data_soil)=="Sample_F01_01_M620F.M1041F.RC__Run3"):which(names(manag_18s_data_soil)=="Sample_F03_05_M620F.M1041F.RC__Run3"))]),
                                         c(rep("034", 6),
                                           rep("050", 6),
                                           rep("052", 6),
                                           rep("081", 6),
                                           rep("083", 6),
                                           rep("095", 6),
                                           rep("132", 6),
                                           rep("F01", 6),
                                           rep("F02", 6),
                                           rep("F03", 5)),
                                         c(rep("Reduced", 18),
                                           rep("Intermediate", 18),
                                           rep("Enhanced", 6),
                                           rep("Forest", 17)))) %>%
  rename("Sample" = "X1",
         "Plot" = "X2",
         "Understory" = "X3") %>%
  mutate(Plot=factor(Plot,
                     levels=c("034", "050", "052",
                              "081", "083", "095",
                              "132", "F01", "F02", "F03")),
         Understory=factor(Understory,
                           levels=c("Reduced",
                                    "Intermediate",
                                    "Enhanced",
                                    "Forest"))) %>%
  mutate(Biomass=c(read.csv("path/to/Fieldwork metadata.csv")$Biomass.mean))%>%
  mutate(NMDS_centroid=NA)%>%
  mutate(NMDS_centroid=ifelse(Understory=="Reduced" | Understory =="Intermediate",
                              "Y",
                              "N"))

for (e in c(which(names(manag_18s_data_soil)=="Sample_034_01_M620F.M1041F.RC__Run2"):which(names(manag_18s_data_soil)=="Sample_132_06_M620F.M1041F.RC__Run2"))){
  manag_18s_data_soil[,e]<-c(manag_18s_data_soil[,e]-manag_18s_data_soil$Sample_Neg18s_M620F.M1041F.RC__Run2)
}
for (e in c(which(names(manag_18s_data_soil)=="Sample_F01_01_M620F.M1041F.RC__Run3"):which(names(manag_18s_data_soil)=="Sample_F03_05_M620F.M1041F.RC__Run3"))){
  manag_18s_data_soil[,e]<-c(manag_18s_data_soil[,e]-manag_18s_data_soil$Sample_ForNeg_M620F.M1041F.RC__Run3)
}
for (e in c(which(names(manag_18s_data_soil)=="Sample_034_01_M620F.M1041F.RC__Run2"):which(names(manag_18s_data_soil)=="Sample_F03_05_M620F.M1041F.RC__Run3"))){
  manag_18s_data_soil[which(manag_18s_data_soil[,e]<0),e]<-0
}
for (e in c(which(names(manag_18s_data_soil)=="Sample_034_01_M620F.M1041F.RC__Run2"):which(names(manag_18s_data_soil)=="Sample_F03_05_M620F.M1041F.RC__Run3"))){
  manag_18s_data_soil[which(manag_18s_data_soil[,e]==1),e]<-0
}
manag_18s_data_soil<-manag_18s_data_soil[,
                                         -c(which(names(manag_18s_data_soil)=="Sample_Neg18s_M620F.M1041F.RC__Run2"),
                                            which(names(manag_18s_data_soil)=="Sample_ForNeg_M620F.M1041F.RC__Run3"))]

non_occur_ESV<-c()
for (e in c(1:nrow(manag_18s_data_soil))){
  non_occur_ESV<-c(non_occur_ESV,
                   sum(manag_18s_data_soil[e,
                                           c(which(names(manag_18s_data_soil)=="Sample_034_01_M620F.M1041F.RC__Run2"):which(names(manag_18s_data_soil)=="Sample_F03_05_M620F.M1041F.RC__Run3"))]))
}
manag_18s_data_soil<-manag_18s_data_soil[-c(which(non_occur_ESV==0)),]

manag_18s_data_soil<-subset(manag_18s_data_soil,
                            manag_18s_data_soil$kingdom=="Animalia")

### Rarefaction and normalisation -------------------------------------------
manag_18s_soil_rarefaction<-manag_18s_data_soil[,
                                                c(1:which(names(manag_18s_data_soil)=="Sample_F03_05_M620F.M1041F.RC__Run3"))] %>%
  remove_rownames(.) %>%
  column_to_rownames(.,
                     var="ESV") %>%
  t(.) %>%
  as.data.frame(.) %>%
  filter(!row_number() %in% c(1:3))%>%
  mutate(across(everything(), as.numeric))

manag_18s_soil_normal<-
  DESeqDataSetFromMatrix(countData = manag_18s_soil_rarefaction%>%
                           t() %>%
                           as.data.frame(),
                         colData = plotsite_metadata_soil %>%
                           column_to_rownames("Sample"),
                         design = ~Understory)%>%
  DESeq2::estimateSizeFactors(.,type="poscounts")%>%
  #DESeq(.) %>% 
  DESeq2::counts(.,
                 normalized=TRUE)%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_18s_data_soil%>%
              select(-c("Sample_034_01_M620F.M1041F.RC__Run2":"Sample_F03_05_M620F.M1041F.RC__Run3")),
            by=c("ESV"))

set.seed(2432)
manag_18s_soil_rarefied<-manag_18s_soil_rarefaction %>%
  rarefy_and_average(.,
                     min(sapply(manag_18s_data_soil[c(which(names(manag_18s_data_soil)=="Sample_034_01_M620F.M1041F.RC__Run2"):which(names(manag_18s_data_soil)=="Sample_F03_05_M620F.M1041F.RC__Run3"))], 
                                sum)),
                     500) %>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_18s_data_soil%>%
              select(-c("Sample_034_01_M620F.M1041F.RC__Run2":"Sample_F03_05_M620F.M1041F.RC__Run3")),
            by=c("ESV"))

### Making taxonomic assignment user-friendly -------------------------------
manag_18s_soil_normal$phylum[which(is.na(manag_18s_soil_normal$phylum))]<-"unclassified_Animalia"
manag_18s_soil_normal$phylum[which(grepl("unclassified",
                                         manag_18s_soil_normal$phylum))]<-"unclassified_Animalia"
manag_18s_soil_normal$class[which(is.na(manag_18s_soil_normal$class))]<-paste("unclassified",
                                                                              manag_18s_soil_normal$phylum[which(is.na(manag_18s_soil_normal$class))],
                                                                              sep="_")
manag_18s_soil_normal$class[which(grepl("unclassified_unclassified",
                                        manag_18s_soil_normal$class))]<-manag_18s_soil_normal$phylum[c(grep("unclassified_unclassified", manag_18s_soil_normal$class))]
manag_18s_soil_normal$order[which(is.na(manag_18s_soil_normal$order))]<-paste("unclassified",
                                                                              manag_18s_soil_normal$class[which(is.na(manag_18s_soil_normal$order))],
                                                                              sep="_")
manag_18s_soil_normal$order[which(grepl("unclassified_unclassified",
                                        manag_18s_soil_normal$order))]<-manag_18s_soil_normal$class[c(grep("unclassified_unclassified", manag_18s_soil_normal$order))]
manag_18s_soil_normal$family[which(is.na(manag_18s_soil_normal$family))]<-paste("unclassified",
                                                                                manag_18s_soil_normal$order[which(is.na(manag_18s_soil_normal$family))],
                                                                                sep="_")
manag_18s_soil_normal$family[which(grepl("unclassified_unclassified",
                                         manag_18s_soil_normal$family))]<-manag_18s_soil_normal$order[c(grep("unclassified_unclassified", manag_18s_soil_normal$family))]
manag_18s_soil_normal$genus[which(is.na(manag_18s_soil_normal$genus))]<-paste("unclassified",
                                                                              manag_18s_soil_normal$family[which(is.na(manag_18s_soil_normal$genus))],
                                                                              sep="_")
manag_18s_soil_normal$genus[which(grepl("unclassified_unclassified",
                                        manag_18s_soil_normal$genus))]<-manag_18s_soil_normal$family[c(grep("unclassified_unclassified", manag_18s_soil_normal$genus))]



manag_18s_data_soil$phylum[which(is.na(manag_18s_data_soil$phylum))]<-"unclassified_Animalia"
manag_18s_data_soil$phylum[which(grepl("unclassified",
                                       manag_18s_data_soil$phylum))]<-"unclassified_Animalia"
manag_18s_data_soil$class[which(is.na(manag_18s_data_soil$class))]<-paste("unclassified",
                                                                          manag_18s_data_soil$phylum[which(is.na(manag_18s_data_soil$class))],
                                                                          sep="_")
manag_18s_data_soil$class[which(grepl("unclassified_unclassified",
                                      manag_18s_data_soil$class))]<-manag_18s_data_soil$phylum[c(grep("unclassified_unclassified", manag_18s_data_soil$class))]
manag_18s_data_soil$order[which(is.na(manag_18s_data_soil$order))]<-paste("unclassified",
                                                                          manag_18s_data_soil$class[which(is.na(manag_18s_data_soil$order))],
                                                                          sep="_")
manag_18s_data_soil$order[which(grepl("unclassified_unclassified",
                                      manag_18s_data_soil$order))]<-manag_18s_data_soil$class[c(grep("unclassified_unclassified", manag_18s_data_soil$order))]
manag_18s_data_soil$family[which(is.na(manag_18s_data_soil$family))]<-paste("unclassified",
                                                                            manag_18s_data_soil$order[which(is.na(manag_18s_data_soil$family))],
                                                                            sep="_")
manag_18s_data_soil$family[which(grepl("unclassified_unclassified",
                                       manag_18s_data_soil$family))]<-manag_18s_data_soil$order[c(grep("unclassified_unclassified", manag_18s_data_soil$family))]
manag_18s_data_soil$genus[which(is.na(manag_18s_data_soil$genus))]<-paste("unclassified",
                                                                          manag_18s_data_soil$family[which(is.na(manag_18s_data_soil$genus))],
                                                                          sep="_")
manag_18s_data_soil$genus[which(grepl("unclassified_unclassified",
                                      manag_18s_data_soil$genus))]<-manag_18s_data_soil$family[c(grep("unclassified_unclassified", manag_18s_data_soil$genus))]

### Ordination of community composition -------------------------------------
manag_18s_soil_NMDS<-manag_18s_soil_normal %>%
  select("ESV", "Sample_034_01_M620F.M1041F.RC__Run2":"Sample_F03_05_M620F.M1041F.RC__Run3") %>%
  column_to_rownames("ESV") %>%
  t()%>%
  as.data.frame()


set.seed(4242)
manag_18s_soil_NMDS_plot<-as.matrix(manag_18s_soil_NMDS) %>%
  vegdist(., method = "bray") %>%
  metaMDS(., trymax=500) %>%
  scores(.) %>%
  as.data.frame() %>%
  mutate("Understory" = plotsite_metadata_soil$Understory,
         "Plot" = plotsite_metadata_soil$Plot,
         .before="NMDS1") %>%
  mutate(Understory = factor(Understory,
                             levels=c("Reduced",
                                      "Intermediate",
                                      "Enhanced",
                                      "Forest")),
         Plot = factor(Plot,
                       levels=c("034", "050", "052",
                                "081", "083", "095",
                                "132", "F01", "F02",
                                "F03")))

polygon_data_18s_soil<-manag_18s_soil_NMDS_plot[manag_18s_soil_NMDS_plot %>%
                                                  group_by(Understory) %>%
                                                  reframe(chull(NMDS1,NMDS2)) %>%
                                                  rename("Rows" = "chull(NMDS1, NMDS2)") %>%
                                                  mutate(Rows=ifelse(Understory=="Intermediate",
                                                                     Rows+18,
                                                                     Rows),
                                                         Rows=ifelse(Understory=="Enhanced",
                                                                     Rows+36,
                                                                     Rows),
                                                         Rows=ifelse(Understory=="Forest",
                                                                     Rows+42,
                                                                     Rows)) %>%
                                                  .$Rows,]

set.seed(42342)
manag_18s_soil_NMDS_envfit<-envfit(manag_18s_soil_NMDS_plot[,c("NMDS1",
                                                               "NMDS2")], 
                                   plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                   permutations=10000, 
                                   na.rm=T) %>%
  scores(., "vectors") %>%
  as.data.frame(.) %>%
  mutate("PVal"=c(envfit(manag_18s_soil_NMDS_plot[,c("NMDS1",
                                                     "NMDS2")], 
                         plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                         permutations=10000, 
                         na.rm=T)$vectors$pval),
         "RSq"=c(envfit(manag_18s_soil_NMDS_plot[,c("NMDS1",
                                                    "NMDS2")], 
                        plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                        permutations=10000, 
                        na.rm=T)$vectors$r))%>%
  filter(PVal<0.05)%>%
  mutate(NMDS1start=envfit(manag_18s_soil_NMDS_plot[,c("NMDS1",
                                                       "NMDS2")], 
                           plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                           permutations=10000, 
                           na.rm=T) %>%
           scores(., "factors") %>%
           as.data.frame()%>%
           select(NMDS1)%>%
           as.numeric(.),
         NMDS2start=envfit(manag_18s_soil_NMDS_plot[,c("NMDS1",
                                                       "NMDS2")], 
                           plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                           permutations=10000, 
                           na.rm=T) %>%
           scores(., "factors") %>%
           as.data.frame()%>%
           select(NMDS2)%>%
           as.numeric(.),
         NMDS1end=NMDS1+(envfit(manag_18s_soil_NMDS_plot[,c("NMDS1",
                                                            "NMDS2")], 
                                plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                permutations=10000, 
                                na.rm=T) %>%
                           scores(., "factors") %>%
                           as.data.frame()%>%
                           select(NMDS1)%>%
                           as.numeric(.)),
         NMDS2end=NMDS2+(envfit(manag_18s_soil_NMDS_plot[,c("NMDS1",
                                                            "NMDS2")], 
                                plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                permutations=10000, 
                                na.rm=T) %>%
                           scores(., "factors") %>%
                           as.data.frame()%>%
                           select(NMDS2)%>%
                           as.numeric(.)))%>%
  mutate(NMDS1end=NMDS1end/50,
         NMDS2end=NMDS2end/50)

ggplot(data=manag_18s_soil_NMDS_plot)+
  geom_point(size=5, shape=21,colour="black",
             aes(x=NMDS1, y=NMDS2, fill=Understory))+
  geom_segment(data = manag_18s_soil_NMDS_envfit, 
               aes(x = NMDS1start, y = NMDS2start,
                   xend = NMDS1end, yend = NMDS2end), 
               arrow = arrow(length = unit(0.1, "inches")),
               col="grey30", linewidth=1, alpha=1,) +
  geom_text_repel(data = manag_18s_soil_NMDS_envfit, 
                  aes(label = rownames(manag_18s_soil_NMDS_envfit),
                      x=NMDS1end, y=NMDS2end), 
                  box.padding = 1,       
                  point.padding = 0.5,     
                  segment.color = 'grey50',
                  max.overlaps = Inf,
                  size=5) + 
  geom_polygon(data = polygon_data_18s_soil, 
               aes(x = NMDS1, y = NMDS2,
                   fill=Understory),
               alpha=0.25,
               show.legend = F)+
  theme_bw()+
  theme(axis.text = element_text(size=),
        axis.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_colour_manual(name = "Treatment", values= c("Orange", "#DC267F",
                                                    "steelblue1", "#2E8B57"))+
  scale_fill_manual(name = "Treatment", values= c("Orange", "#DC267F",
                                                  "steelblue1", "#2E8B57")) 


## Combined soil ordination ------------------------------------------------
manag_soil_NMDS<-manag_ITS_soil_NMDS%>%
  rownames_to_column("X")%>%
  mutate(X=substr(X, 1, nchar(X) - 17))%>%
  column_to_rownames("X")%>%
  add_column(manag_16s_soil_NMDS %>%
               rownames_to_column("X") %>%
               mutate(X=substr(X, 1, nchar(X) - 16))%>%
               column_to_rownames("X"))%>%
  add_column(manag_18s_soil_NMDS %>%
               rownames_to_column("X") %>%
               mutate(X=substr(X, 1, nchar(X) - 22))%>%
               column_to_rownames("X"))%>%
  rename_with(~paste("ESV", seq_along(.), sep = "_"))

set.seed(193)
manag_soil_NMDS_plot<-as.matrix(manag_soil_NMDS) %>%
  vegdist(., method = "bray") %>%
  metaMDS(.) %>%
  scores(.) %>%
  as.data.frame() %>%
  mutate("Understory" = plotsite_metadata_soil$Understory,
         "Plot" = plotsite_metadata_soil$Plot,
         .before="NMDS1") %>%
  mutate(Understory = factor(Understory,
                             levels=c("Reduced",
                                      "Intermediate",
                                      "Enhanced",
                                      "Forest")),
         Plot = factor(Plot,
                       levels=c("034", "050", "052",
                                "081", "083", "095",
                                "132", "F01", "F02",
                                "F03")))

polygon_data_soil<-manag_soil_NMDS_plot[manag_soil_NMDS_plot %>%
                                          group_by(Understory) %>%
                                          reframe(chull(NMDS1,NMDS2)) %>%
                                          rename("Rows" = "chull(NMDS1, NMDS2)") %>%
                                          mutate(Rows=ifelse(Understory=="Intermediate",
                                                             Rows+18,
                                                             Rows),
                                                 Rows=ifelse(Understory=="Enhanced",
                                                             Rows+36,
                                                             Rows),
                                                 Rows=ifelse(Understory=="Forest",
                                                             Rows+42,
                                                             Rows)) %>%
                                          .$Rows,]

set.seed(42342)
manag_soil_NMDS_envfit<-envfit(manag_soil_NMDS_plot[,c("NMDS1",
                                                       "NMDS2")], 
                               plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                               permutations=10000, 
                               na.rm=T) %>%
  scores(., "vectors") %>%
  as.data.frame(.) %>%
  mutate("PVal"=c(envfit(manag_soil_NMDS_plot[,c("NMDS1",
                                                 "NMDS2")], 
                         plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                         permutations=10000, 
                         na.rm=T)$vectors$pval),
         "RSq"=c(envfit(manag_soil_NMDS_plot[,c("NMDS1",
                                                "NMDS2")], 
                        plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                        permutations=10000, 
                        na.rm=T)$vectors$r))%>%
  filter(PVal<0.05)%>%
  mutate(NMDS1start=envfit(manag_soil_NMDS_plot[,c("NMDS1",
                                                   "NMDS2")], 
                           plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                           permutations=10000, 
                           na.rm=T) %>%
           scores(., "factors") %>%
           as.data.frame()%>%
           select(NMDS1)%>%
           as.numeric(.),
         NMDS2start=envfit(manag_soil_NMDS_plot[,c("NMDS1",
                                                   "NMDS2")], 
                           plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                           permutations=10000, 
                           na.rm=T) %>%
           scores(., "factors") %>%
           as.data.frame()%>%
           select(NMDS2)%>%
           as.numeric(.),
         NMDS1end=NMDS1+(envfit(manag_soil_NMDS_plot[,c("NMDS1",
                                                        "NMDS2")], 
                                plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                permutations=10000, 
                                na.rm=T) %>%
                           scores(., "factors") %>%
                           as.data.frame()%>%
                           select(NMDS1)%>%
                           as.numeric(.)),
         NMDS2end=NMDS2+(envfit(manag_soil_NMDS_plot[,c("NMDS1",
                                                        "NMDS2")], 
                                plotsite_metadata_soil[,c(4:ncol(plotsite_metadata_soil))], 
                                permutations=10000, 
                                na.rm=T) %>%
                           scores(., "factors") %>%
                           as.data.frame()%>%
                           select(NMDS2)%>%
                           as.numeric(.)))

ggplot(data=manag_soil_NMDS_plot)+
  geom_point(size=5, shape=21,colour="black",
             aes(x=NMDS1, y=NMDS2, fill=Understory))+
  #geom_segment(data = manag_soil_NMDS_envfit, 
  #             aes(x = NMDS1start, y = NMDS2start, 
  #                 xend = NMDS1end, yend = NMDS2end), 
  #             arrow = arrow(length = unit(0.1, "inches")),
  #             col="grey30", linewidth=1, alpha=1,) +
  #geom_text_repel(data = manag_soil_NMDS_envfit, 
  #                aes(label = rownames(manag_soil_NMDS_envfit),
  #                    x=NMDS1end, y=NMDS2end), 
  #                box.padding = 1,       
  #                point.padding = 0.5,     
  #                segment.color = 'grey50',
  #                max.overlaps = Inf,
  #                size=5,
  #                vjust=-2,
  #                hjust=-1) + 
  geom_polygon(data = polygon_data_soil, 
               aes(x = NMDS1, y = NMDS2,
                   fill=Understory),
               alpha=0.25,
               show.legend = F)+
  theme_bw()+
  theme(axis.text = element_text(size=),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_colour_manual(name = "Treatment", values= c("Orange", "#DC267F",
                                                    "steelblue1", "#2E8B57"))+
  scale_fill_manual(name = "Treatment", values= c("Orange", "#DC267F",
                                                  "steelblue1", "#2E8B57"))


# Integrating Soil and Root Fungi -----------------------------------------------
## Dataframe setup ---------------------------------------------------------

manag_ITS_data_soil<-read.csv("path/to/fungal_dataset.csv")%>%
  select(-c(grep("Run1", names(.), value = TRUE),
            grep("Run3", names(.), value = TRUE))) #re-run relative to above code for the fungal soil dataset to remove forest samples (Run3), which didn't have an associated coconut root.
manag_ITS_data_root<-read.csv("path/to/fungal_dataset.csv")%>%
  select(-c(grep("Run2", names(.), value = TRUE),
            grep("Run3", names(.), value = TRUE)))

plotsite_metadata_soil<-data.frame(cbind(c(names(manag_ITS_data_soil)[which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2")]),
                                         c(rep("034", 6),
                                           rep("050", 6),
                                           rep("052", 6),
                                           rep("081", 6),
                                           rep("083", 6),
                                           rep("095", 6),
                                           rep("132", 6)),
                                         c(rep("Reduced", 18),
                                           rep("Intermediate", 18),
                                           rep("Enhanced", 6)))) %>%
  rename("Sample" = "X1",
         "Plot" = "X2",
         "Understory" = "X3") %>%
  mutate(Plot=factor(Plot,
                     levels=c("034", "050", "052",
                              "081", "083", "095",
                              "132")),
         Understory=factor(Understory,
                           levels=c("Reduced",
                                    "Intermediate",
                                    "Enhanced"))) %>%
  mutate(Biomass=c(read.csv("path/to/Fieldwork metadata.csv")%>%
                     filter(Understory!="Forest")%>%
                     .$Biomass.mean))%>%
  mutate(NMDS_centroid=NA)%>%
  mutate(NMDS_centroid=ifelse(Understory=="Reduced" | Understory =="Intermediate",
                              "Y",
                              "N"))

plotsite_metadata_root<-data.frame(cbind(c(names(manag_ITS_data_root)[which(names(manag_ITS_data_root)=="Sample_034_01_gITS7.ITS4__Run1"):which(names(manag_ITS_data_root)=="Sample_132_06_gITS7.ITS4__Run1")]),
                                         c(rep("034", 6),
                                           rep("050", 6),
                                           rep("052", 6),
                                           rep("081", 6),
                                           rep("083", 6),
                                           rep("095", 6),
                                           rep("132", 6)),
                                         c(rep("Reduced", 18),
                                           rep("Intermediate", 18),
                                           rep("Enhanced", 6)))) %>%
  rename("Sample" = "X1",
         "Plot" = "X2",
         "Understory" = "X3") %>%
  mutate(Plot=factor(Plot,
                     levels=c("034", "050", "052",
                              "081", "083", "095",
                              "132")),
         Understory=factor(Understory,
                           levels=c("Reduced",
                                    "Intermediate",
                                    "Enhanced"))) %>%
  mutate(Biomass=c(read.csv("path/to/Fieldwork metadata.csv")%>%
                     filter(Understory!="Forest")%>%
                     .$Biomass.mean))%>%
  mutate(NMDS_centroid=NA)%>%
  mutate(NMDS_centroid=ifelse(Understory=="Reduced" | Understory =="Intermediate",
                              "Y",
                              "N"))

for (e in c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2"))){
  manag_ITS_data_soil[,e]<-c(manag_ITS_data_soil[,e]-manag_ITS_data_soil$Sample_NegITS_gITS7.ITS4__Run2)
}
for (e in c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2"))){
  manag_ITS_data_soil[which(manag_ITS_data_soil[,e]<0),e]<-0
}
for (e in c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2"))){
  manag_ITS_data_soil[which(manag_ITS_data_soil[,e]==1),e]<-0
}
manag_ITS_data_soil<-manag_ITS_data_soil[,
                                         -which(names(manag_ITS_data_soil)=="Sample_NegITS_gITS7.ITS4__Run2")]
non_occur_ESV<-c()
for (e in c(1:nrow(manag_ITS_data_soil))){
  non_occur_ESV<-c(non_occur_ESV,
                   sum(manag_ITS_data_soil[e,
                                           c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2"))]))
}
manag_ITS_data_soil<-manag_ITS_data_soil[-c(which(non_occur_ESV==0)),]


for (e in c(which(names(manag_ITS_data_root)=="Sample_034_01_gITS7.ITS4__Run1"):which(names(manag_ITS_data_root)=="Sample_132_06_gITS7.ITS4__Run1"))){
  manag_ITS_data_root[,e]<-c(manag_ITS_data_root[,e]-manag_ITS_data_root$Sample_NegManagITS_gITS7.ITS4__Run1)
}
for (e in c(which(names(manag_ITS_data_root)=="Sample_034_01_gITS7.ITS4__Run1"):which(names(manag_ITS_data_root)=="Sample_132_06_gITS7.ITS4__Run1"))){
  manag_ITS_data_root[which(manag_ITS_data_root[,e]<0),e]<-0
}
for (e in c(which(names(manag_ITS_data_root)=="Sample_034_01_gITS7.ITS4__Run1"):which(names(manag_ITS_data_root)=="Sample_132_06_gITS7.ITS4__Run1"))){
  manag_ITS_data_root[which(manag_ITS_data_root[,e]==1),e]<-0
}
manag_ITS_data_root<-manag_ITS_data_root[,
                                         -which(names(manag_ITS_data_root)=="Sample_NegManagITS_gITS7.ITS4__Run1")]
non_occur_ESV<-c()
for (e in c(1:nrow(manag_ITS_data_root))){
  non_occur_ESV<-c(non_occur_ESV,
                   sum(manag_ITS_data_root[e,
                                           c(which(names(manag_ITS_data_root)=="Sample_034_01_gITS7.ITS4__Run1"):which(names(manag_ITS_data_root)=="Sample_132_06_gITS7.ITS4__Run1"))]))
}
manag_ITS_data_root<-manag_ITS_data_root[-c(which(non_occur_ESV==0)),]

manag_ITS_data_root<-subset(manag_ITS_data_root,
                            manag_ITS_data_root$kingdom=="Fungi")
manag_ITS_data_soil<-subset(manag_ITS_data_soil,
                            manag_ITS_data_soil$kingdom=="Fungi")


## Rarefaction and normalisation -------------------------------------------
manag_ITS_root_rarefaction<-manag_ITS_data_root[,
                                                c(1:which(names(manag_ITS_data_root)=="Sample_132_06_gITS7.ITS4__Run1"))] %>%
  remove_rownames(.) %>%
  column_to_rownames(.,
                     var="ESV") %>%
  t(.) %>%
  as.data.frame(.) %>%
  filter(!row_number() %in% c(1:3))%>%
  mutate(across(everything(), as.numeric))

manag_ITS_soil_rarefaction<-manag_ITS_data_soil[,
                                                c(1:which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2"))] %>%
  remove_rownames(.) %>%
  column_to_rownames(.,
                     var="ESV") %>%
  t(.) %>%
  as.data.frame(.) %>%
  filter(!row_number() %in% c(1:3))%>%
  mutate(across(everything(), as.numeric))

manag_ITS_soil_normal<-
  DESeqDataSetFromMatrix(countData = manag_ITS_soil_rarefaction %>%
                           t() %>%
                           as.data.frame(),
                         colData = plotsite_metadata_soil %>%
                           column_to_rownames("Sample"),
                         design = ~Understory) %>%
  DESeq2::estimateSizeFactors(.,type="poscounts")%>%
  #DESeq(.) %>% 
  DESeq2::counts(.,
                 normalized=TRUE)%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_ITS_data_soil %>%
              select(-c("Sample_034_01_gITS7.ITS4__Run2":"Sample_132_06_gITS7.ITS4__Run2")),
            by=c("ESV"))

manag_ITS_root_normal<-
  DESeqDataSetFromMatrix(countData = manag_ITS_root_rarefaction %>%
                           t() %>%
                           as.data.frame(),
                         colData = plotsite_metadata_root %>%
                           column_to_rownames("Sample"),
                         design = ~Understory) %>%
  DESeq2::estimateSizeFactors(.,type="poscounts")%>%
  #DESeq(.) %>% 
  DESeq2::counts(.,
                 normalized=TRUE)%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_ITS_data_root%>%
              select(-c("Sample_034_01_gITS7.ITS4__Run1":"Sample_132_06_gITS7.ITS4__Run1")),
            by=c("ESV"))

set.seed(23492)
manag_ITS_root_rarefied<-manag_ITS_root_rarefaction %>%
  rarefy_and_average(.,
                     min(sapply(manag_ITS_data_soil[c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2"))], 
                                sum),
                         sapply(manag_ITS_data_root[c(which(names(manag_ITS_data_root)=="Sample_034_01_gITS7.ITS4__Run1"):which(names(manag_ITS_data_root)=="Sample_132_06_gITS7.ITS4__Run1"))], 
                                sum)),
                     500) %>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_ITS_data_root%>%
              select(-c("Sample_034_01_gITS7.ITS4__Run1":"Sample_132_06_gITS7.ITS4__Run1")),
            by=c("ESV"))

set.seed(23492)
manag_ITS_soil_rarefied<-manag_ITS_soil_rarefaction %>%
  rarefy_and_average(.,
                     min(sapply(manag_ITS_data_soil[c(which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2"))], 
                                sum),
                         sapply(manag_ITS_data_root[c(which(names(manag_ITS_data_root)=="Sample_034_01_gITS7.ITS4__Run1"):which(names(manag_ITS_data_root)=="Sample_132_06_gITS7.ITS4__Run1"))], 
                                sum)),
                     500) %>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_ITS_data_soil%>%
              select(-c("Sample_034_01_gITS7.ITS4__Run2":"Sample_132_06_gITS7.ITS4__Run2")),
            by=c("ESV"))


## Making taxonomic assignment user-friendly -------------------------------
manag_ITS_root_normal$phylum[which(manag_ITS_root_normal$phylum=="unidentified")]<-"unclassified_Fungi"
manag_ITS_root_normal$class[which(manag_ITS_root_normal$phylum=="unclassified_Fungi")]<-manag_ITS_root_normal$phylum[which(manag_ITS_root_normal$phylum=="unclassified_Fungi")] 
manag_ITS_root_normal$class[c(grep("unidentified", manag_ITS_root_normal$class))]<-paste("unclassified",
                                                                                         manag_ITS_root_normal$phylum[c(grep("unidentified", manag_ITS_root_normal$class))],
                                                                                         sep="_")
manag_ITS_root_normal$class[c(grep("Incertae_sedis", manag_ITS_root_normal$class))]<-paste("unclassified",
                                                                                           manag_ITS_root_normal$phylum[c(grep("Incertae_sedis", manag_ITS_root_normal$class))],
                                                                                           sep="_")
manag_ITS_root_normal$order[which(is.na(manag_ITS_root_normal$order))]<-manag_ITS_root_normal$class[which(is.na(manag_ITS_root_normal$order))]
manag_ITS_root_normal$order[c(grep("unidentified", manag_ITS_root_normal$order))]<-paste("unclassified",
                                                                                         manag_ITS_root_normal$class[c(grep("unidentified", manag_ITS_root_normal$order))],
                                                                                         sep="_")
manag_ITS_root_normal$order[c(grep("Incertae_sedis", manag_ITS_root_normal$order))]<-paste("unclassified",
                                                                                           manag_ITS_root_normal$class[c(grep("Incertae_sedis", manag_ITS_root_normal$order))],
                                                                                           sep="_")
manag_ITS_root_normal$order[c(grep("unclassified_unclassified", manag_ITS_root_normal$order))]<-manag_ITS_root_normal$class[c(grep("unclassified_unclassified", manag_ITS_root_normal$order))]
manag_ITS_root_normal$family[which(is.na(manag_ITS_root_normal$family))]<-manag_ITS_root_normal$order[which(is.na(manag_ITS_root_normal$family))]
manag_ITS_root_normal$family[c(grep("unidentified", manag_ITS_root_normal$family))]<-paste("unclassified",
                                                                                           manag_ITS_root_normal$order[c(grep("unidentified", manag_ITS_root_normal$family))],
                                                                                           sep="_")
manag_ITS_root_normal$family[c(grep("Incertae_sedis", manag_ITS_root_normal$family))]<-paste("unclassified",
                                                                                             manag_ITS_root_normal$order[c(grep("Incertae_sedis", manag_ITS_root_normal$family))],
                                                                                             sep="_")
manag_ITS_root_normal$family[c(grep("unclassified_unclassified", manag_ITS_root_normal$family))]<-manag_ITS_root_normal$order[c(grep("unclassified_unclassified", manag_ITS_root_normal$family))]
manag_ITS_root_normal$genus[which(is.na(manag_ITS_root_normal$genus))]<-manag_ITS_root_normal$family[which(is.na(manag_ITS_root_normal$genus))]
manag_ITS_root_normal$genus[c(grep("unidentified", manag_ITS_root_normal$genus))]<-paste("unclassified",
                                                                                         manag_ITS_root_normal$family[c(grep("unidentified", manag_ITS_root_normal$genus))],
                                                                                         sep="_")
manag_ITS_root_normal$genus[c(grep("Incertae_sedis", manag_ITS_root_normal$genus))]<-paste("unclassified",
                                                                                           manag_ITS_root_normal$family[c(grep("Incertae_sedis", manag_ITS_root_normal$genus))],
                                                                                           sep="_")
manag_ITS_root_normal$genus[c(grep("unclassified_unclassified", manag_ITS_root_normal$genus))]<-manag_ITS_root_normal$family[c(grep("unclassified_unclassified", manag_ITS_root_normal$genus))]



manag_ITS_data_root$phylum[which(manag_ITS_data_root$phylum=="unidentified")]<-"unclassified_Fungi"
manag_ITS_data_root$class[which(manag_ITS_data_root$phylum=="unclassified_Fungi")]<-manag_ITS_data_root$phylum[which(manag_ITS_data_root$phylum=="unclassified_Fungi")] 
manag_ITS_data_root$class[c(grep("unidentified", manag_ITS_data_root$class))]<-paste("unclassified",
                                                                                     manag_ITS_data_root$phylum[c(grep("unidentified", manag_ITS_data_root$class))],
                                                                                     sep="_")
manag_ITS_data_root$class[c(grep("Incertae_sedis", manag_ITS_data_root$class))]<-paste("unclassified",
                                                                                       manag_ITS_data_root$phylum[c(grep("Incertae_sedis", manag_ITS_data_root$class))],
                                                                                       sep="_")
manag_ITS_data_root$order[which(is.na(manag_ITS_data_root$order))]<-manag_ITS_data_root$class[which(is.na(manag_ITS_data_root$order))]
manag_ITS_data_root$order[c(grep("unidentified", manag_ITS_data_root$order))]<-paste("unclassified",
                                                                                     manag_ITS_data_root$class[c(grep("unidentified", manag_ITS_data_root$order))],
                                                                                     sep="_")
manag_ITS_data_root$order[c(grep("Incertae_sedis", manag_ITS_data_root$order))]<-paste("unclassified",
                                                                                       manag_ITS_data_root$class[c(grep("Incertae_sedis", manag_ITS_data_root$order))],
                                                                                       sep="_")
manag_ITS_data_root$order[c(grep("unclassified_unclassified", manag_ITS_data_root$order))]<-manag_ITS_data_root$class[c(grep("unclassified_unclassified", manag_ITS_data_root$order))]
manag_ITS_data_root$family[which(is.na(manag_ITS_data_root$family))]<-manag_ITS_data_root$order[which(is.na(manag_ITS_data_root$family))]
manag_ITS_data_root$family[c(grep("unidentified", manag_ITS_data_root$family))]<-paste("unclassified",
                                                                                       manag_ITS_data_root$order[c(grep("unidentified", manag_ITS_data_root$family))],
                                                                                       sep="_")
manag_ITS_data_root$family[c(grep("Incertae_sedis", manag_ITS_data_root$family))]<-paste("unclassified",
                                                                                         manag_ITS_data_root$order[c(grep("Incertae_sedis", manag_ITS_data_root$family))],
                                                                                         sep="_")
manag_ITS_data_root$family[c(grep("unclassified_unclassified", manag_ITS_data_root$family))]<-manag_ITS_data_root$order[c(grep("unclassified_unclassified", manag_ITS_data_root$family))]
manag_ITS_data_root$genus[which(is.na(manag_ITS_data_root$genus))]<-manag_ITS_data_root$family[which(is.na(manag_ITS_data_root$genus))]
manag_ITS_data_root$genus[c(grep("unidentified", manag_ITS_data_root$genus))]<-paste("unclassified",
                                                                                     manag_ITS_data_root$family[c(grep("unidentified", manag_ITS_data_root$genus))],
                                                                                     sep="_")
manag_ITS_data_root$genus[c(grep("Incertae_sedis", manag_ITS_data_root$genus))]<-paste("unclassified",
                                                                                       manag_ITS_data_root$family[c(grep("Incertae_sedis", manag_ITS_data_root$genus))],
                                                                                       sep="_")
manag_ITS_data_root$genus[c(grep("unclassified_unclassified", manag_ITS_data_root$genus))]<-manag_ITS_data_root$family[c(grep("unclassified_unclassified", manag_ITS_data_root$genus))]

manag_ITS_soil_normal$phylum[which(manag_ITS_soil_normal$phylum=="unidentified")]<-"unclassified_Fungi"
manag_ITS_soil_normal$class[which(manag_ITS_soil_normal$phylum=="unclassified_Fungi")]<-manag_ITS_soil_normal$phylum[which(manag_ITS_soil_normal$phylum=="unclassified_Fungi")] 
manag_ITS_soil_normal$class[c(grep("unidentified", manag_ITS_soil_normal$class))]<-paste("unclassified",
                                                                                         manag_ITS_soil_normal$phylum[c(grep("unidentified", manag_ITS_soil_normal$class))],
                                                                                         sep="_")
manag_ITS_soil_normal$class[c(grep("Incertae_sedis", manag_ITS_soil_normal$class))]<-paste("unclassified",
                                                                                           manag_ITS_soil_normal$phylum[c(grep("Incertae_sedis", manag_ITS_soil_normal$class))],
                                                                                           sep="_")
manag_ITS_soil_normal$order[which(is.na(manag_ITS_soil_normal$order))]<-manag_ITS_soil_normal$class[which(is.na(manag_ITS_soil_normal$order))]
manag_ITS_soil_normal$order[c(grep("unidentified", manag_ITS_soil_normal$order))]<-paste("unclassified",
                                                                                         manag_ITS_soil_normal$class[c(grep("unidentified", manag_ITS_soil_normal$order))],
                                                                                         sep="_")
manag_ITS_soil_normal$order[c(grep("Incertae_sedis", manag_ITS_soil_normal$order))]<-paste("unclassified",
                                                                                           manag_ITS_soil_normal$class[c(grep("Incertae_sedis", manag_ITS_soil_normal$order))],
                                                                                           sep="_")
manag_ITS_soil_normal$order[c(grep("unclassified_unclassified", manag_ITS_soil_normal$order))]<-manag_ITS_soil_normal$class[c(grep("unclassified_unclassified", manag_ITS_soil_normal$order))]
manag_ITS_soil_normal$family[which(is.na(manag_ITS_soil_normal$family))]<-manag_ITS_soil_normal$order[which(is.na(manag_ITS_soil_normal$family))]
manag_ITS_soil_normal$family[c(grep("unidentified", manag_ITS_soil_normal$family))]<-paste("unclassified",
                                                                                           manag_ITS_soil_normal$order[c(grep("unidentified", manag_ITS_soil_normal$family))],
                                                                                           sep="_")
manag_ITS_soil_normal$family[c(grep("Incertae_sedis", manag_ITS_soil_normal$family))]<-paste("unclassified",
                                                                                             manag_ITS_soil_normal$order[c(grep("Incertae_sedis", manag_ITS_soil_normal$family))],
                                                                                             sep="_")
manag_ITS_soil_normal$family[c(grep("unclassified_unclassified", manag_ITS_soil_normal$family))]<-manag_ITS_soil_normal$order[c(grep("unclassified_unclassified", manag_ITS_soil_normal$family))]
manag_ITS_soil_normal$genus[which(is.na(manag_ITS_soil_normal$genus))]<-manag_ITS_soil_normal$family[which(is.na(manag_ITS_soil_normal$genus))]
manag_ITS_soil_normal$genus[c(grep("unidentified", manag_ITS_soil_normal$genus))]<-paste("unclassified",
                                                                                         manag_ITS_soil_normal$family[c(grep("unidentified", manag_ITS_soil_normal$genus))],
                                                                                         sep="_")
manag_ITS_soil_normal$genus[c(grep("Incertae_sedis", manag_ITS_soil_normal$genus))]<-paste("unclassified",
                                                                                           manag_ITS_soil_normal$family[c(grep("Incertae_sedis", manag_ITS_soil_normal$genus))],
                                                                                           sep="_")
manag_ITS_soil_normal$genus[c(grep("unclassified_unclassified", manag_ITS_soil_normal$genus))]<-manag_ITS_soil_normal$family[c(grep("unclassified_unclassified", manag_ITS_soil_normal$genus))]



manag_ITS_data_soil$phylum[which(manag_ITS_data_soil$phylum=="unidentified")]<-"unclassified_Fungi"
manag_ITS_data_soil$class[which(manag_ITS_data_soil$phylum=="unclassified_Fungi")]<-manag_ITS_data_soil$phylum[which(manag_ITS_data_soil$phylum=="unclassified_Fungi")] 
manag_ITS_data_soil$class[c(grep("unidentified", manag_ITS_data_soil$class))]<-paste("unclassified",
                                                                                     manag_ITS_data_soil$phylum[c(grep("unidentified", manag_ITS_data_soil$class))],
                                                                                     sep="_")
manag_ITS_data_soil$class[c(grep("Incertae_sedis", manag_ITS_data_soil$class))]<-paste("unclassified",
                                                                                       manag_ITS_data_soil$phylum[c(grep("Incertae_sedis", manag_ITS_data_soil$class))],
                                                                                       sep="_")
manag_ITS_data_soil$order[which(is.na(manag_ITS_data_soil$order))]<-manag_ITS_data_soil$class[which(is.na(manag_ITS_data_soil$order))]
manag_ITS_data_soil$order[c(grep("unidentified", manag_ITS_data_soil$order))]<-paste("unclassified",
                                                                                     manag_ITS_data_soil$class[c(grep("unidentified", manag_ITS_data_soil$order))],
                                                                                     sep="_")
manag_ITS_data_soil$order[c(grep("Incertae_sedis", manag_ITS_data_soil$order))]<-paste("unclassified",
                                                                                       manag_ITS_data_soil$class[c(grep("Incertae_sedis", manag_ITS_data_soil$order))],
                                                                                       sep="_")
manag_ITS_data_soil$order[c(grep("unclassified_unclassified", manag_ITS_data_soil$order))]<-manag_ITS_data_soil$class[c(grep("unclassified_unclassified", manag_ITS_data_soil$order))]
manag_ITS_data_soil$family[which(is.na(manag_ITS_data_soil$family))]<-manag_ITS_data_soil$order[which(is.na(manag_ITS_data_soil$family))]
manag_ITS_data_soil$family[c(grep("unidentified", manag_ITS_data_soil$family))]<-paste("unclassified",
                                                                                       manag_ITS_data_soil$order[c(grep("unidentified", manag_ITS_data_soil$family))],
                                                                                       sep="_")
manag_ITS_data_soil$family[c(grep("Incertae_sedis", manag_ITS_data_soil$family))]<-paste("unclassified",
                                                                                         manag_ITS_data_soil$order[c(grep("Incertae_sedis", manag_ITS_data_soil$family))],
                                                                                         sep="_")
manag_ITS_data_soil$family[c(grep("unclassified_unclassified", manag_ITS_data_soil$family))]<-manag_ITS_data_soil$order[c(grep("unclassified_unclassified", manag_ITS_data_soil$family))]
manag_ITS_data_soil$genus[which(is.na(manag_ITS_data_soil$genus))]<-manag_ITS_data_soil$family[which(is.na(manag_ITS_data_soil$genus))]
manag_ITS_data_soil$genus[c(grep("unidentified", manag_ITS_data_soil$genus))]<-paste("unclassified",
                                                                                     manag_ITS_data_soil$family[c(grep("unidentified", manag_ITS_data_soil$genus))],
                                                                                     sep="_")
manag_ITS_data_soil$genus[c(grep("Incertae_sedis", manag_ITS_data_soil$genus))]<-paste("unclassified",
                                                                                       manag_ITS_data_soil$family[c(grep("Incertae_sedis", manag_ITS_data_soil$genus))],
                                                                                       sep="_")
manag_ITS_data_soil$genus[c(grep("unclassified_unclassified", manag_ITS_data_soil$genus))]<-manag_ITS_data_soil$family[c(grep("unclassified_unclassified", manag_ITS_data_soil$genus))]


## Functional assignment ---------------------------------------------------
Taxonomy<-c()
for (e in c(1:nrow(manag_ITS_data_root))){
  Taxonomy<-c(Taxonomy,
              paste(manag_ITS_data_root$kingdom[e],
                    manag_ITS_data_root$phylum[e],
                    manag_ITS_data_root$class[e],
                    manag_ITS_data_root$order[e],
                    manag_ITS_data_root$family[e],
                    manag_ITS_data_root$genus[e],
                    manag_ITS_data_root$species[e],
                    sep=";"))
}
manag_ITS_root_funguild<-cbind(manag_ITS_data_root[,
                                                   c(which(names(manag_ITS_data_root)=="ESV"),
                                                     which(names(manag_ITS_data_root)=="Sample_034_01_gITS7.ITS4__Run1"):which(names(manag_ITS_data_root)=="Sample_132_06_gITS7.ITS4__Run1"))],
                               Taxonomy) %>%
  funguild_assign(.) %>%
  mutate(trophicMode = if_else(is.na(trophicMode) |
                                 confidenceRanking=="Possible" |
                                 trophicMode=="Pathotroph-Saprotroph" |
                                 trophicMode=="Pathotroph-Saprotroph-Symbiotroph" |
                                 trophicMode=="Saprotroph-Symbiotroph" |
                                 trophicMode=="Pathotroph-Symbiotroph" |
                                 trophicMode==" Pathotroph-Pathotroph-Saprotroph" |
                                 trophicMode==" Saprotroph",
                               "Unknown",
                               trophicMode))%>%  
  mutate(guild = if_else(is.na(guild) |
                           trophicMode=="Unknown",
                         "Unknown",
                         guild)) %>%
  mutate(guild = if_else(guild == "Lichenized-Wood Saprotroph",
                         "Lichenized", 
                         guild)) %>%
  mutate(guild = if_else(guild == "Litter Saprotroph-Plant Pathogen" |
                           guild == "|Plant Pathogen|" |
                           guild == "|Plant Pathogen|-Plant Saprotroph",
                         "Plant Pathogen", 
                         guild)) %>%
  mutate(guild = if_else(guild == "NULL",
                         "Undefined Pathotroph", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Plant Saprotroph|",
                         "Plant Saprotroph", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Wood Saprotroph|",
                         "Wood Saprotroph", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Fungal Parasite|",
                         "Fungal Parasite", 
                         guild)) %>%
  mutate(guild = if_else(guild == "Undefined Saprotroph-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Undefined Saprotroph-Wood Saprotroph" |
                           guild == "Soil Saprotroph-Undefined Saprotroph" |
                           guild == "Plant Saprotroph-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Plant Saprotroph-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Soil Saprotroph" |
                           guild == "|Plant Saprotroph|-Undefined Saprotroph" |
                           guild == "Undefined Saprotroph-|Wood Saprotroph|" |
                           guild == "|Undefined Saprotroph|" |
                           guild == "Litter Saprotroph-Soil Saprotroph-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Fungal Parasite-Undefined Saprotroph",
                         "Undefined Saprotroph",
                         guild))%>%
  left_join(.,
            manag_ITS_data_root %>%
              select(ESV, genus),
            by="ESV")%>%
  mutate(genus=ifelse(is.na(genus),
                      "Unknown",
                      genus),
         guild=ifelse(genus=="Fusarium" | genus=="Ganoderma" | genus=="Corticium" |
                        genus=="Ceratobasidium" | genus=="Lasiodiplodia" | genus=="Phyllosticta" | 
                        genus=="Mycosphaerella" | genus=="Melanconium" | genus=="Cercospora" | 
                        genus=="Alternaria" | genus=="Graphiola" | genus=="Rhizoctonia" | genus=="Graphium" |
                        genus=="Bipolaris" | genus=="Chalara" | genus=="Ceratocystis",
                      "Plant Pathogen",
                      guild),
         trophicMode=ifelse(genus=="Fusarium" | genus=="Ganoderma" | genus=="Corticium" |
                              genus=="Ceratobasidium" | genus=="Lasiodiplodia" | genus=="Phyllosticta" | 
                              genus=="Mycosphaerella" | genus=="Melanconium" | genus=="Cercospora" | 
                              genus=="Alternaria" | genus=="Graphiola" | genus=="Rhizoctonia" | genus=="Graphium" |
                              genus=="Bipolaris" | genus=="Chalara" | genus=="Ceratocystis",
                            "Pathotroph",
                            trophicMode),
         guild=ifelse(genus=="Castanediella",
                      "Unknown",
                      guild),
         trophicMode=ifelse(genus=="Castanediella",
                            "Unknown",
                            trophicMode)) 

Taxonomy<-c()
for (e in c(1:nrow(manag_ITS_data_soil))){
  Taxonomy<-c(Taxonomy,
              paste(manag_ITS_data_soil$kingdom[e],
                    manag_ITS_data_soil$phylum[e],
                    manag_ITS_data_soil$class[e],
                    manag_ITS_data_soil$order[e],
                    manag_ITS_data_soil$family[e],
                    manag_ITS_data_soil$genus[e],
                    manag_ITS_data_soil$species[e],
                    sep=";"))
}
manag_ITS_soil_funguild<-cbind(manag_ITS_data_soil[,
                                                   c(which(names(manag_ITS_data_soil)=="ESV"),
                                                     which(names(manag_ITS_data_soil)=="Sample_034_01_gITS7.ITS4__Run2"):which(names(manag_ITS_data_soil)=="Sample_132_06_gITS7.ITS4__Run2"))],
                               Taxonomy) %>%
  funguild_assign(.) %>%
  mutate(trophicMode = if_else(is.na(trophicMode) |
                                 confidenceRanking=="Possible" |
                                 trophicMode=="Pathotroph-Saprotroph" |
                                 trophicMode=="Pathotroph-Saprotroph-Symbiotroph" |
                                 trophicMode=="Saprotroph-Symbiotroph" |
                                 trophicMode=="Pathotroph-Symbiotroph" |
                                 trophicMode==" Pathotroph-Pathotroph-Saprotroph" |
                                 trophicMode==" Saprotroph",
                               "Unknown",
                               trophicMode))%>%  
  mutate(guild = if_else(is.na(guild) |
                           trophicMode=="Unknown",
                         "Unknown",
                         guild)) %>%
  mutate(guild = if_else(guild == "Lichenized-Wood Saprotroph",
                         "Lichenized", 
                         guild)) %>%
  mutate(guild = if_else(guild == "Litter Saprotroph-Plant Pathogen" |
                           guild == "|Plant Pathogen|" |
                           guild == "|Plant Pathogen|-Plant Saprotroph",
                         "Plant Pathogen", 
                         guild)) %>%
  mutate(guild = if_else(guild == "NULL" |
                           guild == "Plant Pathogen-Undefined Parasite-Undefined Saprotroph" |
                           guild == "|Animal Parasite|-Animal Pathogen-Arthropod Parasite",
                         "Undefined Pathotroph", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Plant Saprotroph|",
                         "Plant Saprotroph", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Lichen Parasite|",
                         "Lichen Parasite", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Wood Saprotroph|",
                         "Wood Saprotroph", 
                         guild)) %>%
  mutate(guild = if_else(guild == "|Fungal Parasite|",
                         "Fungal Parasite", 
                         guild)) %>%
  mutate(guild = if_else(guild == "Undefined Saprotroph-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Undefined Saprotroph-Wood Saprotroph" |
                           guild == "Soil Saprotroph-Undefined Saprotroph" |
                           guild == "Plant Saprotroph-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Plant Saprotroph-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Soil Saprotroph" |
                           guild == "Dung Saprotroph-Plant Saprotroph" |
                           guild == "Fungal Parasite-Wood Saprotroph" |
                           guild == "Plant Saprotroph-|Undefined Saprotroph|" |
                           guild == "Dung Saprotroph-Plant Saprotroph-|Wood Saprotroph|" |
                           guild == "|Plant Saprotroph|-Wood Saprotroph" |
                           guild == "|Plant Saprotroph|-Undefined Saprotroph" |
                           guild == "Undefined Saprotroph-|Wood Saprotroph|" |
                           guild == "|Undefined Saprotroph|" |
                           guild == "Dung Saprotroph-Plant Saprotroph-Soil Saprotroph" |
                           guild == "Litter Saprotroph-Soil Saprotroph-Wood Saprotroph" |
                           guild == "Dung Saprotroph-Fungal Parasite-Undefined Saprotroph",
                         "Undefined Saprotroph",
                         guild))%>%
  left_join(.,
            manag_ITS_data_soil %>%
              select(ESV, genus),
            by="ESV")%>%
  mutate(genus=ifelse(is.na(genus),
                      "Unknown",
                      genus),
         guild=ifelse(genus=="Fusarium" | genus=="Ganoderma" | genus=="Corticium" |
                        genus=="Ceratobasidium" | genus=="Lasiodiplodia" | genus=="Phyllosticta" | 
                        genus=="Mycosphaerella" | genus=="Melanconium" | genus=="Cercospora" | 
                        genus=="Alternaria" | genus=="Graphiola" | genus=="Rhizoctonia" | genus=="Graphium" |
                        genus=="Bipolaris" | genus=="Chalara" | genus=="Ceratocystis",
                      "Plant Pathogen",
                      guild),
         trophicMode=ifelse(genus=="Fusarium" | genus=="Ganoderma" | genus=="Corticium" |
                              genus=="Ceratobasidium" | genus=="Lasiodiplodia" | genus=="Phyllosticta" | 
                              genus=="Mycosphaerella" | genus=="Melanconium" | genus=="Cercospora" | 
                              genus=="Alternaria" | genus=="Graphiola" | genus=="Rhizoctonia" | genus=="Graphium" |
                              genus=="Bipolaris" | genus=="Chalara" | genus=="Ceratocystis",
                            "Pathotroph",
                            trophicMode),
         guild=ifelse(genus=="Castanediella",
                      "Unknown",
                      guild),
         trophicMode=ifelse(genus=="Castanediella",
                            "Unknown",
                            trophicMode))

manag_ITS_root_funguild_normal<-manag_ITS_root_funguild %>%
  filter(trophicMode!="Unknown") %>%
  select(ESV,
         Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1) %>%
  column_to_rownames("ESV")%>%
  DESeqDataSetFromMatrix(countData = .,
                         colData = plotsite_metadata_root %>%
                           column_to_rownames("Sample"),
                         design = ~Understory) %>%
  DESeq2::estimateSizeFactors(.,type="poscounts")%>%
  #DESeq(.) %>% 
  DESeq2::counts(.,
                 normalized=TRUE)%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_ITS_root_funguild %>%
              select(-c(Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1)),
            by="ESV")

manag_ITS_soil_funguild_normal<-manag_ITS_soil_funguild %>%
  filter(trophicMode!="Unknown") %>%
  select(ESV,
         Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2) %>%
  column_to_rownames("ESV")%>%
  DESeqDataSetFromMatrix(countData = .,
                         colData = plotsite_metadata_soil %>%
                           column_to_rownames("Sample"),
                         design = ~Understory) %>%
  DESeq2::estimateSizeFactors(.,type="poscounts")%>%
  #DESeq(.) %>% 
  DESeq2::counts(.,
                 normalized=TRUE)%>%
  as.data.frame()%>%
  rownames_to_column("ESV")%>%
  left_join(.,
            manag_ITS_soil_funguild %>%
              select(-c(Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2)),
            by="ESV")

###Mycorrhizal investigation --------------------------------------------
manag_ITS_root_AMF<-manag_ITS_root_funguild_normal %>%
  filter(guild=="Arbuscular Mycorrhizal") %>%
  select(ESV,
         Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1) %>%
  left_join(.,
            manag_ITS_root_funguild[,c(which(names(manag_ITS_root_funguild)=="ESV"),
                                       which(names(manag_ITS_root_funguild)=="guild"))],
            by="ESV") %>%
  left_join(.,
            manag_ITS_root_normal[,c(which(names(manag_ITS_root_normal)=="ESV"),
                                     which(names(manag_ITS_root_normal)=="class"),
                                     which(names(manag_ITS_root_normal)=="order"),
                                     which(names(manag_ITS_root_normal)=="family"),
                                     which(names(manag_ITS_root_normal)=="genus"))],
            by="ESV") %>%
  select(-guild) %>%
  group_by(family) %>%
  summarise(across(matches("^Sample_"), c(abundance=sum))) %>%
  pivot_longer(cols = c(2:ncol(.)),
               names_to = "Sample",
               values_to = "Abundance")%>%
  mutate(Sample=substr(Sample, 1, nchar(Sample) - 10)) %>%
  left_join(.,
            manag_ITS_root_funguild %>%
              select(-c("Sample_034_01_gITS7.ITS4__Run1":"Sample_132_06_gITS7.ITS4__Run1")) %>%
              left_join(.,
                        manag_ITS_root_rarefied %>%
                          select(ESV,
                                 Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1),
                        by="ESV") %>%
              relocate(ESV, Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1,
                       Taxonomy:citationSource)%>%
              filter(guild=="Arbuscular Mycorrhizal") %>%
              select(ESV,
                     Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1) %>%
              left_join(.,
                        manag_ITS_root_funguild[,c(which(names(manag_ITS_root_funguild)=="ESV"),
                                                   which(names(manag_ITS_root_funguild)=="guild"))],
                        by="ESV") %>%
              left_join(.,
                        manag_ITS_root_normal[,c(which(names(manag_ITS_root_normal)=="ESV"),
                                                 which(names(manag_ITS_root_normal)=="class"),
                                                 which(names(manag_ITS_root_normal)=="order"),
                                                 which(names(manag_ITS_root_normal)=="family"),
                                                 which(names(manag_ITS_root_normal)=="genus"))],
                        by="ESV") %>%
              group_by(family) %>%
              summarise(across(matches("^Sample"),
                               c(richness=~my_shannon(.)))) %>%
              pivot_longer(cols = c(2:ncol(.)),
                           names_to = "Sample",
                           values_to = "Richness")%>%
              mutate(Sample=substr(Sample, 1, nchar(Sample) - 9)),
            by=c("family", "Sample"))%>%
  mutate(Understory=c(rep(plotsite_metadata_root$Understory,
                          length(unique(family)))),
         .before = "Abundance") %>%
  pivot_longer(cols = c(Abundance:Richness),
               names_to = "Metric",
               values_to = "Value")

manag_ITS_soil_AMF<-manag_ITS_soil_funguild_normal %>%
  filter(guild=="Arbuscular Mycorrhizal") %>%
  select(ESV,
         Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2) %>%
  left_join(.,
            manag_ITS_soil_funguild[,c(which(names(manag_ITS_soil_funguild)=="ESV"),
                                       which(names(manag_ITS_soil_funguild)=="guild"))],
            by="ESV") %>%
  left_join(.,
            manag_ITS_soil_normal[,c(which(names(manag_ITS_soil_normal)=="ESV"),
                                     which(names(manag_ITS_soil_normal)=="class"),
                                     which(names(manag_ITS_soil_normal)=="order"),
                                     which(names(manag_ITS_soil_normal)=="family"),
                                     which(names(manag_ITS_soil_normal)=="genus"))],
            by="ESV")  %>%
  select(-guild) %>%
  group_by(family) %>%
  summarise(across(matches("^Sample_"), c(abundance=sum))) %>%
  pivot_longer(cols = c(2:ncol(.)),
               names_to = "Sample",
               values_to = "Abundance")%>%
  mutate(Sample=substr(Sample, 1, nchar(Sample) - 10)) %>%
  left_join(.,
            manag_ITS_soil_funguild %>%
              select(-c("Sample_034_01_gITS7.ITS4__Run2":"Sample_132_06_gITS7.ITS4__Run2")) %>%
              left_join(.,
                        manag_ITS_soil_rarefied %>%
                          select(ESV,
                                 Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2),
                        by="ESV") %>%
              relocate(ESV, Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2,
                       Taxonomy:citationSource)%>%
              filter(guild=="Arbuscular Mycorrhizal") %>%
              select(ESV,
                     Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2) %>%
              left_join(.,
                        manag_ITS_soil_funguild[,c(which(names(manag_ITS_soil_funguild)=="ESV"),
                                                   which(names(manag_ITS_soil_funguild)=="guild"))],
                        by="ESV") %>%
              left_join(.,
                        manag_ITS_soil_normal[,c(which(names(manag_ITS_soil_normal)=="ESV"),
                                                 which(names(manag_ITS_soil_normal)=="class"),
                                                 which(names(manag_ITS_soil_normal)=="order"),
                                                 which(names(manag_ITS_soil_normal)=="family"),
                                                 which(names(manag_ITS_soil_normal)=="genus"))],
                        by="ESV")  %>%
              select(-guild) %>%
              group_by(family) %>%
              summarise(across(matches("^Sample"),
                               c(richness=~my_shannon(.)))) %>%
              pivot_longer(cols = c(2:ncol(.)),
                           names_to = "Sample",
                           values_to = "Richness")%>%
              mutate(Sample=substr(Sample, 1, nchar(Sample) - 9)),
            by=c("family", "Sample"))%>%
  mutate(Understory=c(rep(plotsite_metadata_soil$Understory,
                          length(unique(family)))),
         .before = "Abundance") %>%
  pivot_longer(cols = c(Abundance:Richness),
               names_to = "Metric",
               values_to = "Value")

manag_ITS_root_AMF_colonisation<-manag_ITS_root_funguild %>%
  select(ESV, guild)%>%
  left_join(manag_ITS_root_rarefied %>%
              select(ESV, 
                     Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1),
            by="ESV")%>%
  filter(guild=="Arbuscular Mycorrhizal") %>%
  column_to_rownames("ESV") %>%
  select(Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1) %>%
  rownames_to_column("ESV") %>%
  summarise(across(matches("^Sample_"), list(ESVs = ~sum(. != 0), shannon = ~my_shannon(.)))) %>%
  t() %>%
  as.data.frame %>%
  mutate(Metric=c(rep(c("ESVs", "Shannons"),
                      42))) %>%
  rownames_to_column("Sample")%>%
  mutate(Sample=substr(Sample,1,24)) %>%
  pivot_wider(names_from = "Metric",
              values_from = "V1")%>%
  left_join(.,
            manag_ITS_soil_funguild %>%
              select(ESV, guild)%>%
              left_join(manag_ITS_soil_rarefied %>%
                          select(ESV, 
                                 Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2),
                        by="ESV")%>%
              filter(guild=="Arbuscular Mycorrhizal") %>%
              column_to_rownames("ESV") %>%
              select(Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2) %>%
              rownames_to_column("ESV") %>%
              summarise(across(matches("^Sample_"), list(ESVs = ~sum(. != 0), shannon = ~my_shannon(.)))) %>%
              t() %>%
              as.data.frame %>%
              mutate(Metric=c(rep(c("ESVs", "Shannons"),
                                  42))) %>%
              rownames_to_column("Sample") %>%
              mutate(Sample=substr(Sample,1,24)) %>%
              pivot_wider(names_from = "Metric",
                          values_from = "V1"),
            by="Sample") %>%
  rename("RootESVs" = "ESVs.x",
         "RootShannons" = "Shannons.x",
         "SoilESVs" = "ESVs.y",
         "SoilShannons" = "Shannons.y") %>%
  mutate(Understory=plotsite_metadata_root$Understory,
         .before="RootESVs") %>%
  mutate(ProportionESVs=RootESVs/SoilESVs) %>%
  mutate(ProportionShannon=RootShannons/SoilShannons)%>%
  select(Sample, RootShannons, SoilShannons, ProportionShannon)%>%
  pivot_longer(cols=c(2:ncol(.)),
               names_to = "Marker",
               values_to = "Value")%>%
  left_join(.,
            plotsite_metadata_root%>%
              select(Sample,Understory)%>%
              mutate(Sample=substr(Sample,
                                   1,24)),
            by="Sample")%>%
  mutate(Marker=factor(Marker,
                       levels=c("SoilShannons",
                                "RootShannons",
                                "ProportionShannon")))

ggplot(data=manag_ITS_root_AMF_colonisation%>%
         filter(Marker!="ProportionShannon"),
       aes(y=Understory, x=Value))+
  geom_boxplot(aes(fill=Understory),
               alpha=0.5, outlier.colour = NA)+
  geom_point(alpha=0.25,
             size=3,
             position=position_jitter(w=0,
                                      h=0.2))+
  scale_fill_manual(name="Understory", 
                    values=c("Orange", "#DC267F",
                             "steelblue1"))+
  scale_x_continuous(name="Shannon's Diversity")+
  scale_y_discrete(position="right")+
  facet_wrap(facets=vars(Marker), scales="free_y", ncol=1,
             labeller = labeller(Marker=c("SoilShannons" = "Soil",
                                          "RootShannons" = "Root",
                                          "ProportionShannon" = "Proportion of soil diversity")))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(size=20),
        strip.text = element_text(size=25),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.width = unit(2, "cm"))
ggplot(data=manag_ITS_root_AMF_colonisation%>%
         filter(Marker!="SoilShannons" &
                  Marker!="RootShannons"),
       aes(y=Understory, x=Value))+
  geom_boxplot(aes(fill=Understory),
               alpha=0.5, outlier.colour = NA)+
  geom_point(alpha=0.25,
             size=3,
             position=position_jitter(w=0,
                                      h=0.2))+
  scale_fill_manual(name="Understory", 
                    values=c("Orange", "#DC267F",
                             "steelblue1"))+
  scale_x_continuous(name="Shannon's Diversity",
                     breaks=c(0,1.0, 2.0),
                     limits = c(0,2))+
  scale_y_discrete(position="right")+
  facet_wrap(facets=vars(Marker), scales="free_y", ncol=1,
             labeller = labeller(Marker=c("SoilShannons" = "Soil",
                                          "RootShannons" = "Root",
                                          "ProportionShannon" = "Proportion of soil diversity")))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_text(size=25),
        legend.text = element_text(size=25),
        strip.text = element_text(size=25),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.width = unit(2, "cm"))

manag_ITS_guild_zscore<-manag_ITS_root_funguild_normal%>%
  group_by(guild) %>%
  summarise(across(matches("^Sample_"), c(abundance=sum))) %>%
  ungroup() %>%
  filter(guild=="Arbuscular Mycorrhizal") %>%
  rowwise() %>%
  mutate(across(matches("^Sample_"),
                ~ {
                  values <- c_across(matches("^Sample_"))
                  max_value <- max(values, na.rm = TRUE)
                  # Exclude max value for mean and SD calculation
                  trimmed_values <- values[values != max_value]
                  # If the current value is the max, assign NA, otherwise calculate z-score
                  if (. == max_value) {
                    NA_real_
                  } else {
                    (. - mean(trimmed_values, na.rm = TRUE)) / sd(trimmed_values, na.rm = TRUE)
                  }
                })) %>%
  #mutate(across(matches("^Sample_"),
  #              ~((. - mean(c_across(matches("^Sample_"))))/
  #                  sd(c_across(matches("^Sample_")))))) 
  ungroup() %>%
  pivot_longer(cols=c(2:ncol(.)),
               names_to = "Sample",
               values_to = "ZScore")%>%
  mutate(Understory=c(plotsite_metadata_root$Understory),
         .before=Sample)%>%
  rename("ZScoreRoot"="ZScore")%>%
  mutate(Sample=substr(Sample, 1, nchar(Sample)-17))%>%
  left_join(.,
            manag_ITS_soil_funguild_normal%>%
              group_by(guild) %>%
              summarise(across(matches("^Sample_"), c(abundance=sum))) %>%
              ungroup() %>%
              filter(guild=="Arbuscular Mycorrhizal") %>%
              rowwise() %>%
              mutate(across(matches("^Sample_"),
                            ~ {
                              values <- c_across(matches("^Sample_"))
                              max_value <- max(values, na.rm = TRUE)
                              # Exclude max value for mean and SD calculation
                              trimmed_values <- values[values != max_value]
                              # If the current value is the max, assign NA, otherwise calculate z-score
                              if (. == max_value) {
                                NA_real_
                              } else {
                                (. - mean(trimmed_values, na.rm = TRUE)) / sd(trimmed_values, na.rm = TRUE)
                              }
                            })) %>%
              #mutate(across(matches("^Sample_"),
              #              ~((. - mean(c_across(matches("^Sample_"))))/
              #                  sd(c_across(matches("^Sample_")))))) 
              ungroup() %>%
              pivot_longer(cols=c(2:ncol(.)),
                           names_to = "Sample",
                           values_to = "ZScore")%>%
              mutate(Understory=c(plotsite_metadata_root$Understory),
                     .before=Sample)%>%
              rename("ZScoreSoil"="ZScore")%>%
              mutate(Sample=substr(Sample, 1, nchar(Sample)-17)),
            by=c("guild", "Understory", "Sample"))%>%
  mutate(ZScoreDiff= ZScoreRoot - ZScoreSoil) %>% #proportion obviously doesn't work here
  select(-Sample) %>%
  group_by(guild, Understory) %>%
  summarise(ZScoreSoil = mean(ZScoreSoil, na.rm = TRUE),
            ZScoreRoot = mean(ZScoreRoot, na.rm = TRUE),
            ZScoreDiff = mean(ZScoreDiff, na.rm = TRUE)) %>%
  ungroup()

ggplot(data=manag_ITS_guild_zscore,
       aes(x=Understory, y=guild, fill=ZScoreRoot))+ #Change fill= depending on which heatmap panel you are after.
  geom_tile(colour="black")+
  geom_text(aes(label=c("High", "Intermediate", "Low")),
            size=10)+
  scale_x_discrete(expand=c(0,0),
                   position="bottom")+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradient2(name="Abundance Z-Score",
                       low="#0000FF",
                       mid = "white",
                       high="#FF0000",
                       guide = guide_colorbar(title.position = "top",
                                              title.hjust = 0.5))+
  theme_bw()+
  theme(axis.line = element_blank()) +
  theme_bw()+
  theme(legend.direction = "horizontal",
        #legend.box.background = element_rect(colour = "black"),
        legend.position = "bottom",
        legend.justification = "centre",
        legend.key.size = unit(2.5,"cm"),
        legend.title = element_text(size = 20),
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())+
  coord_fixed(ratio=1)+
  geom_hline(yintercept=c(6.5,
                          3.5),
             linewidth=1)


### Pathogen investigation --------------------------------------------------
manag_ITS_root_pathogen<-manag_ITS_root_funguild_normal %>%
  filter(guild=="Plant Pathogen") %>%
  select(ESV,
         Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1) %>%
  left_join(.,
            manag_ITS_root_funguild[,c(which(names(manag_ITS_root_funguild)=="ESV"),
                                       which(names(manag_ITS_root_funguild)=="guild"))],
            by="ESV") %>%
  left_join(.,
            manag_ITS_root_normal[,c(which(names(manag_ITS_root_normal)=="ESV"),
                                     which(names(manag_ITS_root_normal)=="class"),
                                     which(names(manag_ITS_root_normal)=="order"),
                                     which(names(manag_ITS_root_normal)=="family"),
                                     which(names(manag_ITS_root_normal)=="genus"))],
            by="ESV") %>%
  group_by(genus) %>%
  summarise(across(matches("^Sample_"), c(abundance=sum))) %>%
  pivot_longer(cols = c(2:ncol(.)),
               names_to = "Sample",
               values_to = "Abundance")%>%
  mutate(Sample=substr(Sample, 1, nchar(Sample) - 10)) %>%
  left_join(.,
            manag_ITS_root_funguild %>%
              select(-c("Sample_034_01_gITS7.ITS4__Run1":"Sample_132_06_gITS7.ITS4__Run1")) %>%
              left_join(.,
                        manag_ITS_root_rarefied %>%
                          select(ESV,
                                 Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1),
                        by="ESV") %>%
              relocate(ESV, Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1,
                       Taxonomy:citationSource)%>%
              filter(guild=="Plant Pathogen") %>%
              select(ESV,
                     Sample_034_01_gITS7.ITS4__Run1:Sample_132_06_gITS7.ITS4__Run1) %>%
              left_join(.,
                        manag_ITS_root_funguild[,c(which(names(manag_ITS_root_funguild)=="ESV"),
                                                   which(names(manag_ITS_root_funguild)=="guild"))],
                        by="ESV") %>%
              left_join(.,
                        manag_ITS_root_normal[,c(which(names(manag_ITS_root_normal)=="ESV"),
                                                 which(names(manag_ITS_root_normal)=="class"),
                                                 which(names(manag_ITS_root_normal)=="order"),
                                                 which(names(manag_ITS_root_normal)=="family"),
                                                 which(names(manag_ITS_root_normal)=="genus"))],
                        by="ESV") %>%
              select(-guild) %>%
              group_by(genus) %>%
              summarise(across(matches("^Sample"),
                               c(richness=~my_shannon(.)))) %>%
              pivot_longer(cols = c(2:ncol(.)),
                           names_to = "Sample",
                           values_to = "Richness")%>%
              mutate(Sample=substr(Sample, 1, nchar(Sample) - 9)),
            by=c("genus", "Sample"))%>%
  mutate(Understory=c(rep(plotsite_metadata_root$Understory,
                          length(unique(genus)))),
         .before = "Abundance") %>%
  pivot_longer(cols = c(Abundance:Richness),
               names_to = "Metric",
               values_to = "Value")

manag_ITS_soil_pathogen<-manag_ITS_soil_funguild_normal %>%
  filter(guild=="Plant Pathogen") %>%
  select(ESV,
         Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2) %>%
  left_join(.,
            manag_ITS_soil_funguild[,c(which(names(manag_ITS_soil_funguild)=="ESV"),
                                       which(names(manag_ITS_soil_funguild)=="guild"))],
            by="ESV") %>%
  left_join(.,
            manag_ITS_soil_normal[,c(which(names(manag_ITS_soil_normal)=="ESV"),
                                     which(names(manag_ITS_soil_normal)=="class"),
                                     which(names(manag_ITS_soil_normal)=="order"),
                                     which(names(manag_ITS_soil_normal)=="family"),
                                     which(names(manag_ITS_soil_normal)=="genus"))],
            by="ESV") %>%
  select(-guild) %>%
  group_by(genus) %>%
  summarise(across(matches("^Sample_"), c(abundance=sum))) %>%
  pivot_longer(cols = c(2:ncol(.)),
               names_to = "Sample",
               values_to = "Abundance")%>%
  mutate(Sample=substr(Sample, 1, nchar(Sample) - 10)) %>%
  left_join(.,
            manag_ITS_soil_funguild %>%
              select(-c("Sample_034_01_gITS7.ITS4__Run2":"Sample_132_06_gITS7.ITS4__Run2")) %>%
              left_join(.,
                        manag_ITS_soil_rarefied %>%
                          select(ESV,
                                 Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2),
                        by="ESV") %>%
              relocate(ESV, Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2,
                       Taxonomy:citationSource)%>%
              filter(guild=="Plant Pathogen") %>%
              select(ESV,
                     Sample_034_01_gITS7.ITS4__Run2:Sample_132_06_gITS7.ITS4__Run2) %>%
              left_join(.,
                        manag_ITS_soil_funguild[,c(which(names(manag_ITS_soil_funguild)=="ESV"),
                                                   which(names(manag_ITS_soil_funguild)=="guild"))],
                        by="ESV") %>%
              left_join(.,
                        manag_ITS_soil_normal[,c(which(names(manag_ITS_soil_normal)=="ESV"),
                                                 which(names(manag_ITS_soil_normal)=="class"),
                                                 which(names(manag_ITS_soil_normal)=="order"),
                                                 which(names(manag_ITS_soil_normal)=="family"),
                                                 which(names(manag_ITS_soil_normal)=="genus"))],
                        by="ESV") %>%
              group_by(genus) %>%
              summarise(across(matches("^Sample"),
                               c(richness=~my_shannon(.)))) %>%
              pivot_longer(cols = c(2:ncol(.)),
                           names_to = "Sample",
                           values_to = "Richness")%>%
              mutate(Sample=substr(Sample, 1, nchar(Sample) - 9)),
            by=c("genus", "Sample"))%>%
  mutate(Understory=c(rep(plotsite_metadata_soil$Understory,
                          length(unique(genus)))),
         .before = "Abundance") %>%
  pivot_longer(cols = c(Abundance:Richness),
               names_to = "Metric",
               values_to = "Value")

manag_ITS_root_pathogen%>%
  filter(Metric=="Abundance") %>%
  filter(genus=="Pestalotiopsis" | genus=="Fusarium" | genus=="Lasiodiplodia") %>%
  add_row(manag_ITS_soil_pathogen %>%
            filter(Metric=="Abundance") %>%
            filter(genus=="Pestalotiopsis" | genus=="Fusarium" | genus=="Lasiodiplodia"),
          .) %>%
  mutate(Axis=c(paste(genus, Sample,sep="_"))) %>%
  mutate(Type=c(rep("Soil",
                    length(grep("Run2", Sample))),
                rep("Root",
                    length(grep("Run1", Sample))))) %>%
  mutate(Label=c(substr(Sample, 8, 13)))%>%
  mutate(Axis=factor(Axis,
                     levels=Axis[order(-Value)])) %>%
  mutate(Type=factor(Type,
                     levels=c("Soil",
                              "Root")))%>%
  mutate(genus = factor(genus,
                        levels = c("Pestalotiopsis", "Lasiodiplodia", "Fusarium"),
                        labels = c("italic(Pestalotiopsis)", "italic(Lasiodiplodia)", "italic(Fusarium)")))%>%
  mutate(ManagementIntensity = NA)%>%
  mutate(ManagementIntensity=ifelse(Understory=="Reduced",
                                    "High",
                                    ManagementIntensity),
         ManagementIntensity=ifelse(Understory=="Intermediate",
                                    "Intermediate",
                                    ManagementIntensity),
         ManagementIntensity=ifelse(Understory=="Enhanced",
                                    "Low",
                                    ManagementIntensity))%>%
  mutate(ManagementIntensity=factor(ManagementIntensity,
                                    levels=c("High", "Intermediate",
                                             "Low", "Forest")))%>%
  #add_row(data.frame(
  #  Axis = "Pestalotiopsis_Sample_132_06_gITS7.ITS4__Run1",
  #  Value = 0.1,
  #  Type = "Soil",  # or "Root", doesn't matter
  #  genus = "italic(Pestalotiopsis)",
  #  ManagementIntensity = factor("Forest", levels=c("High", "Intermediate", "Low", "Forest"))
  #))%>% #unhash if you need forest in the legend.
  ggplot(data=.,
         aes(x=Axis, y=Value, fill=ManagementIntensity))+
  geom_bar(data=.%>%
             filter(Value!=0),
           aes(fill = ManagementIntensity, y = 100000),
           width = 1, stat = "identity",
           alpha=0.5) +
  geom_bar(data=.%>%
             filter(Value!=0),
           aes(fill = ManagementIntensity, y = -5000), 
           width = 1, stat = "identity",
           alpha=0.5) +
  geom_bar(data=.%>%
             filter(Value==0),
           aes(y = 100000),
           fill="grey80",
           width = 1, stat = "identity",
           alpha=0.5) +
  geom_bar(data=.%>%
             filter(Value==0),
           aes(y = -5000), 
           fill="grey80",
           width = 1, stat = "identity",
           alpha=0.5) +
  geom_point(shape=21,
             size=3)+
  geom_point(data = . %>% filter(Value == 0),
             aes(x = Axis, y = Value),
             shape = 21, size = 3, fill = "grey80") +
  coord_cartesian(ylim = c(-5, 6000))+
  #scale_y_continuous(name="Relative Abundance",
  #                   breaks = c(0,0.2,0.4,0.6,0.8))+
  scale_x_discrete(name="Ranked Sample")+
  scale_y_continuous(name="Abundance", breaks=c(0,2500,5000))+
  facet_grid2(cols = vars(genus),
              rows = vars(Type),
              scales = "free",
              independent = "x",
              labeller = labeller(
                genus = label_parsed,  # Italicize genus names
                Type = label_value     # Keep Type as plain text
              ))+
  theme(axis.text.x=element_blank())+
  scale_fill_manual(name="Management Intensity",
                    values=c("Orange", "#DC267F",
                             "steelblue1", "#2E8B57"),
                    guide = guide_legend(
                      title.position = "top",  # puts title above items
                      title.hjust = 0.5  ))+
  theme_bw()+
  theme(legend.position = "none",
        legend.text = element_text(size=25),
        legend.title = element_text(size=30,
                                    margin = margin(b = 20)),
        #legend.box.background = element_rect(colour = "black"),
        legend.key.size = unit(2,"cm"),
        legend.key.width = unit(3, "cm"),
        #legend.title.margin = margin(b = 30),
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        axis.title = element_text(size=25),
        strip.text = element_text(size=25),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept=0.025)



# Production investigation ------------------------------------------------

## Production data ---------------------------------------------------------
metadata<-read.csv("path/to/Fieldwork metadata.csv") %>%
  select(c(1:6))%>%
  mutate(TreeID=c(paste(metadata$Plot,
                        metadata$Line,
                        metadata$Tree,
                        sep="_")))
set.seed(3090)
Required_trees<-read.csv("path/to/Required trees.csv")%>%
  mutate(Sample=c(paste(Plot,
                        Subplot,
                        sep = "_0")),
         TreeID=c(paste(Plot,
                        Line,
                        Tree,
                        sep="_")))%>%
  left_join(.,
            read.csv("path/to/dead_trees.csv")%>%
              select(TreeID, STATUTS)%>%
              rename("Status"="STATUTS"),
            by="TreeID")%>%
  mutate(Status=ifelse(is.na(Status),
                       "Alive",
                       "Dead"),
         Status=ifelse(TreeID %in% c("52_30_12",
                                     "52_30_14",
                                     "52_31_12",
                                     "52_32_12",
                                     "81_45_17"),
                       "Problem",
                       Status))%>%
  mutate(Duplicate=ifelse(duplicated(TreeID) | duplicated(TreeID, fromLast = TRUE),
                          "Y",
                          "N"),
         SampledTree=ifelse(TreeID %in% metadata$TreeID,
                            "Y",
                            "N"))%>%
  filter(Status=="Alive")%>%
  filter(Duplicate=="N")%>%
  group_by(Sample) %>%
  mutate(
    Existing_Y = sum(SampledTree == "Y"),
    Additional_Y_needed = max(0, 4 - Existing_Y),
    Total_trees_in_sample = n(),
    SampledTree = case_when(SampledTree == "Y" ~ "Y",
                            Additional_Y_needed > 0 &
                              SampledTree == "N" &
                              row_number() %in% sample(which(SampledTree == "N"),
                                                       size = min(Additional_Y_needed, sum(SampledTree == "N")),
                                                       replace = FALSE) ~ "Y",
                            TRUE ~ "N")) %>%
  ungroup() %>%
  select(-Existing_Y, -Additional_Y_needed, -Total_trees_in_sample)%>%
  filter(SampledTree=="Y")

production_data<-
  read.csv("path/to/production_034_1998.csv")%>%
  bind_rows(.,
            read.csv("path/to/production_050_2002.csv")%>%
              select(-c(X, X.1)),
            read.csv("path/to/production_052_2002.csv")%>%
              select(-Sent_later.),
            read.csv("path/to/production_083_1997.csv"),
            read.csv("path/to/production_081_2002.csv") %>%
              select(-c(X)),
            read.csv("path/to/production_095_2005.csv"))%>%
  mutate(TreeID=paste(field,
                      row,
                      tree,
                      sep="_"))%>%
  filter(TreeID %in% Required_trees$TreeID)%>%
  arrange(TreeID, month)%>%
  group_by(TreeID)%>%
  mutate(Months=n(),
         Month=c(1:n()))%>%
  ungroup()%>%
  mutate(month=factor(month,
                      levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
                               "10", "11", "12")),
         Month=factor(Month,
                      levels=c("1", "2", "3", "4", "5", "6")))%>%
  group_by(TreeID) %>%
  mutate(Coverage = paste(month, collapse = "_"))%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Plot_Coverage = paste(field,Coverage, collapse="_"))%>%
  ungroup()

production_data_mean<-production_data%>%
  group_by(TreeID)%>%
  summarise(AnnualNuts=sum(nb_fr),
            MeanNuts=mean(nb_fr))%>%
  left_join(production_data,
            .,
            by="TreeID")%>%
  rename("Plot" = "field")%>%
  mutate(Plot=paste("0", Plot, sep=""))%>%
  left_join(.,
            Required_trees%>%
              select(TreeID, Sample),
            by="TreeID")%>%
  left_join(.,
            data.frame(Plot=c("034", "050", "052",
                              "081", "083", "095"),
                       Understory=c(rep("Reduced", 3),
                                    rep("Intermediate", 3))),
            by="Plot")%>%
  group_by(Sample) %>%
  mutate(Sample_meannuts = mean(MeanNuts),
         Understory=factor(Understory,
                           levels=c("Reduced", "Intermediate")),
         Sample_meanannual = mean(AnnualNuts))

ggplot(data=production_data_mean,
       aes(x=Understory, y=MeanNuts, fill=Understory))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size=2.5, alpha=0.02)+
  scale_fill_manual(values=c("orange", "#DC267F"))+
  scale_y_continuous(name = "Coconuts produced")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=10))

Sample_level_production<-production_data_mean %>%
  group_by(Plot,Sample, Understory) %>%
  summarise(Sample_meannuts = mean(Sample_meannuts, na.rm = TRUE), .groups = "drop")%>%
  mutate(Sample=paste("0",
                      Sample,
                      sep=""))

Tree_level_production<-production_data_mean %>%
  group_by(TreeID,Plot, Sample, Understory) %>%
  summarise(MeanNuts = mean(MeanNuts, na.rm = TRUE), .groups = "drop")


## Pestalotiopsis-production relationship ----------------------------------
manag_ITS_soil_normal_clean<-manag_ITS_data_soil%>%
  remove_rownames(.)%>%
  select(ESV,Sample_034_01_gITS7.ITS4__Run2:Sample_095_06_gITS7.ITS4__Run2)%>%
  column_to_rownames("ESV")%>%
  DESeqDataSetFromMatrix(countData = .,
                         colData = data.frame(Sample=names(manag_ITS_data_soil%>%
                                                             select(ESV,Sample_034_01_gITS7.ITS4__Run2:Sample_095_06_gITS7.ITS4__Run2)%>%
                                                             remove_rownames()%>%
                                                             column_to_rownames("ESV")),
                                              Plot=c(rep("034", 6),
                                                     rep("050", 6),
                                                     rep("052", 6),
                                                     rep("081", 6),
                                                     rep("083", 6),
                                                     rep("095", 6)),
                                              Understory=c(rep("Reduced", 18),
                                                           rep("Intermediate", 18)))%>%
                           column_to_rownames("Sample")%>%
                           mutate(Understory=factor(Understory,
                                                    levels=c("Reduced",
                                                             "Intermediate"))),
                         design = ~Understory) %>%
  DESeq2::estimateSizeFactors(.,type="poscounts")%>%
  #DESeq(.) %>% 
  DESeq2::counts(.,
                 normalized=TRUE)%>%
  as.data.frame() %>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("Sample")%>%
  mutate(Sample=substr(Sample,
                       8,
                       nchar(Sample)-17))%>%
  column_to_rownames("Sample")%>%
  t()%>%
  as.data.frame()%>%
  filter(rowSums(.)!=0)%>%
  t()%>%
  as.data.frame()

manag_ITS_root_normal_clean<-manag_ITS_data_root%>%
  remove_rownames(.)%>%
  select(ESV,Sample_034_01_gITS7.ITS4__Run1:Sample_095_06_gITS7.ITS4__Run1)%>%
  column_to_rownames("ESV")%>%
  DESeqDataSetFromMatrix(countData = .,
                         colData = data.frame(Sample=names(manag_ITS_data_root%>%
                                                             select(ESV,Sample_034_01_gITS7.ITS4__Run1:Sample_095_06_gITS7.ITS4__Run1)%>%
                                                             remove_rownames()%>%
                                                             column_to_rownames("ESV")),
                                              Plot=c(rep("034", 6),
                                                     rep("050", 6),
                                                     rep("052", 6),
                                                     rep("081", 6),
                                                     rep("083", 6),
                                                     rep("095", 6)),
                                              Understory=c(rep("Reduced", 18),
                                                           rep("Intermediate", 18)))%>%
                           column_to_rownames("Sample")%>%
                           mutate(Understory=factor(Understory,
                                                    levels=c("Reduced",
                                                             "Intermediate"))),
                         design = ~Understory) %>%
  DESeq2::estimateSizeFactors(.,type="poscounts")%>%
  #DESeq(.) %>% 
  DESeq2::counts(.,
                 normalized=TRUE)%>%
  as.data.frame() %>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("Sample")%>%
  mutate(Sample=substr(Sample,
                       8,
                       nchar(Sample)-17))%>%
  column_to_rownames("Sample")%>%
  t()%>%
  as.data.frame()%>%
  filter(rowSums(.)!=0)%>%
  t()%>%
  as.data.frame()

overall_biodiversity_dataframe_normal<-manag_ITS_soil_normal_clean%>%
  add_column(.,
             manag_ITS_root_normal_clean)%>%
  rename_with(~paste("ESV", seq_along(.), sep = "_"))

overall_ESV_info<-manag_ITS_data_soil%>%
  remove_rownames()%>%
  select(ESV, kingdom:genus)%>%
  mutate(Sample="Soil",
         Marker="ITS")%>%
  filter(ESV %in% names(manag_ITS_soil_normal_clean))%>%
  add_row(.,
          manag_ITS_data_root%>%
            remove_rownames()%>%
            select(ESV, kingdom:genus)%>%
            mutate(Sample="Root",
                   Marker="ITS")%>%
            filter(ESV %in% names(manag_ITS_root_normal_clean)))%>%
  t()%>%
  as.data.frame()%>%
  rename_with(~paste("ESV", seq_along(.), sep = "_"))%>%
  t()%>%
  as.data.frame()%>%
  select(-ESV)%>%
  rownames_to_column("ESV")

genus_correlation_investigation<-
  overall_biodiversity_dataframe_normal %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "ESV",
    values_to = "abundance"
  ) %>%
  full_join(
    overall_ESV_info %>%
      select(ESV, genus, Sample),
    by = "ESV"
  ) %>%
  group_by(SampleID, genus, Sample) %>%
  summarise(total_abundance = sum(as.numeric(abundance), na.rm = TRUE), .groups = "drop")%>%
  mutate(Genus=paste(genus, Sample,sep="_"))%>%
  select(-c(genus, Sample))%>%
  pivot_wider(names_from = "Genus", values_from = "total_abundance")%>%
  column_to_rownames("SampleID")

genus_correlation_investigation%>%
  select("Pestalotiopsis_Soil")%>%
  rownames_to_column("Sample")%>%
  full_join(Sample_level_production%>%
              select(Sample, Understory, Sample_meannuts))%>%
  mutate(Understory=ifelse(Understory=="Reduced",
                           "High Intensity",
                           "Intermediate Intensity"))%>%
  mutate(Understory=factor(Understory,
                           levels = c("High Intensity",
                                      "Intermediate Intensity")))%>%
  #filter(Sample!="034_05")%>%
  #filter(Sample!="050_02")%>%
  ggplot(data=.,
         aes(x=log(Pestalotiopsis_Soil),
             y=Sample_meannuts))+
  geom_point(shape=21, size=5, aes(fill=Understory))+
  geom_smooth(method = "lm",
              alpha=0.2, colour="grey40", 
              lty=2, se=T)+
  facet_wrap(facets=vars(Understory), scales="free_x")+
  scale_x_continuous(name=expression(""*italic("Pestalotiopsis")*" Abundance (Log)"))+
  scale_y_continuous(name="Coconuts produced")+
  scale_fill_manual(values=c("orange", "#DC267F"))+
  theme_bw()+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=12.5),
        strip.text = element_text(size=25),
        panel.spacing = unit(1.5, "lines"),
        legend.position = "none")








