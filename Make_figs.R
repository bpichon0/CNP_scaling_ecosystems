rm(list=ls())
source("CNP_lake_functions.R")

Scale_numeric=function(x){
  d[,colnames(dplyr::select_if(d, is.numeric))]=
    apply(d[,colnames(dplyr::select_if(d, is.numeric))],2,
          function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm = T))})
  return(d)
}

Add_name_panel=function(d){
  d$Name_panel=paste0("Allochtonous flows: ",d$alpha_A," N:C, ",d$beta_A," P:C")
  
  d=d%>%dplyr::mutate(., Name_panel=recode_factor(Name_panel,
                                                  "Allochtonous flows: high N:C, high P:C"="Nutrient rich flows",
                                                  "Allochtonous flows: low N:C, high P:C"="Phosphorus rich flows",
                                                  "Allochtonous flows: high N:C, low P:C"="Nitrogen rich flows",
                                                  "Allochtonous flows: low N:C, low P:C"="Carbon rich flows"))
  return(d)
}

Get_data_thresholds=function(d,color_label=c("A","A","B","A")){
  
  tibble(beta_A=c("low","high","low","high"),alpha_A=c("high","low","low","high"),color_label=color_label,
         threshold=c(min(d$ID[which(d$beta_A=="low" & d$alpha_A=="high" & d$Limitation_Decompo %in% c("N","P"))]),
                     min(d$ID[which(d$beta_A=="high" & d$alpha_A=="low" & d$Limitation_Decompo %in% c("N","P"))]),
                     min(d$ID[which(d$beta_A=="low" & d$alpha_A=="low" & d$Limitation_Decompo %in% c("N","P"))]),
                     min(d$ID[which(d$beta_A=="high" & d$alpha_A=="high" & d$Limitation_Decompo %in% c("N","P"))])))%>%
    Add_name_panel(.)
  
}

scaling_vector_x_0_1=function(x){
  return((x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))
}

Add_limitation_label=function(p){
  return(
    
    p+  geom_label(data=tibble(Name_panel=c("Carbon rich flows",
                                            "Phosphorus rich flows",
                                            "Nitrogen rich flows"),
                               x=c(45,45,53),y=c(1.15,1.15,1.15),label=c("N-limited","N-limited","P-limited")),
                   aes(x=x,y=y,label=label))+
      geom_label(data=tibble(Name_panel=c("Carbon rich flows",
                                          "Phosphorus rich flows",
                                          "Nitrogen rich flows",
                                          "Nutrient rich flows"),
                             x=c(10,10,20,32),y=c(1.15,1.15,1.15,1.15),label=rep("C-limited",4)),
                 aes(x=x,y=y,label=label))+ylim(0,1.2)
  )
}


# ------------------ Final figures ----

## DATA ----

d=read.table("./data/Empirical/Stoichio_flows.csv",sep=";")

p=ggplot(d%>%
           mutate(.,Chemical=recode_factor(Chemical,"NP"="N:P of exported subsidies",
                                           "CP"="C:P of exported subsidies",
                                           "CN"="C:N of exported subsidies")),
         aes(x=Exporter_ecosyst,y=Ratio,shape=Mattyp1))+
  geom_boxplot(aes(fill=Exporter_ecosyst,group=Exporter_ecosyst),width=0.4,outlier.shape = NA,color="white",
               position = position_dodge(width = 1),alpha=.4)+ 
  geom_point(aes(color=Exporter_ecosyst),position=position_jitterdodge(dodge.width=.5,jitter.width=.03),size=2)+
  geom_hline(data=tibble(x=1:3,
                         Chemical=rep(c("N:P of exported subsidies","C:P of exported subsidies","C:N of exported subsidies"),each=1),y=c(116/16,116,16),
                         Exporter_ecosyst=rep(c(""),3),
                         Mattyp1=""),
             aes(x=x,y=y,yintercept = y,group=Exporter_ecosyst))+
  geom_text(data=tibble(x=rep(1:2,3),
                        Chemical=rep(c("N:P of exported subsidies","C:P of exported subsidies","C:N of exported subsidies"),each=2),y=1,
                        Exporter_ecosyst=rep(c("Forest","Grassland"),3),
                        label=paste0("n = ",c(46,6,65,9,104,19)),
                        Mattyp1=""),
            aes(x=x,y=y,label=label,group=Exporter_ecosyst),size=3.5)+
  labs(x="",fill="",y="",shape="",color="",shape="")+scale_y_log10()+
  scale_fill_manual(values=c("#C4DC93","#3D7543"))+
  scale_color_manual(values=c("#C4DC93","#3D7543"))+
  facet_wrap(.~Chemical,scales = "free",switch = "y",nrow = 1)+
  the_theme2+theme(axis.title.x = element_blank())+
  guides(color=F,fill=F)+
  theme(strip.placement = "outside")

ggsave("./Figures/Data_empirical.pdf",p,width = 9,height = 4)

## Along gradient of ID: densities, feedbacks, niche, CNP seston ----

d=read.table("./data/Simulations/Simulation_allochtonous.csv",sep=";")%>%
  mutate(., 
         alpha_A=ifelse(alpha_A==.02,"low","high"),
         beta_A=ifelse(beta_A==.002,"low","high"),
         Simulation_ID=as.factor(Simulation_ID)
  )%>%
  dplyr::filter(., Simulation_ID==3)%>%
  Add_name_panel(.)%>%
  dplyr::group_by(., Name_panel)%>%
  dplyr::mutate(., 
                Decomposers_C_scaled=scaling_vector_x_0_1(Decomposers_C),
                Non_fixers_C_scaled=scaling_vector_x_0_1(Non_fixers_C),
                PC_detritus_scaled=scaling_vector_x_0_1(PC_detritus),
                NC_detritus_scaled=scaling_vector_x_0_1(NC_detritus),
                Fixers_C_scaled=scaling_vector_x_0_1(Fixers_C),
                Frac_decomp_scaled=scaling_vector_x_0_1(Fixers_C),
                NP_threshold_NF_scaled=scaling_vector_x_0_1(NP_threshold_NF),
                PC_diff_scaled=scaling_vector_x_0_1(beta_B-PC_detritus),
                NC_diff_scaled=scaling_vector_x_0_1(alpha_B-NC_detritus),
  )

color_graph=c("C"="#D2B96F","N"="#8EBAEF")

table("alpha"=d$alpha_A,"beta"=d$beta_A,d$Limitation_Decompo)

p1=ggplot(NULL)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = threshold,xmax=65,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = 0,xmax=threshold,ymin=-Inf,ymax = Inf),lwd=1,alpha=.3,fill="#FFF4D5")+
  geom_line(data=d%>%melt(.,measure.vars=c("Decomposers_C_scaled","Non_fixers_C_scaled","Fixers_C_scaled"),
                          value.name = "Species_C",variable.name = "Species_name"),
            aes(x=ID,y=Species_C,color=Species_name,group=interaction(Species_name,beta_A)),lwd=1)+
  geom_point(data=d%>%melt(.,measure.vars=c("Decomposers_C_scaled","Non_fixers_C_scaled","Fixers_C_scaled"),
                           value.name = "Species_C",variable.name = "Species_name")%>%dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)]),
             aes(x=ID,y=Species_C,color=Species_name),size=3,shape=21)+
  scale_shape_manual(values = c(23,21,22))+
  facet_wrap(.~Name_panel,scales = "free")+#,
  the_theme2+
  labs(x="Allochtonous detritus inputs (ID)",y="Density of organisms (scaled)",color="")+
  guides(fill="none")+
  scale_fill_manual(values=c("A"="#A7C3D6","B"="#AA8DB5","C"="transparent"))+
  scale_color_manual(values=c("#FFA963","#80BFB5","#5F9203"),labels=c("Decomposers","Non-fixers","Fixers"))


ggsave("./Figures/Along_ID_gradient_density.pdf",Add_limitation_label(p1),width = 8,height = 7)


p1=ggplot(NULL)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = threshold,xmax=65,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = 0,xmax=threshold,ymin=-Inf,ymax = Inf),lwd=1,alpha=.3,fill="#FFF4D5")+
  geom_line(data=d%>%melt(.,measure.vars=c("Non_fixers_C_scaled","Fixers_C_scaled","NC_diff_scaled","PC_diff_scaled","NP_threshold_NF_scaled"),
                          value.name = "Species_C",variable.name = "Species_name"),
            aes(x=ID,y=Species_C,color=Species_name,group=interaction(Species_name)),lwd=1)+
  geom_point(data=d%>%melt(.,measure.vars=c("Non_fixers_C_scaled","Fixers_C_scaled","NC_diff_scaled","PC_diff_scaled","NP_threshold_NF_scaled"),
                           value.name = "Species_C",variable.name = "Species_name")%>%dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)]),
             aes(x=ID,y=Species_C,color=Species_name,shape=alpha_A),size=3)+
  scale_shape_manual(values = c(23,23))+
  facet_wrap(.~Name_panel,scales = "free")+#,
  the_theme2+
  labs(x="Allochtonous detritus inputs (ID)",y="Density of organisms (scaled)",color="")+
  guides(fill="none",shape="none")+
  scale_fill_manual(values=c("A"="#A7C3D6","B"="#AA8DB5","C"="transparent"))+
  scale_color_manual(values=c("#80BFB5","#5F9203","#FFBA80","#863886","#84412A"),
                     labels=c("Non-fixers","Fixers",
                              expression(paste(alpha[B]," - ",alpha[D])),
                              expression(paste(beta[B]," - ",beta[D])),
                              "Nitrogen/Phosphorus"))

ggsave("./Figures/SI/Along_ID_gradient_mecanisms.pdf",Add_limitation_label(p1),width = 8,height = 7)


# Indirect effects

color_graph=c("C"="#D2B96F","N"="#8EBAEF")

p1=ggplot(NULL)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = threshold,xmax=65,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = 0,xmax=threshold,ymin=-Inf,ymax = Inf),lwd=1,alpha=.3,fill="#FFF4D5")+
  geom_line(data=d%>%melt(.,measure.vars=c(colnames(d)[grep("indirect",colnames(d))])[c(1:2,5:6)],
                          value.name = "Species_C",variable.name = "Species_name"),
            aes(x=ID,y=Species_C,color=Species_name,group=interaction(Species_name,beta_A)),lwd=1)+
  geom_point(data=d%>%melt(.,measure.vars=c(colnames(d)[grep("indirect",colnames(d))])[c(1:2,5:6)],
                           value.name = "Species_C",variable.name = "Species_name")%>%dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)]),
             aes(x=ID,y=Species_C,color=Species_name),size=3,shape=21)+
  scale_shape_manual(values = c(23,23))+
  facet_wrap(.~Name_panel,scales = "free")+#,
  the_theme2+
  geom_hline(yintercept = 0)+
  labs(x="Allochtonous detritus inputs (ID)",y="Indirect effects between \n functional groups",color="")+
  guides(fill="none",shape="none")+
  scale_fill_manual(values=c("A"="#A7C3D6","B"="#AA8DB5","C"="transparent"))+
  scale_color_manual(values=c("#FFBA80","#80BFB5",
                              "#A95F00","#216156"),
                     label=c("Decomposers on fixers",
                             "Fixers on decomposers",
                             "Decomposers on non-fixers",
                             "Non-fixers on decomposers"))+
  guides(color=guide_legend(ncol=2))+
  geom_label(data=tibble(Name_panel=c("Carbon rich flows",
                                      "Phosphorus rich flows",
                                      "Nitrogen rich flows"),
                         x=c(45,45,53),y=c(.2,.2,.2),label=c("N-limited","N-limited","P-limited")),
             aes(x=x,y=y,label=label))+
  geom_label(data=tibble(Name_panel=c("Carbon rich flows",
                                      "Phosphorus rich flows",
                                      "Nitrogen rich flows",
                                      "Nutrient rich flows"),
                         x=c(10,10,20,32),y=c(.2,.2,.2,.2),label=rep("C-limited",4)),
             aes(x=x,y=y,label=label))



ggsave("./Figures/Along_ID_gradient_indirect.pdf",p1,width = 8,height = 7)

p1=ggplot(NULL)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = threshold,xmax=65,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = 0,xmax=threshold,ymin=-Inf,ymax = Inf),lwd=1,alpha=.3,fill="#FFF4D5")+
  geom_line(data=d%>%melt(.,measure.vars=c("NC_detritus_scaled","PC_detritus_scaled"),
                          value.name = "Species_C",variable.name = "Species_name"),
            aes(x=ID,y=Species_C,color=Species_name,group=interaction(Species_name,beta_A)),lwd=1)+
  geom_point(data=d%>%melt(.,measure.vars=c("NC_detritus_scaled","PC_detritus_scaled"),
                           value.name = "Species_C",variable.name = "Species_name")%>%dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)]),
             aes(x=ID,y=Species_C,color=Species_name),size=3,shape=21)+
  scale_shape_manual(values = c(23,21,22))+
  facet_wrap(.~Name_panel,scales = "free")+#,
  #labeller = label_bquote(cols=N:C~allochtonous~low~(alpha[A])~P:C~allochtonous~.(as.character(beta_A))~(beta[A])))+
  the_theme2+
  labs(x="Allochtonous detritus inputs (ID)",y="Elemental ratios (scaled)",color="")+
  guides(fill="none")+
  scale_fill_manual(values=c("A"="#A7C3D6","B"="#AA8DB5","C"="transparent"))+
  scale_color_manual(values=c("black","blue"),labels=c("N:C detritus","P:C detritus"))

ggsave("./Figures/SI/Along_ID_gradient_NC_PC.pdf",Add_limitation_label(p1),width = 7,height = 7)

#Computing changes in species niche' as a proxy of indirect effects

d=read.table("./data/Simulations/Simulation_allochtonous.csv",sep=";")%>%
  mutate(., 
         alpha_A=ifelse(alpha_A==.02,"low","high"),
         beta_A=ifelse(beta_A==.002,"low","high"),
         Simulation_ID=as.factor(Simulation_ID)
  )%>%
  dplyr::filter(., Simulation_ID==3)

table(d$Limitation_Decompo)
d2=read.table("./data/Simulations/Simulation_allochtonous_species_alone.csv",sep=";")%>%
  mutate(., 
         alpha_A=ifelse(alpha_A==.02,"low","high"),
         beta_A=ifelse(beta_A==.002,"low","high"),
         Simulation_ID=as.factor(Simulation_ID)
  )%>%
  dplyr::filter(., Simulation_ID==3)

d$Name_panel=paste0("Allochtonous flows: ",d$alpha_A," N:C, ",d$beta_A," P:C")

for (k in unique(d2$Which_species_alone)[1:3]){
  d2_fil=dplyr::filter(d2,Which_species_alone==k)
  d$NIntA=2*(d[,paste0(k,"_C")]-d2_fil[,paste0(k,"_C")])/(d2_fil[,paste0(k,"_C")]+abs((d[,paste0(k,"_C")]-d2_fil[,paste0(k,"_C")])))
  colnames(d)[ncol(d)]= paste0("NIntA_",k)
}

d=d%>%Add_name_panel(.)%>%
  dplyr::group_by(., Name_panel)

p1=ggplot(NULL)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = threshold,xmax=65,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = 0,xmax=threshold,ymin=-Inf,ymax = Inf),lwd=1,alpha=.3,fill="#FFF4D5")+
  geom_line(data=d%>%melt(.,measure.vars=colnames(d)[grep("NInt",colnames(d))],
                          value.name = "Species_C",variable.name = "Species_name"),
            aes(x=ID,y=Species_C,color=Species_name),lwd=1)+
  geom_point(data=d%>%melt(.,measure.vars=colnames(d)[grep("NInt",colnames(d))],
                           value.name = "Species_C",variable.name = "Species_name")%>%
               dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)]),
             aes(x=ID,y=Species_C,color=Species_name),size=3,shape=21)+
  scale_pattern_manual(values = c("C" = "none", "N" = "stripe"))+
  facet_wrap(.~Name_panel,scales = "free")+#,
  #labeller = label_bquote(cols=N:C~allochtonous~low~(alpha[A])~P:C~allochtonous~.(as.character(beta_A))~(beta[A])))+
  the_theme2+
  geom_hline(yintercept = 0)+
  labs(x="Allochtonous detritus inputs (ID)",y="Change in organism niche",color="")+
  guides(fill="none",shape="none")+
  scale_fill_manual(values=c("A"="#A7C3D6","B"="#AA8DB5","C"="transparent"))+
  scale_color_manual(values=c("#FFA963","#80BFB5","#5F9203"),labels=c("Decomposers","Non-fixers","Fixers"))

ggsave("./Figures/SI/Along_ID_gradient_niche.pdf",Add_limitation_label(p1)+ylim(-.5,2),width = 8,height = 7)



d=read.table("./data/Simulations/Simulation_allochtonous.csv",sep=";")%>%
  dplyr::mutate(., 
                beta_diff=.$alpha_B-.$PC_detritus,
                alpha_diff=.$beta_B-.$NC_detritus)%>%
  mutate(., 
         alpha_A=ifelse(alpha_A==.02,"low","high"),
         beta_A=ifelse(beta_A==.002,"low","high"),
         Simulation_ID=as.factor(Simulation_ID)
  )%>%
  Add_name_panel(.)%>%
  dplyr::group_by(., Name_panel)%>%
  dplyr::filter(.,Simulation_ID==3)%>%
  dplyr::mutate(., 
                Decomposers_C_scaled=scaling_vector_x_0_1(Decomposers_C),
                Non_fixers_C_scaled=scaling_vector_x_0_1(Non_fixers_C),
                beta_diff_scaled=scaling_vector_x_0_1(beta_diff),
                alpha_diff_scaled=scaling_vector_x_0_1(alpha_diff),
                Fixers_C_scaled=scaling_vector_x_0_1(Fixers_C),
                PC_diff_scaled=scaling_vector_x_0_1(beta_B-PC_detritus),
                NC_diff_scaled=scaling_vector_x_0_1(alpha_B-NC_detritus),
                NC_seston_scaled=scaling_vector_x_0_1(1/CN_seston),
                PC_seston_scaled=scaling_vector_x_0_1(1/CP_seston),
                NP_seston_scaled=scaling_vector_x_0_1(NP_seston),
                PC_seston=1/CP_seston,
                NC_seston=1/CN_seston)

color_graph=c("C"="#D2B96F","N"="#8EBAEF")

Rename_seston=function(d){
  return(d%>%mutate(., Species_name=recode_factor(Species_name,
                                                  "CN_seston"="Seston C:N",
                                                  "NC_seston"="Seston N:C",
                                                  "PC_seston"="Seston P:C",
                                                  "NP_seston"="Seston N:P",
                                                  "CP_seston"="Seston C:P",
                                                  "Deviation_Redfield"="Deviation from Redfield ratio")))
}

for (k in 1:4){
  
  d_plot=dplyr::filter(d,alpha_A==c("high","low","low","high")[k],
                       beta_A==c("low","high","low","high")[k])
  
  assign(paste0("p",k),
         ggplot(NULL)+
           geom_rect(data=Get_data_thresholds(d,c("B","A","A","C"))[k,],
                     aes(xmin = threshold,xmax=65,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
           geom_rect(data=Get_data_thresholds(d,c("B","A","A","C"))[k,],
                     aes(xmin = 0,xmax=threshold,ymin=-Inf,ymax = Inf),lwd=1,alpha=.3,fill="#FFF4D5")+
           geom_line(data=d_plot%>%melt(.,measure.vars=c("PC_seston","NC_seston"),
                                   value.name = "Species_C",variable.name = "Species_name")%>%Rename_seston(.),
                     aes(x=ID,y=Species_C,color=interaction(beta_A,alpha_A),group=interaction(beta_A,alpha_A)),lwd=1)+
           geom_point(data=d_plot%>%melt(.,measure.vars=c("NC_seston",
                                                     "PC_seston"),
                                    value.name = "Species_C",variable.name = "Species_name")%>%
                        dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)])%>%Rename_seston(.),
                      aes(x=ID,y=Species_C,color=interaction(beta_A,alpha_A),shape=beta_A),size=3)+
           scale_shape_manual(values = c(23,23))+
           facet_wrap(Species_name~.,scales = "free",switch = "y")+
           geom_hline(data=tibble(Species_name=c("Seston N:C","Seston P:C"),
                                  Ratio=c(1/16,1/106)),aes(yintercept = Ratio))+
           the_theme2+
           labs(x="Allochtonous detritus inputs (ID)",y="",color="")+
           guides(fill="none",shape="none",color="none")+
           scale_fill_manual(values=c("A"="#A7C3D6","B"="#AA8DB5","C"="transparent"))+
           theme(strip.text.x = element_blank(),strip.placement = "outside")+
           scale_color_manual(values=c("black","blue"))
         )
}

p_tot=ggarrange(p4+ggtitle(" Nitrogen rich flows"),
                p1+ggtitle(" Phosphorus rich flows"),
                p2+ggtitle(" Carbon rich flows"),
                p3+ggtitle(" Nutrient rich flows"),ncol=2,nrow=2)

ggsave("./Figures/Along_ID_gradient_CNP.pdf",p_tot,width = 12,height = 7)


d2=dplyr::filter(d)%>%
  melt(., measure.vars=c("Fixers_C","Decomposers_C","Non_fixers_C","Detritus_C"))%>%
  dplyr::group_by(.,ID,Name_panel)%>%
  mutate(., value=100*value/sum(value))

p1=ggplot(d2)+
  geom_area(aes(x=ID,y=value,fill=variable),alpha=.5)+
  # geom_vline(data=Get_data_thresholds(d),aes(xintercept = threshold))+
  geom_segment(data=Get_data_thresholds(d),
               aes(x = threshold, y = 20, xend = threshold, yend = 0),
               arrow = arrow(length = unit(0.2, "cm")),lwd=1)+
  the_theme2+
  facet_wrap(.~factor(Name_panel,c("Carbon rich flows","Nitrogen rich flows",
                                   "Nutrient rich flows","Phosphorus rich flows")),scales = "free",nrow = 2,ncol = 2)+
  labs(x="Allochtonous detritus inputs (ID)",y="Origin of the carbon in the seston (%)",color="",fill="")+
  guides(shape="none")+
  scale_fill_manual(values=c("#80BFB5","#FFA963","#5F9203","brown"),
                    labels=c("Fixers","Decomposers","Non-fixers","Detritus"))

ggsave("./Figures/Along_ID_gradient_stacked.pdf",p1,width = 8,height = 7)



## Homeostasis ----
Compute_slope=function(CP,beta){
  return(unlist(as.numeric(coef(lm(1/CP~beta))[2])))
}

d=rbind(read.table("./data/Simulations/Simulation_allochtonous_alpha_beta_A.csv",sep=";")%>%
          dplyr::mutate(., 
                        beta_diff=.015-.$PC_detritus,
                        alpha_diff=.15-.$NC_detritus)%>%
          Add_name_panel(.)%>%
          dplyr::filter(., ID %in% unique(.$ID)[seq(1,40,by=2)],
                        Simulation_ID==1)%>%
          dplyr::group_by(., ID,Limitation_Decompo)%>%
          dplyr::do(., slope=Compute_slope(.$CP_seston,.$beta_A))%>%
          add_column(., Ratio="P:C ratio"),
        read.table("./data/Simulations/Simulation_allochtonous_alpha_beta_A.csv",sep=";")%>%
          dplyr::mutate(., 
                        beta_diff=.015-.$PC_detritus,
                        alpha_diff=.15-.$NC_detritus)%>%
          Add_name_panel(.)%>%
          dplyr::filter(.,ID %in% unique(.$ID)[seq(1,40,by=2)],
                        Simulation_ID==2)%>%
          dplyr::group_by(., ID,Limitation_Decompo)%>%
          dplyr::do(., slope=Compute_slope(.$CN_seston,.$alpha_A))%>%
          add_column(., Ratio="N:C ratio")
  
)


d$slope=unlist(d$slope)
d$ID=as.numeric(d$ID)
p=ggplot(d)+
  geom_point(aes(x=ID,y=1/slope,shape=Limitation_Decompo,fill=as.numeric(ID)),size=2,color="black")+
  scale_fill_stepsn(colours = viridis(5,option = "G"))+
  scale_shape_manual(values=c(21,22,23))+
  the_theme2+
  facet_wrap(.~Ratio,scales="free",nrow=2)+
  labs(x="Allochtonous detritus \n inputs (ID)",y="Homeostasis (1/slope)",shape="Decomposer' limitation   ",
       fill="Allochtonous   \n detritus inputs   \n   (ID)  ")+
  geom_hline(yintercept = 1)+
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=13),
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        strip.text.x = element_text(size=13))

p_legend=get_legend(p+theme(legend.box = "vertical"))


d=read.table("./data/Simulations/Simulation_allochtonous_alpha_beta_A.csv",sep=";")%>%
  dplyr::mutate(., 
                beta_diff=.015-.$PC_detritus,
                alpha_diff=.15-.$NC_detritus)%>%
  Add_name_panel(.)

p11=ggplot(d%>%
             mutate(., 
                    alpha_A=recode_factor(alpha_A,"0.02"="low","0.14"="high"))%>%
             dplyr::filter(.,
                           ID %in% unique(.$ID)[seq(1,40,by=8)],
                           alpha_A=="low",
                           Simulation_ID==1))+
  geom_line(aes(x=beta_A,y=1/CP_seston,color=ID,group=interaction(ID)),lwd=1)+
  geom_point(data=d%>%
               mutate(., 
                      alpha_A=recode_factor(alpha_A,"0.02"="low","0.14"="high"))%>%
               dplyr::filter(.,
                             Simulation_ID==1,
                             alpha_A=="low",
                             ID %in% unique(.$ID)[seq(1,40,by=8)],
                             beta_A %in% unique(.$beta_A)[seq(1,40,by=5)]),
             aes(x=beta_A,y=1/CP_seston,fill=ID,shape=Limitation_Decompo),size=3,color="black")+
  scale_fill_stepsn(colours = viridis(5,option = "G"))+
  scale_color_stepsn(colours = viridis(5,option = "G"))+
  geom_hline(yintercept = 1/106,color="pink",linetype=2,lwd=1)+
  scale_shape_manual(values = c(21,22,23))+
  labs(x=expression(paste("P:C of allochtonous flows (",beta[A],")")),
       y=expression(paste("Seston P:C ratio (",beta[seston],")")),
       color="Detritus inflow (ID)",shape="Decomposer' limitation")+
  the_theme2+
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=13),
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        strip.text.x = element_text(size=13))


p21=ggplot(d%>%
             mutate(., 
                    beta_A=recode_factor(beta_A,"0.002"="low","0.012"="high"))%>%
             dplyr::filter(.,
                           ID %in% unique(.$ID)[seq(1,40,by=8)],
                           beta_A=="low",
                           Simulation_ID==2))+
  geom_line(aes(x=alpha_A,y=1/CN_seston,color=ID,group=interaction(ID,IN,IP)),lwd=1)+
  geom_point(data=d%>%
               mutate(., 
                      beta_A=recode_factor(beta_A,"0.002"="low","0.012"="high"))%>%
               dplyr::filter(.,
                             Simulation_ID==2,
                             beta_A=="low",
                             ID %in% unique(.$ID)[seq(1,40,by=8)],
                             alpha_A %in% unique(.$alpha_A)[seq(1,40,by=5)]),
             aes(x=alpha_A,y=1/CN_seston,fill=ID,shape=Limitation_Decompo),size=3,color="black")+
  geom_hline(yintercept = 16/106,color="pink",linetype=2,lwd=1)+
  scale_fill_stepsn(colours = viridis(5,option = "G"))+
  scale_color_stepsn(colours = viridis(5,option = "G"))+
  scale_shape_manual(values = c(21,22,23))+
  labs(x=expression(paste("N:C of allochtonous flows (",alpha[A],")")),
       y=expression(paste("Seston N:C ratio (",alpha[seston],")")),
       color="Detritus inflow (ID)",shape="Decomposer' limitation")+
  the_theme2+
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=13),
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        strip.text.x = element_text(size=13))



ggsave("./Figures/Variation_NC_seston_input_all_species.pdf",
       ggarrange(
         ggarrange(
           ggarrange(
             ggplot()+theme_void(),p+theme(legend.position = "none"),widths = c(1,.45),labels = c("a","b")),
                   ggarrange(p11,p21,ncol=2,labels = c("c","d"),common.legend = T,legend="none",heights = c(1,1.05)),
                   nrow=2,heights = c(1.3,1),labels = c("","")),p_legend,heights = c(1,.15),nrow=2),
       width = 8,height = 9)




## Other indirect effects ----


d=read.table("./data/Simulations/Simulation_allochtonous_I_to_A.csv",sep=";")%>%
  mutate(., 
         alpha_A=ifelse(alpha_A==.02,"low","high"),
         beta_A=ifelse(beta_A==.002,"low","high"),
         Simulation_ID=as.factor(Simulation_ID)
  )%>%
  dplyr::filter(., Simulation_ID==3)%>%
  Add_name_panel(.)%>%
  dplyr::group_by(., Name_panel)%>%
  dplyr::mutate(., 
                Decomposers_C_scaled=scaling_vector_x_0_1(Decomposers_C),
                Non_fixers_C_scaled=scaling_vector_x_0_1(Non_fixers_C),
                PC_detritus_scaled=scaling_vector_x_0_1(PC_detritus),
                NC_detritus_scaled=scaling_vector_x_0_1(NC_detritus),
                Fixers_C_scaled=scaling_vector_x_0_1(Fixers_C),
                Frac_decomp_scaled=scaling_vector_x_0_1(Fixers_C),
                NP_threshold_NF_scaled=scaling_vector_x_0_1(NP_threshold_NF),
                PC_diff_scaled=scaling_vector_x_0_1(beta_B-PC_detritus),
                NC_diff_scaled=scaling_vector_x_0_1(alpha_B-NC_detritus),
  )


color_graph=c("C"="#D2B96F","N"="#8EBAEF")

p1=ggplot(NULL)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = threshold,xmax=65,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = 0,xmax=threshold,ymin=-Inf,ymax = Inf),lwd=1,alpha=.3,fill="#FFF4D5")+
  geom_line(data=d%>%melt(.,measure.vars=c(colnames(d)[grep("indirect",colnames(d))])[c(1:2,5:6)],
                          value.name = "Species_C",variable.name = "Species_name"),
            aes(x=ID,y=Species_C,color=Species_name,group=interaction(Species_name,beta_A)),lwd=1)+
  geom_point(data=d%>%melt(.,measure.vars=c(colnames(d)[grep("indirect",colnames(d))])[c(1:2,5:6)],
                           value.name = "Species_C",variable.name = "Species_name")%>%dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)]),
             aes(x=ID,y=Species_C,color=Species_name),size=3,shape=21)+
  scale_shape_manual(values = c(23,23))+
  facet_wrap(.~Name_panel,scales = "free")+#,
  the_theme2+
  geom_hline(yintercept = 0)+
  labs(x="Allochtonous detritus inputs (ID)",y="Indirect effects between \n functional groups",color="")+
  guides(fill="none",shape="none")+
  scale_fill_manual(values=c("A"="#A7C3D6","B"="#AA8DB5","C"="transparent"))+
  scale_color_manual(values=c("#FFBA80","#80BFB5",
                              "#A95F00","#216156"),
                     label=c("Decomposers on fixers",
                             "Fixers on decomposers",
                             "Decomposers on non-fixers",
                             "Non-fixers on decomposers"))+
  guides(color=guide_legend(ncol=2))+
  geom_label(data=tibble(Name_panel=c("Carbon rich flows",
                                      "Phosphorus rich flows",
                                      "Nitrogen rich flows"),
                         x=c(45,45,53),y=c(.2,.2,.2),label=c("N-limited","N-limited","P-limited")),
             aes(x=x,y=y,label=label))+
  geom_label(data=tibble(Name_panel=c("Carbon rich flows",
                                      "Phosphorus rich flows",
                                      "Nitrogen rich flows",
                                      "Nutrient rich flows"),
                         x=c(10,10,20,32),y=c(.2,.2,.2,.2),label=rep("C-limited",4)),
             aes(x=x,y=y,label=label))+ylim(-.28,.22)+xlim(0,65)

ggsave("./Figures/SI/Along_ID_gradient_indirect_other.pdf",p1,width = 8,height = 7)


# ------------------ 2 species: decomposers-fixers ----

d=read.table("./data/Simulations/Simulation_fixers_and_decomposers_2D_quality_quantity.csv",sep=";")%>%
  dplyr::filter(., Simulation_ID=="Beta_vary",beta_A%in%c(.002,.012))%>%
  mutate(., 
         alpha_A=ifelse(alpha_A==.02,"low","high"),
         beta_A=ifelse(beta_A==.002,"low","high"),
         Simulation_ID=as.factor(Simulation_ID)
  )%>%
  Add_name_panel(.)%>%
  dplyr::group_by(., Name_panel)%>%
  dplyr::mutate(., 
                Decomposers_C_scaled=scale(Decomposers_C)[,1],
                Non_fixers_C_scaled=scale(Non_fixers_C)[,1],
                PC_detritus_scaled=scale(PC_detritus)[,1],
                NC_detritus_scaled=scale(NC_detritus)[,1],
                Fixers_C_scaled=scale(Fixers_C)[,1],
                Frac_decomp_scaled=scale(Fixers_C)[,1],
                NP_threshold_NF_scaled=scale(NP_threshold_NF)[,1]
  )


p1=ggplot(NULL)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = threshold,xmax=50,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
  geom_line(data=d%>%melt(.,measure.vars=c("Decomposers_C_scaled","Fixers_C_scaled"),
                          value.name = "Species_C",variable.name = "Species_name"),
            aes(x=ID,y=Species_C,color=Species_name,group=interaction(Species_name,beta_A)),lwd=1)+
  geom_point(data=d%>%melt(.,measure.vars=c("Decomposers_C_scaled","Fixers_C_scaled"),
                           value.name = "Species_C",variable.name = "Species_name")%>%dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)]),
             aes(x=ID,y=Species_C,color=Species_name),size=3,shape=21)+
  scale_shape_manual(values = c(23,21,22))+
  facet_wrap(.~Name_panel,scales = "free")+#,
  #labeller = label_bquote(cols=N:C~allochtonous~low~(alpha[A])~P:C~allochtonous~.(as.character(beta_A))~(beta[A])))+
  the_theme2+
  labs(x="Allochtonous detritus inputs (ID)",y="Density of organisms (scaled)",color="")+
  guides(fill="none")+
  scale_fill_manual(values=c("A"="#9BBDFF","B"="violet","C"="transparent"))+
  scale_color_manual(values=c("#FFA963","#5F9203"),labels=c("Decomposers","Fixers"))



ggsave("./Figures/SI/Decomposers_fixers_PC_detritus.pdf",p1,width =  8,height = 7)


d2=dplyr::filter(d)%>%
  melt(., measure.vars=c("Fixers_C","Decomposers_C","Detritus_C"))%>%
  dplyr::group_by(.,ID,Name_panel)%>%
  mutate(., value=100*value/sum(value))

p1=ggplot(d2)+
  geom_area(aes(x=ID,y=value,fill=variable),alpha=.5)+
  # geom_vline(data=Get_data_thresholds(d),aes(xintercept = threshold))+
  geom_segment(data=Get_data_thresholds(d),
               aes(x = threshold, y = 20, xend = threshold, yend = 0),
               arrow = arrow(length = unit(0.2, "cm")),lwd=1)+
  the_theme2+
  facet_wrap(.~Name_panel,scales = "free",nrow = 2,ncol = 2)+
  labs(x="Allochtonous detritus inputs (ID)",y="Origin of the carbon in the seston (%)",color="",fill="")+
  guides(shape="none")+
  scale_fill_manual(values=c("#80BFB5","#FFA963","brown"),
                    labels=c("Fixers","Decomposers","Detritus"))

ggsave("./Figures/SI/Along_ID_gradient_stacked_fixers_and_decomposers.pdf",p1,width = 8,height = 7)






p1=ggplot(NULL)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = threshold,xmax=50,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
  geom_line(data=d%>%melt(.,measure.vars=c(colnames(d)[grep("indirect",colnames(d))])[c(1:2)],
                          value.name = "Species_C",variable.name = "Species_name"),
            aes(x=ID,y=Species_C,color=Species_name,group=interaction(Species_name,beta_A)),lwd=1)+
  geom_point(data=d%>%melt(.,measure.vars=c(colnames(d)[grep("indirect",colnames(d))])[c(1:2)],
                           value.name = "Species_C",variable.name = "Species_name")%>%dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)]),
             aes(x=ID,y=Species_C,color=Species_name),size=3,shape=21)+
  scale_shape_manual(values = c(23,23))+
  facet_wrap(.~Name_panel,scales = "free")+#,
  #labeller = label_bquote(cols=N:C~allochtonous~low~(alpha[A])~P:C~allochtonous~.(as.character(beta_A))~(beta[A])))+
  the_theme2+
  geom_hline(yintercept = 0)+
  labs(x="Allochtonous detritus inputs (ID)",y="Indirect effects between \n functional groups",color="")+
  guides(fill="none",shape="none")+
  scale_fill_manual(values=c("A"="#9BBDFF","B"="violet","C"="transparent"))+
  scale_color_manual(values=c("#FFBA80","#80BFB5"),#"#5F9203",
                              #"#000467",
                              #"#A95F00","#216156"),
                     label=c("Decomposers on fixers",
                             "Fixers on decomposers"#,
                             #"Non-fixers on fixers",
                             # "Fixers on non-fixers",
                             #"Decomposers on non-fixers",
                             #"Non-fixers on decomposers"
                             ))+
  guides(color=guide_legend(ncol=2))


ggsave("./Figures/SI/Along_ID_gradient_indirect_fix_decompo.pdf",p1,width = 8,height = 7)


d=read.table("./data/Simulations/Simulation_fixers_and_decomposers_2D_quality_quantity.csv",sep=";")%>%
  dplyr::mutate(., 
                beta_diff=.015-.$PC_detritus,
                alpha_diff=.15-.$NC_detritus)

p11=ggplot(d%>%
             mutate(., 
                    alpha_A=recode_factor(alpha_A,"0.02"="low","0.14"="high"))%>%
             dplyr::filter(.,
                           Simulation_ID=="Beta_vary",
                           ID %in% unique(.$ID)[seq(1,40,by=8)]))+
  geom_line(aes(x=beta_A,y=1/CP_seston,color=ID,group=ID),lwd=1)+
  geom_point(data=d%>%
               mutate(., 
                      alpha_A=recode_factor(alpha_A,"0.02"="low","0.14"="high"))%>%
               dplyr::filter(.,
                             Simulation_ID=="Beta_vary",
                             ID %in% unique(.$ID)[seq(1,40,by=8)],
                             beta_A %in% unique(.$beta_A)[seq(1,40,by=5)]),
             aes(x=beta_A,y=1/CP_seston,color=ID,shape=Limitation_Decompo),size=3,fill="white")+
  scale_colour_stepsn(colours = viridis(5,option = "G"))+
  facet_wrap(.~alpha_A,scales = "free",
             labeller = label_bquote(cols=alpha[A]==.(as.character(alpha_A))))+
  scale_shape_manual(values = c(21,23))+
  labs(x=expression(paste("P:C of allochtonous flows (",beta[A],")")),
       y=expression(paste("Seston P:C ratio (",beta[seston],")")))+
  the_theme2

p21=ggplot(d%>%
             mutate(., 
                    beta_A=recode_factor(beta_A,"0.002"="low","0.012"="high"))%>%
             dplyr::filter(.,
                           Simulation_ID=="Alpha_vary",
                           ID %in% unique(.$ID)[seq(1,40,by=8)]))+
  geom_line(aes(x=alpha_A,y=1/CN_seston,color=ID,group=ID),lwd=1)+
  geom_point(data=d%>%
               mutate(., 
                      beta_A=recode_factor(beta_A,"0.002"="low","0.012"="high"))%>%
               dplyr::filter(.,
                             Simulation_ID=="Alpha_vary",
                             ID %in% unique(.$ID)[seq(1,40,by=8)],
                             alpha_A %in% unique(.$alpha_A)[seq(1,40,by=5)]),
             aes(x=alpha_A,y=1/CN_seston,color=ID,shape=Limitation_Decompo),size=3,fill="white")+
  scale_colour_stepsn(colours = viridis(5,option = "G"))+
  facet_wrap(.~beta_A,scales = "free",
             labeller = label_bquote(cols=beta[A]==.(as.character(beta_A))))+
  scale_shape_manual(values = c(21,23))+
  labs(x=expression(paste("N:C of allochtonous flows (",alpha[A],")")),
       y=expression(paste("Seston N:C ratio (",alpha[seston],")")))+
  the_theme2

ggsave("./Figures/SI/Variation_NC_seston_input_fixer_and_decomposers.pdf",
       ggarrange(p11,p21,nrow=2,labels = letters[1:2],common.legend = T,legend="bottom",heights = c(1,1.05)),
       width = 7,height = 8)



# ------------------ 2 species: decomposers-non.fixers ----

d=read.table("./data/Simulations/Simulation_Non_fixers_and_decomposers_2D_quality_quantity.csv",sep=";")%>%
  dplyr::filter(., Simulation_ID=="Beta_vary",beta_A%in%c(.002,.012))%>%
  mutate(., 
         alpha_A=ifelse(alpha_A==.02,"low","high"),
         beta_A=ifelse(beta_A==.002,"low","high"),
         Simulation_ID=as.factor(Simulation_ID)
  )%>%
  Add_name_panel(.)%>%
  dplyr::group_by(., Name_panel)%>%
  dplyr::mutate(., 
                Decomposers_C_scaled=scale(Decomposers_C)[,1],
                Non_fixers_C_scaled=scale(Non_fixers_C)[,1],
                PC_detritus_scaled=scale(PC_detritus)[,1],
                NC_detritus_scaled=scale(NC_detritus)[,1],
                Fixers_C_scaled=scale(Fixers_C)[,1],
                Frac_decomp_scaled=scale(Fixers_C)[,1],
                NP_threshold_NF_scaled=scale(NP_threshold_NF)[,1]
  )



p1=ggplot(NULL)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = threshold,xmax=50,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
  geom_line(data=d%>%melt(.,measure.vars=c("Decomposers_C_scaled","Non_fixers_C_scaled"),
                          value.name = "Species_C",variable.name = "Species_name"),
            aes(x=ID,y=Species_C,color=Species_name,group=interaction(Species_name,beta_A)),lwd=1)+
  geom_point(data=d%>%melt(.,measure.vars=c("Decomposers_C_scaled","Non_fixers_C_scaled"),
                           value.name = "Species_C",variable.name = "Species_name")%>%dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)]),
             aes(x=ID,y=Species_C,color=Species_name),size=3,shape=21)+
  scale_shape_manual(values = c(23,21,22))+
  facet_wrap(.~Name_panel,scales = "free")+#,
  #labeller = label_bquote(cols=N:C~allochtonous~low~(alpha[A])~P:C~allochtonous~.(as.character(beta_A))~(beta[A])))+
  the_theme2+
  labs(x="Allochtonous detritus inputs (ID)",y="Density of organisms (scaled)",color="")+
  guides(fill="none")+
  scale_fill_manual(values=c("A"="#9BBDFF","B"="violet","C"="transparent"))+
  scale_color_manual(values=c("#FFA963","#80BFB5"),labels=c("Decomposers","Non-fixers"))


ggsave("./Figures/SI/Decomposers_non_fixers_PC_NC_detritus.pdf",p1,width =  8,height = 7)

p1=ggplot(NULL)+
  geom_rect(data=Get_data_thresholds(d,c("B","A","A","C")),
            aes(xmin = threshold,xmax=50,ymin=-Inf,ymax = Inf,fill=color_label),lwd=1,alpha=.3)+
  geom_line(data=d%>%melt(.,measure.vars=c(colnames(d)[grep("indirect",colnames(d))])[c(5:6)],
                          value.name = "Species_C",variable.name = "Species_name"),
            aes(x=ID,y=Species_C,color=Species_name,group=interaction(Species_name,beta_A)),lwd=1)+
  geom_point(data=d%>%melt(.,measure.vars=c(colnames(d)[grep("indirect",colnames(d))])[c(5:6)],
                           value.name = "Species_C",variable.name = "Species_name")%>%dplyr::filter(., ID %in% unique(.$ID)[seq(1,100,by=10)]),
             aes(x=ID,y=Species_C,color=Species_name),size=3,shape=21)+
  scale_shape_manual(values = c(23,23))+
  facet_wrap(.~Name_panel,scales = "free")+#,
  #labeller = label_bquote(cols=N:C~allochtonous~low~(alpha[A])~P:C~allochtonous~.(as.character(beta_A))~(beta[A])))+
  the_theme2+
  geom_hline(yintercept = 0)+
  labs(x="Allochtonous detritus inputs (ID)",y="Indirect effects between \n functional groups",color="")+
  guides(fill="none",shape="none")+
  scale_fill_manual(values=c("A"="#9BBDFF","B"="violet","C"="transparent"))+
  scale_color_manual(values=c("#A95F00","#216156"),
                     label=c("Decomposers on non-fixers",
                             "Non-fixers on decomposers"
                     ))+
  guides(color=guide_legend(ncol=2))


ggsave("./Figures/SI/Along_ID_gradient_indirect_nonfix_decompo.pdf",p1,width = 8,height = 7)



d2=dplyr::filter(d)%>%
  melt(., measure.vars=c("Decomposers_C","Non_fixers_C","Detritus_C"))%>%
  dplyr::group_by(.,ID,Name_panel)%>%
  mutate(., value=100*value/sum(value))

p1=ggplot(d2)+
  geom_area(aes(x=ID,y=value,fill=variable),alpha=.5)+
  # geom_vline(data=Get_data_thresholds(d),aes(xintercept = threshold))+
  geom_segment(data=Get_data_thresholds(d),
               aes(x = threshold, y = 20, xend = threshold, yend = 0),
               arrow = arrow(length = unit(0.2, "cm")),lwd=1)+
  the_theme2+
  facet_wrap(.~Name_panel,scales = "free",nrow = 2,ncol = 2)+
  labs(x="Allochtonous detritus inputs (ID)",y="Origin of the carbon in the seston (%)",color="",fill="")+
  guides(shape="none")+
  scale_fill_manual(values=c("#FFA963","#5F9203","brown"),
                    labels=c("Decomposers","Non-fixers","Detritus"))

ggsave("./Figures/SI/Along_ID_gradient_stacked_non_fixers_and_decomposers.pdf",p1,width = 8,height = 7)



d=read.table("./data/Simulations/Simulation_Non_fixers_and_decomposers_2D_quality_quantity.csv",sep=";")%>%
  dplyr::mutate(., 
                beta_diff=.015-.$PC_detritus,
                alpha_diff=.15-.$NC_detritus)

p11=ggplot(d%>%
             mutate(., 
                    alpha_A=recode_factor(alpha_A,"0.02"="low","0.14"="high"))%>%
             dplyr::filter(.,
                           Simulation_ID=="Beta_vary",
                           ID %in% unique(.$ID)[seq(1,40,by=8)]))+
  geom_line(aes(x=beta_A,y=1/CP_seston,color=ID,group=ID),lwd=1)+
  geom_point(data=d%>%
               mutate(., 
                      alpha_A=recode_factor(alpha_A,"0.02"="low","0.14"="high"))%>%
               dplyr::filter(.,
                             Simulation_ID=="Beta_vary",
                             ID %in% unique(.$ID)[seq(1,40,by=8)],
                             beta_A %in% unique(.$beta_A)[seq(1,40,by=5)]),
             aes(x=beta_A,y=1/CP_seston,color=ID,shape=Limitation_Decompo),size=3,fill="white")+
  scale_colour_stepsn(colours = viridis(5,option = "G"))+
  facet_wrap(.~alpha_A,scales = "free",
             labeller = label_bquote(cols=alpha[A]==.(as.character(alpha_A))))+
  scale_shape_manual(values = c(21,23))+
  labs(x=expression(paste("P:C of allochtonous flows (",beta[A],")")),
       y=expression(paste("Seston P:C ratio (",beta[seston],")")))+
  the_theme2

p21=ggplot(d%>%
             mutate(., 
                    beta_A=recode_factor(beta_A,"0.002"="low","0.012"="high"))%>%
             dplyr::filter(.,
                           Simulation_ID=="Alpha_vary",
                           ID %in% unique(.$ID)[seq(1,40,by=8)]))+
  geom_line(aes(x=alpha_A,y=1/CN_seston,color=ID,group=ID),lwd=1)+
  geom_point(data=d%>%
               mutate(., 
                      beta_A=recode_factor(beta_A,"0.002"="low","0.012"="high"))%>%
               dplyr::filter(.,
                             Simulation_ID=="Alpha_vary",
                             ID %in% unique(.$ID)[seq(1,40,by=8)],
                             alpha_A %in% unique(.$alpha_A)[seq(1,40,by=5)]),
             aes(x=alpha_A,y=1/CN_seston,color=ID,shape=Limitation_Decompo),size=3,fill="white")+
  scale_colour_stepsn(colours = viridis(5,option = "G"))+
  facet_wrap(.~beta_A,scales = "free",
             labeller = label_bquote(cols=beta[A]==.(as.character(beta_A))))+
  scale_shape_manual(values = c(21,23))+
  labs(x=expression(paste("N:C of allochtonous flows (",alpha[A],")")),
       y=expression(paste("Seston N:C ratio (",alpha[seston],")")))+
  the_theme2


ggsave("./Figures/SI/Variation_NC_seston_input_nonfixer_and_decomposers.pdf",
       ggarrange(p11,p21,nrow=2,labels = letters[1:2],common.legend = T,legend="bottom",heights = c(1,1.05)),
       width = 7,height = 8)



# ------------------ Transient increase in stoichiometry ----

ode_lake_CNP_R_flexible = function(t,y,param){ 
  
  y[y < 10^-3] = 0 # prevent numerical problems
  
  iP = param$pulse_P(t)
  iN = param$pulse_N(t)
  
  with(as.list(c(y, param)),{
    
    
    
    if (DC_){ #donnor controlled
      
      gO_NP = (muO * min(P, N))             # growth rate non-fixer
      gF_P  = (muF * P)                     # growth rate fixer
      uptake_D = (eB * aD * DC)             # uptake detritus
      uptake_N = (aN * N)                   # uptake nitrogen
      uptake_P = (aP * P)                   # uptake phosphorous
      
    } else{
      if (functional_response_phyto==1){
        
        gO_NP = muO * min(P, N) * OC     # growth rate non-fixer
        gF_P  = muF * P * FC             # growth rate fixer
        
      }else{ #type 2
        
        gO_NP = muO * min(P/(kP+P),N/(kN+N)) * OC # growth rate non-fixer
        gF_P  = (muF * P)/(kP+P) * FC             # growth rate fixer
      }
      uptake_D = eB * aD * DC * BC             # uptake detritus
      uptake_N = aN * N * BC                   # uptake nitrogen
      uptake_P = aP * P * BC                   # uptake phosphorous
    }
    
    beta_D  = DP/DC                     # P:C ratio detritus allocht
    alpha_D = DN/DC                     # N:C ratio detritus
    
    if (same_stoichio){
      alpha_allo = alpha_D
      beta_allo  = beta_D
    }
    
    # we below write phi_i_j the decomposer function (immobilization/mineralization or decomposition) under the limitation j
    # phi_P_C would for instance be the immobilization/mineralization of P under C-limitation of decomposers
    
    #immobilization/mineralization P
    phi_P_C_lim = (uptake_D * (beta_B - beta_D)) / beta_B
    phi_P_N_lim = ((beta_B - beta_D) * alpha_B / ( (alpha_B - alpha_D) * beta_B)) * uptake_N
    phi_P_P_lim = uptake_P
    
    phi_P =  min(c(phi_P_C_lim,phi_P_N_lim,phi_P_P_lim))      # immobilization/mineralization of P
    phi_N =  (( (alpha_B - alpha_D) * beta_B) / ((beta_B - beta_D) * alpha_B)) * phi_P      # immobilization/mineralization of N
    phi_D =  beta_B  * phi_P / (beta_B  - beta_D)      # decomposition of detritus
    
    #Now the dynamics :
    du=rep(NA,14)
    
    #DECOMPOSERS     
    du[1] = phi_D                             -           dB * BC           - m * BC  -           sB * BC * BC  #BC
    du[2] = phi_D * alpha_D + phi_N * alpha_B - alpha_B * dB * BC - alpha_B * m * BC  - alpha_B * sB * BC * BC  #BN
    du[3] = phi_D * beta_D  + phi_P * beta_B  - beta_B  * dB * BC - beta_B  * m * BC  - beta_B  * sB * BC * BC  #BP
    
    #PLANKTON
    du[4] =            gF_P  - dF * FC    -            sF * FC * FC     #FC
    du[5] = alpha_F * (gF_P  - dF * FC)   - alpha_F  * sF * FC * FC     #FN
    du[6] = beta_F  * (gF_P  - dF * FC)   - beta_F   * sF * FC * FC     #FP
    du[7] =            gO_NP - dO * OC    -            sO * OC * OC     #OC
    du[8] = alpha_O * (gO_NP - dO * OC)   - alpha_O  * sO * OC * OC     #ON
    du[9] = beta_O  * (gO_NP - dO * OC)   - beta_O   * sO * OC * OC     #OP
    
    #DETRITUS 
    du[10] = ID - lD * DC + dB * BC + dF * FC + dO * OC - phi_D                #DC  
    
    du[11] = ID * (alpha_allo*iN) - lD * DN + alpha_B * dB * BC +                   #DN
      alpha_F * dF * FC + alpha_O * dO * OC - phi_D * alpha_D                                    
    
    du[12] = ID * (beta_allo*iP)  - lD * DP + beta_B  * dB * BC +                   #DP
      beta_F  * dF * FC + beta_O  * dO * OC - phi_D * beta_D                                    
    
    #NITROGEN 
    du[13] = IN - lN * N - alpha_O * gO_NP  - #- alpha_F * gF_P * FC       #N                           
      alpha_B * phi_N + alpha_B * m * BC
    
    #PHOSPHOROUS 
    du[14] = IP - lP * P - beta_O  * gO_NP   - beta_F  * gF_P  -        #P                      
      beta_B * phi_P  + beta_B  * m * BC
    
    #sum(du)
    #IP+ID+IN+ID * alpha_allo+ID * beta_allo-m * BC-lN * N-lP * P-lD * DC-lD * DN-lD * DP
    
    list(du)
    
  })
}


Get_transient_dynamics=function(nutrient="P",ID_=10,beta_A_=.002,alpha_A_=.02,
                                value_increase=2,sin=T){
  
  # reaching equilibrium
  p=Get_classical_param_lake()
  ini=Get_initial_values(p)
  
  p$functional_response_phyto=1
  p[c("sB","sF","sO")]=0.1
  p$IN=5
  p$ID=ID_
  p$IP=5
  p$beta_allo=beta_A_
  p$alpha_allo=alpha_A_
  
  
  p$pulse_P <- function(t_) {
    i_M_t = ifelse(t_ > 10 & t_ <= (20), 0, 0)
    return(i_M_t)
  }
  p$pulse_N <- function(t_) {
    i_M_t = ifelse(t_ > 10 & t_ <= (20), 0, 0)
    return(i_M_t)
  }
  
  d=Compute_ode(ini,p,julia = F,n_time = 4000,model = "pulse")
  Eq1=Get_equilibrium(d,p) #Equilibrium
  
  # perturbing the spatial flows
  ini[1:14]=as.numeric(Eq1[1:14])
  

  if (sin){
    
    if (nutrient =="P"){
      p$pulse_P <- function(t_) {
        #i_M_t = ifelse(t_ > 100 & t_ <= (200), value_increase, 0)
        i_M_t = (value_increase * (1+sin(10*t_)))
        return(i_M_t)
      }
      p$pulse_N <- function(t_) {
        i_M_t = ifelse(t_ > 100 & t_ <= (200), 0, 0)
        return(i_M_t)
      }
      
      
    }else {
      
      p$pulse_P <- function(t_) {
        i_M_t = ifelse(t_ > 100 & t_ <= (200), 0, 0)
        return(i_M_t)
      }
      p$pulse_N <- function(t_) {
        i_M_t = (value_increase * (1+sin(10*t_)))
        return(i_M_t)
      }
      
    }
    
  }else {
    if (nutrient =="P"){
      p$pulse_P <- function(t_) {
        i_M_t = ifelse(t_ > 100 & t_ <= (200), value_increase, 0)
        return(i_M_t)
      }
      p$pulse_N <- function(t_) {
        i_M_t = ifelse(t_ > 100 & t_ <= (200), 0, 0)
        return(i_M_t)
      }
      
      
    }else {
      
      p$pulse_P <- function(t_) {
        i_M_t = ifelse(t_ > 100 & t_ <= (200), 0, 0)
        return(i_M_t)
      }
      p$pulse_N <- function(t_) {
        i_M_t = ifelse(t_ > 100 & t_ <= (200), value_increase, 0)
        return(i_M_t)
      }
      
    }
  }
  d2=Compute_ode(ini,p,julia = F,n_time = 300,model = "pulse")
  d2$Time=c(4001:(4001+nrow(d2)-1))
  
  return(rbind(d,d2))
}



param_for_plot=expand.grid(
  nutrient=c("P"),
  value_increase = c(.5,2),
  ID=c(50),
  beta_A_=c(.012),
  alpha_A_=c(.14),
  sin=c(T,F))


all_sim=list()
for (i in 1:nrow(param_for_plot)){
  print(i)
  all_sim[[i]]=Get_transient_dynamics(
    nutrient=param_for_plot$nutrient[i],
    ID=param_for_plot$ID[i],value_increase = param_for_plot$value_increase[i],
    beta_A_=param_for_plot$beta_A_[i],
    alpha_A_=param_for_plot$alpha_A_[i],
    sin=param_for_plot$sin[i]
  )
}



p2=ggplot(all_sim[[2]]%>%
            add_column(.,
                       C_seston = .$Decomposers_C + .$Fixers_C + .$Non_fixers_C + .$Detritus_C,
                       N_seston = .$Decomposers_N + .$Fixers_N + .$Non_fixers_N + .$Detritus_N,
                       P_seston = .$Decomposers_P + .$Fixers_P + .$Non_fixers_P + .$Detritus_P)%>%
            add_column(.,
                       NC_seston = .$N_seston/.$C_seston,
                       PC_seston = .$P_seston/.$C_seston)%>%dplyr::mutate(., Time=Time-3950)%>%
            melt(., measure.vars=c("PC_seston","NC_seston","Fixers_C","Decomposers_C","Non_fixers_C"))%>%
            mutate(., variable=recode_factor(variable,
                                             "PC_seston"="P:C seston","NC_seston"="N:C seston",
                                             "Decomposers_C"="Decomposers","Fixers_C"="Fixers","Non_fixers_C"="Non-fixers")))+
  geom_line(aes(x=Time,y=value),lwd=1)+
  # geom_hline(data = tibble(variable=c("P:C seston","N:C seston"),value=c(1/106,1/16)),aes(yintercept = value))+
  facet_wrap(.~variable,scales="free")+the_theme2+
  xlim(1,150)+
  theme(strip.text.x = element_text(size=12),legend.text = element_text(size=15))+
  guides(color = guide_legend(override.aes = list(size = 10)))



p3=ggplot(all_sim[[3]]%>%
            add_column(.,
                       C_seston = .$Decomposers_C + .$Fixers_C + .$Non_fixers_C + .$Detritus_C,
                       N_seston = .$Decomposers_N + .$Fixers_N + .$Non_fixers_N + .$Detritus_N,
                       P_seston = .$Decomposers_P + .$Fixers_P + .$Non_fixers_P + .$Detritus_P)%>%
            add_column(.,
                       NC_seston = .$N_seston/.$C_seston,
                       PC_seston = .$P_seston/.$C_seston)%>%dplyr::mutate(., Time=Time-3950)%>%
            melt(., measure.vars=c("PC_seston","NC_seston","Fixers_C","Decomposers_C","Non_fixers_C"))%>%
            mutate(., variable=recode_factor(variable,
                                             "PC_seston"="P:C seston","NC_seston"="N:C seston",
                                             "Decomposers_C"="Decomposers","Fixers_C"="Fixers","Non_fixers_C"="Non-fixers")))+
  geom_line(aes(x=Time,y=value),lwd=1)+
  # geom_hline(data = tibble(variable=c("P:C seston","N:C seston"),value=c(1/106,1/16)),aes(yintercept = value))+
  geom_rect(data=tibble(xmin=150,xmax=250,ymin=-Inf,ymax=Inf),
            aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax),fill="grey",alpha=.5)+
  facet_wrap(.~variable,scales="free")+the_theme2+
  theme(strip.text.x = element_text(size=12),legend.text = element_text(size=15))+
  guides(color = guide_legend(override.aes = list(size = 10)))+  xlim(0,350)



p4=ggplot(all_sim[[4]]%>%
            add_column(.,
                       C_seston = .$Decomposers_C + .$Fixers_C + .$Non_fixers_C + .$Detritus_C,
                       N_seston = .$Decomposers_N + .$Fixers_N + .$Non_fixers_N + .$Detritus_N,
                       P_seston = .$Decomposers_P + .$Fixers_P + .$Non_fixers_P + .$Detritus_P)%>%
            add_column(.,
                       NC_seston = .$N_seston/.$C_seston,
                       PC_seston = .$P_seston/.$C_seston)%>%dplyr::mutate(., Time=Time-3950)%>%
            melt(., measure.vars=c("PC_seston","NC_seston","Fixers_C","Decomposers_C","Non_fixers_C"))%>%
            mutate(., variable=recode_factor(variable,
                                             "PC_seston"="P:C seston","NC_seston"="N:C seston",
                                             "Decomposers_C"="Decomposers","Fixers_C"="Fixers","Non_fixers_C"="Non-fixers")))+
  geom_line(aes(x=Time,y=value),lwd=1)+
  # geom_hline(data = tibble(variable=c("P:C seston","N:C seston"),value=c(1/106,1/16)),aes(yintercept = value))+
  geom_rect(data=tibble(xmin=150,xmax=250,ymin=-Inf,ymax=Inf),
            aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax),fill="grey",alpha=.5)+
  facet_wrap(.~variable,scales="free")+the_theme2+
  theme(strip.text.x = element_text(size=12),legend.text = element_text(size=15))+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  xlim(0,350)


p_tot=ggarrange(p2+labs(y="")+ggtitle("Seasonal change of allochtonous P:C"),
                p4+labs(y="")+ggtitle("Transient increase of allochtonous P:C"),nrow = 2)


ggsave("./Figures/SI/Seasonal_change.pdf",p_tot,width = 9,height = 10)




