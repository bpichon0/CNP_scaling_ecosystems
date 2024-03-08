rm(list=ls())

source("CNP_lake_functions.R")

# Step 1: Allochtonous flows ----

## >> Along gradient ----

d=read.table("../Table/Simulation_allochtonous.csv",sep=";")
table(d$Limitation_Decompo)
table(d$Limitation_NF)

list_plots=list()
id=1

for (x2 in colnames(d)[c(1,4,7,27,13,14)]){
  
  for (x1 in unique(d$Simulation_ID)){
    
    d2=Plot_with_limitation(d%>%filter(., Simulation_ID==x1),
                         c("IN","IP","ID")[x1],x2)
    
    list_plots[[id]]=d2$p+
             geom_line(data=d2$data,
                       aes(x=value_driver,y=value_y))+
             the_theme+
             labs(x=c("IN","IP","ID")[x1],y=x2)+theme(legend.position = "none")
    
    id=id+1
    
  }
}

p_tot=ggarrange(plotlist=list_plots,ncol = 3,nrow = 6,common.legend = T,legend="none")
ggsave("../Figures/Allochtonous_flows1.pdf",p_tot,width = 10,height = 12)


d=read.table("../Table/Simulation_allochtonous.csv",sep=";")
table(d$Limitation_Decompo)
table(d$Limitation_NF)

list_plots=list()
id=1

for (x2 in colnames(d)[c(21:26)]){
  
  for (x1 in unique(d$Simulation_ID)){
    
    d2=Plot_with_limitation(d%>%filter(., Simulation_ID==x1),
                            c("IN","IP","ID")[x1],x2)
    
    list_plots[[id]]=d2$p+
      geom_line(data=d2$data,
                aes(x=value_driver,y=value_y))+
      the_theme+
      labs(x=c("IN","IP","ID")[x1],y=x2)+theme(legend.position = "none")
    
    id=id+1
    
  }
}

p_tot=ggarrange(plotlist=list_plots,ncol = 3,nrow = 6,common.legend = T,legend="none")
ggsave("../Figures/Allochtonous_flows2.pdf",p_tot,width = 10,height = 12)


## >> Along gradient, indirect effects ----

pdf("../Figures/Allochtonous_flows_indirect.pdf",width = 8,height = 8)
d=read.table("../Table/Simulation_allochtonous_indirect.csv",sep=";")%>%
  filter(., Type=="A_to_I")

list_plots=list()
id=1
par(mfrow=c(3,3))
for (x1 in unique(d$Simulation_ID)){
  
  d_fil=d%>%filter(., Simulation_ID==x1)%>%
    add_column(., Driver=.[,c("IN","IP","ID")[x1]])%>%
    filter(., Driver %in% unique(.$Driver)[c(ifelse(x1==2,11,5),15,100)])
  
  Plot_net_effects(d_fil[1,],param="",compute_effects = F,Indirect = d_fil[1,])
  mtext(paste0("Low ",c("IN","IP","ID")[x1]))
  Plot_net_effects(d_fil[2,],param="",compute_effects = F,Indirect = d_fil[2,])
  mtext(paste0("Medium ",c("IN","IP","ID")[x1]))
  Plot_net_effects(d_fil[3,],param="",compute_effects = F,Indirect = d_fil[3,])
  mtext(paste0("High ",c("IN","IP","ID")[x1]))
  
  id=id+1
  
}
dev.off()


## >> 2D gradient CNP seston ----

d=read.table("../Table/Simulation_allochtonous_2D.csv",sep=";")

d%>%
  filter(., d$Limitation_Decompo=="P" & converged_==T)%>%
  ggplot(.)+
  geom_tile(aes(x=IN,y=ID,fill=CP_seston))+
  scale_fill_viridis_c(option = "A")+
  facet_wrap(.~IP,scales = "free")+
  the_theme

d%>%
  filter(., d$Limitation_Decompo=="P" & converged_==T)%>%
  ggplot(.)+
  geom_tile(aes(x=IN,y=ID,fill=CN_seston))+
  scale_fill_viridis_c(option = "A")+
  facet_wrap(.~IP,scales = "free")+
  the_theme

d%>%
  filter(., d$Limitation_Decompo=="C" & converged_==T)%>%
  ggplot(.)+
  geom_tile(aes(x=IN,y=ID,fill=NP_seston))+
  scale_fill_viridis_c(option = "A")+
  facet_wrap(.~IP,scales = "free")+
  the_theme


# Step 2: Stoichio allochtonous ----


d=read.table("../Table/Simulation_allochtonous_ratios.csv",sep=";")

list_plots=list()
id=1

for (x2 in colnames(d)[c(1,4,7,27,13,14,21:26)]){
  
  for (x1 in unique(d$Ratio_changed)){
    
    d2=Plot_with_limitation(d%>%filter(., Ratio_changed==x1,ID==10),
                            c("beta_allo","alpha_allo")[ifelse(x1=="PC",1,2)],x2)
    
    list_plots[[id]]=d2$p+
      geom_line(data=d2$data,
                aes(x=value_driver,y=value_y))+
      the_theme+
      labs(x=c("P:C allochtonous detritus","N:C allochtonous detritus")[ifelse(x1=="PC",1,2)],
           y=x2)+theme(legend.position = "none")
    
    id=id+1
    
  }
}

p_tot=ggarrange(plotlist=list_plots,ncol = 2,nrow = 10,common.legend = T,legend="none")
ggsave("../Figures/Allochtonous_flows_ratio.pdf",p_tot,width = 7,height = 20)



# Step 3: Self-regulation ----

d=read.table("../Table/Simulation_allochtonous_self_regulation.csv",sep=";")

list_plots=list()
id=1

for (x2 in colnames(d)[c(1,4,7,13,14)]){
  
  for (IN_ in unique(d$IN)){
    for (ID_ in unique(d$ID)){
      
      d2=Plot_with_limitation(d%>%filter(., ID==ID_,IN==IN_,converged_==T),
                              "self",x2)
      
      list_plots[[id]]=d2$p+
        geom_line(data=d2$data,
                  aes(x=value_driver,y=value_y))+
        the_theme+
        labs(x="Self-regulation",y=x2)+theme(legend.position = "none")
      
      id=id+1
    }
    
  }
}

p_tot=ggarrange(plotlist=list_plots,ncol = 4,nrow = 5,common.legend = T,legend="none")
ggsave("../Figures/Allochtonous_flows_self_1.pdf",p_tot,width = 12.5,height = 10)


d=read.table("../Table/Simulation_allochtonous_self_regulation.csv",sep=";")

list_plots=list()
id=1

for (x2 in colnames(d)[21:26]){
  
  for (IN_ in unique(d$IN)){
    for (ID_ in unique(d$ID)){
      
      d2=Plot_with_limitation(d%>%filter(., ID==ID_,IN==IN_,converged_==T),
                              "self",x2)
      
      list_plots[[id]]=d2$p+
        geom_line(data=d2$data,
                  aes(x=value_driver,y=value_y))+
        the_theme+
        labs(x="Self-regulation",y=x2)+theme(legend.position = "none")
      
      id=id+1
    }
    
  }
}
p_tot=ggarrange(plotlist=list_plots,ncol = 4,nrow = 6,common.legend = T,legend="none")
ggsave("../Figures/Allochtonous_flows_self_2.pdf",p_tot,width = 12.5,height = 12)


# Step 4: Self-regulation, each ----
  
for (orga_ in c("B","F","O")){  
  
  d=read.table("../Table/Simulation_allochtonous_self_regulation_each.csv",sep=";")%>%
    filter(., organism==orga_)

  list_plots=list()
  id=1
  
  for (x2 in colnames(d)[c(1,4,7,13,14)]){
    
    for (IN_ in unique(d$IN)){
      for (ID_ in unique(d$ID)){
        
        d2=Plot_with_limitation(d%>%filter(., ID==ID_,IN==IN_,converged_==T),
                                "self",x2)
        
        list_plots[[id]]=d2$p+
          geom_line(data=d2$data,
                    aes(x=value_driver,y=value_y))+
          the_theme+
          labs(x="Self-regulation",y=x2)+theme(legend.position = "none")
        
        id=id+1
      }
    }
  }
  
  p_tot=ggarrange(plotlist=list_plots,ncol = 4,nrow = 5,common.legend = T,legend="none")
  ggsave(paste0("../Figures/Allochtonous_flows_self_",orga_,"_1.pdf"),
         p_tot,width = 12.5,height = 10)
  

  list_plots=list()
  id=1
  
  for (x2 in colnames(d)[21:26]){
    
    for (IN_ in unique(d$IN)){
      for (ID_ in unique(d$ID)){
        
        d2=Plot_with_limitation(d%>%filter(., ID==ID_,IN==IN_,converged_==T),
                                "self",x2)
        
        list_plots[[id]]=d2$p+
          geom_line(data=d2$data,
                    aes(x=value_driver,y=value_y))+
          the_theme+
          labs(x="Self-regulation",y=x2)+theme(legend.position = "none")
        
        id=id+1
      }
    }
  }
  p_tot=ggarrange(plotlist=list_plots,ncol = 4,nrow = 6,common.legend = T,legend="none")
  ggsave(paste0("../Figures/Allochtonous_flows_self_",orga_,"_2.pdf"),
         p_tot,width = 12.5,height = 12)
}

# Step 5: Stoichiometric ratio decomposers & self-regulation decomposers----

d=read.table("../Table/Simulation_traits_decomposers.csv",sep=";")
d$Deviation_Redfield=sapply(1:nrow(d),function(x){
  return(sqrt(sum((c(d$CN_seston[x],d$CP_seston[x],d$NP_seston[x])-
                     c(106/16,106,16))**2)))
})

list_plots=list()
id=1

for (x2 in colnames(d)[21:26]){
  for (IN_ in unique(d$IN)){
    for (ID_ in unique(d$ID)){
      
      list_plots[[id]]=ggplot(d%>%filter(., ID==ID_,IN==IN_,sB>min(.$sB)))+
        geom_tile(aes(x=alphaB,y=sB,fill=Deviation_Redfield))+
        scale_fill_viridis_c(option = "A")+
        the_theme+
        labs(x="alpha_B",y="sB",fill=x2)
      
      id=id+1
    }
    
  }
}

p_tot=ggarrange(plotlist=list_plots,ncol = 4,nrow = 6,common.legend = T,legend="none")
ggsave("../Figures/Allochtonous_flows_decomposers_traits.pdf",p_tot,width = 12.5,height = 10)



# Step 6: Random parameters simulations ----

d=read.table("../Table/Simulation_random_parameters.csv",sep=";")
d$Limitation = paste0("Decomp: ",d$Limitation_Decompo,", Non-fix. ",d$Limitation_NF)
d$Deviation_Redfield=sapply(1:nrow(d),function(x){
  return(sqrt(sum((c(d$CN_seston[x],d$CP_seston[x],d$NP_seston[x])-
                     c(106/16,106,16))**2)))
})

p=ggplot(d%>%
           melt(., measure.vars = c("Deviation_Redfield","CN_seston",
                                    "CP_seston","NP_seston","Decomposers_C",
                                    "Fixers_C","Non_fixers_C","NC_detritus","Frac_decomp")))+
  geom_boxplot_pattern(aes(y=value,x=Limitation,fill=Limitation_Decompo,pattern=Limitation_NF),
                       pattern_density = 0.01)+
  facet_wrap(variable~.,scales = "free",nrow = 3)+
  scale_pattern_manual(values=c("N"="none","P"="stripe"))+
  labs(x="",y="",fill="Limitation Decomp.",pattern="Limitation Non-Fix.")+
  scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
  the_theme+
  theme(axis.text.x = element_blank())

ggsave("../Figures/Random/Boxplot_limitation.pdf",p,width = 10,height = 7)

p=ggarrange(
  ggarrange(
    ggplot(d)+
      geom_point(aes(x=CP_seston,y=CN_seston,color=alpha_B))+
      labs(x="CP seston",y="CN seston",fill=TeX("$\\alpha_B$"),pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=CP_seston,y=NP_seston,color=alpha_B))+
      labs(x="CN seston",y="NP seston",fill=TeX("$\\alpha_B$"),pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=CN_seston,y=NP_seston,color=alpha_B))+
      labs(x="CN seston",y="NP seston",fill=TeX("$\\alpha_B$"),pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ncol=3
  ),
  ggarrange(
    ggplot(d)+
      geom_point(aes(x=CP_seston,y=CN_seston,color=beta_B))+
      labs(x="CP seston",y="CN seston",fill=TeX("$\\beta_B$"),pattern="")+
      scale_color_viridis_c()+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=CP_seston,y=NP_seston,color=beta_B))+
      labs(x="CN seston",y="NP seston",fill=TeX("$\\beta_B$"),pattern="")+
      scale_color_viridis_c()+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=CN_seston,y=NP_seston,color=beta_B))+
      labs(x="CN seston",y="NP seston",fill=TeX("$\\beta_B$"),pattern="")+
      scale_color_viridis_c()+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ncol=3
  ),ggarrange(
    ggplot(d)+
      geom_point(aes(x=CP_seston,y=CN_seston,color=sB))+
      labs(x="CP seston",y="CN seston",color=TeX("$\\s_B$"),pattern="")+
      scale_color_viridis_c(option = "A")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=CP_seston,y=NP_seston,color=sB))+
      labs(x="CN seston",y="NP seston",color=TeX("$\\s_B$"),pattern="")+
      scale_color_viridis_c(option = "A")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=CN_seston,y=NP_seston,color=sB))+
      labs(x="CN seston",y="NP seston",color=TeX("$\\s_B$"),pattern="")+
      scale_color_viridis_c(option = "A")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ncol=3
  ),nrow=3,labels = LETTERS[1:3]
)
  

ggsave("../Figures/Random/Large_influence_stoichio_decompo_ratio.pdf",p,width = 10,height = 10)

p=ggarrange(
  ggarrange(
    ggplot(d)+
      geom_point(aes(x=C_seston,y=N_seston,color=IP))+
      geom_line(data=tibble(x=c(0,max(d$C_seston)),y=c(0,16*max(d$C_seston)/106)),aes(x=x,y=y),lwd=1)+
      labs(x="C seston",y="N seston",fill="Inflow P",pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=C_seston,y=P_seston,color=IP))+
      geom_line(data=tibble(x=c(0,max(d$C_seston)),y=c(0,max(d$C_seston)/106)),aes(x=x,y=y),lwd=1)+
      labs(x="C seston",y="P seston",fill="Inflow P",pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=N_seston,y=P_seston,color=IP))+
      geom_line(data=tibble(x=c(0,max(d$N_seston)),y=c(0,max(d$N_seston)/16)),aes(x=x,y=y),lwd=1)+
      labs(x="N seston",y="P seston",fill="Inflow P",pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ncol=3
  ),
  ggarrange(
    ggplot(d)+
      geom_point(aes(x=C_seston,y=N_seston,color=Frac_decomp))+
      geom_line(data=tibble(x=c(0,max(d$C_seston)),y=c(0,16*max(d$C_seston)/106)),aes(x=x,y=y),lwd=1)+
      labs(x="C seston",y="N seston",color="Level of heterotrophy \n (fraction decompo)",pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=C_seston,y=P_seston,color=Frac_decomp))+
      geom_line(data=tibble(x=c(0,max(d$C_seston)),y=c(0,max(d$C_seston)/106)),aes(x=x,y=y),lwd=1)+
      labs(x="C seston",y="P seston",color="Level of heterotrophy \n (fraction decompo)",pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=N_seston,y=P_seston,color=Frac_decomp))+
      geom_line(data=tibble(x=c(0,max(d$N_seston)),y=c(0,max(d$N_seston)/16)),aes(x=x,y=y),lwd=1)+
      labs(x="N seston",y="P seston",color="Level of heterotrophy \n (fraction decompo)",pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ncol=3
  ),ggarrange(
    ggplot(d)+
      geom_point(aes(x=C_seston,y=N_seston,color=sF))+
      geom_line(data=tibble(x=c(0,max(d$C_seston)),y=c(0,16*max(d$C_seston)/106)),aes(x=x,y=y),lwd=1)+
      labs(x="C seston",y="N seston",color="sF",pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=C_seston,y=P_seston,color=sF))+
      geom_line(data=tibble(x=c(0,max(d$C_seston)),y=c(0,max(d$C_seston)/106)),aes(x=x,y=y),lwd=1)+
      labs(x="C seston",y="P seston",color="sF",pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ggplot(d)+
      geom_point(aes(x=N_seston,y=P_seston,color=sF))+
      geom_line(data=tibble(x=c(0,max(d$N_seston)),y=c(0,max(d$N_seston)/16)),aes(x=x,y=y),lwd=1)+
      labs(x="N seston",y="P seston",color="sF",pattern="")+
      scale_color_viridis_c(option = "F")+
      #scale_fill_manual(values=c("#D2B96F","#8EBAEF","#D7B2F9"))+
      the_theme,
    ncol=3
  ),nrow=3,labels = LETTERS[1:3]
)

ggsave("../Figures/Random/Large_influence_stoichio_decompo_elements.pdf",p,width = 10,height = 10)


# PCA
metric_PCA=colnames(d)[c(21:26,38)]#[c(1,2,3,4,7,21:26,38)]
res.pca=PCA(d[,which(colnames(d) %in% metric_PCA)], scale.unit = T, ncp = 3,  graph=F)
axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

for (metric_color in c("ID","IN","IP","sB","sO","sF","alpha_B","beta_B")){
  
  for (i in 1:3){
    assign(paste0("p",i),
           d%>%
             melt(., measure.vars = metric_color)%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = value,fill=value),alpha=.5)+
             scale_fill_viridis_c()+
             scale_color_viridis_c()+
             # geom_point(data=centroids,aes(x=centroids[,axes_for_plot$x[i]+1],y=centroids[,axes_for_plot$y[i]+1]),shape=24,fill="black",color="white")+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="",fill="")+
             ggtitle("")+guides(shape="none")+
             theme_classic()+theme(legend.position = "bottom")+
             guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")  
    )
  }
  
  p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                        p2+theme(legend.position = "none"),
                        p3+theme(legend.position = "none"),
                        ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
              nrow=2,heights = c(1,.1))
  ggsave(paste0("../Figures/Random/PCA_colored_by_",metric_color,".pdf"),p, width=9,height = 4)
}

