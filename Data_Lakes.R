rm(list=ls())
source("./CNP_lake_functions.R")


# 1) Global lake database ----

library(LAGOSNE)


data_lake=lagosne_load(fpath = "./Data/data_1.087.1.rds")

# Pull lake nutrient data for variance components paper  
# Criteria
# 1: Include lakes from 1990 - 2011
# 2: Use epi data from summer stratified period (15 June - 15 September)
# 3: Take median nutrient value from each lake within a season (single value/season)
# 4: If TN data don't exits, but TKH and N02N03 do - then calculate TN for those lakes as TKN + N02N03
# 5: Response variables (4): TP, TN, nitrate, Secchi

# Turn epi.nutr data frame into data table
d_CNP=data.frame(data_lake$epi_nutr)%>%
  filter(., sampleyear >=1990 & sampleyear <= 2011)%>%
  add_column(., TN_computed = .$tkn + .$no2no3)%>%
  dplyr::select(., eventida1087,lagoslakeid,programname,programtype,lagosversion,
                sampledate,chla,colora,colort,dkn,doc,nh4,
                no2,no2no3,srp,tdn,tn,toc,ton,tp,secchi,sampleyear,
                samplemonth,TN_computed)

d_CNP$TN_new=  ifelse(is.na(d_CNP$tn), d_CNP$TN_computed, d_CNP$tn)
d_CNP$sampledate=ymd(d_CNP$sampledate) #cleaning dates
d_CNP$day=day(d_CNP$sampledate) #day 

# Use epi data from summer stratified period (15 June - 15 September)
d_CNP = d_CNP[d_CNP$samplemonth == 6 & (d_CNP$day >= 15) | d_CNP$samplemonth==7 | 
                d_CNP$samplemonth==8 | d_CNP$samplemonth== 9 & (d_CNP$day <= 15), ]
d_CNP$Area=sapply(1:nrow(d_CNP),function(x){
  return(data_lake$locus$lake_area_ha[which(d_CNP$lagoslakeid[x]==data_lake$locus$lagoslakeid)])
})






# 2) Hessen data n1: Comsat = OK ----

d1=readxl::read_xlsx("./Data/COMSAT.xlsx")%>%
  mutate(., 
         Temperature = as.numeric(.$Temperature),
         Depth = log(as.numeric(.$Depth)),
         Secchi = log(as.numeric(.$Secchi)),
         Area = log(as.numeric(.$Area)))%>%
  add_column(., 
             TOC_TN=.$TOC/.$TN,
             TP_TN=.$TP/.$TN,
             TOC_TP=.$TOC/.$TP
  )%>%
  mutate(., 
         TOC_TP = log(TOC_TP),
         TP_TN = sqrt(TP_TN))
  
colnames(d1)

d_slope=d_residuals=tibble()
for (stoichio_var in c("TOC_TN","TOC_TP","TP_TN")){
  
  d_mod=d1%>%melt(., measure.vars = stoichio_var)
  d_mod[,c("Latitude","Longitude","Area","Altitude","Depth","Temperature","Secchi","pH","value")]=
    apply(d_mod[,c("Latitude","Longitude","Area","Altitude","Depth","Temperature","Secchi","pH","value")],2,scale)
  
  mod_stoichio=(lm(value~Latitude+Longitude+Area+Altitude+Depth+Temperature+Secchi+pH,data=d_mod,))
  res=visreg::visreg(mod_stoichio,xvar="Secchi",plot=F)
  d_residuals=rbind(d_residuals,tibble(Residuals=res$res$visregRes,xvar=res$res$Secchi,
                                       Stoichio_ratio=stoichio_var,Name_xvar="Secchi"))
  res=visreg::visreg(mod_stoichio,xvar="Area",plot=F)
  d_residuals=rbind(d_residuals,tibble(Residuals=res$res$visregRes,xvar=res$res$Area,
                                       Stoichio_ratio=stoichio_var,Name_xvar="Area"))
  
  # 
  # 
  # d_slope=rbind(d_slope,
  #               tibble(slope=summary(mod_stoichio)$coefficient[2:9,1],
  #                      Predictor=rownames(summary(mod_stoichio)$coefficient)[2:9],
  #                      Low_int=confint(mod_stoichio)[2:9,1],
  #                      High_int=confint(mod_stoichio)[2:9,2],
  #                      With_secchi=T,
  #                      Stoichio_ratio=stoichio_var
  #               ))

}


p=ggplot(d_residuals%>%
         mutate(.,Stoichio_ratio=recode_factor(Stoichio_ratio,
                                               "TOC_TN"="TOC:TN",
                                               "TOC_TP"="TOC:TP",
                                               "TP_TN"="TP:TN") ))+
  geom_point(aes(x=xvar,y=Residuals),color="grey",alpha=.5)+
  geom_smooth(aes(x=xvar,y=Residuals),method = "lm")+
  facet_grid(Stoichio_ratio~Name_xvar,scales = "free")+
  the_theme2+
  labs(x="",y="Partial residuals stoichiometric ratios")
ggsave("../Figures/Empirical_data/Data_Comsat.pdf",p,width = 6,height = 6)


# ggplot(d_slope%>%filter(., With_secchi==T))+
#   geom_pointrange(aes(x=slope,xmin=Low_int,xmax=High_int,y=Predictor),color="lightblue")+
#   the_theme+
#   facet_wrap(.~Stoichio_ratio,scales = "free")+
#   geom_vline(xintercept = 0)


# 3) Hessen data n2: Hessen et al 2009 ----

d=as_tibble(apply(readxl::read_xlsx("./Data/1000Lakes.xlsx"),2,as.numeric))%>%
  dplyr::select(.,-starts_with("t."),-starts_with("SDEP"),-starts_with("STD"),-starts_with("TRS"))%>%
  add_column(., 
             TOC_TN=.$TOC/.$N.tot,
             TOC_TON=.$TOC/.$TON,
             TN_TP=.$N.tot/.$P.tot,
             TOC_TP=.$TOC/.$P.tot)%>%
  dplyr::rename(., 
                Lat=LATITUDE,Long=LONGITUDE, 
                Precipitation=PRECIP,Temperature=TEMP...19,
                Arable=ARABLE,Allochot_flux=TOCFLUX,Altitude=ALTI,
                Area=LAKE.AREA,Bog=BOG,N_deposition=NDEP,
                Tree_cover=Treecov_AVG,Herb_cover=Herbcov_AVG,
                Runoff=RUNOFF,Modis=Modis_AVG,Drainage_area=DRAIN_AREA,
                Slope=SLOPE,Lake_shape=LAKE)%>%
  add_column(., Long_sin=sin(.$Long),Long_cos=cos(.$Long))

d=d%>%
  Normalize_data(.)%>%
  Standardize_data(.)

d_slope=d_residuals=tibble()
for (stoichio_var in c("TOC_TN","TOC_TP","TN_TP")){
  
  d_mod=d%>%melt(., measure.vars = stoichio_var)
    
  mod_stoichio=(lm(value~Lat+Long_sin+Long_cos+Altitude+Temperature+pH+Area+Slope+N_deposition+NDVI+
                     Herb_cover+Tree_cover+Bog+Allochot_flux,data=d_mod))
  
  res=visreg::visreg(mod_stoichio,"Area",plot = F)
  d_residuals=rbind(d_residuals,tibble(Residuals=res$res$visregRes,xvar=res$res$Area,
                                       Stoichio_ratio=stoichio_var,Name_xvar="Area"))
  
  res=visreg::visreg(mod_stoichio,"Tree_cover",plot = F)
  d_residuals=rbind(d_residuals,tibble(Residuals=res$res$visregRes,xvar=res$res$Tree_cover,
                                       Stoichio_ratio=stoichio_var,Name_xvar="Tree_cover"))
  
  res=visreg::visreg(mod_stoichio,"NDVI",plot = F)
  d_residuals=rbind(d_residuals,tibble(Residuals=res$res$visregRes,xvar=res$res$NDVI,
                                       Stoichio_ratio=stoichio_var,Name_xvar="NDVI"))
  
  res=visreg::visreg(mod_stoichio,"N_deposition",plot = F)
  d_residuals=rbind(d_residuals,tibble(Residuals=res$res$visregRes,xvar=res$res$N_deposition,
                                       Stoichio_ratio=stoichio_var,Name_xvar="N_deposition"))
  
  res=visreg::visreg(mod_stoichio,"Allochot_flux",plot = F)
  d_residuals=rbind(d_residuals,tibble(Residuals=res$res$visregRes,xvar=res$res$Allochot_flux,
                                       Stoichio_ratio=stoichio_var,Name_xvar="Allochot_flux"))
  
}


p=ggplot(d_residuals%>%
           mutate(.,Stoichio_ratio=recode_factor(Stoichio_ratio,
                                                 "TOC_TN"="TOC:TN",
                                                 "TOC_TP"="TOC:TP",
                                                 "TN_TP"="TN:TP") )%>%filter(., Name_xvar=="Area"))+
  geom_point(aes(x=xvar,y=Residuals),color="grey",alpha=.5)+
  geom_smooth(aes(x=xvar,y=Residuals),method = "lm")+
  facet_grid(Stoichio_ratio~Name_xvar,scales = "free")+
  the_theme2+
  labs(x="",y="Partial residuals stoichiometric ratios")
ggsave("../Figures/Empirical_data/Data_Hessen.pdf",p,width = 8,height = 6)



mod_TCTN=(lm(TOC_TN~Lat+Long_sin+Long_cos+Altitude+Temperature+pH+Area+Slope+N_deposition+NDVI+#Drainage_area+
          Herb_cover+Tree_cover+Bog+Allochot_flux,data=d,na.action = "na.omit"))


mod_TNTP=(lm(TN_TP~Lat+Long_sin+Long_cos+Altitude+Temperature+pH+Area+Slope+N_deposition+NDVI+#Drainage_area+
          Herb_cover+Tree_cover+Bog+Allochot_flux,data=d,na.action = "na.omit"))

visreg::visreg(mod_TNTP,"Area")
visreg::visreg(mod_TNTP,"Bog")
visreg::visreg(mod_TNTP,"Tree_cover")
visreg::visreg(mod_TNTP,"Herb_cover")
visreg::visreg(mod_TNTP,"NDVI")
visreg::visreg(mod_TNTP,"Allochot_flux")
visreg::visreg(mod_TNTP,"N_deposition")


#SEM
model_lm=lm("value ~ Long_sin+Long_cos + Lat + Altitude + Slope + pH + Temperature",
              data = d%>%melt(., measure.vars="TOC_TN")%>%filter(., !is.na(value)),
              na.action = na.omit)

resid_model=residuals(model_lm) #extract residuals of the distance to the tipping point after controlling for all covariates

d_sem=d[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
  add_column(., Resid_mod=resid_model)

SEM_TCTN=psem(
  lm(Resid_mod ~ Area+N_deposition+Allochot_flux, d_sem),
  lm(Area ~ Allochot_flux+N_deposition, d_sem),
  lm(Allochot_flux ~ Tree_cover+Herb_cover+NDVI+Bog, d_sem),
  lm(N_deposition ~ Tree_cover+Herb_cover+NDVI+Bog, d_sem))

summary(SEM_TCTN)




