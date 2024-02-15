rm(list=ls())

source("CNP_lake_functions.R")

the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90),
                                strip.text.x = element_text(size = 8),axis.text = element_text(size=11),axis.title = element_text(size=13),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))


pal=colorRampPalette(c("#5C7ECC","#32C4E2","#7DCE9A","#B4CE7D","#EFDD35","#FFB700"))
pal_aqua=colorRampPalette((c("white","#CFDCDE","#A4DEE6","#49CADC","#3B7FE0","#193C82","#060D61")))
pal_terr=colorRampPalette((c("white","#D8ECCD","#BBE0A7","#86D45C","#3BA23B","#066F06","#033E03")))


the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90),
                                strip.text.x = element_text(size = 8),axis.text = element_text(size=11),axis.title = element_text(size=13),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))

