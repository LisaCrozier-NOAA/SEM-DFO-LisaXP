source("LisaXP/functions/sem_plot_fxn.r")
load(file="LisaXP/outputs_4/DAG1.topmodels.rdata",verbose=T)

sem_output_summary_fxn(DAG1Blong.topmodel)
sem_output_summary_fxn(DAG1Bshort.topmodel.reduced)
sem_output_summary_fxn(DAG1A_short_topmodel)
sem_output_summary_fxn(DAG1Ashort.topmodel.reduced)

#notes: DAG1Along & DAGDAG1A_short_topmodel_reduced#notes: DAG1Along & DAG1Ashort & DAG1Bshort: X12_copepodBiomass_WGoA
#notes: DAG1Blong: X12_DFA_biomassEuphShelfSum
# model_id    PreyNCCindNames       PredNCCindNames             PreyAKindNames          PredAKindNames             
# 1 DAG1A_long  05.ForageFishNCC_DFA1 09.PredFishNCC_b_DFA1       12.ZooPreyAK_0_smoothed 10.PredMammalNCC_1_smoothed
# 2 DAG1A_short 05.ForageFishNCC_DFA1 08.PredBirdNCC_b_0_smoothed 12.ZooPreyAK_0_smoothed 10.PredMammalNCC_0_smoothed
# 3 DAG1B_long  05.ForageFishNCC_DFA1 09.PredFishNCC_b_DFA1       12.ZooPreyAK_DFA1       10.PredMammalNCC_1_smoothed
# 4 DAG1B_short 05.ForageFishNCC_DFA1 08.PredBirdNCC_b_0_smoothed 12.ZooPreyAK_0_smoothed 10.PredMammalNCC_0_smoothed

# 1. Define the exact names as they appear in your model string
# (You can also map these to shorter nicknames if you use the rename_with trick)

var_x01a    <- "X01.ZooPreyNCC_JSOES_DFA_sumPreyOfPrey_planktonJun"
var_x01b    <- "X01.ZooPreyNCC_JSOES_Cassin_s_aukl_WS"
var_x03    <- "X03.FishPreyNCC_b_HakeAge1"
var_x04    <- "X04.CompFishNCC_smoothed_marketsquid_GAM"
var_x05      <- "X05_DFA_abundSardine"
var_x08      <- "X08_commonMurre_JSOES"
var_x08b      <- "X08.PredBirdNCC_hump_Common_murre_WS"
var_x09      <- "X09_DFA_HakeAge5Plus"
var_x10      <- "X10_Harbour_s_2yrLead_WS"
var_x12      <- "X12_DFA_biomassEuphShelfSum"
var_x12a      <- "X12_copepodBiomass_WGoA"
var_igf       <- "X06_DFA_IGF_mu"
var_cond      <- "X07_DFA_cpue_IntSprJunHW"
var_sar       <- "X16_SAR"

preyNCC   <-"X05_DFA_abundSardine"
predNCC_long   <-"X09_DFA_HakeAge5Plus"
predNCC_short   <-"X08_commonMurre_JSOES"
preyAK <-  "X12_copepodBiomass_WGoA"
preyAK_1Blong <-  "X12_DFA_biomassEuphShelfSum"
predAK_long   <-"X10_Harbour_s_2yrLead_WS"
predAK_short   <-"X10_Harbor_seal_CR_2yrLead"
preyNCC   <-"X05_DFA_abundSardine"
predNCC_long   <-"X09_DFA_HakeAge5Plus"
predNCC_short   <-"X08_commonMurre_JSOES"
preyAK <-  "X12_copepodBiomass_WGoA"
preyAK_1Blong <-  "X12_DFA_biomassEuphShelfSum"
predAK_long   <-"X10_Harbour_s_2yrLead_WS"
predAK_short   <-"X10_Harbor_seal_CR_2yrLead"




preyNCC_long_1A   <-"X05_DFA_abundSardine"
preyNCC_short_1A   <-"X05_DFA_abundSardine"
preyNCC_long_1B   <-"X05_DFA_abundSardine"
preyNCC_short_1B   <-"X05_DFA_abundSardine"
preyNCC_long_1C   <-"X05_DFA_abundSardine"
preyNCC_short_1C   <-"X05_DFA_abundSardine"
predNCC_long_1A   <-"X09_DFA_HakeAge5Plus"
predNCC_short_1A   <-"X08_commonMurre_JSOES"
...
preyAK_long_1A   <-"X05_DFA_abundSardine"
preyAK_short_1A   <-"X05_DFA_abundSardine"
preyAK_long_1B   <-"X05_DFA_abundSardine"
preyAK_short_1B   <-"X05_DFA_abundSardine"
preyAK_long_1C   <-"X05_DFA_abundSardine"
preyAK_short_1C   <-"X05_DFA_abundSardine"
predAK_long_1A   <-"X09_DFA_HakeAge5Plus"
predAK_short_1A   <-"X08_commonMurre_JSOES"
...
growth_long_1A
...
var_cond      <- "X07_DFA_cpue_IntSprJunHW"
var_sar       <- "X16_SAR"


#for DAG1C
sem_output_summary_fxn(DAG1C_long_topmodel)
sem_output_summary_fxn(DAG1C_short_topmodel)

var_growth_short       <- "X06_StomFull_May"
var_growth_long       <- "X06_Lmu_IntSprJunH"
var_cond      <- "X07_DFA_cpue_IntSprJunHW"
var_sar       <- "X16_SAR"

preyNCC_long   <-"X04_marketsquid_GAM"
preyNCC_short   <-"X01_habCompInd"
predNCC_long   <-"X09_DFA_HakeAge5Plus"
predNCC_short   <-"X08_commonMurre_JSOES"
preyAK_long <-  "X12_copepodCom_EGoA"
preyAK_short <-  "X14_pinkSalmon"
predAK_long   <-"X10_Northern_f_s_2yrLead_WS"
predAK_short   <-"X10_Harbor_seal_CR_2yrLead"


#-------------
my_pretty_names <- c(
  "X01_habCompInd"  = "Habitat Compress",
  "X01.ZooPreyNCC_JSOES_DFA_sumPreyOfPrey_planktonJun"  = "JSOES Plankton",
  "X01.ZooPreyNCC_JSOES_Cassin_s_aukl_WS"  = "Cassin's auklet",
  "X03.FishPreyNCC_b_HakeAge1"           = "FishPrey_HakeAge1",
  "X04_marketsquid_GAM"  = "Market squid_GAM",
  "X05_DFA_abundSardine"  = "Sardine DFA",
  "X08_commonMurre_JSOES" = "Murre JSOES",
  "X08.PredBirdNCC_hump_Common_murre_WS" = "Common Murre WS",
  "X09_DFA_HakeAge5Plus"  = "Hake/Mackerel",
  "X10_Harbour_s_2yrLead_WS"  = "Harbour_seal_2yrLead_WS",
  "X10_Harbor_seal_CR_2yrLead"  = "Harbour_seal_2yrLead_CR",
  "X10_Northern_f_s_2yrLead_WS"  = "No_fur_seal_WS",
  "X12_DFA_biomassEuphShelfSum"  = "WCVI Euph/Amphipods",
  "X12_copepodBiomass_WGoA"  = "WGOA copepod Biomass",
  "X12_copepodCom_EGoA"  = "EGOA copepod Biomass",
  "X14_pinkSalmon"  = "Pink Salmon",
  
  
  "X06_Lmu_IntSprJunH"             = "Lmu_JunH",
  "X06_StomFull_May"             = "StomFull_May",
  "X06_DFA_IGF_mu"             = "Salmon Growth",
  "X07_DFA_cpue_IntSprJunHW"   = "Salmon CPUE",
  "X16_SAR"                              = "Salmon SAR"
)


#notes: mod2/4= a: JSOES_DFA_sumPreyOfPrey & PredBirdNCC_0_commonMurre_JSOES   


# All Path layouts--------- 
        path_layout_DAG1C_long <- get_layout(
              
          preyNCC_long      ,"", "",    "", predNCC_long      ,  
              ""      ,"","", var_cond      ,  ""      ,
              "" , ""      , var_growth_long     ,"", ""      ,
          preyAK_long    , "", ""      , ""      , predAK_long     , 
              # ""    , "", ""      , ""      , ""      , 
              ""    , ""   , var_sar          ,   ""   ,  ""  , 
              rows = 5
            )

        path_layout_DAG1C_short <- get_layout(
          
              preyNCC_short      ,"", "",    "", predNCC_short      ,  
              ""      ,var_growth_short,""     , var_cond ,  ""      ,
              "" , ""      , ""     ,"", ""      ,
              preyAK_short    , "", ""      , ""      , predAK_short     , 
              # ""    , "", ""      , ""      , ""      , 
              ""    , ""   , var_sar          ,   ""   ,  ""  , 
              rows = 5
            )


          path_layout_DAG1A_long <- get_layout(
            
                    preyNCC      ,"", "",    "", predNCC_long      ,  
                   ""      ,var_igf,"", var_cond      ,  ""      ,
                    "" , ""      , ""     ,"", ""      ,
                   preyAK    , "", ""      , ""      , predAK_long     , 
                  # ""    , "", ""      , ""      , ""      , 
                  ""    , ""   , var_sar          ,   ""   ,  ""  , 
                    rows = 5
                  )

          path_layout_DAG1B_long <- get_layout(
                  preyNCC      ,"", "",    "", predNCC_long      ,  
                  ""      ,var_igf,"", var_cond      ,  ""      ,
                  "" , ""      , ""     ,"", ""      ,
                  preyAK_1Blong    , "", ""      , ""      , predAK_long     , 
                  # ""    , "", ""      , ""      , ""      , 
                  ""    , ""   , var_sar          ,   ""   ,  ""  , 
                  rows = 5
                )

          path_layout_DAG1A_short <- get_layout(
                  preyNCC      ,"", "",    "", predNCC_short      ,  
                  ""      ,var_igf,"", var_cond      ,  ""      ,
                  "" , ""      , ""     ,"", ""      ,
                  preyAK    , "", ""      , ""      , predAK_short     , 
                  # ""    , "", ""      , ""      , ""      , 
                  ""    , ""   , var_sar          ,   ""   ,  ""  , 
                  rows = 5
                )
          
          path_layout_DAG1B_short <- get_layout(
            preyNCC      ,"", "",    "", predNCC_short      ,  
            ""      ,var_igf,"", var_cond      ,  ""      ,
            "" , ""      , ""     ,"", ""      ,
            preyAK    , "", ""      , ""      , predAK_short     , 
            # ""    , "", ""      , ""      , ""      , 
            ""    , ""   , var_sar          ,   ""   ,  ""  , 
            rows = 5
          )
          
          
          
          
#create plots -----------          
model="DAG1A_long"
          path_layout<-get(paste0("path_layout_",model))
          finalmodel<-get(paste(model,"topmodel",sep="_"))
          finalmodel_reduced<-get(paste(model,"topmodel","reduced",sep="_"))

          p<- Lisa_sem_graph.v2(finalmodel, 
                  layout = path_layout, 
                  title = paste(model,"topmodel",sep="_"),
                  node_labels = my_pretty_names, 
                  r2_nodes = c("X06_DFA_IGF_mu", 
                               "X07_DFA_cpue_IntSprJunHW",
                               "X16_SAR"),
                  width=10,height=6, 
                  angle = 180, rect_width=3,rect_height=0.8,
                  ellipses_width=2,ellipses_height=1.2,
                  variance_diameter=.7,text_size=4.5, curvature=10,
                  node_width = 10,  
                  node_height = 1.5,
                   save = TRUE,
                  savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel",sep="_"),".png"))
          
        p3<- plot(p)+ labs(title = paste(model, "topmodel", sep = "_"))+
          theme(plot.title = element_text(hjust = 0.5))
        

          assign(paste(model,"Dougtopmodel","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","plot",sep="_"))
          plot(DAG1A_long_Dougtopmodel_plot)

          p.reduced<- Lisa_sem_graph.v2(finalmodel_reduced, 
                                layout = path_layout, 
                                title = paste(model,"topmodel","reduced",sep="_"),
                                node_labels = my_pretty_names, 
                                r2_nodes = c("X06_DFA_IGF_mu", 
                                             "X07_DFA_cpue_IntSprJunHW",
                                             "X16_SAR"),
                                width=10,height=6, 
                                angle = 180, rect_width=3,rect_height=0.8,
                                ellipses_width=2,ellipses_height=1.2,
                                variance_diameter=.7,text_size=4.5, curvature=10,
                                node_width = 10,  
                                node_height = 1.5,
                                save = TRUE,
                                savename=paste0("LisaXP/outputs_4/",paste(model,"topmodel","reduced",sep="_"),".png"))
          
          p3<- plot(p.reduced)+ labs(title = paste(model, "topmodel","reduced", sep = "_"))+
            theme(plot.title = element_text(hjust = 0.5))
          
          
          assign(paste(model,"Dougtopmodel","reduced","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","reduced","plot",sep="_"))
          
          
          
          
          
                    
model="DAG1A_short"
          path_layout<-get(paste0("path_layout_",model))
          finalmodel<-get(paste(model,"topmodel",sep="_"))
          finalmodel_reduced<-get(paste(model,"topmodel","reduced",sep="_"))
          
          p<- Lisa_sem_graph.v2(finalmodel, 
                                layout = path_layout, 
                                title = paste(model,"topmodel",sep="_"),
                                node_labels = my_pretty_names, 
                                r2_nodes = c("X06_DFA_IGF_mu", 
                                             "X07_DFA_cpue_IntSprJunHW",
                                             "X16_SAR"),
                                width=10,height=6, 
                                angle = 180, rect_width=3,rect_height=0.8,
                                ellipses_width=2,ellipses_height=1.2,
                                variance_diameter=.7,text_size=4.5, curvature=10,
                                node_width = 10,  
                                node_height = 1.5,
                                save = TRUE,
                                savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel",sep="_"),".png"))
          
          p3<- plot(p)+ labs(title = paste(model, "topmodel", sep = "_"))+
            theme(plot.title = element_text(hjust = 0.5))
          
          assign(paste(model,"Dougtopmodel","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","plot",sep="_"))
          
          p.reduced<- Lisa_sem_graph.v2(finalmodel_reduced, 
                                        layout = path_layout, 
                                        title = paste(model,"topmodel","reduced",sep="_"),
                                        node_labels = my_pretty_names, 
                                        r2_nodes = c("X06_DFA_IGF_mu", 
                                                     "X07_DFA_cpue_IntSprJunHW",
                                                     "X16_SAR"),
                                        width=10,height=6, 
                                        angle = 180, rect_width=3,rect_height=0.8,
                                        ellipses_width=2,ellipses_height=1.2,
                                        variance_diameter=.7,text_size=4.5, curvature=10,
                                        node_width = 10,  
                                        node_height = 1.5,
                                        save = TRUE,
                                        savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel","reduced",sep="_"),".png"))
          
          p3<- plot(p.reduced)+ labs(title = paste(model, "topmodel","reduced", sep = "_"))+
            theme(plot.title = element_text(hjust = 0.5))
          
          assign(paste(model,"Dougtopmodel","reduced","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","reduced","plot",sep="_"))
          
model="DAG1B_long"
          path_layout<-get(paste0("path_layout_",model))
          finalmodel<-get(paste(model,"topmodel",sep="_"))
          finalmodel_reduced<-get(paste(model,"topmodel","reduced",sep="_"))
          
          p<- Lisa_sem_graph.v2(finalmodel, 
                                layout = path_layout, 
                                title = paste(model,"topmodel",sep="_"),
                                node_labels = my_pretty_names, 
                                r2_nodes = c("X06_DFA_IGF_mu", 
                                             "X07_DFA_cpue_IntSprJunHW",
                                             "X16_SAR"),
                                width=10,height=6, 
                                angle = 180, rect_width=3,rect_height=0.8,
                                ellipses_width=2,ellipses_height=1.2,
                                variance_diameter=.7,text_size=4.5, curvature=10,
                                node_width = 10,  
                                node_height = 1.5,
                                save = TRUE,
                                savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel",sep="_"),".png"))
          
          p3<- plot(p)+ labs(title = paste(model, "topmodel", sep = "_"))+
            theme(plot.title = element_text(hjust = 0.5))
          
          assign(paste(model,"Dougtopmodel","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","plot",sep="_"))
          
          p.reduced<- Lisa_sem_graph.v2(finalmodel_reduced, 
                                        layout = path_layout, 
                                        title = paste(model,"topmodel","reduced",sep="_"),
                                        node_labels = my_pretty_names, 
                                        r2_nodes = c("X06_DFA_IGF_mu", 
                                                     "X07_DFA_cpue_IntSprJunHW",
                                                     "X16_SAR"),
                                        width=10,height=6, 
                                        angle = 180, rect_width=3,rect_height=0.8,
                                        ellipses_width=2,ellipses_height=1.2,
                                        variance_diameter=.7,text_size=4.5, curvature=10,
                                        node_width = 10,  
                                        node_height = 1.5,
                                        save = TRUE,
                                        savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel","reduced",sep="_"),".png"))
          
          p3<- plot(p.reduced)+ labs(title = paste(model, "topmodel","reduced", sep = "_"))+
            theme(plot.title = element_text(hjust = 0.5))
          
          assign(paste(model,"Dougtopmodel","reduced","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","reduced","plot",sep="_"))
          
          
model="DAG1B_short"
          path_layout<-get(paste0("path_layout_",model))
          finalmodel<-get(paste(model,"topmodel",sep="_"))
          finalmodel_reduced<-get(paste(model,"topmodel","reduced",sep="_"))
          
          p<- Lisa_sem_graph.v2(finalmodel, 
                                layout = path_layout, 
                                title = paste(model,"topmodel",sep="_"),
                                node_labels = my_pretty_names, 
                                r2_nodes = c("X06_DFA_IGF_mu", 
                                             "X07_DFA_cpue_IntSprJunHW",
                                             "X16_SAR"),
                                width=10,height=6, 
                                angle = 180, rect_width=3,rect_height=0.8,
                                ellipses_width=2,ellipses_height=1.2,
                                variance_diameter=.7,text_size=4.5, curvature=10,
                                node_width = 10,  
                                node_height = 1.5,
                                save = TRUE,
                                savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel",sep="_"),".png"))
          p3<- plot(p)+ labs(title = paste(model, "topmodel", sep = "_"))+
            theme(plot.title = element_text(hjust = 0.5))
          
          assign(paste(model,"Dougtopmodel","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","plot",sep="_"))
          
          
          p.reduced<- Lisa_sem_graph.v2(finalmodel_reduced, 
                                        layout = path_layout, 
                                        title = paste(model,"topmodel","reduced",sep="_"),
                                        node_labels = my_pretty_names, 
                                        r2_nodes = c("X06_DFA_IGF_mu", 
                                                     "X07_DFA_cpue_IntSprJunHW",
                                                     "X16_SAR"),
                                        width=10,height=6, 
                                        angle = 180, rect_width=3,rect_height=0.8,
                                        ellipses_width=2,ellipses_height=1.2,
                                        variance_diameter=.7,text_size=4.5, curvature=10,
                                        node_width = 10,  
                                        node_height = 1.5,
                                        save = TRUE,
                                        savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel","reduced",sep="_"),".png"))
         
          p3<- plot(p.reduced)+ labs(title = paste(model, "topmodel","reduced", sep = "_"))+
            theme(plot.title = element_text(hjust = 0.5))
          
          assign(paste(model,"Dougtopmodel","reduced","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","reduced","plot",sep="_"))
          
#DAG1C==========          
          model="DAG1C_long"
          path_layout<-get(paste0("path_layout_",model))
          finalmodel<-get(paste(model,"topmodel",sep="_"))
          finalmodel_reduced<-get(paste(model,"topmodel","reduced",sep="_"))
          
          p<- Lisa_sem_graph.v2(finalmodel, 
                                layout = path_layout, 
                                title = paste(model,"topmodel",sep="_"),
                                node_labels = my_pretty_names, 
                                r2_nodes = c("X06_DFA_IGF_mu", 
                                             "X07_DFA_cpue_IntSprJunHW",
                                             "X16_SAR"),
                                width=10,height=6, 
                                angle = 180, rect_width=2.2,rect_height=0.8,
                                ellipses_width=2,ellipses_height=1.2,
                                variance_diameter=.7,text_size=4.5, curvature=10,
                                node_width = 10,  
                                node_height = 1.5,
                                save = TRUE,
                                savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel",sep="_"),".png"))
          
          p3<- plot(p)+ labs(title = paste(model, "topmodel", sep = "_"))+
            theme(plot.title = element_text(hjust = 0.5))
          
          # ggsave(filename = paste0("LisaXP/outputs_3/", paste(model, "topmodel", sep="_"), ".png"), 
          #       plot = p, width = 10, height = 6)
          
          assign(paste(model,"Dougtopmodel","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","plot",sep="_"))
          plot(DAG1A_long_Dougtopmodel_plot)
          
          p.reduced<- Lisa_sem_graph.v2(finalmodel_reduced, 
                                        layout = path_layout, 
                                        title = paste(model,"topmodel","reduced",sep="_"),
                                        node_labels = my_pretty_names, 
                                        r2_nodes = c("X06_DFA_IGF_mu", 
                                                     "X07_DFA_cpue_IntSprJunHW",
                                                     "X16_SAR"),
                                        width=10,height=6, 
                                        angle = 180, rect_width=3,rect_height=0.8,
                                        ellipses_width=2,ellipses_height=1.2,
                                        variance_diameter=.7,text_size=4.5, curvature=10,
                                        node_width = 10,  
                                        node_height = 1.5,
                                        save = TRUE,
                                        savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel","reduced",sep="_"),".png"))
          
          p3<- plot(p.reduced)
            # + labs(title = paste(model, "topmodel","reduced", sep = "_"))+
            # theme(plot.title = element_text(hjust = 0.5))
          
          
          assign(paste(model,"Dougtopmodel","reduced","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","reduced","plot",sep="_"))
          
          
          
          
          
          
          model="DAG1C_short"
          path_layout<-get(paste0("path_layout_",model))
          finalmodel<-get(paste(model,"topmodel",sep="_"))
          finalmodel_reduced<-get(paste(model,"topmodel","reduced",sep="_"))
          
          p<- Lisa_sem_graph.v2(finalmodel, 
                                layout = path_layout, 
                                title = paste(model,"topmodel",sep="_"),
                                node_labels = my_pretty_names, 
                                r2_nodes = c("X06_DFA_IGF_mu", 
                                             "X07_DFA_cpue_IntSprJunHW",
                                             "X16_SAR"),
                                width=10,height=6, 
                                angle = 180, rect_width=2,rect_height=0.8,
                                ellipses_width=2,ellipses_height=1.2,
                                variance_diameter=.7,text_size=4.5, curvature=10,
                                node_width = 10,  
                                node_height = 1.5,
                                save = TRUE,
                                savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel",sep="_"),".png"))
          
          p3<- plot(p)
          # + labs(title = paste(model, "topmodel", sep = "_"))+
          #   theme(plot.title = element_text(hjust = 0.5))
          
          assign(paste(model,"Dougtopmodel","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","plot",sep="_"))
          
          p.reduced<- Lisa_sem_graph.v2(finalmodel_reduced, 
                                        layout = path_layout, 
                                        title = paste(model,"topmodel","reduced",sep="_"),
                                        node_labels = my_pretty_names, 
                                        r2_nodes = c("X06_DFA_IGF_mu", 
                                                     "X07_DFA_cpue_IntSprJunHW",
                                                     "X16_SAR"),
                                        width=10,height=6, 
                                        angle = 180, rect_width=1.5,rect_height=0.8,
                                        ellipses_width=2,ellipses_height=1.2,
                                        variance_diameter=.7,text_size=4.5, curvature=10,
                                        node_width = 10,  
                                        node_height = 1.5,
                                        save = TRUE,
                                        savename=paste0("LisaXP/outputs_3/",paste(model,"topmodel","reduced",sep="_"),".png"))
          
          p3<- plot(p.reduced)
          # + labs(title = paste(model, "topmodel","reduced", sep = "_"))+
          #   theme(plot.title = element_text(hjust = 0.5))
          
          assign(paste(model,"Dougtopmodel","reduced","plot",sep="_"),p3,envir = .GlobalEnv)
          print(paste(model,"Dougtopmodel","reduced","plot",sep="_"))
          
          
#make plots
          plot(DAG1A_long_Dougtopmodel_plot) 
          
          DAG1C_short_Dougtopmodel_reduced_plot 
