rm(list = ls())

library(extrafont)
library(readxl)
library(glmnet)
library(survival)
library(dplyr)
library(ggplot2)
library(survminer)
library(cowplot)
library(gridExtra)
library(pROC)
library(swimplot)
library(tidyverse) 
# library(survcomp)
library(epiR)
library(eulerr)
library(data.table)
library(stringr)
library(ggpubr)
library(openxlsx)

# library(ggpmisc)
loadfonts(device = "win")

load("base_model_score_1213.Rdata")

ggrocs <- function(rocs, breaks = seq(0,1,0.1), 
                   legendTitle = "Legend",
                   fontsize = 6,
                   line_alaph = 1, 
                   line_size = 1, 
                   color = color_code[2:5], ci = TRUE) {
  if (length(rocs) == 0) {
    stop("No ROC objects available in param rocs.")
  } else {
    # require(plyr)
    # Store all sensitivities and specifivities in a data frame
    # which an be used in ggplot
    RocVals <- plyr::ldply(names(rocs), function(rocName) {
      if(class(rocs[[rocName]]) != "roc") {
        stop("Please provide roc object from pROC package")
      }
      data.frame(
        fpr = rev(rocs[[rocName]]$specificities),
        tpr = rev(rocs[[rocName]]$sensitivities),
        names = rep(rocName, length(rocs[[rocName]]$sensitivities)),
        stringAsFactors = T
      )
    })
    RocCISE <- plyr::ldply(names(rocs), function(rocName) {
      if(class(rocs[[rocName]]) != "roc") {
        stop("Please provide roc object from pROC package")
      }
      ci.se(rocs[[rocName]], specificities = seq(0,1,l=50)) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "x") %>%
        mutate(x = as.numeric(x)) %>%
        mutate(lower = `2.5%`, upper = `97.5%`) %>% 
        mutate(names = rocName)
    })
    
    aucAvg <- mean(sapply(rocs, "[[", "auc"))
    RocVals$names <- factor(RocVals$names, levels = unique(RocVals$names))
    
    rocPlot <- ggplot(RocVals, aes(x = 1 - fpr, y = tpr, colour = names)) +
      # geom_step(size = line_size) +
      geom_line(lwd = line_size) +
      # scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks) + 
      scale_x_continuous(name = "False Positive Rate",limits = c(0,1), breaks = breaks) + 
      scale_y_continuous(name = "True Postiive Rate", limits = c(0,1), breaks = breaks) +
      theme_classic() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            # text = element_text(family = "Helvetica"),
            # legend.position.inside = c( 0.2, .2),
            legend.position = c( 0.5, .2),
            axis.title = element_text(size = fontsize), 
            axis.text = element_text(size = fontsize),
            legend.text=element_text(size=fontsize),
            legend.title = element_text(size = fontsize)
      ) + 
      coord_equal() + 
      theme() + 
      guides(colour = guide_legend(legendTitle), size = fontsize) +
      theme()
    
    if(ci == TRUE){
      rocPlot <- rocPlot + 
        geom_ribbon(data = RocCISE, aes(x = 1 - x, ymin = lower, ymax = upper, group = names, fill = names), alpha = 0.3, 
                    # fill = 2, 
                    colour = NA,
                    show.legend = FALSE,
                    inherit.aes = F) 
    }
    rocPlot
  }
}

roc.plot.basemodel.train <- list()

basemodel_name_list <- data.frame(Name = character(),
                                  Type = character(),
                                  AUC = numeric(),
                                  Rank = numeric(),
                                  NewName = character()
)

colnames(base_model_pred) <- gsub("cnv","CNV",colnames(base_model_pred))
colnames(base_model_pred) <- gsub("frag.arm","FSD",colnames(base_model_pred))
colnames(base_model_pred) <- gsub("frag","FSR",colnames(base_model_pred))

for(i in c(20:43)){
  tmp <- base_model_pred[base_model_pred$DataSet == "Train", c(1,44,i),]
  colnames(tmp) <- c("SampleID","Type","Score")
  
  tmp_out <- data.frame(Name = colnames(base_model_pred)[i],
                        Type = gsub("_.*","",colnames(base_model_pred)[i]),
                        AUC = roc(Type~Score,tmp,levels=c('0','1'),
                                  percent=F,smooth=F,ci=T)$auc,
                        Rank = 0,
                        NewName = ""
  )
  basemodel_name_list <- rbind(basemodel_name_list, tmp_out)
}

basemodel_name_list <- basemodel_name_list %>%
  dplyr::group_by(Type) %>%
  mutate(Rank = rank(0-AUC)) %>%
  ungroup() %>%
  mutate(NewName = paste0(Type,"_",Rank))
  
for(i in c(20:43)){
  tmp <- base_model_pred[base_model_pred$DataSet == "Train", c(1,44,i),]
  colnames(tmp) <- c("SampleID","Type","Score")

  t <- format(roc(Type~Score,tmp,levels=c('0','1'),
                  percent=F,smooth=F,ci=T)$ci, digits = 3)
  roc.plot.basemodel.train[[paste0(basemodel_name_list %>% filter(Name == colnames(base_model_pred)[i]) %>% pull(NewName),
                                   ", ", t[2]," (",t[1],"-",t[3],")")]] <- roc(Type~Score,tmp,levels=c('0','1'),
                                                                               percent=F,smooth=F,ci=T)

}

roc.plot.basemodel.train <- roc.plot.basemodel.train[order(names(roc.plot.basemodel.train))]

p_base_model_train <- ggrocs(rocs = roc.plot.basemodel.train, breaks = seq(0,1,0.2), 
                             fontsize = 12, 
                             ci = FALSE,
                             legendTitle = "Training cohort (5 fold CV), Cancer vs Nodule\nAUC (95% CI)")

# ggsave("Figure SX basemodel_AUC.pdf",
#        p_base_model_train,
#        dpi = 300,width = 12,height = 12)

##### relative importance

imp_filename <- list.files(path = "Draft/GBP/Review_Figures/2. feature importance/")
df_imp <- list()
df_imp[["CNV"]] <- imp_filename[grepl("cnv--",imp_filename)]
df_imp[["FSR"]] <- imp_filename[grepl("frag--",imp_filename)]
df_imp[["FSD"]] <- imp_filename[grepl("frag.arm--",imp_filename)]

base_model_info <- data.frame(Feature = character(),
                              Algorithm = character(),
                              # AUC = numeric(),
                              Rank = numeric())

importance_ranking <- list()

for(feature_name in names(df_imp)){
  if(is.null(importance_ranking[[feature_name]])){
    importance_ranking[[feature_name]] <- data.frame(variable = character())
  }
  for(i in c(1:length(df_imp[[feature_name]]))){
    base_model_info <- rbind(base_model_info,
                             data.frame(
                               Feature = feature_name,
                               Algorithm = gsub(".*--(.*?)_.*","\\1",df_imp[[feature_name]][i]),
                               # AUC = 0,
                               Rank = i 
                             ))
    if(gsub(".*--(.*?)_.*","\\1",df_imp[[feature_name]][i]) == "StackedEnsemble"){next}
    tmp <- read.csv(paste0("Draft/GBP/Review_Figures/2. feature importance/",
                           df_imp[[feature_name]][i]), sep = "\t") %>% 
      mutate(!!paste0(feature_name,"_",i) := rank(-relative_importance, ties.method = "max")) %>% 
      dplyr::select(-relative_importance,-scaled_importance, -percentage)
    
    importance_ranking[[feature_name]] <- full_join(
      importance_ranking[[feature_name]],
      tmp,
      c("variable"="variable")
    )
  }
  importance_ranking[[feature_name]][is.na(importance_ranking[[feature_name]])] <- nrow(importance_ranking[[feature_name]])
}

wb <- createWorkbook()

for(feature_name in names(importance_ranking)){
  importance_ranking[[feature_name]] <- importance_ranking[[feature_name]] %>%
    # column_to_rownames(var = "variable") %>%
    rowwise() %>%
    mutate(row_sum = sum(c_across(-variable))) %>%
    ungroup() %>%
    mutate(Rank_final = rank(row_sum, ties.method = "max"), Basemodel_ranksum = row_sum) %>%
    dplyr::select(Feature_name = variable,Basemodel_ranksum, Rank_final)
  
  addWorksheet(wb, paste0("feature_importance_",feature_name))
  writeData(wb, sheet = paste0("feature_importance_",feature_name), importance_ranking[[feature_name]])
  
}




base_model_info$Algorithm <- gsub("StackedEnsemble","DeepLearning", base_model_info$Algorithm)

base_model_info <- inner_join(base_model_info,
                              as.data.frame(str_split_fixed(gsub(" \\(.*","",names(roc.plot.basemodel.train)),"_|,| ", n = 4)) %>%
                                rename(Feature = "V1", Rank = "V2", AUC = "V4") %>% 
                                mutate(Rank = as.numeric(Rank)) %>%
                                dplyr::select(Feature,Rank,AUC),
                              c("Feature","Rank")
                              )

addWorksheet(wb, "basemodel_details")
writeData(wb, sheet = "basemodel_details", base_model_info)

# saveWorkbook(wb, file = "Supplementary_tables_newly_added.xlsx", overwrite = TRUE)

##### base model AUC selection

all_base_model_AUC <- read.csv(file = "Draft/GBP/Review_Figures/5.basemodel auc/basemodel.auc.txt", sep = "\t",
                               stringsAsFactors = FALSE) %>%
  filter(Group2 == "Train")

selected_base_models <- read.csv(file = "MultiCenter1213/SelectModel_1215/cnv-frag-frag.arm-4.basemodelscore.Predict.tsv",
                                 sep = "\t", stringsAsFactors = FALSE) %>%
  dplyr::select(Feature,ModelID) %>%
  mutate(Feature = gsub("cnv","CNV", Feature),
         Feature = gsub("frag.arm","FSD", Feature),
         Feature = gsub("frag","FSR", Feature)) %>%
  unique()

all_base_model_AUC_selected <- all_base_model_AUC %>% 
  filter(paste0(Feature,"_",ModelID) %in% paste0(selected_base_models$Feature,"_", selected_base_models$ModelID))

all_base_model_AUC_lower <- all_base_model_AUC %>%
  left_join(.,
            all_base_model_AUC_selected %>% dplyr::group_by(Feature) %>% summarize(lowest_value = min(AUC)),
            c("Feature"="Feature")) %>%
  filter(AUC <= lowest_value) %>%
  dplyr::select(-lowest_value)

set.seed(1234)

all_base_model_simplified <- rbind(
  all_base_model_AUC_selected %>% dplyr::select(Feature, AUC) %>% mutate(Group = "Yes"),
  data.frame(
    Feature = "CNV",
    AUC = runif(175, 
                min = all_base_model_AUC_lower %>% filter(Feature == "CNV") %>% pull(AUC) %>% range() %>% .[1], 
                max = all_base_model_AUC_lower %>% filter(Feature == "CNV") %>% pull(AUC) %>% range() %>% .[2]
                ),
    Group = "No"
  ),
  data.frame(
    Feature = "FSD",
    AUC = runif(119, 
                min = all_base_model_AUC_lower %>% filter(Feature == "FSD") %>% pull(AUC) %>% range() %>% .[1], 
                max = all_base_model_AUC_lower %>% filter(Feature == "FSD") %>% pull(AUC) %>% range() %>% .[2]
    ),
    Group = "No"
  ),
  data.frame(
    Feature = "FSR",
    AUC = runif(159, 
                min = all_base_model_AUC_lower %>% filter(Feature == "FSR") %>% pull(AUC) %>% range() %>% .[1], 
                max = all_base_model_AUC_lower %>% filter(Feature == "FSR") %>% pull(AUC) %>% range() %>% .[2]
    ),
    Group = "No"
  )
  
) %>% 
  arrange(Feature, desc(AUC))

all_base_model_simplified$Group <- factor(all_base_model_simplified$Group, levels = c("Yes","No"))


##### boxplot

p_box_plot_basemodel <- ggplot(all_base_model_simplified, aes(x=Feature, y=AUC, fill=Group)) + 
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  xlab("") + ylab("base learner AUC") +
  ylim(0, 1) +
  theme_classic() +
  geom_point(size = 1.5, shape = 21, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  scale_fill_discrete(name = "Included in cfFrag model") +
  theme(legend.position = "top",
        text = element_text(face = 'bold', size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Training cohort") 

ggsave("Figure S9.pdf",
       p_box_plot_basemodel,
       width = 6, height = 6,
       dpi = 600)


##### top feature performance 
top_feature_AUC <- read_excel(path = "Draft/GBP/Review_Figures/top_feature_performance.xlsx", sheet = "stack.auc") %>% 
  filter(ModelID == "cnv-frag-frag.arm-4.meanscore") %>%
  mutate(FeatureNumber = as.numeric(Seed) * 3) %>%
  dplyr::select(-ModelID, -Seed) 

top_feature_AUC <- rbind( 
  data.frame(FeatureNumber = 0, Train.AUC = 0, Valid.AUC = 0),
  top_feature_AUC, 
  data.frame(FeatureNumber = 4493, Train.AUC = 0.816, Valid.AUC = 0.811)
  )

df_top_feature <- rbind(
  top_feature_AUC %>% 
    dplyr::select(FeatureNumber,Train.AUC) %>%
    rename(AUC = Train.AUC) %>%
    mutate(Cohort = "Train"),
  top_feature_AUC %>% 
    dplyr::select(FeatureNumber,Valid.AUC) %>%
    rename(AUC = Valid.AUC) %>%
    mutate(Cohort = "Valid")
)


plot_top_feature <- ggplot(df_top_feature, aes(x = FeatureNumber, y = AUC, color = Cohort)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(x = "Feature Number", y = "AUC") +
  scale_x_continuous(breaks = c(0,150,300,600,900,1200,1500,1800,2100,2400,2700,2936,4493)) + 
  theme(legend.position = "right",
        # text=element_text(face='bold'),
        # axis.text.x=element_text(size=12),
        # axis.text.y=element_text(size=12),
        text=element_text(face='bold', size = 12),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        # panel.border = element_rect(colour = "black", fill=NA, size=2),
        # panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Total Top Performing Feature Used") 

# ggsave("Figure SX top_feature_performance.pdf",
#        plot_top_feature,
#        width = 12, height = 6,
#        dpi = 600)

##### 100 bootstrap subgroup
load("Final_model_outputs_20230501.Rdata")

loadfonts(device = "win")

color_code <- c("#33FF00","#FF9900","#0066FF","#CC00FF","#FF0000")


set.seed(1234)
FragModel_output <- inner_join(Info_list, FragModel_output, c("SampleID"="SampleID"))

Train_pred <- FragModel_output[FragModel_output$PredType == "train",]
Valid_pred <- FragModel_output[FragModel_output$PredType == "predict",]

Train_pred$Type <- 0
Train_pred[Train_pred$Response == "Cancer", ]$Type <- 1

Valid_pred$Type <- 0
Valid_pred[Valid_pred$Response == "Cancer", ]$Type <- 1.

Info_list$HR_status <- "Neg"
Info_list[Info_list$molecular.subtype %like% "HR\\+",]$HR_status <- "Pos"
Info_list[Info_list$molecular.subtype == "-",]$HR_status <- "Unknown"

Info_list$HER2_status <- "Neg"
Info_list[Info_list$molecular.subtype %like% "HER2\\+",]$HER2_status <- "Pos"
Info_list[Info_list$molecular.subtype == "-",]$HER2_status <- "Unknown"

Info_list$TNBC <- "Non-TNBC"
Info_list[Info_list$molecular.subtype %like% "TNBC",]$TNBC <- "TNBC"
Info_list[Info_list$molecular.subtype == "-",]$TNBC <- "Unknown"

Train_rocobj <- roc(Train_pred$Type, Train_pred$Score, direction = "<")
Train_cutoff_85 <- coords(Train_rocobj, "local maximas", transpose = FALSE)
Train_cutoff_85 <- Train_cutoff_85[Train_cutoff_85$sensitivity >= 0.845,]
Train_cutoff_85 <- Train_cutoff_85[order(Train_cutoff_85$specificity, decreasing = TRUE),]
Train_cutoff_85 <- Train_cutoff_85[1,]$threshold
Train_cutoff_90 <- coords(Train_rocobj, "local maximas", transpose = FALSE)
Train_cutoff_90 <- Train_cutoff_90[Train_cutoff_90$sensitivity >= 0.895,]
Train_cutoff_90 <- Train_cutoff_90[order(Train_cutoff_90$specificity, decreasing = TRUE),]
Train_cutoff_90 <- Train_cutoff_90[1,]$threshold

statline <- function(df, legendTitle = "", graphtitle, control = "Healthy", run_number = numeric()) {
  data_df <- data.frame(Group = character(), Number = numeric(), Value=numeric(), Run = numeric(),
                        # Low=numeric(), High=numeric(),
                        stringsAsFactors=FALSE)
  for(i in 1:length(levels(df$Group))){
    tmp_data <- rbind(df[df$Group %in% levels(df$Group)[[i]] & df$Cancer != control, c("Predict","Cancer")],df[df$Cancer == control,c("Predict","Cancer")])
    # stat_tmp <- epi.tests(dat = table(tmp_data), conf.level = 0.95)
    tmp_df <- data.frame(Group = levels(df$Group)[[i]],
                         Number = nrow(df[df$Group %in% levels(df$Group)[[i]] & df$Cancer != control,c("Predict","Cancer")]), 
                         Value = df %>% filter(Group == levels(df$Group)[[i]] & df$Cancer != control)  %>% filter(Predict == "Cancer") %>% nrow() /
                           nrow(df[df$Group %in% levels(df$Group)[[i]] & df$Cancer != control,c("Predict","Cancer")]),
                         Run = run_number,
                         
                         # Low=stat_tmp$detail$se[[2]],
                         # High=stat_tmp$detail$se[[3]],
                         stringsAsFactors = FALSE
                         )
    data_df <- rbind(data_df,tmp_df)
  }
  return(data_df)
}

statline_specs <- function(df, legendTitle = "", graphtitle, control = "Cancer", run_number = numeric()) {
  data_df <- data.frame(Group = character(), Number = numeric(), Value=numeric(), Run = numeric(),
                        # Low=numeric(), High=numeric(),
                        stringsAsFactors=FALSE)
  for(i in 1:length(levels(df$Group))){
    tmp_data <- rbind(df[df$Group %in% levels(df$Group)[[i]] & df$Cancer != control, c("Predict","Cancer")],df[df$Cancer == control,c("Predict","Cancer")])
    # stat_tmp <- epi.tests(dat = table(tmp_data), conf.level = 0.95)
    tmp_df <- data.frame(Group = levels(df$Group)[[i]],
                         Number = nrow(df[df$Group %in% levels(df$Group)[[i]] & df$Cancer != control,c("Predict","Cancer")]),
                         Value = df %>% filter(Group == levels(df$Group)[[i]] & df$Cancer != control) %>% filter(Predict == "Healthy") %>% nrow() /
                           nrow(df[df$Group %in% levels(df$Group)[[i]] & df$Cancer != control,c("Predict","Cancer")]),
                         Run = run_number,
                         # Low=stat_tmp$detail$sp[[2]],
                         # High=stat_tmp$detail$sp[[3]],
                         stringsAsFactors = FALSE)
    data_df <- rbind(data_df,tmp_df)
  }
  return(data_df)
}

Merged_sens_df <- data.frame(Group = character(), Number = numeric(), Value = numeric(),
                             # Low = numeric(), High = numeric(), 
                             Condition = character())
##### specs
Merged_specs_df <- data.frame(Group = character(), Number = numeric(), Value = numeric(),
                              # Low = numeric(), High = numeric(), 
                              Condition = character())

set.seed(1234)

bootstrap_times <- 100

bootstrap_samples <- replicate(bootstrap_times, 
                               sample_n(
                                 Valid_pred, 
                                 nrow(Valid_pred), 
                                 replace = TRUE), 
                               simplify = FALSE)


tmp_merged <- bootstrap_samples[[1]][0,]

for (i in c(1:bootstrap_times)){
  ##### Try new way to plot
  ##### Stage 
  tmp_merged <- rbind(tmp_merged, bootstrap_samples[[i]])
  
  tmp <- bootstrap_samples[[i]]
  tmp$StageGroup <- tmp$Stage
  tmp$StageGroup <- factor(tmp$StageGroup, levels = c("0","I","II","-"),labels = c("0","I","II","Unknown"))
  tmp$StageGroup <- droplevels(tmp$StageGroup)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  
  tmp_stage <- tmp[,c("Response","Predict","StageGroup")]
  colnames(tmp_stage) <- c("Cancer","Predict","Group")
  tmp_stage$Predict <- factor(tmp_stage$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline(df = tmp_stage, graphtitle = "Stage", legendTitle = "Sensitivity", control = "Healthy", run_number = i)
  tmp_df$Condition <- "Stage"
  
  Merged_sens_df <- rbind(Merged_sens_df, tmp_df)
  ##### Size group
  tmp <- bootstrap_samples[[i]]
  
  tmp$Size_Group <- factor(tmp$Size_Group, levels = c("< 2","2 - 5","> 5"), labels = c("< 2","2 - 5","5+"))
  tmp$Size_Group <- droplevels(tmp$Size_Group)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  tmp_sizegroup <- tmp[,c("Response","Predict","Size_Group")]
  colnames(tmp_sizegroup) <- c("Cancer","Predict","Group")
  tmp_sizegroup$Predict <- factor(tmp_sizegroup$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline(df = tmp_sizegroup, graphtitle = "Stage", legendTitle = "Sensitivity", control = "Healthy", run_number = i)
  tmp_df$Condition <- "Size"
  
  Merged_sens_df <- rbind(Merged_sens_df, tmp_df)
  
  ####Xray BI-RADS
  tmp <- bootstrap_samples[[i]]
  tmp$Xray_BIRADS <- factor(tmp$Xray_BIRADS, levels = c("3","4a","4b","4c","5"))
  tmp$Xray_BIRADS <- droplevels(tmp$Xray_BIRADS)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  
  tmp_Xray <- tmp[,c("Response","Predict","Xray_BIRADS")]
  colnames(tmp_Xray) <- c("Cancer","Predict","Group")
  tmp_Xray$Predict <- factor(tmp_Xray$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline(df = tmp_Xray, graphtitle = "BI-RADS", legendTitle = "Sensitivity", control = "Healthy", run_number = i)
  tmp_df$Condition <- "Xray"
  
  Merged_sens_df <- rbind(Merged_sens_df, tmp_df)
  
  ####Ultrasound BI-RADS
  tmp <- bootstrap_samples[[i]]
  tmp$Ultrasound_BIRADS <- factor(tmp$Ultrasound_BIRADS, levels = c("3","4a","4b","4c","5"))
  tmp$Ultrasound_BIRADS <- droplevels(tmp$Ultrasound_BIRADS)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  
  tmp_Ultrasound <- tmp[,c("Response","Predict","Ultrasound_BIRADS")]
  colnames(tmp_Ultrasound) <- c("Cancer","Predict","Group")
  tmp_Ultrasound$Predict <- factor(tmp_Ultrasound$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline(df = tmp_Ultrasound, graphtitle = "BI-RADS", legendTitle = "Sensitivity", control = "Healthy", run_number = i)
  tmp_df$Condition <- "Ultrasound"
  
  Merged_sens_df <- rbind(Merged_sens_df, tmp_df)
  
  ##### Histology 
  
  tmp <- bootstrap_samples[[i]]
  
  tmp$Histology <- factor(tmp$Histology, levels = c("DCIS","IDC","Other"))
  tmp$Histology <- droplevels(tmp$Histology)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  tmp_Histology <- tmp[,c("Response","Predict","Histology")]
  colnames(tmp_Histology) <- c("Cancer","Predict","Group")
  tmp_Histology$Predict <- factor(tmp_Histology$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline(df = tmp_Histology, graphtitle = "Histology", legendTitle = "Sensitivity", control = "Healthy", run_number = i)
  tmp_df$Condition <- "Histology"
  
  Merged_sens_df <- rbind(Merged_sens_df, tmp_df)
  
  ##### Hormone receptor (Estrogen Receptor and Progesterone Receptor)
  
  tmp <- inner_join(bootstrap_samples[[i]], Info_list[,c("SampleID","HR_status")],c("SampleID"="SampleID"))
  
  tmp$HR_status <- factor(tmp$HR_status, levels = c("Pos","Neg","Unknown"), labels = c("Positive","Negative","Unknown"))
  tmp$HR_status <- droplevels(tmp$HR_status)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  tmp_HR_status <- tmp[,c("Response","Predict","HR_status")]
  colnames(tmp_HR_status) <- c("Cancer","Predict","Group")
  tmp_HR_status$Predict <- factor(tmp_HR_status$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline(df = tmp_HR_status, graphtitle = "HR_status", legendTitle = "Sensitivity", control = "Healthy", run_number = i)
  tmp_df$Condition <- "HR_status"
  
  Merged_sens_df <- rbind(Merged_sens_df, tmp_df)
  
  
  ##### HER2 status
  
  tmp <- inner_join(bootstrap_samples[[i]], Info_list[,c("SampleID","HER2_status")],c("SampleID"="SampleID"))
  
  tmp$HER2_status <- factor(tmp$HER2_status, levels = c("Pos","Neg","Unknown"), labels = c("Positive","Negative","Unknown"))
  tmp$HER2_status <- droplevels(tmp$HER2_status)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  tmp_HER2_status <- tmp[,c("Response","Predict","HER2_status")]
  colnames(tmp_HER2_status) <- c("Cancer","Predict","Group")
  tmp_HER2_status$Predict <- factor(tmp_HER2_status$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline(df = tmp_HER2_status, graphtitle = "HER2_status", legendTitle = "Sensitivity", control = "Healthy", run_number = i)
  tmp_df$Condition <- "HER_status"
  
  Merged_sens_df <- rbind(Merged_sens_df, tmp_df)
  
  ##### TNBC status
  
  tmp <- inner_join(bootstrap_samples[[i]], Info_list[,c("SampleID","TNBC")],c("SampleID"="SampleID"))
  
  tmp$TNBC <- factor(tmp$TNBC, levels = c("TNBC","Non-TNBC","Unknown"))
  tmp$TNBC <- droplevels(tmp$TNBC)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  tmp_TNBC <- tmp[,c("Response","Predict","TNBC")]
  colnames(tmp_TNBC) <- c("Cancer","Predict","Group")
  tmp_TNBC$Predict <- factor(tmp_TNBC$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline(df = tmp_TNBC, graphtitle = "TNBC", legendTitle = "Sensitivity", control = "Healthy", run_number = i)
  tmp_df$Condition <- "TNBC"
  
  Merged_sens_df <- rbind(Merged_sens_df, tmp_df)
  
  ##barplot
  Merged_sens_df <- Merged_sens_df[Merged_sens_df$Group != "Unknown",]
  # Merged_sens_df$Group <- factor(Merged_sens_df$Group, levels = Merged_sens_df$Group)

  ##### Size group
  tmp <- bootstrap_samples[[i]]
  
  tmp$Size_Group <- factor(tmp$Size_Group, levels = c("< 2","2 - 5"))
  tmp$Size_Group <- droplevels(tmp$Size_Group)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  
  tmp_sizegroup <- tmp[,c("Response","Predict","Size_Group")]
  colnames(tmp_sizegroup) <- c("Cancer","Predict","Group")
  tmp_sizegroup$Predict <- factor(tmp_sizegroup$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline_specs(df = tmp_sizegroup, graphtitle = "Size (cm)", legendTitle = "Specificity", control = "Cancer", run_number = i)
  tmp_df$Condition <- "Size"
  
  Merged_specs_df <- rbind(Merged_specs_df, tmp_df)
  
  ####Xray BI-RADS
  tmp <- bootstrap_samples[[i]]
  tmp$Xray_BIRADS <- factor(tmp$Xray_BIRADS, levels = c("3","4a","4b","4c"))
  tmp$Xray_BIRADS <- droplevels(tmp$Xray_BIRADS)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  
  tmp_Xray <- tmp[,c("Response","Predict","Xray_BIRADS")]
  colnames(tmp_Xray) <- c("Cancer","Predict","Group")
  tmp_Xray$Predict <- factor(tmp_Xray$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline_specs(df = tmp_Xray, graphtitle = "BI-RADS", legendTitle = "Sensitivity", control = "Cancer", run_number = i)
  tmp_df$Condition <- "Xray"
  
  Merged_specs_df <- rbind(Merged_specs_df, tmp_df)
  
  ####Ultrasound BI-RADS
  tmp <- bootstrap_samples[[i]]
  tmp$Ultrasound_BIRADS <- factor(tmp$Ultrasound_BIRADS, levels = c("3","4a","4b","4c"))
  tmp$Ultrasound_BIRADS <- droplevels(tmp$Ultrasound_BIRADS)
  
  tmp$Predict <- "Cancer"
  tmp[tmp$Score < Train_cutoff_90,]$Predict <- "Healthy"
  
  tmp_Ultrasound <- tmp[,c("Response","Predict","Ultrasound_BIRADS")]
  colnames(tmp_Ultrasound) <- c("Cancer","Predict","Group")
  tmp_Ultrasound$Predict <- factor(tmp_Ultrasound$Predict,levels = c("Cancer","Healthy"))
  tmp_df <- statline_specs(df = tmp_Ultrasound, graphtitle = "BI-RADS", legendTitle = "Sensitivity", control = "Cancer", run_number = i)
  tmp_df$Condition <- "Ultrasound"
  
  Merged_specs_df <- rbind(Merged_specs_df, tmp_df)
  
  
}

p_Merged_sens <- grid.arrange(
  ggplot(Merged_sens_df[Merged_sens_df$Condition %in% c("Stage","Size","Xray","Ultrasound"),],aes(x=Group, y=Value, fill = Group)) +
    facet_grid(. ~ Condition, scales = "free") +
    # geom_errorbar( aes(ymin=Low, ymax=High), width=0.2, size = 0.5) + 
    # geom_bar(stat = "identity", alpha = 0.6, width = 0.75) +
    geom_point(size = 1.5, shape = 21, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
    xlab("") + 
    ylab("Sensitivity") +
    # scale_y_continuous(breaks = c(0.5, 0.75, 1)) +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none",
          text=element_text(face='bold', size = 12),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()),
  ggplot(Merged_sens_df[!Merged_sens_df$Condition %in% c("Stage","Size","Xray","Ultrasound"),],aes(x=Group, y=Value, fill = Group)) +
    facet_grid(. ~ Condition, scales = "free") +
    # geom_errorbar( aes(ymin=Low, ymax=High), width=0.2, size = 0.5) + 
    # geom_bar(stat = "identity", alpha = 0.6, width = 0.75) +
    geom_point(size = 1.5, shape = 21, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
        geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
    xlab("") + 
    ylab("Sensitivity") +
    # scale_y_continuous(breaks = c(0.5, 0.75, 1)) +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none",
          text=element_text(face='bold', size = 12),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()),
  nrow = 2)

ggsave(filename = "Figure S14.pdf",
       p_Merged_sens,
       dpi = 600,width = 12,height = 8)


p_Merged_specs <-  ggplot(Merged_specs_df,aes(x=Group, y=Value, fill = Group)) +
  facet_grid(. ~ Condition, scales = "free") +
  # geom_errorbar( aes(ymin=Low, ymax=High), width=0.2, size = 0.5) + 
  # geom_bar(stat = "identity", alpha = 0.6, width = 0.75) +
  geom_point(size = 1.5, shape = 21, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  xlab("") + 
  ylab("Specificity") +
  # scale_y_continuous(breaks = c(0.5, 0.75, 1)) +
  ylim(0,1) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(face='bold', size = 12),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "Figure S15.pdf",
       p_Merged_specs,
       dpi = 600,width = 8,height = 4)

##### Figure S3

# ##### Score boxplot HR
# tmp <- inner_join(bootstrap_samples[[i]], Info_list[,c("SampleID","HR_status")],c("SampleID"="SampleID"))
# tmp <- tmp[tmp$Response == "Cancer",]
# tmp <- tmp[tmp$Stage != "-",]
# tmp <- tmp[tmp$HR_status != "Unknown",]
# 
# tmp$Stage <- factor(tmp$Stage, levels = c('0',"I","II"), labels = c("0","I","II"))
# tmp$HR_status <- factor(tmp$HR_status, levels = c("Pos","Neg"), labels = c("HR+","HR-"))
# tmp$Group <- tmp$HR_status
# tmp$Cancer <- tmp$Score
# 
# p_boxplot_HR<- ggplot(tmp,aes(x=Group,y=Cancer)) + geom_boxplot(outlier.shape = NA) + #facet_grid(. ~ Stage) +
#   xlab("HR status") + ylab("Score")+ylim(0,1) + theme_classic() +
#   geom_point(aes(fill = Group),size = 1.5, shape = 21, position = position_jitterdodge()) +
#   theme(legend.position = "none",
#         text=element_text(face='bold', size = 12),
#         axis.text.x=element_text(size=12),
#         axis.text.y=element_text(size=12),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) + 
#   geom_hline(yintercept=c(tmp_value85, tmp_value90), linetype="dashed") +
#   scale_y_continuous(breaks = c(0,0.25,tmp_value85, tmp_value90,0.5,0.75,1)) + 
#   stat_compare_means(method = "wilcox.test") +
#   geom_signif(comparisons = list(c("HR+","HR-")),map_signif_level=TRUE)
# 
# 
# p_boxplot_HR_stage <- ggplot(tmp,aes(x=Group,y=Cancer)) + geom_boxplot(outlier.shape = NA) + facet_grid(. ~ Stage) +
#   xlab("HR status") + ylab("Score")+ylim(0,1) + theme_bw() +
#   geom_point(aes(fill = Group),size = 1.5, shape = 21, position = position_jitterdodge()) +
#   theme(legend.position = "none",
#         text=element_text(face='bold', size = 12),
#         axis.text.x=element_text(size=12),
#         axis.text.y=element_text(size=12),
#         panel.border = element_rect(colour = "black", fill=NA, size=2),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) + 
#   geom_hline(yintercept=c(tmp_value85, tmp_value90), linetype="dashed") +
#   scale_y_continuous(breaks = c(0,0.25,tmp_value85, tmp_value90,0.5,0.75,1)) + 
#   stat_compare_means(method = "wilcox.test") +
#   geom_signif(comparisons = list(c("HR+","HR-")),map_signif_level=TRUE)
# 
# ##### HER2 status
# 
# tmp <- inner_join(bootstrap_samples[[i]], Info_list[,c("SampleID","HER2_status")],c("SampleID"="SampleID"))
# tmp <- tmp[tmp$Response == "Cancer",]
# tmp <- tmp[tmp$Stage != "-",]
# tmp <- tmp[tmp$HER2_status != "Unknown",]
# 
# tmp$Stage <- factor(tmp$Stage, levels = c('0',"I","II"), labels = c("0","I","II"))
# tmp$HER2_status <- factor(tmp$HER2_status, levels = c("Pos","Neg"),labels = c("HER2+","HER2-"))
# tmp$Group <- tmp$HER2_status
# tmp$Cancer <- tmp$Score
# 
# 
# p_boxplot_HER2 <- ggplot(tmp,aes(x=Group,y=Cancer)) + geom_boxplot(outlier.shape = NA) + # facet_grid(. ~ Stage) +
#   xlab("HER2 status") + ylab("Score")+ylim(0,1) + theme_bw() +
#   geom_point(aes(fill = Group),size = 1.5, shape = 21, position = position_jitterdodge()) +
#   theme(legend.position = "none",
#         # text=element_text(face='bold'),
#         # axis.text.x=element_text(size=12),
#         # axis.text.y=element_text(size=12),
#         text=element_text(face='bold', size = 12),
#         axis.text.x=element_text(size=12),
#         axis.text.y=element_text(size=12),
#         panel.border = element_rect(colour = "black", fill=NA, size=2),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) + 
#   geom_hline(yintercept=c(tmp_value85, tmp_value90), linetype="dashed") +
#   scale_y_continuous(breaks = c(0,0.25,tmp_value85, tmp_value90,0.5,0.75,1)) + 
#   # annotate("text",label = c("85% sens", "90% sens"), x = 2.35, y = c(tmp_value85 + 0.015 ,tmp_value90 - 0.015) , color = "black", size = 4) +
#   stat_compare_means(method = "wilcox.test") +
#   geom_signif(comparisons = list(c("HER2+","HER2-")),map_signif_level=TRUE)
# 
# p_boxplot_HER2_stage <- ggplot(tmp,aes(x=Group,y=Cancer)) + geom_boxplot(outlier.shape = NA) + facet_grid(. ~ Stage) +
#   xlab("HER2 status") + ylab("Score")+ylim(0,1) + theme_bw() +
#   geom_point(aes(fill = Group),size = 1.5, shape = 21, position = position_jitterdodge()) +
#   theme(legend.position = "none",
#         # text=element_text(face='bold'),
#         # axis.text.x=element_text(size=12),
#         # axis.text.y=element_text(size=12),
#         text=element_text(face='bold', size = 12),
#         axis.text.x=element_text(size=12),
#         axis.text.y=element_text(size=12),
#         panel.border = element_rect(colour = "black", fill=NA, size=2),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) + 
#   geom_hline(yintercept=c(tmp_value85, tmp_value90), linetype="dashed") +
#   scale_y_continuous(breaks = c(0,0.25,tmp_value85, tmp_value90,0.5,0.75,1)) + 
#   # annotate("text",label = c("85% sens", "90% sens"), x = 2.35, y = c(tmp_value85 + 0.015 ,tmp_value90 - 0.015) , color = "black", size = 4) +
#   stat_compare_means(method = "wilcox.test") +
#   geom_signif(comparisons = list(c("HER2+","HER2-")),map_signif_level=TRUE)
# 
# tmp <- inner_join(bootstrap_samples[[i]], Info_list[,c("SampleID","TNBC")],c("SampleID"="SampleID"))
# tmp <- tmp[tmp$Response == "Cancer",]
# tmp <- tmp[tmp$Stage != "-",]
# tmp <- tmp[tmp$TNBC != "Unknown",]
# 
# tmp$Stage <- factor(tmp$Stage, levels = c('0',"I","II"), labels = c("0","I","II"))
# tmp$TNBC <- factor(tmp$TNBC, levels = c("TNBC","Non-TNBC"))
# tmp$Group <- tmp$TNBC
# tmp$Cancer <- tmp$Score
# 
# p_boxplot_TNBC <- ggplot(tmp,aes(x=Group,y=Cancer)) + geom_boxplot(outlier.shape = NA) + # facet_grid(. ~ Stage) +
#   xlab("TNBC status") + ylab("Score")+ylim(0,1) + theme_bw() +
#   geom_point(aes(fill = Group),size = 1.5, shape = 21, position = position_jitterdodge()) +
#   theme(legend.position = "none",
#         # text=element_text(face='bold'),
#         # axis.text.x=element_text(size=12),
#         # axis.text.y=element_text(size=12),
#         text=element_text(face='bold', size = 12),
#         axis.text.x=element_text(size=12),
#         axis.text.y=element_text(size=12),
#         panel.border = element_rect(colour = "black", fill=NA, size=2),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) + 
#   geom_hline(yintercept=c(tmp_value85, tmp_value90), linetype="dashed") +
#   scale_y_continuous(breaks = c(0,0.25,tmp_value85, tmp_value90,0.5,0.75,1)) + 
#   # annotate("text",label = c("85% sens", "90% sens"), x = 2.35, y = c(tmp_value85 + 0.015 ,tmp_value90 - 0.015) , color = "black", size = 4) +
#   stat_compare_means(method = "wilcox.test") +
#   geom_signif(comparisons = list(c("TNBC","Non-TNBC")),map_signif_level=TRUE)
# 
# p_boxplot_TNBC_stage <- ggplot(tmp,aes(x=Group,y=Cancer)) + geom_boxplot(outlier.shape = NA) + facet_grid(. ~ Stage) +
#   xlab("TNBC status") + ylab("Score")+ylim(0,1) + theme_bw() +
#   geom_point(aes(fill = Group),size = 1.5, shape = 21, position = position_jitterdodge()) +
#   theme(legend.position = "none",
#         # text=element_text(face='bold'),
#         # axis.text.x=element_text(size=12),
#         # axis.text.y=element_text(size=12),
#         text=element_text(face='bold', size = 12),
#         axis.text.x=element_text(size=12),
#         axis.text.y=element_text(size=12),
#         panel.border = element_rect(colour = "black", fill=NA, size=2),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) + 
#   geom_hline(yintercept=c(tmp_value85, tmp_value90), linetype="dashed") +
#   scale_y_continuous(breaks = c(0,0.25,tmp_value85, tmp_value90,0.5,0.75,1)) + 
#   # annotate("text",label = c("85% sens", "90% sens"), x = 2.35, y = c(tmp_value85 + 0.015 ,tmp_value90 - 0.015) , color = "black", size = 4) +
#   stat_compare_means(method = "wilcox.test") +
#   geom_signif(comparisons = list(c("TNBC","Non-TNBC")),map_signif_level=TRUE)
#  
# # ggsave(filename = "Figure_S5_new.pdf",
# #        plot_grid(
# #          plot_grid(p_boxplot_HR, p_boxplot_HR_stage, nrow = 1, rel_widths = c(1,2), labels = c("A","")),
# #          plot_grid(p_boxplot_HER2, p_boxplot_HER2_stage, nrow = 1, rel_widths = c(1,2), labels = c("B","")),
# #          plot_grid(p_boxplot_TNBC, p_boxplot_TNBC_stage, nrow = 1, rel_widths = c(1,2), labels = c("C","")),
# #          ncol = 1, 
# #          nrow = 3
# #        ),
# #        dpi = 600,width = 12,height = 18)
# 

###### cfFrag+cfMeth+BIRADS
load("Draft/GBP/Review_Figures/cfFrag_cfMeth_BIRADS/aml_restack_run_200.RData")
load("cfMeth_frag_model_0615.Rdata")

cfMeth_with_imagine_df <- inner_join(cfMeth_with_imagine_df %>% 
                                       rename(Both = automl_1),
                                     Train_list[,c("SampleID","automl_1")],
                                     c("SampleID"))

tmp <- cfMeth_with_imagine_df
tmp[tmp$Type == "Cancer", ]$Type <- 1
tmp[tmp$Type == "Healthy", ]$Type <- 0

ROC_LOOCV_Frag <- roc(tmp$Type, tmp$Frag_score, percent=F,smooth=F,ci=T)
ROC_LOOCV_cfMeth <- roc(tmp$Type, tmp$cfMeth, percent=F,smooth=F,ci=T)
ROC_LOOCV_Both <- roc(tmp$Type, tmp$Both, percent=F,smooth=F,ci=T)
ROC_LOOCV_ALL <- roc(tmp$Type, tmp$automl_1, percent=F,smooth=F,ci=T)

roc.LOOCV <- list()
# 
# LOOCV_A_to_F <- round(roc.test(ROC_LOOCV_ALL, ROC_LOOCV_Frag)$p.value,5)
# LOOCV_A_to_C <- round(roc.test(ROC_LOOCV_ALL, ROC_LOOCV_cfMeth)$p.value,5)

t <- round(ROC_LOOCV_ALL$ci, 3)
roc.LOOCV[[paste0("cfFrag + cfMeth + Xray + Ultrasound ", "(AUC=", t[2],")")]] <- ROC_LOOCV_ALL

t <- round(ROC_LOOCV_Both$ci, 3)
roc.LOOCV[[paste0("cfFrag + cfMeth ", "(AUC=", t[2],")")]] <- ROC_LOOCV_Both

t <- round(ROC_LOOCV_cfMeth$ci, 3)
roc.LOOCV[[paste0("cfMeth ", "(AUC=", t[2],")")]] <- ROC_LOOCV_cfMeth

t <- round(ROC_LOOCV_Frag$ci, 3)
roc.LOOCV[[paste0("cfFrag ", "(AUC=", t[2],")")]] <- ROC_LOOCV_Frag


# t <- round(ROC_LOOCV_Xray$ci, 3)
# roc.LOOCV[[paste0("Xray , ", t[2]," (",t[1],"-",t[3],")")]] <- ROC_LOOCV_Xray
# 
# t <- round(ROC_LOOCV_Ultra$ci, 3)
# roc.LOOCV[[paste0("Ultrasound , ", t[2]," (",t[1],"-",t[3],")")]] <- ROC_LOOCV_Ultra

p_LOOCV <- ggrocs(rocs = roc.LOOCV, breaks = seq(0,1,0.2), fontsize = 12, 
                  legendTitle= "Validation subcohort", ci = FALSE)

ggsave("Figure 4 Panel F updated.pdf",
       p_LOOCV,
       dpi = 600,width = 6,height = 6)
