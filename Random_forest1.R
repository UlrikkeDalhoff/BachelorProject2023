# Preparing of data
library(ape)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(GenomicRanges)
library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(ranger)
library(rsample)
library(pROC)

#NB requires BiocManager 
# install.packages("BiocManager")
list.of.packages <- c("ape", "BSgenome", "BSgenome.Athaliana.TAIR.TAIR9", "GenomicRanges", "tidyverse", "ggplot2", "rtracklayer", "BiocParallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

#Importing the iCLIP data
my_data_iCLIP <- rtracklayer::import("D:/Data_bachelor/iCLIP.bed")


#Importing the m6a data
my_data_m6a <- rtracklayer::import("D:/Data_bachelor/m6a_nanopore.bed")

#Make the datasets comparable 
seqlevels(my_data_m6a, pruning.mode = "coarse") <- c("Chr1","Chr2","Chr3","Chr4","Chr5")


#Importing the araport data 
my_data_Araport <- rtracklayer::import("D:/Data_bachelor/Araport11.gtf") 

#Finding the unique genes for iclip by using the araport data
araport_iclip <- my_data_Araport[subjectHits(findOverlaps(my_data_iCLIP,my_data_Araport))]

length(unique(araport_iclip$gene_id))



#Finding unique genes for m6a by using araport data
araport_m6a <- my_data_Araport[subjectHits(findOverlaps(my_data_m6a, my_data_Araport))]

length(unique(araport_m6a$gene_id))

#Import the background data
background <- rtracklayer::import("D:/Data_bachelor/background.bed")

#make the background comparable 
seqlevels(background, pruning.mode= "coarse") <- c("Chr1","Chr2","Chr3","Chr4","Chr5")


#making the groups for my 5 fold cross validation
genes_iclip <- my_data_Araport[nearest(my_data_iCLIP, my_data_Araport)]$gene_id

genes_background <- my_data_Araport[nearest(background, my_data_Araport)]$gene_id


genes_all <- c(genes_iclip, genes_background)

sample_vec <- sample(1:6, length(unique(genes_all)), replace = TRUE)

unique_genes <- unique(genes_all)
names(sample_vec) <- unique_genes
my_group_vector <- sample_vec[genes_all]
as.numeric(my_group_vector)




# we import the data needed. 

# Import of 50 nt window full significant set
matrix_50_full <- read.csv("D:/new_data_bachelor/matrix_50_full_new.csv") #done
matrix_50_full_background <- read.csv("D:/new_data_bachelor/matrix_50_full_background_new.csv") #done


#Import of 50 nt window with 50 most significant hits
matrix_50_50 <- read.csv("D:/new_data_bachelor/matrix_50_50_new.csv") #done
matrix_50_50_background <- read.csv("D:/new_data_bachelor/matrix_50_50_background_new.csv") #done


#import of 50 nt window with 10 most significant hits
matrix_50_10 <- read.csv("D:/new_data_bachelor/matrix_50_10_new.csv") #done
matrix_50_10_background <- read.csv("D:/new_data_bachelor/matrix_50_10_background_new.csv") #done 


#Import of 25 nt window full significant set
matrix_25_full <- read.csv("D:/new_data_bachelor/matrix_25_full_new.csv") #done
matrix_25_full_background <- read.csv("D:/new_data_bachelor/matrix_25_full_background_new.csv") #done


#Import of 25 nt window with 50 most significant hits
matrix_25_50 <- read.csv("D:/new_data_bachelor/matrix_25_50_new.csv") #done
matrix_25_50_background <- read.csv("D:/new_data_bachelor/matrix_25_50_background_new.csv") #done


#Import of 25 nt window with 10 most significant hits
matrix_25_10 <- read.csv("D:/new_data_bachelor/matrix_25_10_new.csv") #done
matrix_25_10_background <- read.csv("D:/new_data_bachelor/matrix_25_10_background_new.csv") #done


#Import of 10 nt window full significant set
matrix_10_full <- read.csv("D:/new_data_bachelor/matrix_10_full_new.csv") #need to be run
matrix_10_full_background <- read.csv("D:/new_data_bachelor/matrix_10_full_background_new.csv") #need to be run


#Import of 10 nt window with 50 most significant hits
matrix_10_50 <- read.csv("D:/new_data_bachelor/matrix_10_50_new.csv") #done
matrix_10_50_background <- read.csv("D:/new_data_bachelor/matrix_10_50_background_new.csv") #done


#Import of 10 nt window with 10 most significant hits 
matrix_10_10 <- read.csv("D:/new_data_bachelor/matrix_10_10_new.csv") #done
matrix_10_10_background <- read.csv("D:/new_data_bachelor/matrix_10_10_background_new.csv") #done


?group_vfold_cv

set.seed(2305)
?ranger
# Define the model as a function, to make fine-tuning easier. 
model_rf <- function(data_set = seq_matrix, data_background = bachground, n_size = 10, n_trees = 500, importance = "impurity",
                  CV = FALSE, m_try = NULL, split_rule = "gini"){

  #how many kmers we want. #Combine into data matrix and label vector
  print(length(data_set))
  data_matrix <- data_set[,2:(n_size+1)]
  background_matrix <- data_background[,2:(n_size+1)]
  data_combined <- rbind(data_matrix, background_matrix)
  
  n1 <- nrow(data_matrix)
  n2 <- nrow(background_matrix)
  
  vec1 <- rep(1, n1)
  vec2 <- rep(0, n2)
  
  labels <- c(vec1, vec2)
  
  #We combine the sub_labels (Our answers) into the sub_combined matrix
  final <- cbind(labels, data_combined)
  
  #final <- final %>% mutate(split = my_group_vector)
  
  #this part if CV = TRUE:
  # 5fold cross validation
  #Create the 5 folds
  if (CV == TRUE){
    final <- final %>% mutate(split = my_group_vector)
    data_split <- group_initial_split(final, prop = 0.80, group = split)
    
    training_data <- data_split %>% training()
    testing_data <- data_split %>% testing()
    

    #we need to correct this, so the folds do not include the same genes 
    cv_split <- group_vfold_cv(training_data, group = split, v = 5)
    
    #Building the model 
    cv_data <- cv_split %>%
      mutate(
        train = map(splits, ~.x %>% training),
        validate = map(splits, ~.x %>% testing)
      )
    
    
    randomforest_function <- ~ranger(formula = labels~.,
                                     data = .x, 
                                     num.trees = n_trees,
                                     classification = TRUE,
                                     probability = TRUE,
                                     mtry = m_try,
                                     splitrule = split_rule)
    cv_models_rf <- cv_data %>% 
      mutate(model = map(train, randomforest_function))
    
    #making prediction
    prediction_function <- ~predict(.x, .y)$predictions
    
    identity_fun <- function(x){
      return(x)
    }
    
    cv_prep_rf <- cv_models_rf %>%
      mutate(
        actual_site = map(validate, ~.x %>% pull(labels)),
        predicted_site = map2(.x = model, .y = validate, prediction_function),
        models = map(.x = model, identity_fun))
    
    #print(cv_prep_rf$models)
    
   # cv_prep_rf$predicted_site <- cv_prep_rf %>% 
    #  mutate(
     #   predicted_site = if_else(predicted_site$predictions >= 0.5, 1, 0))
    
  
    #Evaluating our results on the validation set 
    prediction_df <- cv_prep_rf %>% select(id, actual_site, predicted_site) %>% unnest(c(actual_site, predicted_site))
    

    
    my_accuracies <- prediction_df %>%
      group_by(id) %>%
      summarise("accuracy" = sum(actual_site == if_else(predicted_site[,2] >= 0.5, 1, 0))/n())
    
    print(my_accuracies)
    
    
    
    #mae = mean absolut error, so the smaller mae = better model.
    #mae_df <- prediction_df %>% group_by(id) %>% summarise("mae" = mean(abs(actual_site - if_else(predicted_site[,2] >= 0.5, 1, 0))))
    
    
    
    
    #ROC_df <- prediction_df %>% group_by(id) %>% summarise("AUC" =as.numeric(roc(labels ~ predicted_site[,2], plot = FALSE, print.auc = TRUE)$auc))
    
    
    
    
    #Test and evaluate the final model 
    model_final <- ranger(formula = labels~.,
                            data = training_data, 
                            num.trees = n_trees,
                            classification = TRUE,
                            probability = TRUE,
                            mtry = m_try,
                            splitrule = split_rule)
    
    prediction_df1 <- tibble("actual_site" = testing_data %>% pull(labels),
                             "predicted_site" = predict(model_final, testing_data)$predictions)
    
    
    #pred <- predict(model_final, testing_data)$predictions
 
    
    return(list(prediction_df1 %>% summarise("AUC" = as.numeric(roc(actual_site ~ predicted_site[,2], quiet = TRUE, plot = FALSE, print.auc = TRUE)$auc)),
           my_accuracies %>% summarise("accuracy" = mean(accuracy))))
    
  }
  
  
  rf <- ranger(formula = labels~.,
               data = final, 
               num.trees = n_trees,
               classification = TRUE,
               probability = TRUE,
               mtry = m_try,
               splitrule = split_rule,
               importance = importance)
  return(rf)
}



#grid search
tune_func <- function(data_set = seq_matrix, data_background = bachground, n_size = 10, p1, p2){
  res_df <- data.frame(Number_trees = vector("numeric", length(p1)*length(p2)),
                       M_try = vector("numeric", length(p1)*length(p2)),
                       Accuracy = vector("numeric", length(p1)*length(p2)),
                       AUC = vector("numeric", length(p1)*length(p2)))
  for (i in 1:length(p1)){
    for (j in 1:length(p2)){
      index <- i+(j-1)*length(p1)
      model <- model_rf(data_set = data_set, data_background = data_background, n_size = n_size, n_trees = p1[i], m_try = p2[j], CV = TRUE)
      res_df$Number_trees[index] <- (p1[i])
      res_df$M_try[index] <- (p2[j])
      res_df$Accuracy[index] <- model[2]
      res_df$AUC[index] <- model[1]
      gc()
    } 
  }
  return(res_df)
}

#tuning_50_full <- tune_func(data_set = matrix_50_full, data_background = matrix_50_full_background, n_size = 729, p1 = c(500, 1000, 2000), p2 = c(10, 17, 27))

#tuning_10_full <- tune_func(data_set = matrix_10_full, data_background = matrix_10_full_background, n_size = 602, p1 = c(500, 1000, 2000), p2 = c(10, 14, 24))

#tuning_50_50 <- tune_func(data_set = matrix_50_50, data_background = matrix_50_50_background, n_size = 50, p1 = c(500, 1000, 2000), p2 = c(2, 5, 7))

#tuning_50_10 <- tune_func(data_set = matrix_50_10, data_background = matrix_50_10_background, n_size = 10, p1 = c(500, 1000, 2000), p2 = c(1, 2, 3))

#tuning_25_full <- tune_func(data_set = matrix_25_full, data_background = matrix_25_full_background, n_size = 693, p1 = c(500, 1000, 2000), p2 = c(10, 16, 26))

#tuning_10_full <- tune_func(data_set = matrix_10_full, data_background = matrix_10_full_background, n_size = 602, p1 = c(500, 1000, 2000), p2 = c(10, 14, 24))

#tuning_25_50 <- tune_func(data_set = matrix_25_50, data_background = matrix_25_50_background, n_size = 50, p1 = c(500, 1000, 2000), p2 = c(2, 5, 7))

#tuning_10_50 <- tune_func(data_set = matrix_10_50, data_background = matrix_10_50_background, n_size = 50, p1 = c(500, 1000, 2000), p2 = c(2, 5, 7))

tuning_25_10 <- tune_func(data_set = matrix_25_10, data_background = matrix_25_10_background, n_size = 10, p1 = c(500, 1000, 2000), p2 = c(1, 2, 3))

tuning_25_10

tuning_25_50

tuning_50_10

tuning_50_10 <- tune_func(data_set = matrix_50_10, data_background = matrix_50_10_background, n_size = 50, p1 = c(500, 1000, 2000), p2 = c(1, 2, 3)) 

tuning_10_50 

tuning_10_10

tuning_10_full 

tuning_25_10

tuning_25_50

tuning_25_full

tuning_50_full

tuning_50_50

tuning_50_10 #done


Hypertribe_for_rf <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf_test.csv")

#hypertribe_finalform <- Hypertribe_for_rf[, colnames(Hypertribe_for_rf) %in% colnames(matrix_25_full)]
#hypertribe_finalform


set.seed(2305)
#Creating the final model 
finalmodel_RF_10 <- model_rf(data_set = matrix_10_full, data_background = matrix_10_full_background, n_size = 602, n_trees = 1000, m_try = 14)
save(finalmodel_RF_10, file ="final_RF_model_10.RData")

sorted_importance <- sort(finalmodel_RF_10$variable.importance, decreasing = TRUE)

the_most_important <- sorted_importance[1:20]

# install.packages("Polychrome")
library(Polychrome)

# build-in color palette
Glasbey = glasbey.colors(32)

important_df
#load araport data:
load("D:/new_data_bachelor/final_RF_model_10.RData")

# Convert named list to a data frame
df <- data.frame(Name = names(the_most_important), Value = unlist(the_most_important))

# Sort the data frame by Value in descending order
df <- df[order(-df$Value), ]

# Reorder the levels of Name based on the sorted order
df$Name <- factor(df$Name, levels = df$Name)


# Create the barplot using ggplot
ggplot(df, aes(x = Name, y = Value, fill = Name)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Motifs (variables)", y = "Variable Importance", title = "20 most important variables for RF10") + 
  theme_classic() + 
  coord_flip()

#Doing all of the above just with the RF 50 
load("D:/new_data_bachelor/final_RF_model.RData")
sorted_importance_50 <- sort(finalmodel_RF_50$variable.importance, decreasing = TRUE)

the_most_important_50 <- sorted_importance_50[1:20]

df_50 <- data.frame(Name = names(the_most_important_50), Value = unlist(the_most_important_50))

# Sort the data frame by Value in descending order
df_50 <- df_50[order(-df_50$Value), ]

# Reorder the levels of Name based on the sorted order
df_50$Name <- factor(df_50$Name, levels = df_50$Name)

# Create the barplot using ggplot
ggplot(df_50, aes(x = Name, y = Value, fill = Name)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Motifs (variables)", y = "Variable Importance", title = "20 most important variables for RF50") + 
  theme_classic() + 
  coord_flip()



#load araport data:
hyper_gene1 <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf10_gene1.csv")


#predicting exiting

pred.hyper <- predict(finalmodel_RF_10, hyper_gene1)

length(pred.hyper$predictions[,2])

pred.hyper$predictions[,2]

xes <- 1:127
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT1G01020 from position 6788-6914 Chr1") +
  ylim(0,1)





#Prediction on gene 2

# Load rfmodel:

#load hypertribe data:
Hypertribe_gene2 <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf10_gene2.csv")



#predicting exiting

pred.hyper <- predict(finalmodel_RF_10, Hypertribe_gene2)

length(pred.hyper$predictions[,2])

xes <- 1:533
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))

                                 
ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT1G01080 from position 44970-45502 Chr1") +
  ylim(0,1)






#Prediction on gene 3

# Load rfmodel:
load("D:/new_data_bachelor/final_RF_model.RData")

#load hypertribe data:
Hypertribe_gene3 <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf10_gene3.csv")



#predicting exiting

pred.hyper <- predict(finalmodel_RF_10, Hypertribe_gene3)

length(pred.hyper$predictions[,2])

xes <- 1:248
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT2G01060 from position 73208-73455 Chr2") +
  ylim(0,1)




#Prediction on gene 4

# Load rfmodel:
load("D:/new_data_bachelor/final_RF_model.RData")

#load hypertribe data:
Hypertribe_gene4 <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf10_gene4.csv")



#predicting exiting

pred.hyper <- predict(finalmodel_RF_10, Hypertribe_gene4)

length(pred.hyper$predictions[,2])

xes <- 1:337
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT2G40400 from position 16872570-16872906 Chr2") +
  ylim(0,1)





#Prediction on gene 5

# Load rfmodel:
load("D:/new_data_bachelor/final_RF_model.RData")

#load hypertribe data:
Hypertribe_gene5 <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf10_gene5.csv")



#predicting exiting

pred.hyper <- predict(finalmodel_RF_10, Hypertribe_gene5)

length(pred.hyper$predictions[,2])

xes <- 1:631
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT3G02460 from position 508093 - 508723 Chr3") +
  ylim(0,1)





#Prediction on gene 6

# Load rfmodel:
load("D:/new_data_bachelor/final_RF_model.RData")

#load hypertribe data:
Hypertribe_gene6 <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf10_gene6.csv")



#predicting exiting

pred.hyper <- predict(finalmodel_RF_10, Hypertribe_gene6)

length(pred.hyper$predictions[,2])

xes <- 1:304
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT3G22220 from position 7839504-7839807 Chr3") +
  ylim(0,1)





#Prediction on gene 7

# Load rfmodel:
load("D:/new_data_bachelor/final_RF_model.RData")

#load hypertribe data:
Hypertribe_gene7 <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf10_gene7.csv")



#predicting exiting

pred.hyper <- predict(finalmodel_RF_10, Hypertribe_gene7)

length(pred.hyper$predictions[,2])

xes <- 1:229
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT4G11850 from position 7129123-7129351 Chr4") +
  ylim(0,1)





#Prediction on gene 8

# Load rfmodel:
load("D:/new_data_bachelor/final_RF_model.RData")

#load hypertribe data:
Hypertribe_gene8 <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf10_gene8.csv")



#predicting exiting

pred.hyper <- predict(finalmodel_RF_10, Hypertribe_gene8)

length(pred.hyper$predictions[,2])

xes <- 1:444
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT4G39850 from position 18496763-18497206 Chr4") +
  ylim(0,1)




#Prediction on gene 9

# Load rfmodel:
load("D:/new_data_bachelor/final_RF_model.RData")

#load hypertribe data:
Hypertribe_gene9 <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf10_gene9.csv")



#predicting exiting

pred.hyper <- predict(finalmodel_RF_10, Hypertribe_gene9)

length(pred.hyper$predictions[,2])

xes <- 1:197
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT5G50310 from position 20479282-20479478 Chr5") +
  ylim(0,1)



#Prediction on gene 10

# Load rfmodel:
load("D:/new_data_bachelor/final_RF_model.RData")

#load hypertribe data:
Hypertribe_gene10 <- read.csv("D:/new_data_bachelor/Hypertribe_for_rf10_gene10.csv")



#predicting exiting

pred.hyper <- predict(finalmodel_RF_10, Hypertribe_gene10)

length(pred.hyper$predictions[,2])

xes <- 1:30
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT5G67640 from position 26969516-26969545 Chr5") +
  ylim(0,1)










#Prediction on gene test1 

# Load rfmodel:
load("D:/new_data_bachelor/final_RF_model.RData")

#load m6a data:
araport_gene_test1 <- read.csv("D:/new_data_bachelor/araport_for_rf_gene_test1.csv")



#predicting exiting

pred.araport <- predict(finalmodel_RF_50, araport_gene_test1)

length(pred.araport$predictions[,2])

pred.araport$predictions[,2]

xes <- 1:165
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.araport$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "probabilities from the Araport data of gene AT1G01225 from position 97243-97407 Chr1") +
  ylim(0,1)









# 10 window test

set.seed(2305)

finalmodel_RF_10 <- model_rf(data_set = matrix_10_full, data_background = matrix_10_full_background, n_size = 602, n_trees = 1000, m_try = 14)

#load araport data:
load("D:/new_data_bachelor/final_RF_model_10.RData")
araport_gene_test1_10 <- read.csv("D:/new_data_bachelor/araport_for_rf_gene_test1_10.csv")



#predicting exiting

pred.araport <- predict(finalmodel_RF_10, araport_gene_test1_10)

length(pred.araport$predictions[,2])

pred.araport$predictions[,2]

xes <- 1:366
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.araport$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities from the Araport data of gene AT1G12070 from position 4078547-4078912 Chr1") +
  ylim(0,1)




# load araport data
araport_gene_test2_10 <- read.csv("D:/new_data_bachelor/araport_for_rf_gene_test2_10.csv")



#predicting exiting

pred.araport <- predict(finalmodel_RF_10, araport_gene_test2_10)

length(pred.araport$predictions[,2])

pred.araport$predictions[,2]

xes <- 1:149
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.araport$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities from the Araport data of gene AT1G80810 from position 30368899-30369047 Chr1") +
  ylim(0,1)





#creating the same plots, but with the 50 model 

#load araport data:
load("D:/new_data_bachelor/final_RF_model.RData")
araport_gene_test1_50 <- read.csv("D:/new_data_bachelor/araport_for_rf_gene_test1_50.csv")



#predicting exiting

pred.araport <- predict(finalmodel_RF_50, araport_gene_test1_50)

length(pred.araport$predictions[,2])

pred.araport$predictions[,2]

xes <- 1:366
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.araport$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF50 probabilities from the Araport data of gene AT1G12070 from position 4078547-4078912 Chr1") +
  ylim(0,1)



#load araport data:
araport_gene_test2_50 <- read.csv("D:/new_data_bachelor/araport_for_rf_gene_test2_50.csv")

#predicting exiting

pred.araport <- predict(finalmodel_RF_50, araport_gene_test2_50)

length(pred.araport$predictions[,2])

pred.araport$predictions[,2]

xes <- 1:149
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.araport$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF50 probabilities from the Araport data of gene AT1G80810 from position 30368899-30369047 Chr1") +
  ylim(0,1)




#we do the plots again but with a cut-off maximum on 75

ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF50 probabilities from the Araport data of gene AT1G12070 from position 4078547-4078912 Chr1.") +
  ylim(0,1) +
  geom_hline(yintercept = 0.75, linetype = "dashed")





#And then for our model50 hypertribe

Hypertribe_genecutoff <- read.csv("D:/new_data_bachelor/Hypertribe_for_RF50_gene2cut.csv")



#predicting exiting

pred.hyper <- predict(finalmodel_RF_50, Hypertribe_genecutoff)

length(pred.hyper$predictions[,2])

xes <- 1:533
data_frame_results <- data.frame(xes = c(xes),
                                 predictions = c(pred.hyper$predictions[,2]))


ggplot(data = data_frame_results) + 
  geom_line(aes(x = xes, y = predictions)) + theme_classic() + 
  labs( x = "Positions", y = "probabilities", title = "RF10 probabilities for HyperTRIBE of gene AT1G01080 from position 44970-45502 Chr1") +
  ylim(0,1) +
  geom_hline(yintercept = 0.75, linetype = "dashed")





# Then I make a barplot on my gridsearch AUC results. 
AUC_df <- data.frame(Tuning_model = c("Tuning 10 full", "Tuning 25 full", "Tuning 50 full"),
                     AUC = c(0.695637, 0.724108, 0.742877))

ggplot(data = AUC_df, aes(x=Tuning_model, y = AUC, fill = Tuning_model)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  ggtitle("AUC values of the full models") +
  theme_classic() +
  scale_fill_brewer(palette = "Set3")
  
