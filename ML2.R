#installing the used packages
install.packages(c("tidyr", "magrittr", "dplyr", "GGally", "caret","glmnet", "pROC", "corrplot", "BootValidation"," ElemStatLearn", "skimr","reshape2"))

# the used libraries
library(tidyr)
library(magrittr)
library(dplyr)
library(GGally)
library(caret)
library(glmnet)
library(pROC)
library(corrplot)
library(BootValidation)
library(ElemStatLearn)
library(skimr)
library(reshape2)

#********************************************************************************************************

#download all the data
dataset <- read.csv("./GDS504.clean.csv",header=TRUE)

#subset certain genes 1500:2500
data <- data.frame(dataset[1500:2500,])

#combine Probe_ID with gene_ID to have distinct Genes identifier
data$ID <- paste(data$PROBE_ID_REF,data$GENE_IDENTIFIER)

#number of genes we have is 1001
nrow(distinct(data,ID))

# convert ID to factor with 1001 level
data$ID <- factor(data$ID)


#change coloums names of pateints to numbers 1:20
colnames(data)[3:22] <- 1:20
rownames(data) <- 1:1001


#re-organise the data and store the samples in one column and the expression profile in another coloum
data_PAH<-  data[,-c(1,2)] %>% gather(patient_status, value, -ID)%>% arrange(ID)

#rearrange coloums sequence
data_PAH <- data_PAH[,c(2,1,3)]

#spread data to have each gene as coloum
data_genes <- data_PAH %>% group_by(patient_status)%>% distinct(ID, .keep_all = TRUE) %>% spread(ID , value)

# convert data_genes tibble to dataframe
data_genes <- data.frame(data_genes) %>% arrange(patient_status)

#data_genes <- data.frame(data_genes,stringsAsFactors = FALSE)
#convert patient status to factor with 20 levels
data_genes$patient_status <- factor(data_genes$patient_status)

#combine levels to have two levels hypertensive (first 1:14) =0 and healthy(15:20)=1
levels(data_genes$patient_status) <- list("hyper"=as.character(c(1:14)) , "healthy"= as.character(c(15:20)))
#data_genes$patient_status <- factor(data_genes$patient_status,levels=c(1,0),labels = c("healthy","hyper"))

#**************************************************************************************************************************************************************
#Exploratory Analaysis


# first 5 observations of the data after cleaning
sample_n(data_genes,5)[1:5]

#summary of data before data cleaning
kable(skim(dataset))

#provide part of the correlation matrix of genes
cor(as.matrix(data_genes[, 2:1002]))[cor(as.matrix(data_genes[, 2:1002])) > 0.7][1:10]


#plot part of the correlation matrix with higher value than 0.8
index <- which(cor(as.matrix(data_genes[, 2:1002])) > 0.8)
corrplot(cor(as.matrix(data_genes[,2:1002]))[,index[1:10]],method="number",type="upper")

#partition data in training data %70 and testing data %30
set.seed(345)
set.seed(54321)
train_index <- data_genes$patient_status %>% createDataPartition(p = 0.7, list = FALSE)
train_data  <- data_genes[train_index, ]
test_data <- data_genes[-train_index, ]

#dimension of training data
dim(train_data)
#dimension of testing data
dim(test_data)

#standarise all variables (except for the response variable)
train_data[ ,2:1002] <- apply(train_data[,2:1002] , 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
test_data[ ,2:1002] <- apply(test_data[,2:1002] , 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
#data_genes[ ,2:1002] <- apply(data_genes[,2:1002] , 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))                           

#fit the data using training data predictors
x <- model.matrix(patient_status ~ .-1,data=train_data)
#response
y <- train_data$patient_status

#test data and test respone
x.test <- model.matrix(patient_status ~.-1,data= test_data)
y.test <- test_data$patient_status

n <- dim(x)[1]

##****************************************************************************************************************************************
##
#Ridge regression

#fit ridge regression for differenr values of lambda
set.seed(54321)
fit.ridge <- glmnet(x,y,family= "binomial",alpha=0)
plot(fit.ridge , xvar="lambda" , label=TRUE)

#k-fold cross validation to find the best lambda with minimum error
cv.ridge <- cv.glmnet(x,y,alpha=0,family= "binomial",nfolds=5, type.measure="class")
plot(cv.ridge)

#lambda value associated with the minimum deviance
cv.ridge$lambda.min

#lambda value within one standard error of the one with the minimum deviance
cv.ridge$lambda.1se

#fit the ridge model with the optimum lambda (lambda_1se)
fit.ridge.1se.lambda <- glmnet(x,y,family = "binomial" , alpha=0 , lambda=cv.ridge$lambda.1se )
#coefficients when lambda value within one standard error of the one with the minimum deviance
coef.ridge <- coef(fit.ridge.1se.lambda)

#number of coefficients equal to zero : non of the coef is zero
table(as.data.frame(as.matrix(coef.ridge== 0)))

#fit the ridge model with the optimum lambda (lambda_min)
fit.ridge.min.lambda <- glmnet(x,y,family = "binomial" , alpha=0 , lambda=cv.ridge$lambda.min)
#coefficients when lambda value gives the  minimum deviance
#coef(fit.ridge.min.lambda)


# CV model accuracy using testing data for model with lambda_1se
pred.class.ridge.1se <- fit.ridge.1se.lambda %>% predict(x , type="class")
confusionMatrix(factor(pred.class.ridge.1se , levels= c('hyper','healthy')), reference = y, positive = "healthy")


#model accuracy using testing data for model with lambda_1se
pred.class.ridge.1se <- fit.ridge.1se.lambda %>% predict(x.test , type="class")
confusionMatrix(factor(pred.class.ridge.1se , levels= c('hyper','healthy')), reference = y.test , positive = "healthy")


# Validate glmnet using bootstrap.
vboot(fit.ridge.1se.lambda, x, y, nfolds = 5, B = 200, s= cv.ridge$lambda.1se, cv_replicates = 20, n_cores = 3)


#****************************************************************************************************************
#************************************************************************************************************
#fit LASSO

#fit lasso regression for different values of lambda
fit.lasso <- glmnet(x,y,family= "binomial",alpha=1)
#plot the result 
plot(fit.lasso, xvar="lambda" , label=TRUE )
plot(fit.lasso, xvar="dev" , label=TRUE)

#5-fold cross validation to find the best lambda with minimum error
set.seed(54321)
cv.lasso<- cv.glmnet(x,y,alpha=1 , family= "binomial",nfolds=5, type.measure = "class")
#plot the result showing lambda with minimum error
plot(cv.lasso)

#lambda value associated with the minimum deviance
cv.lasso$lambda.min

#lambda value within one standard error of the one with the minimum deviance
cv.lasso$lambda.1se

#lasso model with lambda_min
fit.lasso.min.lambda <- glmnet(x,y,family = "binomial" , alpha=1 , lambda=cv.lasso$lambda.min)
#coefficients when lambda value gives the minimum deviance
#coef(fit.lasso.min.lambda)

#lasso model with lambda_1se
fit.lasso.1se.lambda <- glmnet(x,y,family = "binomial" , alpha=1 , lambda=cv.lasso$lambda.1se )
#coefficients when lambda value within one standard error of the one with the minimum deviance
coef.lasso<-coef(fit.lasso.1se.lambda)

#model accuracy using testing data for the model with lambda_1se
pred.class.lasso.1se <- fit.lasso.1se.lambda %>% predict(x , type="class")
confusionMatrix(factor(pred.class.lasso.1se , levels= c('hyper','healthy')) , reference = y, positive = "healthy")


#CV model accuracy using testing data for model with lambda value within one standard error of the one with the minimum deviance
pred.class.lasso.1se <- fit.lasso.1se.lambda %>% predict(x.test, type="class")
confusionMatrix(factor(pred.class.lasso.1se , levels= c('hyper','healthy')) , reference = y.test, positive = "healthy")

#number of coefficients equal to zero 
table(as.data.frame(as.matrix(coef.lasso== 0)))

#find non_zero genes names
non_zero_genes_lasso<-rownames(as.data.frame(as.matrix(coef.lasso)))[as.data.frame(as.matrix(coef.lasso)) !=0]

# Validate glmnet using bootstrap.
vboot(fit.lasso.1se.lambda, x, y, nfolds = 5,B = 100, s= cv.lasso$lambda.1se, cv_replicates = 20, n_cores = 3)


#***********************************************************************************************************
#************************************************************************************************************

#Elastic net


#fit Elastic net and tune the model to identify the best alpha and lambda
set.seed(54321)
model.net <- train(patient_status ~.-1, data = train_data, method = "glmnet", trControl = trainControl("CV",number=5),tuneLength = 50 ,metric="Accuracy")

#check coeffieient for the final model with best lambda and alpha
coef.net<-coef(model.net$finalModel,model.net$bestTune$lambda)


#number of zero and non zero genes coeff
table(as.data.frame(as.matrix(coef.net)== 0))

#find non_zero genes names
non_zero_genes_net<-rownames(as.data.frame(as.matrix(coef.net)))[as.data.frame(as.matrix(coef.net)) !=0]

# CV model accuracy using testing data for model with the best lambda and alpha.
pred.class.elastic <- model.net %>% predict(x , type="raw")
confusionMatrix(factor(pred.class.elastic , levels= c('hyper','healthy')) , reference = y , positive = "healthy")

#model accuracy using testing data for model with the best lambda and alpha.
pred.class.elastic <- model.net %>% predict(x.test , type="raw")
confusionMatrix(factor(pred.class.elastic , levels= c('hyper','healthy')) , reference = y.test, positive = "healthy")


#fit the elastic net model with best alpha and different values of lambda
model.net1 <- glmnet(x,y,family='binomial',alpha= model.net$bestTune[1])
plot(model.net1,xvar="lambda")

#with best value of alpha fit the elastic model
cv.net<- cv.glmnet(x,y,family='binomial',alpha= model.net$bestTune[1],type.measure = "class")
plot(cv.net)

# plot of binomial deviance as a function of lambda
plot(model.net$glmnet.fit,xvar="dev")

# Validate glmnet using bootstrap.
vboot(model.net1, x, y, nfolds = 5,B = 100, s= cv.lasso$lambda.1se, cv_replicates = 20, n_cores = 3)

#correlation between selected genes by elastic model
cor(data_genes[,non_zero_genes_net[-1]])


#plot the correlation matrix between selected genes
corrplot(cor(data_genes[,c(non_zero_genes_net[-1],"M21302_at.SPRR2D")]),type="upper",method="circle")
corrplot(cor(data_genes[,c(non_zero_genes_net[-1],"M21302_at.SPRR2D")]),type="upper",method="number")

#******************************************************************************************************************
#the results of cross-validation using 1000 di???erent realizations of patient status
CV_lasso=1:1000
CV_elastic=1:1000
CV_ridge=1:1000
for (i in 1:1000)
{set.seed(i)
  
  # train_index <- data_genes$patient_status %>% createDataPartition(p = 0.7, list = FALSE)
  # train_data  <- data_genes[train_index, ]
  # test_data <- data_genes[-train_index, ]
  
  training_index <- sample(seq_len(nrow(data_genes)), size = 20)
  data_genes <- data_genes[training_index,]
  
  # #dimension of training data
  # dim(train_data)
  # #dimension of testing data
  # dim(test_data)
  
  # #standarise all variables (except for the response variable)
  # train_data[ ,2:1002] <- apply(train_data[,2:1002] , 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
  # test_data[ ,2:1002] <- apply(test_data[,2:1002] , 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
  data_genes[ ,2:1002] <- apply(data_genes[,2:1002] , 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))                           
  # =map(perms$map,~ model.matrix(patient_status ~ .-1,data=data_genes))
  #fit the data using training data
  #predictors
  x <- model.matrix(patient_status ~ .-1,data= data_genes)
  #response
  y <- data_genes$patient_status
  
  # x.test <- model.matrix(patient_status ~.-1,data= test_data)
  # y.test <- test_data$patient_status
  
  set.seed(i)
  cv.lasso<- cv.glmnet(x,y,alpha=1 , family= "binomial",nfolds=5, type.measure = "dev")
  cv.ridge <- cv.glmnet(x,y,alpha=0,family= "binomial",nfolds=5, type.measure="class")
  model.net <- train(patient_status ~.-1, data = data_genes, method = "glmnet", trControl = trainControl("CV",number=5),tuneLength = 50 ,metric="Accuracy")
  #coefficients when lambda value within one standard error of the one with the minimum deviance
  fit.lasso.1se.lambda <- glmnet(x,y,family = "binomial" , alpha=1 , lambda=cv.lasso$lambda.1se )
  fit.ridge.1se.lambda <- glmnet(x,y,family = "binomial" , alpha=0 , lambda=cv.ridge$lambda.1se )
  
  CV_lasso[i] <- sum(y ==predict(fit.lasso.1se.lambda,x,type="class"))
  CV_ridge[i] <- sum(y ==predict( fit.ridge.1se.lambda,x,type="class"))
  CV_elastic[i] <- sum(y== predict(model.net,x , type="raw"))
}
#count the error in each iterations
Cv_lasso_error <- 20 - CV_lasso
CV_elastic_error <- 20 - CV_elastic
CV_ridge_error <- 20 -CV_ridge

#frequency table of the CV error of each model

table(Cv_lasso_error)
mean(CV_lasso_error)
table(CV_ridge_error)
mean(CV_ridge_error)
table(CV_elastic_error)
mean(CV_elastic_error)

ggplot(data = melt(CV_Error), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable)) +xlab("Models")+ylab("CV error (out of 20)") +theme(plot.title = element_text(size=12))+theme(legend.position = "none")


##**************************************************************************************************************
#check for model instability
CV_lasso=1:1000
CV_elastic=1:1000
CV_ridge=1:1000
for (i in 1:1000)
{set.seed(i)
  
  # train_index <- data_genes$patient_status %>% createDataPartition(p = 0.7, list = FALSE)
  # train_data  <- data_genes[train_index, ]
  # test_data <- data_genes[-train_index, ]
  
  training_index <- sample(seq_len(nrow(data_genes)), size = 20)
  data_genes <- data_genes[training_index,]
  
  # #dimension of training data
  # dim(train_data)
  # #dimension of testing data
  # dim(test_data)
  
  # #standarise all variables (except for the response variable)
  # train_data[ ,2:1002] <- apply(train_data[,2:1002] , 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
  # test_data[ ,2:1002] <- apply(test_data[,2:1002] , 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
  data_genes[ ,2:1002] <- apply(data_genes[,2:1002] , 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))                           
  # =map(perms$map,~ model.matrix(patient_status ~ .-1,data=data_genes))
  #fit the data using training data
  #predictors
  x <- model.matrix(patient_status ~ .-1,data= data_genes)
  #response
  y <- data_genes$patient_status
  
  # x.test <- model.matrix(patient_status ~.-1,data= test_data)
  # y.test <- test_data$patient_status
  
  set.seed(i)
  #cv.lasso<- cv.glmnet(x,y,alpha=1 , family= "binomial",nfolds=5, type.measure = "dev")
  #cv.ridge <- cv.glmnet(x,y,alpha=0,family= "binomial",nfolds=5, type.measure="class")
  #model.net <- train(patient_status ~.-1, data = data_genes, method = "glmnet", trControl = trainControl("CV",number=5),tuneLength = 50 ,metric="Accuracy")
  #coefficients when lambda value within one standard error of the one with the minimum deviance
  fit.lasso.1se.lambda <- glmnet(x,y,family = "binomial" , alpha=1 , lambda=cv.lasso$lambda.1se )
  fit.ridge.1se.lambda <- glmnet(x,y,family = "binomial" , alpha=0 , lambda=cv.ridge$lambda.1se )
  model.net1 <- glmnet(x,y,family = "binomial" , alpha=model.net$bestTune[1] , lambda=model.net$bestTune[2] )
  CV_lasso[i] <- sum(y ==predict(fit.lasso.1se.lambda,x,type="class"))
  CV_ridge[i] <- sum(y ==predict( fit.ridge.1se.lambda,x,type="class"))
  CV_elastic[i] <- sum(y== predict(model.net,x , type="raw"))
}

