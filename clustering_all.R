### execute only once: 
install.packages('dsBaseClient', repos = c('https://cran.obiba.org', 'https://stat.ethz.ch/CRAN'))
devtools::install_github('sib-swiss/dsSwissKnifeClient')
devtools::install_github('IulianD/dsHelpersClient', ref = 'DSOpal')
#####
library(dsSwissKnifeClient)
library(dsHelpersClient)
library(dsBaseClient)
library(magrittr)
library(reshape2)
library(ggplot2)
library(scales)
library(grid)
source('./cluster_functions.R')
cohorts <- c(  'abos', 'accelerate', 'andis','dcs', 'godarts')
logindata <- read.delim('./logindata_DSOpal.txt')
logindata <- logindata[logindata$server %in% cohorts,,drop=FALSE] 
opals <- datashield.login(logindata)
# load all the tables:
datashield.assign(opals, 'lb', 'rhapsody.LB')
datashield.assign(opals, 'vs', 'rhapsody.VS')
datashield.assign(opals, 'dm', 'rhapsody.DM')
try(datashield.assign(opals,'mh', 'rhapsody.MH'), silent = TRUE) # not present everywhere

# prepare the kmeans input dataframe
prepare_cluster_godarts(with.hdl = TRUE)
prepare_cluster_andis(with.hdl = TRUE)
prepare_cluster_dcs(with.hdl = TRUE)
prepare_cluster_accelerate(with.hdl = TRUE)
prepare_cluster_abos(with.hdl = TRUE)

# iter.max = 60, nstart = 60+ for significant results - takes a while
clusters_all <- do_clustering(cohorts, centers = 5, iter.max = 10, nstart = 3)


datashield.symbols(opals)

# create the input for random forests
ds.summary('kmeans_input_scaled_km_clust5')
# add the classification column to the input:
ds.cbind(c('kmeans_input_scaled','kmeans_input_scaled_km_clust5'), newobj = 'kmeans_output');
ds.summary('kmeans_output')
#give it a nicer name:
dssColNames('kmeans_output', value = 'cluster_number', to.replace = 'kmeans_input_scaled_km_clust5')
dssShowFactors('kmeans_output')
ds.class('kmeans_output$cluster_number')
# round 1, train on 90% of the data (picked randomly):
dssSubset("forests_train", "kmeans_output", row.filter = 'sample(nrow(kmeans_output),nrow(kmeans_output)*9/10)', datasources = opals)
# and test on what's left:
dssSubset("forests_test", "kmeans_output", row.filter = '!(rownames(kmeans_output) %in% rownames(forests_train))', datasources = opals)
ds.summary('forests_train')
ds.summary('forests_test')
#train:
forests <- dssRandomForest(train = list(what = 'forests_train', dep_var = 'cluster_number'))
#test:
pred <- dssRandomForest(test = list(forests = forests,testData = 'forests_test'), train = NULL)
#check the prediction
sapply(pred, function(x) x$accuracy)
sapply(pred, function(x) x$confusion_matrix)
#same as above but with only half of the data for training:
dssSubset("forests_train", "kmeans_output", row.filter = 'sample(nrow(kmeans_output),nrow(kmeans_output)*5/10)', datasources = opals)
dssSubset("forests_test", "kmeans_output", row.filter = '!(rownames(kmeans_output) %in% rownames(forests_train))', datasources = opals)
forests2 <- dssRandomForest(train = list(what = 'forests_train', dep_var = 'cluster_number'))
pred2 <- dssRandomForest(test = list(forests = forests2,testData = 'forests_test'), train = NULL)
#and test:
sapply(pred2, function(x) x$accuracy)
sapply(pred2, function(x) x$confusion_matrix)

#PCA:
pca <- dssPrincomp('kmeans_output')
ds.summary('kmeans_output_scores')
dssMean('kmeans_output_scores$Comp.1')
biplot(pca$global,choices = c(1,2), draw.arrows = TRUE, levels = 'kmeans_output_scores$cluster_number')
biplot(pca$global,choices = c(1,3), draw.arrows = TRUE, levels = 'kmeans_output_scores$cluster_number')
biplot(pca$global,choices = c(2,3), draw.arrows = TRUE, levels = 'kmeans_output_scores$cluster_number')


# generate plots(if create.image=true a pdf will be created); 
# the cluster.labels msut be supplied in the proper order to match the actual clusters.
gg <- do_ggplot(clusters_all,cluster.labels = c('3/SIRD', '2/SIDD' ,'5/MARD/low HDL' , '6/MARD/high HDL', '4/MOD'), create.pdf = TRUE, include.split = TRUE)

