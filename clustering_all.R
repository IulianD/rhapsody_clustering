pkgload::load_all('/media/sf_shareddisk/datashield/dsSwissKnifeClient')
pkgload::load_all('/media/sf_shareddisk/datashield/dsSwissKnife')
pkgload::load_all('/media/sf_shareddisk/datashield/dsHelpersClient')
library(dsBaseClient)
library(magrittr)
library(reshape2)
library(ggplot2)
library(scales)
library(grid)
source('./cluster_functions.R')
cohorts <- c(  'abos', 'accelerate', 'andis','dcs', 'godarts')
#cohorts <- c( 'andis','dcs', 'godarts')
#cohorts <- c('accelerate')
logindata <- read.delim('../logindata_DSOpal.txt')
logindata <- logindata[logindata$server %in% cohorts,,drop=FALSE] 
opals <- datashield.login(logindata)
#datashield.aggregate(opals, quote(set.stringsAsFactors(TRUE)))
datashield.assign(opals, 'lb', 'rhapsody.LB')
datashield.assign(opals, 'vs', 'rhapsody.VS')
datashield.assign(opals, 'dm', 'rhapsody.DM')
#datashield.assign(opals, 'lb', as.symbol('convert.units("lb")'))
#datashield.assign(opals[c('dcs', 'godarts')],'mh', 'rhapsody.MH')
try(datashield.assign(opals,'mh', 'rhapsody.MH'), silent = TRUE) # not present everywhere
#datashield.assign(opals, 'vs', as.symbol('convert.units("vs")'))
prepare_cluster_godarts(with.hdl = TRUE)
prepare_cluster_andis(with.hdl = TRUE)
prepare_cluster_dcs(with.hdl = TRUE)
prepare_cluster_accelerate(with.hdl = TRUE)
prepare_cluster_abos(with.hdl = TRUE)


clusters_all <- do_clustering(cohorts, centers = 5, iter.max = 10, nstart = 3)

datashield.symbols(opals)
ds.summary('kmeans_input_scaled_km_clust5')
ds.cbind(c('kmeans_input_scaled','kmeans_input_scaled_km_clust5'), newobj = 'kmeans_output');
ds.summary('kmeans_output')
dssColNames('kmeans_output', value = 'cluster_number', to.replace = 'kmeans_input_scaled_km_clust5')
dssShowFactors('kmeans_output')
ds.class('kmeans_output$cluster_number')
dssSubset("forests_train", "kmeans_output", row.filter = 'sample(nrow(kmeans_output),nrow(kmeans_output)*5/10)', datasources = opals)
dssSubset("forests_test", "kmeans_output", row.filter = '!(rownames(kmeans_output) %in% rownames(forests_train))', datasources = opals)
ds.summary('forests_train')
ds.summary('forests_test')
forests <- dssRandomForest(train = list(what = 'forests_train', dep_var = 'cluster_number'))
pred <- dssRandomForest(test = list(forests = forests,testData = 'forests_test'), train = NULL)
sapply(pred, function(x) x$accuracy)
sapply(pred, function(x) x$confusion_matrix)
dssSubset("forests_train", "kmeans_output", row.filter = 'sample(nrow(kmeans_output),nrow(kmeans_output)*9/10)', datasources = opals)
dssSubset("forests_test", "kmeans_output", row.filter = '!(rownames(kmeans_output) %in% rownames(forests_train))', datasources = opals)
forests2 <- dssRandomForest(train = list(what = 'forests_train', dep_var = 'cluster_number'))
pred2 <- dssRandomForest(test = list(forests = forests2,testData = 'forests_test'), train = NULL)
sapply(pred2, function(x) x$accuracy)
sapply(pred2, function(x) x$confusion_matrix)

pca <- dssPrincomp('kmeans_output')
ds.summary('kmeans_output_scores')
#dssMean('kmeans_output_scores$Comp.1')
biplot(pca$global,choices = c(1,2), draw.arrows = FALSE, levels = 'kmeans_output_scores$cluster_number')
biplot(pca$global,choices = c(1,3), draw.arrows = FALSE, levels = 'kmeans_output_scores$cluster_number')
biplot(pca$global,choices = c(2,3), draw.arrows = FALSE, levels = 'kmeans_output_scores$cluster_number')


#xxx_tst <- data.frame(CPEPTIDE =c(-1.21736769, -1.01701768, -0.63200772, -0.11906342,  0.50552169 , 1.20010862,  1.76977853,  0.02650179 ),
#                      HBA1C=c(-1.121462754, -1.028895100, -0.658624487, -0.195786219,  0.452187354,  1.192728582,  1.840702156, -0.006046347),
#                      HDL = c(1.39329268, -1.22391599, -0.71578591, -0.17378049,  0.53760163,  1.31673442,  1.84857724, -0.01881389),
#                      BMI = c(-1.4763975, -1.1975066, -0.7359768, -0.1080810,  0.5409829,  1.2464616,  1.6981353, -0.0284235),
#                      AGE = c(-1.69083875, -1.28110588, -0.53558678,  0.02117467,  0.72060924,  1.18177489,  1.57078733,  0.01402004 ),
#                      cluster_number  = c(1,2,3,4,5,3,2,4))


#f <- .encode.arg(forests, serialize.it = TRUE)
#forestDSS('test', 'testData',f )
#sil <- Reduce(rbind,clusters_all$cluster.object$global$silhouette)
#mean(sil[,3])
clusters_split <- do_clustering(cohorts, centers = 5, iter.max = 30, nstart = 10, ktype = 'split')


xx <-matrix(nrow=0, ncol=3)
for (i in clusters_split$cluster.object){
xx <- rbind(xx,i$silhouette)
}
mean(xx[,3])

gg <- do_ggplot(clusters_all,cluster.labels = c('3/SIRD', '2/SIDD' ,'5/MARD/low HDL' , '6/MARD/high HDL', '4/MOD'), create.pdf = TRUE)
gg_split <- split_ggplot(clusters_all,cluster.lbls = c('3/SIRD', '2/SIDD' ,'5/MARD/low HDL' , '6/MARD/high HDL', '4/MOD'), create.pdf = TRUE)


do_pca(clusters_all)
#c('3/SIRD','5/MARD/low HDL'  , '6/MARD/high HDL', '4/MOD', '2/SIDD')
#c('3/SIRD', '2/SIDD' ,'5/MARD/low HDL' , '6/MARD/high HDL', '4/MOD')
clusters_godarts <- do_clustering(c('godarts'), centers = 5, iter.max = 40, nstart = 30)
gg_godarts <- do_ggplot(clusters_godarts,cluster.labels = c('3/SIRD','2/SIDD'  ,  '5/MARD/low HDL'   ,  '6/MARD/high HDL' , '4/MOD'), create.pdf = FALSE)
#c('3/SIRD','5/MARD/low HDL'  , '6/MARD/high HDL', '4/MOD', '2/SIDD')
#c('3/SIRD','2/SIDD'  ,  '5/MARD/low HDL'   ,  '6/MARD/high HDL' , '4/MOD')
#do_pca(clusters_godarts)
clusters_dcs <- do_clustering(c('dcs'), centers = 5, iter.max = 40, nstart = 30)
gg_dcs <- do_ggplot(clusters_dcs,cluster.labels = c('6/MARD/high HDL','2/SIDD','5/MARD/low HDL', '4/MOD', '3/SIRD'), create.pdf = FALSE)
#c('3/SIRD'          ,'5/MARD/low HDL'  , '6/MARD/high HDL', '4/MOD', '2/SIDD')
#c( '6/MARD/high HDL','2/SIDD','6/MARD/high HDL', '4/MOD', '3/SIRD')

#do_pca(clusters_godarts)
clusters_andis <- do_clustering(c('andis'), centers = 5, iter.max = 40, nstart = 30)
gg_andis <- do_ggplot(clusters_andis,cluster.labels = c('3/SIRD','2/SIDD','5/MARD/low HDL' , '4/MOD' ,'6/MARD/high HDL'), create.pdf = FALSE)
#c('3/SIRD','2/SIDD','5/MARD/low HDL' , '6/MARD/high HDL','4/MOD')
#c('3/SIRD','2/SIDD','5/MARD/low HDL' , '4/MOD' ,'6/MARD/high HDL')

#do_pca(clusters_godarts)
#c('6/MARD/high HDL','5/MARD/low HDL'  , '3/SIRD', '4/MOD', '2/SIDD')
#c('3/SIRD','5/MARD/low HDL'  , '6/MARD/high HDL', '4/MOD', '2/SIDD')
#c('2/SIDD', '3/SIRD', '4/MOD', '5/MARD/low HDL', '6/MARD/high HDL')
#c('6/MARD/high HDL','5/MARD/low HDL'  , '3/SIRD', '4/MOD', '2/SIDD')
clusters_no_hdl <- do_clustering(c('andis', 'dcs'), centers = 4, iter.max = 40)
#do_ggplot(clusters,cluster.labels =c('3/SIRD','2/SIDD','4/MOD','6/MARD/high HDL', '5/MARD/low HDL'))
do_ggplot(clusters_no_hdl,cluster.labels =c('5/MARD','3/SIRD','4/MOD','2/SIDD'))
#c('2/SIDD', '3/SIRD', '4/MOD', '5/MARD/low HDL', '6/MARD/high HDL')
#c('3/SIRD','2/SIDD','4/MOD','6/MARD/high HDL', '5/MARD/low HDL')
#c('3/SIRD','2/SIDD','4/MOD','5/MARD')
#c('5/MARD','3/SIRD','4/MOD','2/SIDD')
do_pca(clusters_no_hdl)

op <-opal.login('administrator', 'IRELAND-BEAUTY-wrong-MONEY', url='https://rhap-fdb01.vital-it.ch/accelerate')
opal.assign(op, 'lb', 'rhapsody.LB')
lb <- opal.execute(op,'lb')
hdl <- lb[lb$LBTESTCD == 'HDL',]

earlies <- hdl[{
              x <- hdl$SUBJID  
              hdl$LBDTC == unlist(sqldf::sqldf(paste0('select min(LBDTC) from hdl where SUBJID =', x)))
              },]

earlies <- hdl[{
  lapply(unique(hdl$SUBJID), function(x){
      y <- hdl
      hdl$LBDTC == min(unlist(y[y$SUBJID == x, 'LBDTC']))
  })
},]

hdl <- hdl[order(hdl$SUBJID),]

earlies <- hdl[{
  Reduce(c,lapply(unique(hdl$SUBJID), function(x){
    y <- hdl
    this.vec <- unlist(y[y$SUBJID == x, 'LBDTC'])
    if(length(levels(this.vec)) == 0 ){
      rt <- c(TRUE, rep(FALSE, length(this.vec) - 1))
    } else {
      dt <- min(as.vector(this.vec))
      rt <-unlist(lapply(this.vec, function(t){
        if (t==dt){
          return(TRUE)
        } else {
          return(FALSE)
       }
    }))
    }
    return(rt)
  }))
expr <- list(as.symbol('partialKmeans'), 'kmeans_input', dsSwissKnifeClient:::.encode.arg(as.data.frame(my_centers)), NULL, TRUE)
# kms <- datashield.aggregate(datasources, as.symbol(expr), async = async)
kms <- datashield.aggregate(opals['abos'], as.call(expr), async = FALSE)
