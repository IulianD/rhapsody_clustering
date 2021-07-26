pkgload::load_all('/media/sf_shareddisk/datashield/dsSwissKnifeClient')
pkgload::load_all('/media/sf_shareddisk/datashield/dsHelpersClient')
library(dsBaseClient)
library(magrittr)
library(reshape2)
library(ggplot2)
library(scales)
library(grid)
source('./cluster_functions.R')
cohorts <- c(  'accelerate')
logindata <- read.delim('../logindata_DSOpal.txt')
logindata <- logindata[logindata$server %in% cohorts,,drop=FALSE] 
opals <- datashield.login(logindata)
#datashield.aggregate(opals, quote(set.stringsAsFactors(TRUE)))
datashield.assign(opals, 'lb', 'rhapsody.LB')
datashield.assign(opals, 'vs', 'rhapsody.VS')
datashield.assign(opals, 'dm', 'rhapsody.DM')
#datashield.assign(opals, 'lb', as.symbol('convert.units("lb")'))
#datashield.assign(opals[c('dcs', 'godarts')],'mh', 'rhapsody.MH')
datashield.assign(opals,'mh', 'rhapsody.MH')
#datashield.assign(opals, 'vs', as.symbol('convert.units("vs")'))
prepare_cluster_godarts(with.hdl = TRUE)
prepare_cluster_andis(with.hdl = TRUE)
prepare_cluster_dcs(with.hdl = TRUE)

clusters_all <- do_clustering(cohorts, centers = 5, iter.max = 50, nstart = 30)
gg <- do_ggplot(clusters_all,cluster.labels = c('3/SIRD', '2/SIDD' ,'5/MARD/low HDL' , '6/MARD/high HDL', '4/MOD'), create.pdf = TRUE)
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
},]

