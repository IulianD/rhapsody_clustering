my.smooth2d <-  function (x, y, npoints = 128, 
                          colour.pool = grep("\\d", grep("dark|deep", colours(), value = T), 
                                             invert = T, value = T), shades = 8, draw.image = FALSE, 
                          labels = TRUE, axes = TRUE, categories = NULL, emphasize_level = 0, type = "combine", 
                          async = FALSE, wait = TRUE, datasources = NULL) 
{
  if (is.null(datasources)) {
    datasources <- dsBaseClient:::findLoginObjects()
  }
  bandwidth <- list(x = .my.bw.args(x, datasources), y = .my.bw.args(y, 
                                                                     datasources))
  lims <- dssRange(x, y, type = type, datasources = datasources)
  holder <- NULL
  if (!is.null(categories)) {
    keep.cols <- NULL
    holder <- "tempSmooth"
    dssSubsetByClass(c(x, y), subsets = holder, variables = categories, 
                      keep.cols = keep.cols, async = async, wait = wait, 
                      datasources = datasources)
  }
  bandwidth.arg <- dsCDISCclient:::.encode.arg(bandwidth)
  lims.arg <- dsCDISCclient:::.encode.arg(lims)
  holder <- dsCDISCclient:::.encode.arg(holder)
  expr <- paste0("partial.kde2d('", holder, "', '", x, "', '", 
                 y, "','", bandwidth.arg, "','", lims.arg, "',", npoints, 
                 ")")
  res <- opal::datashield.aggregate(datasources, as.symbol(expr), 
                                    async = FALSE, wait = wait)
  res <- res[!sapply(res, is.null)]
  if (length(res) == 0) {
    return(NULL)
  }
  if (draw.image) {
    xlab <- ""
    ylab <- ""
    if (labels) {
      xlab <- x
      ylab <- y
    }
  }
  sw <- dssSwapKeys(res)
  if (type == "combine") {
    img <- list(Combined = sapply(sw, function(this.level) {
      list(x = this.level[[1]]$x, y = this.level[[1]]$y, 
           z = Reduce("+", lapply(this.level, function(a) {
             a$len * a$z
           }))/Reduce("+", lapply(this.level, function(a) a$len)))
    }, simplify = FALSE))
  }
  else {
    img <- sapply(res, function(this.node) {
      sapply(this.node, function(this.level) {
        this.level[c("x", "y", "z")]
      }, simplify = FALSE)
    }, simplify = FALSE, USE.NAMES = TRUE)
  }
  my.legend <- NULL
  if (draw.image) {
    startcols <- list()
    endcols <- list()
    levs <- names(sw)
    
    st <- grep("\\d", grep("dark|deep", colour.pool, invert = T, 
                           value = T), invert = T, value = T)
    en <- grep("\\d", grep("light|white", colour.pool, invert = T, 
                           value = T), invert = T, value = T)
    startcols[levs] <- sample(colour.pool, length(levs))
    endcols[levs] <- sample(setdiff(colour.pool, startcols[levs]), 
                            length(levs))
    startcols[levs] <- rep('white', length(levs))
   # endcols[levs] <- c("#4DB3E6","#00B399","#E69900","#CC80B3","#8B1A4F")
    endcols[levs] <- c("red", "green", "darkcyan", "orange", "blue")[1:length(levs)]
    my.legend <- sapply(img, function(this.set) {
      ret <- list()
      categories <- names(this.set)
      if(emphasize_level >0){
        categories <- c(categories[-emphasize_level], categories[emphasize_level]) # move the important one in the upper layer
      }
      sapply(categories, function(this.name) {
        ret[[this.name]] <<- c(startcols[[this.name]], 
                               endcols[[this.name]])
        if(emphasize_level > 0){
          if(levs[emphasize_level] == this.name){
            this.alpha <- 1
          } else {
            this.alpha <- 0.3
  
          }

          image(this.set[[this.name]], col = dsCDISCclient:::.add.alpha(colorRampPalette(c(startcols[[this.name]], 
                                                                                           endcols[[this.name]]))(shades), c(0,rep(this.alpha, length.out = shades-1 ))), 
                xlab = xlab, ylab = ylab, axes = axes)
          
        } else {
          image(this.set[[this.name]], col = dsCDISCclient:::.add.alpha(colorRampPalette(c(startcols[[this.name]], 
                                                                                           endcols[[this.name]]))(shades), c(0, seq(0.2, 
                                                                                                                                    0.8, length.out = shades - 1))), xlab = xlab, 
                ylab = ylab, axes = axes)
        }
        par(new = T)
      })
      par(new = F)
      return(ret)
    }, simplify = FALSE)
  }
  ret <- list(lims = lims, img = img)
  if (!is.null(my.legend)) {
    ret[["legend"]] <- unlist(unname(my.legend), recursive = FALSE)
  }
  invisible(ret)
}





.my.bw.args <- function (x, datasources = NULL) 
{
  if (is.null(datasources)) {
    datasources <- dsBaseClient:::findLoginObjects()
  }
  r <- (ds.quantileMean(x, type = "combine", datasources)[c("25%", 
                                                            "75%")])
  v <- ds.var(x, type = "combine", datasources)[[1]]
  l <- ds.length(x, datasources, type = "combine")[[1]]
  list(quarts = r, var = v, len = l)
}

prepare_cluster_andis <-  function(with.hdl = TRUE){
  this.cohort <- c('andis')
  rowfilter_andis <- "((LBTESTCD=='HBA1C' & LBORRESU=='mmol/mol') | (LBTESTCD== 'CPEP' ) | LBTESTCD == 'HDL')  & (VISIT %in% c('DIAGNOSIS', 'RECRUITMENT') ) "
  dssSubset('lb2', 'lb', row.filter = rowfilter_andis, async = FALSE, datasources = opals[this.cohort] )
  #LB2
  lb_with_visit = c('CPEP', 'HBA1C')
  if(with.hdl){
    lb_with_visit = c(lb_with_visit, 'HDL')
  }
  for (m in lb_with_visit){
    pre <- paste0('pre_', m)
    by_visit <- paste0('by_visit_', m)
    row_filter <- paste0("LBTESTCD ==  '", m, "' & !is.na(LBORRES)")
    dssSubset(pre, 'lb2', row.filter = row_filter , col.filter = "c('SUBJID', 'LBORRES', 'VISIT')", async = FALSE, datasources = opals[this.cohort])
    dssSubsetByClass(pre, by_visit, variables = paste0(pre, '$VISIT'), keep.cols = c('SUBJID', 'LBORRES'), async = FALSE, datasources = opals[this.cohort])
    #keep most from diagnosis, only supplementary ones from recruitment
    rec <- paste0('rec_', m)
    row_fil <- paste0('!(SUBJID %in% ', by_visit, '$DIAGNOSIS$SUBJID)')
    dssSubset( rec, paste0(by_visit, '$RECRUITMENT'), row.filter = row_fil, async = FALSE, datasources = opals[this.cohort])
    final <- paste0('lb_',m)
    dssRbind(final, paste0(by_visit, '$DIAGNOSIS'), rec, async = FALSE, datasources = opals[this.cohort] )
    dssColNames(final, value = c('SUBJID', m), async = FALSE, datasources = opals[this.cohort])
  }
  
  
  # grab HDL and HBA1C without visit information
  lb_without_visit = c('HBA1C')
  if(with.hdl){
    lb_without_visit =c('HBA1C', 'HDL')
  }
  for (m in lb_without_visit){
    novisit <- paste0('no_visit_', m)
    final <- paste0('lb_',m)
    row_filter <- paste0(" LBTESTCD == '", m ,"' & !is.na(LBORRES) & !(SUBJID %in% ", final, "$SUBJID)")
    dssSubset(novisit, 'lb', row.filter = row_filter, col.filter = "c('SUBJID', 'LBORRES', 'LBDTC')", async = FALSE, datasources = opals[this.cohort] )
    ordered_name <- paste0('ordered_',m)
    dssSubset(ordered_name, novisit, row.filter = 'ordered(SUBJID)', async = FALSE, datasources = opals[this.cohort])
    row_filter_block <- paste0('{
                               
                               Reduce(c,lapply(unique(SUBJID), function(my_id){
                               my_dates <- ', ordered_name,  '[', ordered_name, '$SUBJID == my_id, "LBDTC"]; 
                               my_dates <- as.Date(as.vector(my_dates));
                               diag_date <- lb[lb$SUBJID == my_id & lb$VISIT == "DIAGNOSIS" & lb$LBTESTCD == "',m,'", "LBDTC"];
                               diag_date <- as.Date(as.vector(diag_date));
                               if(length(diag_date) == 0 || is.na(diag_date)){
                               #keep the first date
                               out <- c(TRUE, rep(FALSE, length(my_dates) - 1))
                               } else { 
                               out <- (my_dates == my_dates[which.min(abs(my_dates - diag_date))]);
                               #if we have only FALSE and NA, make the first NA TRUE:
                               #if(any(is.na(out)) && all(out[!is.na(out)] == FALSE)){
                               #    true_index <- which(is.na(out))[1];
                               #    out <- rep(FALSE, length(out));
                               #    out[true_index] <- TRUE;
                               #}
                               }
                               
                               out
                               }))
                               
  }')
    #cat(row_filter_block)
    #dssSubset(ordered_name, ordered_name, row.filter = row_filter_block, col.filter = "c('SUBJID', 'LBORRES')",async = FALSE, datasources = opals[this.cohort])
    dssSubset(ordered_name, ordered_name, row.filter = '!duplicated(SUBJID)',col.filter = "c('SUBJID', 'LBORRES')", async = FALSE, datasources = opals[this.cohort])
    dssColNames(ordered_name, value = c('SUBJID', m), async = FALSE, datasources = opals[this.cohort])
    dssRbind(ordered_name, final, ordered_name, async = FALSE, datasources = opals[this.cohort])
}
  
  join_input <- c('ordered_HBA1C', 'lb_CPEP')
  if(with.hdl){
    join_input <- c(join_input, 'ordered_HDL')
  }
  
  dssJoin(join_input, symbol = 'comp_lb', by = 'SUBJID', join.type = 'inner' , async = FALSE, datasources = opals[this.cohort])
  
  #VS
  dssSubset('comp_vs', 'vs', row.filter = "VSTESTCD == 'BMI' & (SUBJID %in% comp_lb$SUBJID) & !is.na(VSORRES)" , col.filter = "c('SUBJID', 'VSORRES')", async = FALSE, datasources = opals[this.cohort])
  dssColNames('comp_vs', value = c('SUBJID', 'BMI'), async = FALSE, datasources = opals[this.cohort])
  
  #DM
  dssSubset('comp_dm', 'dm',  col.filter = "c('SUBJID', 'AGE', 'SEX')", async = FALSE, datasources = opals[this.cohort])
  dssSubset('comp_dm', 'comp_dm', row.filter = '!is.na(AGE)  & (SUBJID %in% comp_vs$SUBJID)', async = FALSE, datasources = opals[this.cohort])
  
  #join
  dssJoin(c('comp_lb', 'comp_vs', 'comp_dm'), 'join_all', by = 'SUBJID', join.type = 'inner', async = FALSE, datasources = opals[this.cohort])
  
  
  if(with.hdl){
    output_cols <- 'c("AGE", "BMI", "CPEP", "HBA1C", "HDL")'
  } else {
    output_cols <- 'c("AGE", "BMI", "CPEP", "HBA1C")'
  }
  dssSubset('kmeans_input', 'join_all', col.filter = output_cols, async = FALSE, datasources = opals[this.cohort])
  #dssSubsetByClass('join_all', subsets = 'kmeans_by_sex', variables = 'join_all$SEX', keep.cols = eval(parse(text = output_cols)), async = FALSE, datasources = opals[this.cohort])
  #dssSubset('kmeans_M', 'kmeans_by_sex$M', row.filter = '1==1', async = FALSE, datasources = opals[this.cohort])
  #dssSubset('kmeans_F', 'kmeans_by_sex$F', row.filter = '1==1', async = FALSE, datasources = opals[this.cohort])
  
  vars <- sapply(ds.colnames('kmeans_input', datasources = opals[this.cohort])[[this.cohort]], function(x){
    this.col <- paste0('kmeans_input$',x)
    ds.var(this.col, datasources = opals[this.cohort])
  }, simplify = FALSE)
  
  fivesds <- sapply(vars, function(x){
    5*sqrt(x$Variance.by.Study[1,1])
    #x$Variance.by.Study[1,1]
  }, simplify = FALSE) %>% unlist
  
  colm <- dssColMeans('kmeans_input', datasources = opals[this.cohort])[[this.cohort]]$means
  
  names(fivesds) <- names(colm)
  
  sapply(names(colm), function(x){
    maxx <- colm[x] + fivesds[x]
    minx <- colm[x] - fivesds[x]
    my.filter <- paste0(x, ' > ', minx, ' & ', x, ' < ', maxx)
    dssSubset('kmeans_input', 'kmeans_input', row.filter = my.filter, async = FALSE, datasources = opals[this.cohort])
  })
  
  
  }

prepare_cluster_dcs <- function(with.hdl = TRUE){
  this.cohort <- c( 'dcs')
  #LB2
  if(with.hdl){
    needed_cols <- "c('CPEPTIDE', 'HBA1C','HDL')"
  } else {
    needed_cols <- "c('CPEPTIDE', 'HBA1C')"
  }
  dssSubset('lb2', 'lb', row.filter = paste0("LBTESTCD %in% ", needed_cols) , col.filter = "c('SUBJID', 'LBORRES', 'LBTESTCD', 'VISIT','LBORRESU', 'LBDTC')", async = FALSE, datasources = opals[this.cohort])
  dssSubset('mh2', 'mh', row.filter = "MHTERM == 'TYPE 2 DIABETES' & MHOCCUR == 'Y' & VISIT == 'DIAGNOSIS'", async = FALSE, datasources = opals[this.cohort])
  #mh DIAGNOSIS is duplicated for some reason
  dssSubset('mh2', 'mh2', row.filter = '!duplicated(mh2)', async = FALSE, datasources = opals[this.cohort])
  dssSubsetByClass('lb2', subsets = 'by_test', variables = 'lb2$LBTESTCD', async = FALSE, datasources = opals[this.cohort])
  available_tests <- ds.names('by_test', datasources = opals[this.cohort])[[this.cohort]]
  final_name <- 'tests_around_diag'
  dssSubset(final_name, paste0('by_test$',available_tests[1]), row.filter = '1==2', col.filter = 'c("SUBJID","LBORRES", "LBORRESU", "LBTESTCD", "VISIT")', async = FALSE, datasources = opals[this.cohort])
  sapply(available_tests, function(this.test){
    ordered_name <- paste0('by_test_', this.test)
    out_name <- paste0(this.test, '_around_diag')
    dssSubset(ordered_name, paste0('by_test$', this.test), row.filter = 'order(SUBJID)', async = FALSE, datasources = opals[this.cohort])
    
    
    row_filter_block <- paste0('{
                               
                               Reduce(c,lapply(unique(SUBJID), function(my_id){
                               my_dates <- ', ordered_name,  '[', ordered_name, '$SUBJID == my_id, "LBDTC"]; 
                               my_dates <- as.Date(as.vector(my_dates));
                               diag_date <- mh2[mh2$SUBJID == my_id, "MHDTC"];
                               diag_date <- as.Date(as.vector(diag_date));
                               if(length(diag_date) == 0 || is.na(diag_date)){
                               out <- c(TRUE, rep(FALSE, length(my_dates) - 1))
                               } else { 
                               out <- (my_dates == my_dates[which.min(abs(my_dates - diag_date))]);
                               }
                               out
                               }))
                               
  }')
  
    dssSubset(out_name, ordered_name, row.filter = row_filter_block, col.filter = 'c("SUBJID","LBORRES", "LBORRESU", "LBTESTCD")', async = FALSE, datasources = opals[this.cohort]) 
    dssRbind(final_name, final_name, out_name, async = FALSE, datasources = opals[this.cohort])
    
})
  
  
  dssPivot('wide_tests', 'tests_around_diag', value.var = 'LBORRES', formula = 'SUBJID ~ LBTESTCD + LBORRESU', async = FALSE, datasources = opals[this.cohort])
  dssSubset('comp_lb', 'wide_tests', row.filter = 'complete.cases(wide_tests)', async = FALSE, datasources = opals[this.cohort])
  dssColNames('comp_lb', value = c('SUBJID', eval(parse(text = needed_cols))), async = FALSE, datasources = opals[this.cohort])
  
  
  
  
  #VS2
  dssSubset('vs2', 'vs', row.filter = "VSTESTCD == 'BMI'" , col.filter = "c('SUBJID', 'VSORRES', 'VSDTC')", async = FALSE, datasources = opals[this.cohort])
  dssSubset('vs2', 'vs2', row.filter = '!is.na(VSORRES)', async = FALSE, datasources = opals[this.cohort])
  ordered_name <- 'by_test_BMI'
  out_name <- paste0('BMI_around_diag')
  dssSubset(ordered_name, 'vs2', row.filter = 'order(SUBJID)', async = FALSE, datasources = opals[this.cohort])
  
  
  row_filter_block <- paste0('{
                             
                             Reduce(c,lapply(unique(SUBJID), function(my_id){
                             my_dates <- ', ordered_name,  '[', ordered_name, '$SUBJID == my_id, "VSDTC"]; 
                             my_dates <- as.Date(as.vector(my_dates));
                             diag_date <- mh2[mh2$SUBJID == my_id, "MHDTC"];
                             diag_date <- as.Date(as.vector(diag_date));
                             if(length(diag_date) == 0 || is.na(diag_date)){
                             out <- c(TRUE, rep(FALSE, length(my_dates) - 1))
                             } else { 
                             out <- (my_dates == my_dates[which.min(abs(my_dates - diag_date))]);
                             }
                             return(out)
                             }))
                             
                             }')
  
  dssSubset(out_name, ordered_name, row.filter = row_filter_block, col.filter = 'c("SUBJID","VSORRES")', async = FALSE, datasources = opals[this.cohort]) 
  dssSubset('comp_vs', out_name, row.filter = paste0('complete.cases(',out_name,')'), async = FALSE, datasources = opals[this.cohort])
  dssColNames('comp_vs', value = c('SUBJID', 'BMI'), async = FALSE, datasources = opals[this.cohort])
  #DM
  
  dssSubset('dm2', 'dm',  col.filter = "c('SUBJID', 'AGE', 'SEX', 'DMDTC')", async = FALSE, datasources = opals[this.cohort])
  dssSubset('dm2', 'dm2', row.filter = '!is.na(AGE) & !is.na(SEX)', async = FALSE, datasources = opals[this.cohort])
  ordered_name <- paste0('by_test_dm')
  out_name <- 'dm_around_diag'
  dssSubset(ordered_name, 'dm2', row.filter = 'order(SUBJID)', async = FALSE, datasources = opals[this.cohort])
  
  
  row_filter_block <- paste0('{
                             
                             Reduce(c,lapply(unique(SUBJID), function(my_id){
                             my_dates <- ', ordered_name,  '[', ordered_name, '$SUBJID == my_id, "DMDTC"]; 
                             my_dates <- as.Date(as.vector(my_dates));
                             diag_date <- mh2[mh2$SUBJID == my_id, "MHDTC"];
                             diag_date <- as.Date(as.vector(diag_date));
                             if(length(diag_date) == 0 || is.na(diag_date)){
                             out <- c(TRUE, rep(FALSE, length(my_dates) - 1))
                             } else { 
                             out <- (my_dates == my_dates[which.min(abs(my_dates - diag_date))]);
                             }
                             out
                             }))
                             
                             }')
  
  dssSubset(out_name, ordered_name, row.filter = row_filter_block, col.filter = 'c("SUBJID","AGE", "SEX" )', async = FALSE, datasources = opals[this.cohort]) 
  
  dssSubset('comp_dm_all', out_name, row.filter = paste0('complete.cases(',out_name,')'), async = FALSE, datasources = opals[this.cohort])
  #age > 35
#  dssSubset('comp_dm', 'comp_dm_all', row.filter = 'AGE > 35', async = FALSE, datasources = opals[this.cohort])
  
  # join
  dssJoin(c('comp_lb', 'comp_vs', 'comp_dm_all'), 'join_all', by = 'SUBJID', join.type = 'inner', async = FALSE, datasources = opals[this.cohort])
  
  
  
  dssSubset('kmeans_input', 'join_all', col.filter = '!(colnames(join_all) %in% c("SUBJID", "SEX"))', async = FALSE, datasources = opals[this.cohort])
#  dssSubsetByClass('join_all', subsets = 'kmeans_by_sex', variables = 'join_all$SEX', keep.cols = eval(parse(text = needed_cols)), async = FALSE, datasources = opals[this.cohort])
#  dssSubset('kmeans_M', 'kmeans_by_sex$M', row.filter = '1==1', async = FALSE, datasources = opals[this.cohort])
#  dssSubset('kmeans_F', 'kmeans_by_sex$F', row.filter = '1==1', async = FALSE, datasources = opals[this.cohort])
  
  vars <- sapply(ds.colnames('kmeans_input', datasources = opals[this.cohort])[[this.cohort]], function(x){
    this.col <- paste0('kmeans_input$',x)
    ds.var(this.col, datasources = opals[this.cohort])
  }, simplify = FALSE)
  
  fivesds <- sapply(vars, function(x){
    5*sqrt(x$Variance.by.Study[1,1])
    #x$Variance.by.Study[1,1]
  }, simplify = FALSE) %>% unlist
  
  colm <- dssColMeans('kmeans_input', datasources = opals[this.cohort])[[this.cohort]]$means
  
  names(fivesds) <- names(colm)
  
  sapply(names(colm), function(x){
    maxx <- colm[x] + fivesds[x]
    minx <- colm[x] - fivesds[x]
    my.filter <- paste0(x, ' > ', minx, ' & ', x, ' < ', maxx)
    
    dssSubset('kmeans_input', 'kmeans_input', row.filter = my.filter, async = FALSE, datasources = opals[this.cohort])
  })
  
}


prepare_cluster_godarts <- function(with.hdl = TRUE){
  this.cohort <- 'godarts'
  if(with.hdl){
    lb_filter_suffix <- "c('CPEPTIDE','HDL')"
  } else {
    lb_filter_suffix <-  "c('CPEPTIDE')"
  }
  dssSubset('mh2', 'mh', row.filter = "(MHTERM %in% c('TYPE 2 DIABETES')) & MHOCCUR == 'Y' ", async = FALSE, datasources = opals[this.cohort])
  rowfilter_godarts <- paste0("(SUBJID %in% mh2$SUBJID) & (LBTESTCD=='HBA1C' & LBMETHOD=='IFCC' & LBORRESU=='mmol/mol') | (LBTESTCD %in% ", lb_filter_suffix, ")")
  dssSubset('lb2', 'lb', row.filter = rowfilter_godarts, async = FALSE, datasources = opals[this.cohort] )
  dssSubsetByClass('lb2', subsets = 'by_test', variables = 'lb2$LBTESTCD', async = FALSE, datasources = opals[this.cohort])
 available_tests <- ds.names('by_test', datasources = opals[this.cohort])[[this.cohort]]
  #final_name <- 'tests_around_diag'
  #dssSubset(final_name, paste0('by_test$',available_tests[1]), row.filter = '1==2', col.filter = 'c("SUBJID","LBORRES", "LBORRESU", "LBTESTCD", "VISIT")', async = FALSE, datasources = opals[this.cohort])
  sapply(available_tests, function(this.test){
    ordered_name <- paste0('by_test_', this.test)
    out_name <- paste0(this.test, '_around_diag')
    dssSubset(ordered_name, paste0('by_test$', this.test), row.filter = 'order(SUBJID)', async = FALSE, datasources = opals[this.cohort])
    
    
    row_filter_block <- paste0('{
                               
                               Reduce(c,lapply(unique(SUBJID), function(my_id){
                               my_dates <- ', ordered_name,  '[', ordered_name, '$SUBJID == my_id, "LBDTC"]; 
                               my_dates <- as.Date(as.vector(my_dates));
                               diag_date <- mh2[mh2$SUBJID == my_id, "MHDTC"];
                               diag_date <- as.Date(as.vector(diag_date));
                               if(length(diag_date) == 0 || is.na(diag_date)){
                               out <- c(TRUE, rep(FALSE, length(my_dates) - 1))
                               } else { 
                               out <- (my_dates == my_dates[which.min(abs(my_dates - diag_date))]);
                               }
                               out
                               }))
                               
  }')
  
   # dssSubset(out_name, ordered_name, row.filter = row_filter_block, col.filter = 'c("SUBJID","LBORRES", "LBORRESU", "LBTESTCD")', async = FALSE, datasources = opals[this.cohort]) 
    dssSubset(out_name, ordered_name, row.filter = row_filter_block, col.filter = 'c("SUBJID","LBORRES")', async = FALSE, datasources = opals[this.cohort]) 
    dssColNames(out_name, value = c('SUBJID', this.test), async = FALSE, datasources = opals[this.cohort])
    dssSubset(out_name, out_name, row.filter = '!duplicated(SUBJID)', datasources = opals[this.cohort])
    #dssRbind(final_name, final_name, out_name, async = FALSE, datasources = opals[this.cohort])
    
})
  
  dssJoin(paste(available_tests, '_around_diag', sep  = ''), 'comp_lb', join.type = 'inner', async = FALSE, datasources = opals[this.cohort])
  
  #VS2
  dssSubset('vs2', 'vs', row.filter = "VSTESTCD == 'BMI'" , col.filter = "c('SUBJID', 'VSORRES', 'VSDTC')", async = FALSE, datasources = opals[this.cohort])
  dssSubset('vs2', 'vs2', row.filter = '!is.na(VSORRES)', async = FALSE, datasources = opals[this.cohort])
  ordered_name <- 'by_test_BMI'
  out_name <- paste0('BMI_around_diag')
  dssSubset(ordered_name, 'vs2', row.filter = 'order(SUBJID)', async = FALSE, datasources = opals[this.cohort])
  
  
  row_filter_block <- paste0('{
                             
                             Reduce(c,lapply(unique(SUBJID), function(my_id){
                             my_dates <- ', ordered_name,  '[', ordered_name, '$SUBJID == my_id, "VSDTC"]; 
                             my_dates <- as.Date(as.vector(my_dates));
                             diag_date <- mh2[mh2$SUBJID == my_id, "MHDTC"];
                             diag_date <- as.Date(as.vector(diag_date));
                             if(length(diag_date) == 0 || is.na(diag_date)){
                             out <- c(TRUE, rep(FALSE, length(my_dates) - 1))
                             } else { 
                             out <- (my_dates == my_dates[which.min(abs(my_dates - diag_date))]);
                             }
                             return(out)
                             }))
                             
                             }')
  
  dssSubset(out_name, ordered_name, row.filter = row_filter_block, col.filter = 'c("SUBJID","VSORRES")', async = FALSE, datasources = opals[this.cohort]) 
  dssSubset('comp_vs', out_name, row.filter = paste0('complete.cases(',out_name,')'), async = FALSE, datasources = opals[this.cohort])
  dssSubset('comp_vs', 'comp_vs', row.filter ='SUBJID %in% comp_lb$SUBJID', async = FALSE, datasources = opals[this.cohort])
  dssSubset('comp_vs', 'comp_vs', row.filter = '!duplicated(SUBJID)', async = FALSE, datasources = opals[this.cohort])
  
  dssColNames('comp_vs', value = c('SUBJID', 'BMI'), async = FALSE, datasources = opals[this.cohort])
  
  #DM
  
  dssSubset('dm2', 'dm', col.filter = "c('SUBJID', 'AGE', 'SEX')", async = FALSE, datasources = opals[this.cohort])
  dssSubset('comp_dm', 'dm2', row.filter = "SUBJID %in% comp_vs$SUBJID", async = FALSE, datasources = opals[this.cohort])
  # join
  dssJoin(c('comp_lb', 'comp_vs', 'comp_dm'), 'join_all', by = 'SUBJID', join.type = 'inner', async = FALSE, datasources = opals[this.cohort])
  
  dssSubset('kmeans_input', 'join_all', col.filter = '!(colnames(join_all) %in% c("SUBJID", "SEX"))', async = FALSE, datasources = opals[this.cohort])
  #  dssSubsetByClass('join_all', subsets = 'kmeans_by_sex', variables = 'join_all$SEX', keep.cols = eval(parse(text = needed_cols)), async = FALSE, datasources = opals[this.cohort])
  #  dssSubset('kmeans_M', 'kmeans_by_sex$M', row.filter = '1==1', async = FALSE, datasources = opals[this.cohort])
  #  dssSubset('kmeans_F', 'kmeans_by_sex$F', row.filter = '1==1', async = FALSE, datasources = opals[this.cohort])
  
  vars <- sapply(ds.colnames('kmeans_input', datasources = opals[this.cohort])[[this.cohort]], function(x){
    this.col <- paste0('kmeans_input$',x)
    ds.var(this.col, datasources = opals[this.cohort])
  }, simplify = FALSE)
  
  fivesds <- sapply(vars, function(x){
    5*sqrt(x$Variance.by.Study[1,1])
    #x$Variance.by.Study[1,1]
  }, simplify = FALSE) %>% unlist
  
  colm <- dssColMeans('kmeans_input', datasources = opals[this.cohort])[[this.cohort]]$means
  
  names(fivesds) <- names(colm)
  
  sapply(names(colm), function(x){
    maxx <- colm[x] + fivesds[x]
    minx <- colm[x] - fivesds[x]
    my.filter <- paste0(x, ' > ', minx, ' & ', x, ' < ', maxx)
    print(my.filter)
    dssSubset('kmeans_input', 'kmeans_input', row.filter = my.filter, async = FALSE, datasources = opals[this.cohort])
  })
  
}


prepare_cluster_abos <- function(with.hdl = TRUE){
  this.cohort <- 'abos'
  if(with.hdl){
    lb_filter_suffix <- "c('CPEPTIDE','HDL')"
  } else {
    lb_filter_suffix <-  "c('CPEPTIDE')"
  }
  
  
#  dssShowFactors('lb')
#  dssSubset('hba', 'lb', row.filter = "LBTESTCD=='HBA1C'")
#  dssShowFactors('hba')
#  dssSubsetByClass('hba', 'hba_vis', 'hba$VISIT')
#  ds.summary('hba_vis$BASELINE')
#  dssIsUnique('hba_vis$BASELINE$SUBJID') ## YESSS
  
  
#  dssSubset('hdl', 'lb', row.filter = "LBTESTCD=='HDL'")
#  dssShowFactors('hdl')
#  dssSubsetByClass('hdl', 'hdl_vis', 'hdl$VISIT')
#  ds.summary('hdl_vis$BASELINE')
#  dssIsUnique('hdl_vis$BASELINE$SUBJID') ## yess
  
  
  
#  dssSubset('cpep', 'lb', row.filter = "LBTESTCD=='CPEPTIDE'")
#  dssShowFactors('cpep')
#  
#  dssSubsetByClass('cpep', 'cpep_vis', 'cpep$VISIT')
#  ds.summary('cpep_vis$BASELINE')
#  dssSubset('cpep_0', 'cpep_vis$BASELINE', row.filter = 'LBTPT ==  "0 min"', async = FALSE)
#  ds.summary('cpep_0')
#  dssIsUnique('cpep_0$SUBJID') #yes
  
  if(with.hdl){
    rowfilter_abos <- "(LBTESTCD=='HBA1C' | (LBTESTCD == 'CPEPTIDE' & LBTPT ==  '0 min') | LBTESTCD == 'HDL') & VISIT == 'BASELINE'"
  } else {
    rowfilter_abos <- "(LBTESTCD=='HBA1C' | (LBTESTCD == 'CPEPTIDE' & LBTPT ==  '0 min')) & VISIT == 'BASELINE' "
  }
  
  ### LB
  dssSubset('lb2', 'lb', row.filter = rowfilter_abos, col.filter = "c('SUBJID', 'LBORRES', 'LBTESTCD')",  async = FALSE, datasources = opals[this.cohort] )
  dssSubsetByClass('lb2', subsets = 'by_test', variables = 'lb2$LBTESTCD', async = FALSE, datasources = opals[this.cohort])
  
  available_tests <- ds.names('by_test', datasources = opals[this.cohort])[[this.cohort]]
  join_source = paste('by_test$', available_tests, sep  = '')
  sapply(available_tests, function(x){
    
    dssColNames(paste0('by_test$', x), c("SUBJID",x , "LBTESTCD"), datasources = opals[this.cohort])
  })
  dssShowFactors('by_test$CPEPTIDE')

  
  dssJoin(join_source, 'comp_lb', join.type = 'inner', async = FALSE, datasources = opals[this.cohort])
  dssSubset('comp_lb', 'comp_lb', col.filter = "c('SUBJID', 'CPEPTIDE', 'HBA1C', 'HDL')",datasources = opals[this.cohort])

  #VS2
  dssSubset('comp_vs', 'vs', row.filter = "VSTESTCD == 'BMI' & VISIT == 'BASELINE'" , col.filter = "c('SUBJID', 'VSORRES')", async = FALSE, datasources = opals[this.cohort])
  
  dssSubset('comp_vs', 'comp_vs', row.filter ='SUBJID %in% comp_lb$SUBJID', async = FALSE, datasources = opals[this.cohort])
 
  
  dssColNames('comp_vs', value = c('SUBJID', 'BMI'), async = FALSE, datasources = opals[this.cohort])
  
  #DM
  
  dssSubset('comp_dm', 'dm', col.filter = "c('SUBJID', 'AGE', 'SEX')", async = FALSE, datasources = opals[this.cohort])
  # join
  dssJoin(c('comp_lb', 'comp_vs', 'comp_dm'), 'join_all', by = 'SUBJID', join.type = 'inner', async = FALSE, datasources = opals[this.cohort])
  
  dssSubset('kmeans_input', 'join_all', col.filter = '!(colnames(join_all) %in% c("SUBJID", "SEX"))', async = FALSE, datasources = opals[this.cohort])
  #  dssSubsetByClass('join_all', subsets = 'kmeans_by_sex', variables = 'join_all$SEX', keep.cols = eval(parse(text = needed_cols)), async = FALSE, datasources = opals[this.cohort])
  #  dssSubset('kmeans_M', 'kmeans_by_sex$M', row.filter = '1==1', async = FALSE, datasources = opals[this.cohort])
  #  dssSubset('kmeans_F', 'kmeans_by_sex$F', row.filter = '1==1', async = FALSE, datasources = opals[this.cohort])
  
  vars <- sapply(ds.colnames('kmeans_input', datasources = opals[this.cohort])[[this.cohort]], function(x){
    this.col <- paste0('kmeans_input$',x)
    ds.var(this.col, datasources = opals[this.cohort])
  }, simplify = FALSE)
  
  fivesds <- sapply(vars, function(x){
    5*sqrt(x$Variance.by.Study[1,1])
    #x$Variance.by.Study[1,1]
  }, simplify = FALSE) %>% unlist
  
  colm <- dssColMeans('kmeans_input', datasources = opals[this.cohort])[[this.cohort]]$means
  
  names(fivesds) <- names(colm)
  
  sapply(names(colm), function(x){
    maxx <- colm[x] + fivesds[x]
    minx <- colm[x] - fivesds[x]
    my.filter <- paste0(x, ' > ', minx, ' & ', x, ' < ', maxx)
    print(my.filter)
    dssSubset('kmeans_input', 'kmeans_input', row.filter = my.filter, async = FALSE, datasources = opals[this.cohort])
  })
  
}
prepare_cluster_accelerate <- function(with.hdl = TRUE){
  this.cohort <- 'accelerate'
  if(with.hdl){
    lb_filter_suffix <- "c('CPEPTIDE','HDL')"
  } else {
    lb_filter_suffix <-  "c('CPEPTIDE')"
  }
  
  
    dssShowFactors('lb')
    dssSubset('hba', 'lb', row.filter = "LBTESTCD=='HBA1C'")
    dssShowFactors('hba')
    
    dssSubsetByClass('hba', 'hba_vis', 'hba$VISIT')
    ds.summary('hba_vis$SCREENI')
    ds.names('hba_vis')
    ds.summary('lb')
  #  dssIsUnique('hba_vis$BASELINE$SUBJID') ## YESSS
  
  
  #  dssSubset('hdl', 'lb', row.filter = "LBTESTCD=='HDL'")
  #  dssShowFactors('hdl')
  #  dssSubsetByClass('hdl', 'hdl_vis', 'hdl$VISIT')
  #  ds.summary('hdl_vis$BASELINE')
  #  dssIsUnique('hdl_vis$BASELINE$SUBJID') ## yess
  
  
  
  #  dssSubset('cpep', 'lb', row.filter = "LBTESTCD=='CPEPTIDE'")
  #  dssShowFactors('cpep')
  #  
  #  dssSubsetByClass('cpep', 'cpep_vis', 'cpep$VISIT')
  #  ds.summary('cpep_vis$BASELINE')
  #  dssSubset('cpep_0', 'cpep_vis$BASELINE', row.filter = 'LBTPT ==  "0 min"', async = FALSE)
  #  ds.summary('cpep_0')
  #  dssIsUnique('cpep_0$SUBJID') #yes
  
  if(with.hdl){
    rowfilter_abos <- "(LBTESTCD=='HBA1C' | (LBTESTCD == 'CPEPTIDE' & LBTPT ==  '0 min') | LBTESTCD == 'HDL') & VISIT == 'BASELINE'"
  } else {
    rowfilter_abos <- "(LBTESTCD=='HBA1C' | (LBTESTCD == 'CPEPTIDE' & LBTPT ==  '0 min')) & VISIT == 'BASELINE' "
  }
  
  ### LB
  dssSubset('lb2', 'lb', row.filter = rowfilter_abos, col.filter = "c('SUBJID', 'LBORRES', 'LBTESTCD')",  async = FALSE, datasources = opals[this.cohort] )
  dssSubsetByClass('lb2', subsets = 'by_test', variables = 'lb2$LBTESTCD', async = FALSE, datasources = opals[this.cohort])
  
  available_tests <- ds.names('by_test', datasources = opals[this.cohort])[[this.cohort]]
  join_source = paste('by_test$', available_tests, sep  = '')
  sapply(available_tests, function(x){
    
    dssColNames(paste0('by_test$', x), c("SUBJID",x , "LBTESTCD"), datasources = opals[this.cohort])
  })
  dssShowFactors('by_test$CPEPTIDE')
  
  
  dssJoin(join_source, 'comp_lb', join.type = 'inner', async = FALSE, datasources = opals[this.cohort])
  dssSubset('comp_lb', 'comp_lb', col.filter = "c('SUBJID', 'CPEPTIDE', 'HBA1C', 'HDL')",datasources = opals[this.cohort])
  
  #VS2
  dssSubset('comp_vs', 'vs', row.filter = "VSTESTCD == 'BMI' & VISIT == 'BASELINE'" , col.filter = "c('SUBJID', 'VSORRES')", async = FALSE, datasources = opals[this.cohort])
  
  dssSubset('comp_vs', 'comp_vs', row.filter ='SUBJID %in% comp_lb$SUBJID', async = FALSE, datasources = opals[this.cohort])
  
  
  dssColNames('comp_vs', value = c('SUBJID', 'BMI'), async = FALSE, datasources = opals[this.cohort])
  
  #DM
  
  dssSubset('comp_dm', 'dm', col.filter = "c('SUBJID', 'AGE', 'SEX')", async = FALSE, datasources = opals[this.cohort])
  # join
  dssJoin(c('comp_lb', 'comp_vs', 'comp_dm'), 'join_all', by = 'SUBJID', join.type = 'inner', async = FALSE, datasources = opals[this.cohort])
  
  dssSubset('kmeans_input', 'join_all', col.filter = '!(colnames(join_all) %in% c("SUBJID", "SEX"))', async = FALSE, datasources = opals[this.cohort])
  #  dssSubsetByClass('join_all', subsets = 'kmeans_by_sex', variables = 'join_all$SEX', keep.cols = eval(parse(text = needed_cols)), async = FALSE, datasources = opals[this.cohort])
  #  dssSubset('kmeans_M', 'kmeans_by_sex$M', row.filter = '1==1', async = FALSE, datasources = opals[this.cohort])
  #  dssSubset('kmeans_F', 'kmeans_by_sex$F', row.filter = '1==1', async = FALSE, datasources = opals[this.cohort])
  
  vars <- sapply(ds.colnames('kmeans_input', datasources = opals[this.cohort])[[this.cohort]], function(x){
    this.col <- paste0('kmeans_input$',x)
    ds.var(this.col, datasources = opals[this.cohort])
  }, simplify = FALSE)
  
  fivesds <- sapply(vars, function(x){
    5*sqrt(x$Variance.by.Study[1,1])
    #x$Variance.by.Study[1,1]
  }, simplify = FALSE) %>% unlist
  
  colm <- dssColMeans('kmeans_input', datasources = opals[this.cohort])[[this.cohort]]$means
  
  names(fivesds) <- names(colm)
  
  sapply(names(colm), function(x){
    maxx <- colm[x] + fivesds[x]
    minx <- colm[x] - fivesds[x]
    my.filter <- paste0(x, ' > ', minx, ' & ', x, ' < ', maxx)
    print(my.filter)
    dssSubset('kmeans_input', 'kmeans_input', row.filter = my.filter, async = FALSE, datasources = opals[this.cohort])
  })
  
}


do_clustering <- function(cohorts, centers,
                          input.name = 'kmeans_input', with.scaling = TRUE, iter.max = 30, nstart =30){

  if(with.scaling){
    
    kmeans.input = paste0(input.name, '_scaled')
    dssScale(kmeans.input, input.name, datasources = opals[cohorts])
  } else {
    kmeans.input <- input.name
  }
  cluster_factor <- paste0(kmeans.input, '_km_clust', centers)
  cluster_measures <- dssColNames(input.name, datasources = opals[cohorts])
  if(length(cohorts) >1){
    if(length(Reduce(setdiff, cluster_measures)) > 0 ){
      stop('The input dataframe has different columns on different nodes.')
    }
  }
  cluster_measures <- cluster_measures[[1]]
  
  ktot <- dssKmeans(kmeans.input, centers = centers, iter.max = iter.max, nstart = nstart, async = TRUE, datasources = opals[cohorts])
  #x <- ktot$good$cluster
  #perc <- round(100*x/sum(x))
  #dssSubsetByClass('kmeans_input', subsets = 'by_cluster', variables = 'kmeans_input_scaled_km_clust5', keep.cols = c('AGE', 'BMI', 'CPEPTIDE', 'HBA1C', 'HDL'))
  #dssSubsetByClass('kmeans_input', subsets = 'by_cluster', variables = 'kmeans_input_scaled_km_clust5', keep.cols = c('AGE', 'BMI', 'CPEPTIDE', 'HBA1C', 'HDL'), async = FALSE, datasources = opals[this.cohort])
  dssSubsetByClass(input.name, subsets = 'by_cluster', variables = cluster_factor, keep.cols = cluster_measures, async = TRUE, datasources = opals[cohorts])
  
 # quants <- sapply(ds.names('by_cluster', datasources = opals[this.cohort])[[this.cohort]], function(x){
#    dfname <- paste0('by_cluster$', x)
#    sapply(ds.colnames(dfname, datasources = opals[this.cohort])[[this.cohort]], function(y){
#      colname <- paste0(dfname, '$', y)
#      ds.quantileMean(colname, datasources = opals[this.cohort])
#    }, simplify = FALSE)
    
#  }, simplify = FALSE)
  
  quants <- sapply(ds.names('by_cluster', datasources = opals[cohorts])[[1]], function(x){
    dfname <- paste0('by_cluster$', x)
    sapply(cluster_measures, function(y){
      colname <- paste0(dfname, '$', y)
      ds.quantileMean(colname, datasources = opals[cohorts])
    }, simplify = FALSE)
    
  }, simplify = FALSE)
  
  
  
  
  
  #values <- c("#4DB3E6","#132B41","#E69900","#CC80B3","#8B1A4F","#00B399")
  #values <- c("#4DB3E6","#CC80B3","#132B41","#00B399","#8B1A4F","#E69900")
  
  list( 
    input = sapply(names(quants), function(this.cluster){
      t(Reduce(cbind, quants[[this.cluster]])) %>% 
        #  data.frame(measurement = names(quants[[this.cluster]]), cluster = rep(this.cluster,5) )
        data.frame(measurement = names(quants[[this.cluster]]), cluster = rep(this.cluster,length(cluster_measures)) )
    }, simplify = FALSE) %>% Reduce(rbind,.),
    cluster.object = ktot, 
    cohorts = cohorts
  )
}

do_ggplot <- function(input.list, cluster.labels =c('2/SIDD', '3/SIRD', '4/MOD', '5/MARD/low HDL', '6/MARD/high HDL'), create.pdf = TRUE){
  
  cohorts <- input.list$cohorts
  pre_gginput <- input.list$input[, c('X5.', 'X25.', 'X50.', 'X75.', 'X95.', 'measurement', 'cluster')] 
  pre_gginput$newclust <- rep(NA,nrow(pre_gginput))
  legend.labels <- cluster.labels
  if(!is.null(cluster.labels)){
    
    for(cl.name in sort(levels(pre_gginput$cluster))){
      pre_gginput[pre_gginput$cluster == cl.name,'newclust'] <- cluster.labels[1]
      cluster.labels <- cluster.labels[-1]
    }
  } else {
    pre_gginput$newclust <- pre_gginput$cluster
  }
  #pre_gginput[pre_gginput$cluster == 'X1','newclust'] <-'5/MARD'
  #pre_gginput[pre_gginput$cluster == 'X3','newclust'] <-'2/SIDD'  
  #pre_gginput[pre_gginput$cluster == 'X2','newclust'] <- '3/SIRD'
  #pre_gginput[pre_gginput$cluster == 'X5','newclust'] <- '5/MARD/low HDL' 
  #pre_gginput[pre_gginput$cluster == 'X4','newclust'] <- '4/MOD'  
  pre_gginput$newclust <- factor(pre_gginput$newclust)
  gginput <- melt(pre_gginput)
  values <- c("#4DB3E6","#00B399","#E69900","#CC80B3","#8B1A4F", "#132B41")
  by_measure <- ggplot(gginput,aes(x = newclust, y = value, fill = newclust)) +
    geom_boxplot() +
    #  facet_wrap(~measurement, scale="free", ncol=6) +
    facet_wrap(~measurement, scale="free", ncol=length(levels(pre_gginput$measurement))) +
    scale_fill_manual(values=values) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))+
    ylab("Level")+
    xlab("Cluster") + ggtitle(paste0("Combined males and females, cohorts: ", paste(cohorts, collapse = ', '), ", N = ", sum(input.list$cluster.object$good$cluster)))
  
  
  xxx <- sapply(levels(pre_gginput$measurement), function(this.level){
    xx <- pre_gginput[pre_gginput$measurement == this.level,]
    r <- c(min(xx[,1]), max(xx[,5]))
    out <- sapply(xx[,1:5],  function(line){
      rescale(line,c(1,100), r)
    })
    data.frame(out, xx[,c(6,7,8)])
  }, simplify = FALSE)

  yyy <- Reduce(rbind, xxx) %>% melt
  n <- input.list$cluster.object$global$size
  names(n) <- paste('X', names(n), sep = '')
  cluster_labeller <- function(val){
  sapply(val, function(a){
      cl <- pre_gginput[pre_gginput$newclust == a, 'cluster'][1]
      paste0(a, ', N = ', n[cl])
    })
  }
  
  by_cluster <- ggplot(yyy,aes(x = measurement, y = value, fill = newclust)) +
    geom_boxplot() +
    facet_wrap(~newclust, scales="fixed", ncol=length(levels(pre_gginput$measurement))) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))+
    theme(legend.position="none") +
    scale_fill_manual(values=values) +
    ggtitle("Combined males and females grouped by cluster")+
    ylab("Level")+
    xlab("Measurement")
  
  dfclust <- as.data.frame(unclass(input.list$cluster.object$global$size))
  names(dfclust) <- 'freq'
  dfclust$clust <- factor(legend.labels)
  dfclust$newy <- (dfclust$freq/2 + c(0, cumsum(dfclust$freq)[-length(dfclust$freq)]))
  dfclust$value <- round(dfclust$freq / sum(dfclust$freq),2)
  p2 <-ggplot(dfclust, aes(x=clust, y = value, fill = clust)) +
    #geom_text(aes(y = freq, label = paste0(freq,"\n", percent(value))), size=4)+
    geom_bar(width = 1, stat='identity')+
    #coord_polar("y", start=0, direction = -1)+
    #theme(axis.text.x=element_blank()) +
    geom_text(aes( label = paste0(freq,"\n", percent(value))), size=4, colour = 'white', hjust = 1.6)+  
    theme_minimal() +
    theme(legend.position="none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
    scale_fill_manual(values=values) + coord_flip()
  nclust <- length(levels(pre_gginput$cluster))
  if(create.pdf){
    pdf(paste0(nclust, ' Clusters combined ', paste(cohorts, collapse = ', '),".pdf"), height=12, width=20)
    pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 6)))
    plot(by_measure, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4))
    plot(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 5:6))
    plot(by_cluster, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:4))
    dev.off()
  }
  invisible(list(by_measure = by_measure, by_cluster = by_cluster, percentages = p2))
}

do_pca <- function( cluster_obj, input_df = 'kmeans_input_scaled'){
  cohorts <- cluster_obj$cohorts
  k <- length(levels(cluster_obj$inp$cluster))
  smooth_categories <- paste0(input_df, '_km_clust', k)
  input_scores <- paste0(input_df, '_scores')
  pca <- dssPrincomp(input_df, type = 'combine', datasources = opals[cohorts])
  par(mfrow = c(1,1))
  plot.eig.percs(pca$global, max.pcs = 5)
  cols <- ds.summary(input_scores, datasources = opals[cohorts])[[1]]$`variables held`

  sapply(combn(cols,2, simplify = FALSE), function(this.combo){
    this.x <- paste0(input_scores, '$', this.combo[1])
    this.y <- paste0(input_scores, '$', this.combo[2])
    par(mfrow = c(3,2))
    for(i in c(1:k)){
      my.smooth2d(this.x, this.y, categories = smooth_categories, 
                  emphasize_level = i, draw.image = TRUE,  type = 'combine', async = FALSE, datasources = opals[cohorts])
    }
  })
  par(mfrow = c(1,1))
}

plot.eig.percs <- function (res, main="Principal Component Scree Plot", sub="", max.pcs=NULL) {
  eig.percs <- res$sdev^2/sum(res$sdev^2)
  cutoff <- 1/length(eig.percs) # theo avg inertia per PC
  eig.percs2 <- c()
  if (!is.null(max.pcs)) {
    if (max.pcs>=length(eig.percs)) max.pcs <- length(eig.percs) - 1
    eig.percs2 <- eig.percs[ (max.pcs+1):length(eig.percs)]
    eig.percs <- eig.percs[ 1:max.pcs ]
    eig.percs2 <- sum(eig.percs2)
    if (eig.percs2>eig.percs[1]) eig.percs2 <- eig.percs[1]
  }
  plot( c(eig.percs, eig.percs2), type="b", ylim=c(0,max(eig.percs, eig.percs2)),
        pch=c(rep(21, length(eig.percs)), ifelse(eig.percs2==eig.percs[1],24,23)), 
        col="black", bg=c(ifelse(eig.percs>cutoff,"blue",NA), NA), 
        xaxt="n", xlab="", ylab="Proportion of variance",
        main=main, sub=sub)
  abline( h=cutoff, lty="dashed", col="purple" )
  labels <- paste("PC",1:length(eig.percs), sep=" ")
  if (length(eig.percs2)>0) {
    labels <- c(labels , paste0("PCs ",max.pcs,"+"))
  } 
  axis(1, at=1:length(c(eig.percs, eig.percs2)), 
       labels = labels,  las=3)

}
#ds2.unmemoise(.my.bw.args)
#ds2.unmemoise(dssRange)

#ds2.memoise(.my.bw.args, debug = TRUE)
#ds2.memoise(dssRange, debug = TRUE)



