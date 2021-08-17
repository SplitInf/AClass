##### AClass #####

## To Do ##
# Fix issue with TYPE. Use name of object instead?
# hardcode Probes rank path and probes_rank (save all in "Training_Results")
# add nano.set.train.settings similar to nano.set.colour so that nano.train can be broken down to smaller
# set subgroup specific threshold?

##### all functions #####

###### initialize.prj ######
# set up paths and ensure all libraries are present

initialize.prj <- function() {
  
  print(paste0("[MSG] Checking for libraries..."))
  # look for all packages #
  
  # if(!require(installr)){ # if Rtools is missing
  #   install.packages("installr")
  #   library(installr)
  # }
  
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  
  if(!require(tcltk)){ # required for choosing folder
    install.packages("tcltk")
    library(tcltk)
  }
  
  if(!require(rlang)){ # ggplot2
    install.packages("rlang")
    library(rlang)
  }
  
  if(!require(vctrs)){ # ggplot2
    install.packages("vctrs")
    library(vctrs)
  }
  
  if(!require(pillar)){ # ggplot2
    install.packages("pillar")
    library(pillar)
  }
  
  if(!require(scales)){ # ggpubr
    install.packages("scales")
    library(scales)
  }
  
  if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)
  }
  
  
  if(!require(gower)){
    install.packages("gower") # caret
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("gower")
    library(gower)
  }
  
  
  if(!require(data.table)){
    #install.packages("data.table") #caret
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("data.table")
    BiocManager::install(c("data.table"))
    library(data.table)
  }
  
  if(!require(ipred)){
    install.packages("ipred") #caret
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("ipred")
    library(ipred)
  }
  
  if(!require(caret)){
    BiocManager::install(c("caret"))
    library(caret)
  }
  
  if(!require(randomForest)){
    BiocManager::install(c("randomForest"))
    library(randomForest)
  }
  
  if(!require(glmnet)){
    BiocManager::install(c("glmnet"))
    library(glmnet)
  }
  
  if(!require(pamr)){
    BiocManager::install(c("pamr"))
    library(pamr)
  }
  
  if(!require(klaR)){
    BiocManager::install(c("klaR"))
    library(klaR)
  }
  
  if(!require(reshape2)){
    install.packages("reshape2")
    library(reshape2)
  }
  
  if(!require(ggrepel)){
    install.packages("ggrepel")
    library(ggrepel)
  }
  
  if(!require(grid)){
    install.packages("grid")
    library(grid)
  }
  
  if(!require(gridExtra)){
    install.packages("gridExtra")
    library(gridExtra)
  }

  if (!require(vsn)){
    BiocManager::install(c("vsn"))
    library("vsn")
  }
  
  if(!require(XML)){ # NanostringNorm
    install.packages("XML", type = "binary")
    library(XML)
  }
  
  if(!require(NanoStringNorm)){
    # install.packages("NanoStringNorm")
    #BiocManager::install(c("NanoStringNorm"))
    install.packages('devtools')
    devtools::install_url('https://cran.r-project.org/src/contrib/Archive/NanoStringNorm/NanoStringNorm_1.2.1.tar.gz')
    library(NanoStringNorm)
  }
  
  if (!require(ggpubr)){
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("ggpubr")
    install.packages("colorspace")
    library(colorspace)
    install.packages("lazyeval")
    library(lazyeval)
    install.packages("ggpubr")
    library("ggpubr")
  }
  
  if(!require(cowplot)){
    install.packages("cowplot")
    library(cowplot)
  }
  
  if(!require(ResourceSelection)){
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("ResourceSelection")
    #install.packages("ResourceSelection")
    BiocManager::install(c("ResourceSelection"))
    library(ResourceSelection)
  }
  
  if(!require(Boruta)){
    #install.packages("Boruta")
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("ResourceSelection")
    BiocManager::install(c("Boruta"))
    library(Boruta)
  }
  if(!require(limma)){
    install.packages("limma")
    library(limma)
  }
  if(!require(plotly)){
    install.packages("plotly")
    library(plotly)
  }
  if(!require(openxlsx)){
    #install.packages("openxlsx")
    source("https://bioconductor.org/biocLite.R")
    biocLite("openxlsx")
    library(openxlsx)
  }

  print(paste0("[MSG] All libraries loaded."))
}

# process.raw #
# running nano.load(), nano.prenorm.qc(), nano.prep(), nano.norm() and nano.MDS()
# work_path - where output directory will be created
# raw_path - path to raw data
# keep_file_path - path to the file that contain the list of samples to be included (one sample per row, no header)
# omit_file_path - path to the file that contain the list of samples not to be included (one sample per row, no header)
# prefix - project prefix that goes to file output(s) *optional*
# SampleContent - how nano.norm() handles normalization. Default is "housekeeping.geo.mean"

process.raw <- function(work_path=getwd(), raw_path=NULL, keep_file_path=NULL, omit_file_path=NULL, prefix=NULL, SampleContent = "housekeeping.geo.mean"){
  
  # check paths #
  if(is.null(raw_path)){stop("[MSG] raw_path missing")}

  if(is.null(omit_file_path)) {omit_file_path=""}
  if(is.null(keep_file_path)) {keep_file_path=""}

  print(paste0("[MSG] work_path: ",work_path, " raw_path: ", raw_path, " keep_file_path: ", keep_file_path, " omit_file_path: ", omit_file_path, " prefix: ",prefix))
  raw.obj <- "" # raw.obj
  
  # create prefix folder #
  if(!is.null(prefix)){
    prj_prefix <- paste(prefix,format(Sys.time(), "%Y%m%d-%H%M"), sep = "_")
  } else {
    prj_prefix <- format(Sys.time(), "%Y%m%d-%H%M")
  }
  out_path <- paste0(work_path,"/",prj_prefix)
  dir.create(out_path, showWarnings = FALSE)

  # [1] load data
  print(paste0("=nano.load="))
  raw.obj <- nano.load(raw_path = raw_path, keep_file_path=keep_file_path, omit_file_path=omit_file_path)
  
  ## check if samples loaded ##
  raw_n <- ncol(raw.obj$raw)-3

  raw.obj$run_info$run_id <- prj_prefix
  
  # [2] Prenorm qc
  if(raw_n > 0){
    print(paste0("=nano.prenorm.qc="))
    raw.obj <- nano.prenorm.qc(data=raw.obj, prefix = prefix, out_path = out_path)
  }
  # [3] Normalization
  if(raw_n > 1){
    print(paste0("=nano.norm="))
    raw.obj <- nano.norm(data=raw.obj, SampleContent = SampleContent) 
  }
  ## [4] prep data
  if(raw_n > 0){
    print(paste0("=nano.prep="))
    raw.obj <- nano.prep(data=raw.obj)# transpose and check for zero's
  }

  return(raw.obj)
  
}

##### Batch Process Raw #####
# batch mode of process.raw() where the end obj comes from combining obj from each data dir
# all parameters matches with process.raw() except:
# raw_dir_path  expects multiple raw folders stored within this path.

batch.process.raw <- function(work_path=getwd(), raw_dir_path=raw_dir_path, keep_file_path=NULL, omit_file_path=NULL, prefix=NULL, SampleContent = "housekeeping.geo.mean"){
  
  # get path to csv's
  raw_dir <- as.data.frame(list.files(path = raw_dir_path, pattern = ".*NormalizedData.*.csv", recursive = TRUE, full.names=TRUE))
  raw_dir <- apply(raw_dir,1,function(x) unlist(strsplit(gsub(paste0(raw_dir_path,"/"),"",x),split = "/"))[1] )
  raw_dir <- unique(raw_dir)
  
  raw.obj <- list()
  for (i in raw_dir) {
    print(paste0("[MGS] Processing dir: ", i))
    raw_path_i <- paste0(raw_dir_path,"/",i)
    train.i <- process.raw(work_path=work_path, raw_path=raw_path_i, keep_file_path=keep_file_path, omit_file_path=omit_file_path, prefix=prefix, SampleContent = SampleContent)
    ### merging results from multiple batches ###

    ## run_info ##
    # csv #
    if(is.null(raw.obj$run_info$csv)) {
      raw.obj$run_info$csv <- train.i$run_info$csv
    } else {
      raw.obj$run_info$csv <- c(raw.obj$run_info$csv,train.i$run_info$csv)
    }
    # samples_found #
    if(is.null(raw.obj$run_info$samples_found)) {
      raw.obj$run_info$samples_found <- train.i$run_info$samples_found
    } else {
      samples_found.i <- as.numeric(unlist(strsplit(train.i$run_info$samples_found,split = " "))[1])
      samples_found.full <- as.numeric(unlist(strsplit(raw.obj$run_info$samples_found,split = " "))[1])
      raw.obj$run_info$samples_found <- paste0(samples_found.i+samples_found.full, " samples found.")
    }
    # samples_loaded #
    if(is.null(raw.obj$run_info$samples_loaded)) {
      raw.obj$run_info$samples_loaded <- train.i$run_info$samples_loaded
    } else {
      samples_loaded.i <- as.numeric(unlist(strsplit(train.i$run_info$samples_loaded,split = " "))[1])
      samples_loaded.full <- as.numeric(unlist(strsplit(raw.obj$run_info$samples_loaded,split = " "))[1])
      raw.obj$run_info$samples_loaded <- paste0(samples_loaded.i+samples_loaded.full, " samples loaded.")
    }
    ## run_id ##
    if(is.null(raw.obj$run_info$run_id)) {
      raw.obj$run_info$run_id <- train.i$run_info$run_id
    } else {
      raw.obj$run_info$run_id <- c(raw.obj$run_info$run_id,train.i$run_info$run_id)
    }
    ## $raw, $prenorm_qc $norm, $norm.t
    merge_df_by_rowname <- function(a,b,df){
      print(a)
      if(df == "raw"){
        header <- c("Code.Class", "Name", "Accession")
        a[[df]] <- merge(a[[df]] , b[[df]], by=header)
      } else if(df != "raw"){
        a[[df]] <- merge(a[[df]] , b[[df]], by="row.names")
        row.names(a[[df]]) <- a[[df]][,1]
        a[[df]] <- a[[df]][,-1]
      }
      return(a[[df]])
    }
    
    if(is.null(raw.obj$raw)){
      raw.obj$raw <- train.i$raw
    } else {
      raw.obj$raw <- merge_df_by_rowname(a=raw.obj, b=train.i, df="raw")
    }
    
    if(is.null(raw.obj$prenorm_qc)){
      raw.obj$prenorm_qc <- train.i$prenorm_qc
    } else {
      raw.obj$prenorm_qc <- rbind(raw.obj$prenorm_qc, train.i$prenorm_qc)
    }
    
    if(is.null(raw.obj$norm)){
      raw.obj$norm <- train.i$norm
    } else {
      raw.obj$norm <- merge_df_by_rowname(a=raw.obj, b=train.i, df="norm")
    }
    
    raw.obj$norm.t <- as.data.frame(t(raw.obj$norm))
    
  }
  
  return(raw.obj)
}

##### classify.data #####
# work_path
classify.data <- function(work_path=getwd(), data, prefix, training_model_obj, alg_list = c("rf","glmnet","pam", "nb", "knn"), keep_file_path = NULL, omit_file_path = NULL, out_path=NULL){
  
  test.obj <- data
  
  # default to use run_id as out_path unless otherwise provided
  if (!is.null(test.obj$run_info$run_id)){
    prj_prefix <- test.obj$run_info$run_id[1]
    out_path <- paste0(work_path,"/",prj_prefix)
  }
  
  if(is.null(out_path)){
    out_path = paste(getwd(),data$run_info$run_id[1], sep = "/")
  }
  
  if (!is.null(keep_file_path)) {
      test.keep <- read.table(file = keep_file_path, stringsAsFactors = FALSE, sep = "\t")
      test.keep <- apply(test.keep, 1, make.names)
      print(paste0("[MSG] Keeping ",length(test.keep), " samples"))
      test$norm.t <- test$norm.t[row.names(test.obj$norm.t) %in% test.keep,,drop=FALSE]
  }
  
  if (!is.null(omit_file_path)) {
      test.omit <- read.table(file = omit_file_path, stringsAsFactors = FALSE, sep = "\t") # files omitted
      test.omit <- apply(test.omit, 1, make.names)
      print(paste0("[MSG] Omitting ",length(test.omit), " samples"))
      test.obj$norm.t <- test$norm.t[!row.names(test.obj$norm.t) %in% test.omit,,drop=FALSE]
  }

  # [6] Test
  # choose algorithms
  test.obj <- nano.test(prefix = prefix, training_model_obj = training_model_obj, data = test.obj, alg_list = alg_list, out_path = out_path) # output text file to out_path
  
  # [7] Consolidate results
  # choose min max range based on model accuracy
  test.obj <- get.nano.test.results(prefix,test.obj, out_path = out_path)
  saveRDS(test.obj, file = paste0(out_path,"/",prefix,"_test.data.tested.RDS"))
  
  # [8] Set Colour Code based on pre-trained models
  # group <- c("Group1","Group2A","Group2B")
  # group <- c("Group1","Group2")
  group <- unique(training_model_obj$train.data.main$Group)
  test.obj <- nano.set.colour(test.obj, group)
  
  # [9] Generate report
  test.obj <- nano.plot(prefix = prefix, data = test.obj, prob= "Avg_Probability", report_type="Summary", print_report = TRUE, thres_avg_prob=0, thres_geomean = 100, out_path=out_path)
  
  return(test.obj)
}


##### Load Nanostring Data #####
# v3 include keep_file_path for choosing files to include
# v2 fixes bug that crashes program when no samples were loaded
#' Looks for Nanostring Data with "NormalizedData.csv" and load as dataframe
#' @param raw_path path of the directory where raw nanostring data is kept. Read recursively.
#' @param omit_file_path path and tab-delimited list of sample names to be omitted. One per row.
#' @return dataframe of read data and number of samples loaded.
#' @example
#' test.omit <- read.table(file = omit_file_path, stringsAsFactors = FALSE, sep = "\t")
#' head(test.omit)
# V1
# 1 COG15-PAUDKK FFPE
# 2 COG16-PAUFJA FFPE
# 3 COG27-PAUPNX FFPE
# 4 COG38-PAVBCB FFPE
# 5 COG41-PAVDHZ FFPE
# 6 COG43-PAVEDH FFPE
#' test.raw <- nano.load(Testing_data_path, test_omit_path)

nano.load <- function(raw_path = getwd(), keep_file_path="", omit_file_path="") {
  raw.summary <- list()
  raw.merge <- data.frame()
  for (nanofile in list.files(path=raw_path, pattern = ".*NormalizedData.*.csv", recursive = FALSE)){

    raw <- read.table(paste(raw_path,nanofile, sep = "/"), sep = ",", skip = 15, stringsAsFactors = FALSE) # skip header
    header <- c("Code.Class", "Name", "Accession")
    info <- read.table(paste(raw_path,nanofile, sep = "/"), sep = ",", skip = 2, stringsAsFactors = FALSE) #add column info in col1-4
    Sample_names <- info[1,4:ncol(info)] # get sample names
    colnames(raw) <- c(header,Sample_names)

    raw.summary[["csv"]][[nanofile]] <- paste0(dim(raw)[1]," features ", dim(raw)[2]-3, " samples.")
    
    print(paste0("[MSG] Loading ",raw.summary[["csv"]][[nanofile]]))
    
    if (ncol(raw.merge)==0) {
      raw.merge <- raw
    } else {
      raw.merge <- merge(raw.merge, raw, by = header)
    }
  }
  
  n_sample <- ifelse(ncol(raw.merge)-3 <0,0,ncol(raw.merge)-3)
  raw.summary[["samples_found"]] <- paste0(n_sample, " samples found.")
  
  print(paste0("[MSG] ",raw.summary[["samples_found"]]))

  ##### Choose Samples #####
  ## fix names ##
  colnames(raw.merge) <- make.names(colnames(raw.merge))

  if (keep_file_path != "") {
    test.keep <- read.table(file = keep_file_path, stringsAsFactors = FALSE, sep = "\t")
    test.keep <- apply(test.keep, 1, make.names)
    test.keep <- c("Code.Class","Name","Accession",test.keep)
    raw.merge <- raw.merge[,colnames(raw.merge) %in% test.keep]
  }
  
  if (omit_file_path != "") {
    test.omit <- read.table(file = omit_file_path, stringsAsFactors = FALSE, sep = "\t") # files omitted
    test.omit <- apply(test.omit, 1, make.names)
    raw.merge <- raw.merge[,!colnames(raw.merge) %in% test.omit]
  }
  
  n_sample_loaded <- ifelse(ncol(raw.merge)-3 <0,0,ncol(raw.merge)-3)
  raw.summary[["samples_loaded"]] <- paste0(n_sample_loaded, " samples loaded.")
  print(paste0("[MSG] ",raw.summary[["samples_loaded"]])) #7
  raw.loaded <- list()

  raw.loaded$run_info <- raw.summary
  raw.loaded$raw <- raw.merge
  if (n_sample_loaded == 0){
    stop("[MSG] No samples were loaded. Check file omit list.")
  }
  return(raw.loaded)
}


##### Prenorm Data QC #####
#' Prenorm Nanostring Data and produce report. Required to calculate HK_geomean
#' @param raw_data Raw dataframe from nano.load.
#' @param code_class Default "Housekeeping".
#' @param prefix project prefix
#' @return Report of the raw data and GeoMean table.
#' @example
#' library(reshape2)
#' library(ggplot2)
#' library(ggrepel)
#' library(gridExtra)
#' library(grid) # dependencies of gridExtra
#' nano.prenorm.qc(test.raw)
nano.prenorm.qc <- function(data, code_class = "Housekeeping", prefix, out_path=NULL){
  raw_data <- data$raw

  if(is.null(out_path)){
    out_path = paste(getwd(),data$run_info$run_id[1], sep = "/")
  }
  
  library(reshape2)
  library(ggplot2)
  library(ggrepel)
  library(gridExtra)
  library(grid) # dependencies of gridExtra
  
  geomean <- function(x){
    exp(mean(log(as.matrix(x))))
  }

  cv <- function(x){
    (sd(x))/mean(x) * 100
  }

  hk <- raw_data[raw_data$Code.Class == code_class,]
  row.names(hk) <- hk$Name
  hk <- subset(hk, select = -c(Code.Class,Name,Accession))

  hk.mean <- colMeans(hk)
  hk.geomean <- apply(hk, 2, geomean)
  hk.cv <- apply(hk, 2, cv)
  hk.prenorm_qc <- data.frame(Sample=colnames(hk),Mean=hk.mean,GeoMean=hk.geomean,CV=hk.cv) # Samples for melt()

  # plot #
  # individual hk genes #
  hk$genes <- row.names(hk)
  hk.melt <- reshape2::melt(hk, variable.names = "genes")
  p.genes <- ggplot(aes(x = variable, y = value, shape = genes),data = hk.melt)
  p.genes <- p.genes + geom_point(aes(color = genes)) +
    ggtitle("Expression") +
    theme(axis.text.x=element_text(angle = 90, hjust = 0),
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.box = "vertical") +
    xlab('Sample') +
    ylab('Values')

  # geometric means of hk genes #
  hk.prenorm_qc.melt <- reshape2::melt(hk.prenorm_qc[c("GeoMean","CV","Mean","Sample")], variable.names = "Sample") # melt requires matrix for rownames.
  #p.qc <- ggplot(aes(x = Sample, y = value), data = hk.prenorm_qc.melt)  # plot both mean and geomean
  p.qc <- ggplot(aes(x = Sample, y = value), data = subset(hk.prenorm_qc.melt, hk.prenorm_qc.melt$variable == "GeoMean"))  # geomean only
  max_val <- max(hk.prenorm_qc.melt[hk.prenorm_qc.melt$variable == "GeoMean",]$value)
  print(paste0("[troubleshoot] max_val:", max_val))
  p.qc <- p.qc + geom_boxplot(aes(color= variable)) +
    ggtitle("Geo_Mean") +
    theme(axis.text.x=element_text(angle = 90, hjust = 0),
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.box = "vertical") +
    xlab('Sample') +
    ylab('Values') +
    scale_y_continuous(breaks = seq(from = 0, to = max_val,by = 500)) +
    coord_cartesian(ylim = c(0, max_val))

  p.cv <- ggplot(aes(x = Sample, y = value), data = subset(hk.prenorm_qc.melt, hk.prenorm_qc.melt$variable == "CV"))  # CV only
  p.cv <- p.cv + geom_point() +
    ggtitle("CV") +
    theme(axis.text.x=element_text(angle = 90, hjust = 0),
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.box = "vertical") +
    xlab('Sample') +
    ylab('Values') +
    scale_y_continuous(breaks = seq(from = 0, to = 200,by = 100)) +
    coord_cartesian(ylim = c(0, 200))

  p.mean_cv <- ggplot(aes(x = Mean, y = GeoMean, col=CV), data = hk.prenorm_qc)  # CV only
  p.mean_cv <- p.mean_cv + geom_point() +
    ggtitle("Geo_Mean vs Mean") +
    #geom_text(aes(label=Sample),hjust=0, vjust=0) +   ## use geom_label_repel
    geom_label_repel(aes(Mean, GeoMean, label = Sample),
                     size = 2,   # control text size via geom_label param
                     box.padding = 0.35, point.padding = 1,
                     segment.color = 'grey50') +
    theme(axis.text.x=element_text(angle = 90, hjust = 0),
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.box = "vertical")

  ### Output ###
  #data_name <- (substitute(raw_data)) # get data name. error 
  data_name <- "raw_data"
  print(paste0("[MSG] Exporting QC report"))
  
  if(!is.null(prefix)){prefix <- paste0(prefix,"_")}
  
  pdf(file = paste0(out_path,"/",prefix,"Prenorm_QC","_",data_name,".pdf"), width= 8,  height = 10.5)
  grid.arrange(p.genes, p.qc, p.cv, p.mean_cv, ncol=2,
               top=textGrob("Prenorm Housekeeping Genes", gp=gpar(fontsize=15,font=8))
  )
  dev.off()
  
  write.table(hk.prenorm_qc, file = paste0(out_path,"/",prefix,"Prenorm_QC","_",data_name,".txt"), sep = "\t", col.names = NA, quote = FALSE)

  data$prenorm_qc <- hk.prenorm_qc
  return(data)
}




##### Pretrain QC #####
#' Transform and check dataframe for missing values
#'
#' @param train.norm Normalized dataframe.
#' @param test.norm Normalized dataframe.
#' @return ??.
#' @example
#' library(caret)
#' test.data <- nano.prep(test.norm)

nano.prep <- function(data){
    # Transpose training and testing samples
    # samples are in rows and features are in columns
    norm <- data$norm
    norm.t <- as.data.frame(t(norm))

    library(caret)
    
    # Return the positions of the variables that are flagged to be problematic.
    norm.t.nzv <- nearZeroVar(norm.t, saveMetrics= TRUE)
    # check for blanks #
    # causes issues for single sample cases and modified to allow such cases and issue warning.

    if(nrow(norm.t.nzv[norm.t.nzv$nzv == "TRUE",]) != 0){
      warning(paste0("[MSG] Zero variance detected in data - proceed with caution. For single sample analysis this is expected and can be ignored."))
      print(norm.t.nzv[norm.t.nzv$nzv == "TRUE",])
      #stopifnot(nrow(norm.t.nzv[norm.t.nzv$nzv == "TRUE",]) == 0)
    }
    
    data$norm.t <- norm.t
    print(paste0("[MSG] Dataset read for testing/training"))
    return(data)

}

##### Training #####

generate.models <- function(data, training_memberships_path, N.train.per, prefix, alg_list = c("rf","glmnet","pam", "nb", "knn"), probes_rank, min_test_features=5, max_test_features=30){

    raw.obj <- train_model <- ""
    
    raw.obj <- nano.trainsplit(data = data,training_memberships_path = training_memberships_path,N.train.per = N.train.per)
    
    training_model_obj <- nano.train(prefix = prefix, data = data , alg_list = alg_list, probes_rank = probes_rank, min_test_features=min_test_features, max_test_features=max_test_features)
    nano.train.report(prefix, training_model_obj=training_model_obj, feature_min, feature_max)
}

##### Train Splitting #####
# v2 fixed bugs
# v3 fix bug where memberships from link were not sorted correctly.
#    retired "have_test_label"
#' Split training data and assign memberships to training data (or testing). Returns index for splitting data.
#' @param train.norm Normalized dataframe.
#' @param training_memberships_path Normalized dataframe.
#' @param N.train.per ## percentage of training samples use for training, the remaining would be used for validation e.g. 0.8
#' @param have_test_label If ground truth labels are avaliable for test data
#' @return train.idx
#' @example
#' train.idx<-  nano.trainsplit(train.data,training_memberships_path, 0.8, FALSE)
#' train.data.main <- train.data[train.idx,,drop=FALSE]
#' train.data.validate <- train.data[-train.idx,,drop=FALSE]

nano.trainsplit <- function(data,training_memberships_path,N.train.per){
    train.data <- data$norm.t
    training_memberships <- read.table(file = training_memberships_path, sep = "\t", row.names = 1, header = FALSE)
    row.names(training_memberships) <- make.names(row.names(training_memberships)) # fix names
    colnames(training_memberships) <- "Group"
    # Create Training_Full labels
    training_memberships <- training_memberships[row.names(training_memberships) %in% row.names(train.data),,drop=FALSE]
    #training_memberships <- training_memberships[match(row.names(training_memberships), row.names(train.data)),,drop=FALSE]
    training_memberships <- training_memberships[match(row.names(train.data), row.names(training_memberships)),,drop=FALSE] #190530
    training_memberships <- training_memberships[!is.na(training_memberships$Group),,drop=FALSE] #190530
    
    train.grp <- training_memberships[row.names(training_memberships) %in% row.names(train.data),"Group",drop=FALSE]
    train.data <- merge(train.data,train.grp,by="row.names") 
    row.names(train.data) <- train.data[,1]
    train.data <- train.data[,-1]

    # Choose subsample from Training_Full
    # createDataPartition
    # 1) takes in account of balance of labels
    # input = samples and labels
    # output = index of selected samples
    train.idx <- createDataPartition(y = train.data$Group, ## the outcome data are needed for random sampling
                                         p = N.train.per,     ## The percentage of data in the training set
                                         list = FALSE)
    train.data.main <- train.data[train.idx,,drop = FALSE]     # main training
    train.data.validate <- train.data[-train.idx,,drop = FALSE]     # training validation

    print(paste0("[MSG] Training: N",nrow(train.data.main),
                     "     Training_validatation: N",nrow(train.data.validate)))

    
    data$train.data.main <- train.data.main
    data$train.data.validate <- train.data.validate
    return(data)
}





########## Training ###########

##### Train Data #####
# v3 read list of genes rather than path 
# v2 tracks cv results and save as training_model_mat_list
#' Prenorm Nanostring Data and produce report
#'
#' @param train.norm Normalized dataframe.
#' @param alg_list list of alg
#' @param probes_rank_path path to probe ranking list. Given list sorted by p value. Looks for probes_list.txt in the working directory if probes_rank_path is missing.
#' @param min_test_features Minimum number of features tested.
#' @param max_test_features Maximum number of features tested.
#' @param out_path output path
#' @return Optimal_Training_Attributes Full_Training_Attributes.txt performance_pdf
#'
#' @example
#' nano.train(raw_data_path, omit_file_path_with_file_name)

nano.train <- function(prefix, data , alg_list = c("rf","glmnet","pam", "nb", "knn"), probes_rank_path=NULL, min_test_features=20, max_test_features=30, c.method = "repeatedcv", c.repeats = 5, c.number = 10, out_path=NULL) {
    
    if(is.null(data$train.data.main) || is.null(data$train.data.validate)){
        stop("[MSG] train.data.main / train.data.validate missing. Did you run nano.trainsplit()?")
    }
    train.data.training_main <- data$train.data.main
    
    if(is.null(probes_rank_path)) {
      if(file.exists("probes_list.txt")){
        probes.list.full <- as.character(unlist(read.table("probes_list.txt",sep = "/")))
      } else  {
        stop("[MSG] Probes rank path missing.")
      }
    } else  {
        probes.list.full <- as.character(unlist(read.table(file = probes_rank_path)))
    }

    if(is.null(out_path)){
      out_path = paste(getwd(),data$run_info$run_id[1], sep = "/")
    }
    
    train.settings <- list(alg_list=alg_list, probes_rank=probes.list.full, min_test_features=min_test_features,max_test_features=max_test_features,c.method =c.method , c.repeats = c.repeats, c.number = c.number )

    training_model_obj <- list()
    ###### Training Settings #####

    # set seeed #
    # set.seed(123)
    # seeds <- vector(mode = "list", length = 101)
    # for(i in 1:50) seeds[[i]] <- sample.int(1000, 100)
    # seeds[[1]] <- sample.int(1000, 1)
    ctrl <- trainControl(method = train.settings$c.method,
                         repeats = train.settings$c.repeats,  # rep20 feature 30-2 CV10, takes about 40 minutes for rf
                         number = train.settings$c.number,  # x fold CV
                         classProbs = TRUE, # save class probability for things like ROC curve
                         #seeds = seeds, # see above
                         savePredictions= TRUE)  # for CaretEnsemble

    ### reset dataframe ##

    training_model_mat_list <- training_model_list <- c()  # stores trained models and results
    training_model_mat_list.full <- train_sub.class.accuracy.full <- data.frame()

    for (alg in alg_list) {

      print(paste0("[MSG] training and optimizing algorithm: ",alg))
      filename = paste0(prefix,"_",alg)

      ## different number of features for loop ##
      for (probe.idx in max_test_features:min_test_features){  ## select features based on list given
        print(paste0("[MSG] ----- ",paste0(alg,"_",probe.idx)," -----"))
        probes.list <- probes.list.full[1:probe.idx]

        # prepare data
        train.data.training_main.selected_probes <- train.data.training_main[,colnames(train.data.training_main) %in% c(probes.list,"Group")]

        ##### Trainig #####
        train_model <- train(subset(train.data.training_main.selected_probes, select = -Group),
                             train.data.training_main.selected_probes$Group ,
                             method = alg,
                             tuneLength = 10,# by default the function will tune through three values of each tuning parameter. Use 10 
                             trControl = ctrl,
                             # metric="ROC",
                             preProc = c("center", "scale"))

        assign(paste0(alg,"_",probe.idx), train_model)

        training_model_list[[paste0(alg,"_",probe.idx)]] <- train_model
        ## save every model and overwrites it ##
        #saveRDS(training_model_list, file=paste0(out_path,"/",prefix,"_Training_Models_List.tmp.RDS")) # replaces original RDS with newer one every loop rather than in the end in case crashing

        ## show variables that were used in the final model
        print(predictors(train_model))

        ##### track testing samples #####

        train_model.names <- data.frame(rowIndex = c(1:nrow(train_model$trainingData)), Sample = row.names(train_model$trainingData))
        train.mat <- merge(train_model$pred, train_model.names, by="rowIndex") 

        bestTuneModel <- train_model$bestTune
        ncol_total <- ncol(train.mat)
        
        ## create a matrix of training results in every Resample rep using if optimal settings were used (may not match final bestTuneModel exactly e.g. glmnet) ##
        if (alg != "glmnet") {
            # number of parameters to tune #
            if (length(bestTuneModel) == 1 ){
                train.mat$grid_index_col <- train.mat[, ncol_total- 2] 
                train.mat.best <- train.mat[train.mat[,"grid_index_col"] == bestTuneModel[[1]],] 
            } else if (length(bestTuneModel) == 2){
                train.mat$grid_index_col <- train.mat[,ncol_total - 3]
                train.mat$grid_index_col2 <- train.mat[,ncol_total - 4]
                train.mat.best <- train.mat[train.mat[,"grid_index_col"] == bestTuneModel[[1]] & train.mat[,"grid_index_col2"] == bestTuneModel[[2]],] 
            }
            train.mat.best <- train.mat.best[with(train.mat.best, order(rowIndex,grid_index_col,Resample)),]
            train.mat.best$Matching <- ifelse(train.mat.best$pred == train.mat.best$obs, 1,0)
            
            train.mat.best.simple <- aggregate(train.mat.best$Matching, by=list(train.mat.best$Sample), FUN=sum)
            colnames(train.mat.best.simple) <- c("Sample",paste0(alg,"_",probe.idx,"_N_Match"))
            
            library(reshape2)
            train.mat.best.cast <- reshape2::dcast(train.mat.best, Sample+pred+obs ~ Resample, value.var="Sample")
            train.mat.best.cast <- merge(train.mat.best.simple, train.mat.best.cast, by="Sample")
    
            if (length(training_model_mat_list.full)==0){
              training_model_mat_list.full <- train.mat.best.simple
            }  else{
              training_model_mat_list.full <- merge(training_model_mat_list.full, train.mat.best.simple, by="Sample")
            }
        } # if using incompatible alg
        
        # confusion matrix #
        train_model.conmat <- confusionMatrix(train_model)
        #training_model_mat_list[[paste0(alg,"_",probe.idx,"_N_match")]] <- train.mat.best.simple
        training_model_mat_list[[paste0(alg,"_",probe.idx,"_matrix")]] <- train.mat.best.cast
        training_model_mat_list[[paste0(alg,"_",probe.idx,"_confmatrix")]] <- train_model.conmat
        training_model_mat_list[["Training_Result_Matrix"]] <- training_model_mat_list.full
        saveRDS(training_model_mat_list, file=paste0(out_path,"/",prefix,"_Training_Models_Mat_List.RDS")) # replaces original RDS with newer one every loop rather than in the end in case crashing
        
      } # probe.idx loop

    } # alg loop
    
    print(paste0("[MSG] ",length(training_model_list)," models created.")) #N119 created
    
    ##### Output #####
    training_model_obj[["training_model_list"]] <- training_model_list
    training_model_obj[["train.data.main"]] <- train.data.training_main
    training_model_obj[["train.settings"]] <- train.settings
    
    saveRDS(training_model_obj, file=paste0(out_path,"/",prefix,"_Training_Models_List.RDS")) 
    return(training_model_obj)
    # pdf(file = paste0(prefix,"_Accuracy_by_Alg.pdf"), width = 8, height = 10.5)
    # gg <- ggline(internal_performance, x = "Num_Features", y = "Accuracy",
    #              linetype = "Alg",
    #              color = "Alg",
    #              shape = "Alg",
    #              size = 0.5,
    #              main = "Accuracy by Alg",
    #              xlab = "Number of probes",
    #              ylab = "Accuracy"
    # ) +
    #   scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.1)) +
    #   coord_cartesian(ylim = c(0, 1))
    # print(gg)
    # dev.off()

}

######### [3b] get training results ##########
#v3 - accepts training_model_obj
#v2 - added Alg performance summary

#' @param prefix - file prefix that goes to algorithms performance report
#' @param training_model_obj - training model object where training performance is pulled out from
#' @param feature_min - controls the number of min plotting range
#' @param feature_max - controls the number of max plotting range
#' @param out_path output path
nano.train.report <- function(prefix, training_model_obj, feature_min, feature_max, out_path=NULL){
  
  if(is.null(out_path)){
    out_path = paste(getwd())
  }
  
  train_list <- training_model_obj[["training_model_list"]]
  full_internal_performance <- internal_performance <- data.frame()
  
  for (model_idx in 1:length(train_list)){
    train_model <- train_list[[model_idx]]
    alg <- as.character(train_model$method)
    probe.idx <- as.numeric(gsub(paste0(alg,"_"),"",names(train_list[model_idx])))
    
    ## print model name ##
    #print(paste0("[MSG] Parsing ",names(train_list[model_idx]))) # comment out for brevity
    
    ########## Evaluate Results ##########
    #print(train_model$results) # Training Results
    if (alg == "svmRadial") {
      Final_Model_Training_Row <- as.numeric(row.names(train_model$bestTune))
    } else if(!(alg %in% c("svmLinear","svmPoly"))){ # ie other than these two
      Final_Model_Training_Row <- as.numeric(rownames(train_model$finalModel$tuneValue)) #Optimal training parameter position
    } else {
      Final_Model_Training_Row <- 1
    } # for alg with no optimization
    Final_Model_Training_Stats <- train_model$results[Final_Model_Training_Row,]
    internal_performance_tmp <-data.frame(Alg = train_model$method, Num_Features = probe.idx, Final_Model_Pos = Final_Model_Training_Row, Accuracy = Final_Model_Training_Stats$Accuracy, AccuracySD = Final_Model_Training_Stats$AccuracySD, Kappa = Final_Model_Training_Stats$Kappa, KappaSD = Final_Model_Training_Stats$KappaSD)
    full_internal_performance_tmp <-data.frame(Alg = train_model$method, Num_Features = probe.idx, Final_Model_Pos = Final_Model_Training_Row, Accuracy = train_model$results$Accuracy, AccuracySD = train_model$results$AccuracySD, Kappa = train_model$results$Kappa, KappaSD = train_model$results$KappaSD)
    internal_performance <- rbind(internal_performance, internal_performance_tmp)
    full_internal_performance <- rbind(full_internal_performance, full_internal_performance_tmp)
    
    ## show variables that were used in the final model
    # print(predictors(train_model)) # comment out for brevity
    
  } # for loop
  
  
  ##### Prepare output #####
  probe_max <- max(internal_performance$Num_Features) # max in trained data
  probe_min <- min(internal_performance$Num_Features) # min in trained data
  stopifnot(feature_max >= probe_min & feature_max <= probe_max)
  stopifnot(feature_min >= probe_min & feature_min <= probe_max)
  
  # calculate avg #
  
  internal_performance$Alg <- as.factor(internal_performance$Alg)
  stats <- stats.i <- data.frame()
  for(i in probe_min:probe_max){
    for(j in probe_min:probe_max){ 
      #print(paste0(i," vs ",j))
      # it is ok for i==j
      if(i<=j){ 
        internal_performance.i <- internal_performance[internal_performance$Num_Features >= i & internal_performance$Num_Features <= j,]
      }else if(i>j){
        internal_performance.i <- internal_performance[internal_performance$Num_Features <= i & internal_performance$Num_Features >= j,]
      }
      internal_performance.i[internal_performance.i$Num_Features > feature_min & internal_performance.i$Num_Features < feature_max ,]
      stats.i <- aggregate.data.frame(x=internal_performance.i$Accuracy, by=list(internal_performance.i$Alg),FUN = mean)
      # stats.i$Num_Features <- paste0(i,"_",j)
      stats.i$Num_Features.i <- i
      stats.i$Num_Features.j <- j
      colnames(stats.i) <- c("Alg","Avg_accuracy","Num_Features.i","Num_Features.j")
      stats <- rbind(stats,stats.i)
    } # j
  } # i
  
  pdf(file = paste0(prefix,"_Accuracy_by_Alg.pdf"), width = 10.5, height = 8)
  ## plot ##
  g_conf_mat.facet <- ggplot(data =stats, aes(Num_Features.i, Num_Features.j,Avg_accuracy))+ 
    geom_tile(aes(fill = Avg_accuracy),colour = "white")   + 
    scale_x_continuous(breaks = seq(from=probe_min, to =probe_max, by=2))+
    scale_y_continuous(breaks = seq(from=probe_min, to =probe_max, by=2))+
    scale_fill_gradientn(colours = c("cyan", "black", "red"))+
    theme_minimal() +  coord_equal(ratio = 1) +
    #facet_grid(Alg~. )
    facet_wrap(~Alg, ncol=2)
  print(g_conf_mat.facet)
  
  g_conf_mat <- ggplot(data =stats, aes(x=Num_Features.i, y=Num_Features.j,z=Avg_accuracy))+ 
    geom_tile(aes(fill = Avg_accuracy),colour = "white")   + 
    scale_x_continuous(breaks = seq(from=probe_min, to =probe_max, by=2))+
    scale_y_continuous(breaks = seq(from=probe_min, to =probe_max, by=2))+
    #scale_fill_gradient2(low = "blue",  high = "red") +
    #scale_fill_gradient2() +
    scale_fill_gradientn(colours = c("cyan", "black", "red"))+
    theme_minimal() + coord_equal(ratio = 1) +  ggtitle("Overall Average") +
    theme(axis.text.x=element_text(angle = 90, hjust = 0))
  print(g_conf_mat)
  
  # g_conf_mat + geom_density_2d(stats, aes(x = Num_Features.i, y = Num_Features.j, z = Avg_accuracy))
  # g_conf_mat + geom_density_2d()
  # g_conf_mat + stat_contour(stats, aes(x = Num_Features.i, y = Num_Features.j, z = Avg_accuracy))
  
  ## line graph ##
  internal_performance.select <- internal_performance[internal_performance$Num_Features >= feature_min & internal_performance$Num_Features <= feature_max,,drop=FALSE]
  internal_performance.select.num_features.agg <- aggregate(internal_performance.select[,"Accuracy",drop=FALSE], by=list(internal_performance.select$Num_Features), FUN=mean)
  colnames(internal_performance.select.num_features.agg) <- c("Num_Features","Avg_accuracy")
  
  internal_performance.select.alg.agg <- aggregate(internal_performance.select[,"Accuracy",drop=FALSE], by=list(internal_performance.select$Alg), FUN=mean)
  colnames(internal_performance.select.alg.agg) <- c("Alg","Avg_accuracy")
  
  gg_line.facet <- ggplot(data=internal_performance, aes(x=Num_Features, y=Accuracy, colour=Alg)) + 
    geom_point(show.legend = FALSE) + 
    facet_grid(Alg~.)+
    scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.02)) +
    scale_x_continuous(breaks = seq(from = feature_min, to = feature_max ,by = 2)) +
    coord_cartesian(xlim = c(feature_min, feature_max)) +
    theme_minimal()+ geom_path(aes(colour = Alg),show.legend = FALSE)  +  
    #stat_summary(fun.y = mean, geom="line")
    #stat_summary(aes(group=Alg), fun.y=mean, geom="line", colour="red")
    #geom_hline(data = internal_performance, aes(yintercept = mean(Accuracy), colour = Alg), color="blue")
    geom_hline(data = internal_performance.select.alg.agg, aes(yintercept=Avg_accuracy, group=Alg), linetype = "dashed", show.legend = FALSE)+
    geom_text(data = internal_performance.select.alg.agg, aes(x=0,y=Avg_accuracy,colour="black", group=Alg, label=paste("avg_acc",round(Avg_accuracy,3))), nudge_x=mean(c(probe_max,probe_min)), nudge_y=0.05, cex=3,show.legend = FALSE)
  print(gg_line.facet)
  
  gg_line.combined <- ggplot(data=internal_performance, aes(x=Num_Features, y=Accuracy, colour=Alg)) + 
    geom_point(show.legend = FALSE) + 
    #facet_grid(Alg~.)+
    scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.02)) +
    scale_x_continuous(breaks = seq(from = feature_min, to = feature_max ,by = 2)) +
    coord_cartesian(xlim = c(feature_min, feature_max)) +
    theme_minimal()+ geom_path(aes(colour = Alg),show.legend = FALSE)  +  
    geom_hline(data = internal_performance.select.alg.agg, aes(yintercept=Avg_accuracy, group=Alg, colour=Alg), linetype = "dashed", show.legend = FALSE)+
    geom_text(data = internal_performance.select.alg.agg, aes(x=0,y=Avg_accuracy,colour="black", group=Alg, label=paste(Alg,"avg_acc",round(Avg_accuracy,3))), nudge_x=probe_min*2, cex=3,show.legend = FALSE)
  print(gg_line.combined)
  
  gg_line <- ggplot(data=internal_performance.select.num_features.agg, aes(x=Num_Features, y=Avg_accuracy)) + 
    geom_point(show.legend = FALSE) + 
    scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.02)) +
    scale_x_continuous(breaks = seq(from = feature_min, to = feature_max ,by = 2)) +
    coord_cartesian(xlim = c(feature_min, feature_max)) +
    theme_minimal()+ geom_path(show.legend = FALSE)  +  
    #stat_summary(fun.y = mean, geom="line")
    #stat_summary(aes(group=Alg), fun.y=mean, geom="line", colour="red")
    #geom_hline(data = internal_performance, aes(yintercept = mean(Accuracy), colour = Alg), color="blue")
    geom_hline(data = internal_performance.select.num_features.agg, aes(yintercept=mean(Avg_accuracy)), linetype = "dashed", show.legend = FALSE)+
    geom_text(data = internal_performance.select.num_features.agg, aes(x=0,y=mean(Avg_accuracy),colour="black", 
                                                                       label=paste("avg_acc",round(mean(Avg_accuracy),3))), 
              nudge_x=mean(c(probe_max,probe_min)), nudge_y=0.05, cex=3,show.legend = FALSE)+  
    ggtitle("Overall Average") +
    theme(axis.text.x=element_text(angle = 90, hjust = 0))
  print(gg_line)
  dev.off()
  
  overivew_internal_performance <- get.training.stats(train_list)
  overivew_internal_performance <- merge(overivew_internal_performance,internal_performance.select.alg.agg, by="Alg")
  
  print(overivew_internal_performance)
  
  ##### Output #####
  print(paste0("[MSG] Check directory for detailed reports."))
  write.table(internal_performance, file = paste0(out_path,"/",prefix,"_Optimal_Training_Attributes.txt"), col.names = NA, sep = "\t")
  write.table(full_internal_performance, file = paste0(out_path,"/",prefix,"_Full_Training_Attributes.txt"), col.names = NA, sep = "\t")
  write.table(overivew_internal_performance,file = paste0(out_path,"/",prefix,"_Overview_Training_Attributes.txt"), col.names = NA, sep = "\t")
}


######### [3c] parse training results ##########
# supporting function used in nano.train.report ()
get.training.stats <- function(training_model_obj=models){
  
  models = models$training_model_list
  
  model.stats.full <- data.frame()
  for(I in 1:length(models)){
    ALG <- models[[I]]$method
    FET <- ncol(subset(models[[I]]$trainingData, select = -.outcome))
    MET <- models[[I]]$control$method
    CV <- models[[I]]$control$number
    REP <- models[[I]]$control$repeats
    model.stats <- data.frame(Alg=ALG, N_Features=FET, Method=MET, CV=CV, Repeats=REP)
    model.stats.full <- rbind(model.stats.full, model.stats)
  }
  Training_N <- aggregate(model.stats.full$Alg, by=list(model.stats.full$Alg), FUN=length)
  Training_Range <- aggregate(model.stats.full$N_Features, by=list(model.stats.full$Alg), FUN=range)
  Training_Met_CV_Rep <- aggregate(model.stats.full[,c(3:5)], by=list(model.stats.full$Alg), FUN=unique)
  Training_summary <- merge(Training_Range,Training_Met_CV_Rep, by="Group.1")
  Training_summary <- data.frame(merge(Training_summary,Training_N, by="Group.1"))
  colnames(Training_summary) <-c("Alg","Feature_Range","Method","CV","Repeats","N_Models")
  return(Training_summary)
  
}

######### [4] norm #########
# perform Nanostring norm using NanoStringNorm. For details refer to help from package
#' @param data test data
#' @param SampleContent what to normalize with
#' @param round.values round values or not
#' @param take.log take log or not
#' @param return.matrix.of.endogenous.probes return matrix of endogenous probes or not
#' @param verbose verbose

nano.norm <- function(data, SampleContent = "housekeeping.geo.mean", round.values = FALSE, take.log = TRUE, return.matrix.of.endogenous.probes = TRUE, verbose = TRUE){

  library(NanoStringNorm)
  
  norm <- NanoStringNorm(x = data$raw, SampleContent = SampleContent,  round.values = round.values,  take.log = take.log,   return.matrix.of.endogenous.probes = return.matrix.of.endogenous.probes, verbose = return.matrix.of.endogenous.probes)
  norm <- as.data.frame(norm)
  if(ncol(norm) == 1) {
      sample_id <- colnames(data$raw)[4] # sample name gets dropped from NanoStringNorm when there's only one sample
      colnames(norm) <- sample_id
  }
    
  data$norm <- as.data.frame(norm)
  return(data)
  
}

######### [5] Test #########
# v2 updated to use train_mod_obj
# perform confusion matrix if "Group" column exist in dataframe
# settings #
#' Prenorm Nanostring Data and produce report
#'
#' @param training_model_list Trained_models.
#' @param alg_list list of alg
#' @param data test data.
#' @param min_test_features Minimum number of features tested.
#' @param max_test_features Maximum number of features tested.
#' @param prefix output prefix e.g. paste(project,format(Sys.time(), "%Y-%m-%d_%H%M"), sep = "_")
#' @param out_path output location
#' @return testing_results_summary_groups_score
#'
#' @example
#' nano.test(raw_data_path, omit_file_path_with_file_name)


nano.test <- function(prefix, training_model_obj, data , alg_list=NULL, min_test_features=NULL, max_test_features=NULL, out_path=NULL) {
    
    if(is.null(data$norm.t)) {
      stop("[MSG] Run nano.prep() first before testing.")
    } else {
      test.df <- data$norm.t
    }
  
    if(is.null(training_model_obj[["training_model_list"]])) {
      stop("[MSG] Run training_model missing. Run nano.train() first.")
    } else {
      training_model_list <- training_model_obj[["training_model_list"]]
    }  
  
    if(is.null(out_path)){
      out_path = paste(getwd(),data$run_info$run_id[1], sep = "/")
    }
  
    ## setting defaults ##
    training_model_list.mat <- data.frame(matrix(unlist(strsplit(names(training_model_list),"_")), ncol=2,byrow = TRUE, dimnames = list(NULL,c("Alg","Model"))), stringsAsFactors = FALSE)

    if(is.null(min_test_features)&is.null(max_test_features)) {
      min_test_features <- min(as.numeric(training_model_list.mat$Model))
      max_test_features <- max(as.numeric(training_model_list.mat$Model))
      print(paste0("[MSG] Using min and max number of features from training model - min: ",min_test_features," and max: ",max_test_features))
    }
    
    if(is.null(alg_list)){
      alg_list <- unique(training_model_list.mat$Alg)
      print(paste0("[MSG] Using algorithm list from training model:"))
      print(paste0(alg_list))
    } else {
      alg_list <- unique(training_model_list.mat$Alg)[unique(training_model_list.mat$Alg) %in% alg_list]
      print(paste0("[MSG] Using algorithm:"))
      print(paste0(alg_list))   
    }
    library(caret)
    
    print(paste0("[MSG] Testing ",nrow(test.df), " samples:"))
    print(row.names(test.df))

    # Initialize Report #
    conf_matrix_list <- c()
    testing_results <- c()
    conf_matrix_results.full <- c() # from conf_matrix
    testing_results_summary_groups_score <- c() # testing class and probability

    for (alg in alg_list) {

      for (i in max_test_features:min_test_features) {

        # Select models #
        testing_model <- training_model_list[[eval(paste0(alg,"_",i))]]

        # Predict #
        #print(paste0("Testing... ",alg,"_",i))

        if (("Group" %in% colnames(test.df))){
          # If Group information is present (ie testing known samples) #
          ## print(paste0("[MGS] Have labels. Contruct Confusion Matrix."))
          if (alg == "svmLinear"){
            preProcValues <- preProcess(subset(train.data.training_main.selected_probes, select = -Group), method= c("center", "scale") )

          }
          Prediction_class <- predict(testing_model, newdata = subset(test.df, select = -Group)) # predicts class or use subset(test.df, select = -Group)
          Prediction_prob <- predict(testing_model, newdata = subset(test.df, select = -Group), type = "prob") # predicts class probability
          class_max_prob <- as.data.frame(apply(Prediction_prob,1,max))
          colnames(class_max_prob) <- "Probability"

          testing_results <- data.frame("Sample" = rownames(test.df), "Alg" = alg, "Num_Features" = i, "Class" = as.character(Prediction_class),"Probability" = as.numeric(class_max_prob$Probability))
          colnames(testing_results) <- c("Sample" , "Alg","Num_Features","Class","Probability")
          testing_results_summary_groups_score <- as.data.frame(rbind(testing_results_summary_groups_score, testing_results))
          colnames(testing_results_summary_groups_score) <- c("Sample","Alg","Num_Features","Class","Probability")

          # Confusion matrix (when actual class is known) #
          # Accuracy(PredictedTest, testTrainData$label)

          conf_matrix<-confusionMatrix(data = Prediction_class, test.df$Group, positive = NULL)
          conf_matrix_list[[eval(paste0(alg,"_",i))]] <- conf_matrix

          conf_matrix_results.current <- as.data.frame(t(c(alg,i,t(as.data.frame(conf_matrix$overall)))))
          colnames(conf_matrix_results.current) <- c("Alg","i","Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull","AccuracyPValue","McnemarPValue")

          conf_matrix_results.full <- rbind(conf_matrix_results.full,conf_matrix_results.current)
          colnames(conf_matrix_results.full) <- c("Alg","i","Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull","AccuracyPValue","McnemarPValue")

        } else {
          # If Group information is absent (most classification cases) #
          ##print("[MSG] Group not detected in dataframe. Skip Confusion Matrix...")
          if (alg == "svmLinear"){
            preProcValues <- preProcess(train.data.training_main.selected_probes, method= c("center", "scale") )

          }
          Prediction_class <- predict(testing_model, newdata = test.df) # predicts class or use subset(test.df, select = -Group)
          Prediction_prob <- predict(testing_model, newdata = test.df, type = "prob") # predicts class probability
          class_max_prob <- as.data.frame(apply(Prediction_prob,1,max))
          colnames(class_max_prob) <- "Probability"

          testing_results <- data.frame("Sample" = rownames(test.df), "Alg" = alg, "Num_Features" = i, "Class" = Prediction_class,"Probability" = as.numeric(class_max_prob$Probability))
          colnames(testing_results) <- c("Sample" , "Alg","Num_Features","Class","Probability")
          testing_results_summary_groups_score <- as.data.frame(rbind(testing_results_summary_groups_score, testing_results))
          colnames(testing_results_summary_groups_score) <- c("Sample","Alg","Num_Features","Class","Probability")

        }
      }
    }
    
    ### Output ###
    print("[MSG] Printing results. Also check conf_matrix for matrix")
    
    if(!is.null(prefix)){prefix <- paste0(prefix,"_")}
    
    write.table(testing_results_summary_groups_score, file=paste0(out_path,"/",prefix,"testing_summary_full.txt"), sep = "\t", col.names = NA)
    
    data$test_results_full <- testing_results_summary_groups_score
    
    data$run_info$test_settings <- c("models" = names(training_model_list), "alg_list" = alg_list,"min_test_features" = min_test_features, "max_test_features" = max_test_features)
    return(data)
    
    # if you have conf_matrix # (need if statement)
    if (("Group" %in% colnames(test.df))){
      write.table(conf_matrix_results.full, file=paste0(out_path,"/",prefix,"_conf_matrix_results_full_Test.txt"), sep = "\t", col.names = NA)
      saveRDS(conf_matrix_list, file=paste0(out_path,"/",prefix,"_Alg_Conf_Matrix_Test.RDS")) # more details
    }

} # end of test



##### Get test results #####
# v2 - added remarks and soft/hard voting distinction
# Plot Test Probability by samples
# perform confusion matrix if "Group" column exist in dataframe
#' @param prefix file prefix (optional)
#' @param testing_results_summary_groups_score
#' @param print_report print or not. default is FALSE
#' @param out_path output path
#' @return Test_aggregate_summary_ _Test_aggregate_summary_MAX_
#' @example
#' library(ggpubr)
#' library(cowplot)
#' nano.test(raw_data_path, omit_file_path_with_file_name)
get.nano.test.results <- function(prefix, data, print_report=FALSE, out_path=NULL) {
    
  if(is.null(out_path)){
    out_path = paste(getwd(),data$run_info$run_id[1], sep = "/")
  }
  
  testing_results_summary_groups_score <- data$test_results_full
  prenorm_qc <- data$prenorm_qc # geomean, mean and cv
  min_test_features <- data$run_info$test_settings[["min_test_features"]]
  max_test_features <- data$run_info$test_settings[["max_test_features"]]
  
  #### Aggregate ####
  
  testing_results_summary_groups_score_i_max_full <- testing_results_summary_groups_score_i_summary_full <- data.frame(SAMPLE=character(0), Num_Features=character(0), Class=character(0), Count=numeric(0), Probability=numeric(0))
  
  ##### Get ensemble results #####
  for (SAMPLE in unique(testing_results_summary_groups_score$Sample)) {
    testing_results_summary_groups_score_i <- testing_results_summary_groups_score[testing_results_summary_groups_score$Sample == SAMPLE,]
    
    ## model summary (by number of features) ##
    for (FEATURES in c(min_test_features:max_test_features,"ALL")){
      # select results for specific # of features #
      if (FEATURES == "ALL"){
        # "ALL" uses all Num_Features #
        testing_results_summary_groups_score_i.feat <- testing_results_summary_groups_score_i
      } else {
        # select individual Num_Features
        testing_results_summary_groups_score_i.feat <- testing_results_summary_groups_score_i[testing_results_summary_groups_score_i$Num_Features == FEATURES,]
      }
    
        ## model summary (combined) ##
        # Count           Number of alg for that Num_Features that chooses the particular class as most probable
        # Num_Features    Number of features used for prediction
        # Class           Prediction Class 
        # N_models        Number of Alg's tested for that Num_Features
        # Agreement       Proportion of models that chooses that particular class as the most probable
        n_models_i <- nrow(testing_results_summary_groups_score_i.feat)
        testing_results_summary_groups_score_i_summary <- data.frame(summary(testing_results_summary_groups_score_i.feat$Class))
        testing_results_summary_groups_score_i_summary$Num_Features <- FEATURES
        testing_results_summary_groups_score_i_summary$Class <- as.character(row.names(testing_results_summary_groups_score_i_summary))
        colnames(testing_results_summary_groups_score_i_summary) <- c("Count", "Num_Features","Class")
        # avg probability from algs with same Num_Features
        avg_prob <- aggregate(Probability ~ Class ,testing_results_summary_groups_score_i.feat, mean)
        
        testing_results_summary_groups_score_i_summary <- merge(testing_results_summary_groups_score_i_summary, avg_prob, by="Class", all = TRUE, sort = TRUE)
        testing_results_summary_groups_score_i_summary <- cbind(SAMPLE,testing_results_summary_groups_score_i_summary, n_models_i)
        testing_results_summary_groups_score_i_summary$Model_Agreement <- testing_results_summary_groups_score_i_summary$Count/testing_results_summary_groups_score_i_summary$n_models_i
        colnames(testing_results_summary_groups_score_i_summary) <- c("Sample","Class","Count","Num_Features" ,"Avg_Probability","N_models","Agreement")
  
        # combine individual results
        testing_results_summary_groups_score_i_summary_full <- rbind(testing_results_summary_groups_score_i_summary_full,testing_results_summary_groups_score_i_summary)
        
        ##### Calculate most probable class ####
        if (FEATURES == "ALL"){
          testing_results_summary_groups_score_i_summary.ALL <- testing_results_summary_groups_score_i_summary[testing_results_summary_groups_score_i_summary$Num_Features == "ALL",]
          # Hard voting (Count priority) - Sort by count then Avg_Probability
          testing_results_summary_groups_score_i_max <- testing_results_summary_groups_score_i_summary.ALL[with(testing_results_summary_groups_score_i_summary.ALL, order(Count,Avg_Probability, decreasing = TRUE)),][1,]
          # Soft voting (Probability priority) - Sort by Avg_Probability then count
          #testing_results_summary_groups_score_i_max <- testing_results_summary_groups_score_i_summary.ALL[with(testing_results_summary_groups_score_i_summary.ALL, order(Avg_Probability, Count, decreasing = TRUE)),][1,]
          
          testing_results_summary_groups_score_i_max_full <- rbind(testing_results_summary_groups_score_i_max_full,testing_results_summary_groups_score_i_max)
        }
    }
  }
  
  #### Output ####
  
  if(!is.null(prefix)){prefix <- paste0(prefix,"_")}
  
  if(print_report == TRUE){
    write.table(testing_results_summary_groups_score_i_summary_full, file=paste0(out_path,"/",prefix,"test_summary_aggregate.txt"), sep = "\t", col.names = NA)
    write.table(testing_results_summary_groups_score_i_max_full, file=paste0(out_path,"/",prefix,"test_summary.txt"), sep = "\t", col.names = NA)
  }
      
  data$test_results_agg <- testing_results_summary_groups_score_i_summary_full
  data$test_results <- testing_results_summary_groups_score_i_max_full
  
  return(data)
} # end of get test results


#### Plot Data ####
# v2 allow avg and avg_cal probabilities
# updated visuals
nano.plot <- function(prefix, data, prob="Avg_Probability", thres_avg_prob=0, thres_geomean = 0, report_type=c("Summary","Detailed"), print_report=FALSE, out_path=NULL){
  
  # Check #
  if(prob == "Avg_Cal_Probability" & is.null(data$test_results_agg$Avg_Cal_Probability)){
    stop("[MSG] Avg_Cal_Probability missing from data. Run nano.calibrate() first to plot results from calibrated probability or proceed with Avg_Probability.")
  } else if(prob == "Avg_Cal_Probability" & is.null(data$test_results_agg$Avg_Cal_Probability)){
    stop("[MSG] Avg_Probability missing from data. Run get.nano.test.results() first.")
  }
  
  if(!is.numeric(thres_avg_prob)){stop("[MSG] thres_avg_prob must be numeric")}
  
  if(is.null(data$colour_code)){
    stop("[MSG] colour_code missing from data. Run nano.set.colour() first.")
  }
  
  if(is.null(out_path)){
    out_path = paste(getwd(),data$run_info$run_id[1], sep = "/")
  }
  
  
  test_results <- data$test_results
  test_results_agg <- data$test_results_agg
  test_results_full <- data$test_results_full
  prenorm_qc <- data$prenorm_qc
  col_code <- data$colour_code
  
  library(gridExtra)
  library(ggpubr)
  
  test_results$prob <- test_results[,prob]
  test_results_agg$prob <- test_results_agg[,prob]

  if(prob == "Avg_Probability"){plot_title <-"Average probability"}
  if(prob == "Avg_Cal_Probability"){plot_title <-"Average calibrated probability"}

  ### Output ###
  
  t.result.list <- list()
  if(!is.null(prefix)){prefix <- paste0(prefix,"_")}
  
  pdf(file = paste0(out_path,"/",prefix,"test_result_",report_type,"_plots.pdf"), width = 10.5, height = 8)
  for (SAMPLE in unique(test_results_full$Sample)) {
    print(SAMPLE)
    
    #### Result stats ####
    result.i <- c(test_results[test_results$Sample == SAMPLE,], prenorm_qc[prenorm_qc$Sample == SAMPLE,])
    
    if (result.i$GeoMean >= thres_geomean){
      remarks1 <- "QC - PASS"
    } else if (result.i$GeoMean < thres_geomean) {
      remarks1 <- "QC - FAIL"
    } 
    
    if (result.i[prob] >= thres_avg_prob && result.i$GeoMean >= thres_geomean){
      remarks2 <- "High Confidence"
    } else if (result.i[prob] > thres_avg_prob && result.i$GeoMean < thres_geomean) {
      remarks2 <- "Caution"
    } else if (result.i[prob] < thres_avg_prob) {
      remarks2 <- "Low Confidence"
    } 
    
    t.result <- t(data.frame("Sample" = result.i$Sample, "Class" = result.i$Class, "Class Prob" = paste0(round(as.numeric(result.i[prob]), 4) * 100, "%"), "Prob Type"=plot_title, "Models Tested" = result.i$N_models, "Models Agreement" =  paste0(round(result.i$Agreement * 100, 4),"%"), "GeoMean" = result.i$GeoMean, "CV" = result.i$CV, "Remarks1" = remarks1, "Remarks2" = remarks2))
    colnames(t.result) <- "Summary Table"
    t.result.list[[SAMPLE]] <- t.result
    ##### Plots ####
    color_chart <- col_code$Group_Colour
    names(color_chart) <- col_code$Group
    ## Format tables ##
    agg.table.i <- test_results_agg[test_results_agg$Sample == SAMPLE,]
    test_results_full.i <- test_results_full[test_results_full$Sample == SAMPLE,]
    
    ## probability by algorithms ##
    # by Num_Features != "ALL"
    agg.i <- agg.table.i[agg.table.i$Num_Features != "ALL",]
    agg.m <- reshape2::melt(agg.i[c("Num_Features",prob,"Class")], variable_name = "Class", id.vars=c("Num_Features",prob,"Class"))
    agg.m[is.na(agg.m)] <- 0 # for cases with no predictions probability
    agg.p <- ggline(agg.m, x = "Num_Features", y = prob,
                    linetype = "Class",
                    color = "Class",
                    palette = color_chart,

                    #facet.by = "Class",
                    nrow = 1,
                    size = 0.5,
                    main = paste0(SAMPLE," - Calibrated probability by features"),
                    xlab = "Number of Genes",
                    ylab = prob
    ) +
      scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.2)) +
      #scale_x_continuous(breaks = seq(from = 0, to = 30 ,by = 1)) +
      coord_cartesian(ylim = c(0, 1))
    
    ## cleveland plot ##
    agg.ALL <- agg.table.i[agg.table.i$Num_Features == "ALL",]
    Groups <- data.frame(Class=col_code$Group)
    agg.ALL <- merge(Groups, agg.ALL,by = "Class", all.x = TRUE) # make sure all three groups are present
    
    # fill NA's with 0 for cases with no predictions. Col 1 and 2 are Class and Sample (factors) #
    agg.ALL.tmp <- agg.ALL[,-c(1,2)]
    agg.ALL.tmp[is.na(agg.ALL.tmp)] <- 0
    agg.ALL <- cbind(agg.ALL[,c(1,2)],agg.ALL.tmp) 

    agg.ALL$Agreement <- as.numeric(round(agg.ALL$Agreement, digits = 4) * 100)
    agg.ALL.p <- ggdotchart(agg.ALL, x= "Class",y = prob,
                    linetype = "Class",
                    color = "Class",
                    palette = color_chart,
                    dot.size = 5,
                    nrow = 1,
                    size = 0.5,
                    legend="",
                    xlab = ""
                  ) +
                  scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.2)) +
                  coord_cartesian(ylim = c(0, 1))+
                  #theme(plot.title = element_text(hjust = 0.0)) +
                  theme_minimal()+
                  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                  #theme(axis.text = element_text(size=10))+
                  geom_hline(yintercept=seq(0,1,0.2), linetype="dashed", colour = "grey70")
    
    # reverse order to order from bottom up#
    agg.ALL$Class <- factor(agg.ALL$Class, levels = rev(levels(agg.ALL$Class)))
    agg.ALL.bar <- ggplot(data = agg.ALL, aes(x= Sample,y=Agreement, fill= Class, label = paste0(Agreement,"%"))) + 
                          geom_bar(data = agg.ALL, aes(x= Sample,y=Agreement, fill= Class),
                                stat="identity", width = 1) + 
                          theme_minimal() +
                          scale_fill_manual(values = color_chart) +
                          geom_text(position = position_stack(vjust = 0.5), size = 4) +  # add percentage
                          theme(legend.position="", # remove legend and axis text
                                #axis.text.y=element_blank(),
                                axis.text.x=element_blank(),
                                axis.title.x=element_blank()
                                )  +
                          scale_y_continuous(expand = c(0,0))  # remove space between axis and values

    
    # add margin (note need ggplot2::margin)
    agg.ALL.bar <- agg.ALL.bar + theme(
      panel.background = element_rect(fill = "white"),
      #plot.margin = ggplot2::margin(60, 10, 70, 20, "pt"), #top, right, bottom, left
      plot.background = element_rect(
        fill = "white",
        colour = "white",
        size = 1
      )
    ) 
    #+ ggtitle(paste0("Fig.2 - Models agreement"))+
    #  theme(plot.title = element_text(hjust = 0.0))
    
    # horiz bar plot #
    all.horbar <- test_results_full.i
    all.horbar$Order <- row.names(all.horbar)
    #all.horbar$Class <- factor(all.horbar$Class, levels = c("Group2B","Group2A","Group1"))
    all.horbar.p <- ggbarplot(all.horbar, x="Order",y="Probability",                    
                              color = "Class",
                              fill = "Class",
                              palette = color_chart,
                              sort.by.groups = TRUE,
                              rotate = TRUE,
                              sort.val = "asc",
                              legend="") + 
                    theme_minimal() +
                    ylab(NULL) +
                    xlab(NULL) +
                    theme(axis.text.y = element_blank(), 
                          axis.ticks.y = element_blank(),
                          axis.text.x = element_blank()
                    ) +
                    theme(legend.position="")

    ## raw line plot ##
    raw.i <- test_results_full[test_results_full$Sample == SAMPLE,]
    raw.m <- reshape2::melt(raw.i[c("Num_Features","Probability","Class","Alg")], variable_name = "Class", id.vars=c("Num_Features","Probability","Class","Alg"))
    raw.p <- ggline(raw.m, x = "Num_Features", y = "Probability",
                    linetype = "Class",
                    color = "Class",
                    palette = color_chart,
                    facet.by = "Alg",
                    nrow = 5,
                    size = 0.5,
                    main = paste0("Fig.3 - ",plot_title, " by algorithms"),
                    xlab = "N_features",
                    ylab = prob
    ) +
      scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.2)) +
      scale_x_continuous(breaks = seq(from = 0, to = 30 ,by = 1)) +
      coord_cartesian(ylim = c(0, 1))+
      theme(plot.title = element_text(hjust = 0.0))
    
   ##### Format Report ####
    
    if (report_type == "Detailed"){
      grid.arrange(
        top = paste0("Test Result Details - ",SAMPLE),
        tableGrob(t.result),
        #agg.p,
        agg.ALL.p,
        agg.ALL.bar,
        raw.p,
        ncol = 4,
        #widths = c(1,2)
        widths = c(1,0.5,0.5,2),
        bottom=paste0("ATRT Classifier - For Research Only")
        
      )
    } else if(report_type == "Summary"){
      #plot_grid(tableGrob(t.result), agg.ALL.p, agg.ALL.bar = c('A', 'B','C'))
      grid.arrange(
        top = paste0("Test Result Summary - ",SAMPLE),
        tableGrob(t.result),
        arrangeGrob(arrangeGrob(agg.ALL.bar,all.horbar.p, ncol=2, widths = c(0.8,1),top=paste0("Fig.1 - Models agreement")),
                    arrangeGrob(agg.ALL.p, widths = 1, top=paste0("Fig.2 - ",plot_title," by class")),
                    ncol=1),
        
        ncol = 2,
        #widths = c(2,1,0.5),
        bottom=paste0("ATRT Classifier - For Research Only")
        #layout_matrix=matrix(c(1,2,1,3), byrow=TRUE,ncol = 2)
        
      )
     }
  } # for loop
  
  dev.off()
  
  ##### Output #####
  
  if(!is.null(prefix)){prefix <- paste0(prefix,"_")}
  
  if(print_report == TRUE){
    write.table(test_results_agg, file=paste0(out_path,"/",prefix,"test_summary_aggregate.txt"), sep = "\t", col.names = NA)
    write.table(test_results, file=paste0(out_path,"/",prefix,"test_summary.txt"), sep = "\t", col.names = NA)
  }
  
  data$test_summary <- t.result.list
  return(data)
} # end of plot


# install.packages("ResourceSelection")
# library(ResourceSelection)
##### Create Calibration Models #####
# v2 - added agreement to the model and glmnet
# Required column: "obs" for ground truth, 
# "Class" for multi-class prediction made from classifier, 
# "Sample" for Names for where the Avg_Probability originated, 
# "Avg_Probability" for avg. probability predicted from the classifier.
#' @param cal_labels.df
#' @return Cal_models list and trained models saved as RDS
#' @example
#' library(ResourceSelection)
#' nano.calibrate(cal_labels.df)
nano.cal_model <- function(prefix, cal_labels.df, method=c("glm","glmnet")) {
  n_data_pt <- nrow(cal_labels.df)
  print(paste0("[MSG] Using ",n_data_pt, " data points for calibration"))
  Cal_models <- list()
  # training a logistic regression model
  #pdf(file = paste0("0_",prefix,"_calibration_model_plots.pdf"), width = 10.5, height = 8)
  
  head(cal_labels.df)
  # Sample   Class Avg_Probability    obs Agreement
  # 1 ATRT115  Group1       0.8669274 Group1         1
  # 2 ATRT115 Group2A              NA Group1         0
  for (c in levels(cal_labels.df$obs)) {
    print(c)
    cal_labels.df.i <- cal_labels.df[cal_labels.df$Class == c, ]
    cal_labels.df.i$obs <- as.character(cal_labels.df.i$obs)
    cal_labels.df.i[is.na(cal_labels.df.i$Avg_Probability),]$Avg_Probability<- 0
    
    if (method == "glm") {
      cal_labels.df.i[cal_labels.df.i$obs != c,]$obs <- 0
      cal_labels.df.i[cal_labels.df.i$obs == c,]$obs <- 1
      
      #head(cal_labels.df.i)
      # Sample  Class Avg_Probability obs Agreement
      # 1  ATRT321 Group1       0.5495129   1       0.8
      # 4  ATRT321 Group1       0.5672007   1       1.0
      
      # 1,0 gives opposite trend,
      #cal_labels.df.i$obs <-factor(cal_labels.df.i$obs, levels = c(1,0))
      # 0,1
      cal_labels.df.i$obs <-factor(cal_labels.df.i$obs, levels = c(0,1)) 
      # Avg_Prob only
      model.i <- glm(obs~Avg_Probability,data = cal_labels.df.i,family = binomial, control =  list(maxit = 100))
      # Avg_Prob and Agreement (turns out not significant for G2 - glm.fit: fitted probabilities numerically 0 or 1 occurred )
      #model.i. <- glm(obs~Avg_Probability + Agreement,data = cal_labels.df.i,family = binomial, control =  list(maxit = 100))
      
      Cal_models[[c]]  <- model.i
      #### Evaluate models ####
      
      # Hosmer-Lemeshow Goodness of Fit
      # How well our model fits depends on the difference between the model and the observed data. 
      #hoslem.test(mtcars$vs, fitted(model.i)) 
      
      plot(model.i) # lots of plots
      summary(model.i)
      # Estimate  Positive or Negative influence
      # Pr        Wald test to check if explanatory variables in a model are significant.
      # null deviance   how well the response variable is predicted by a model that includes only the intercept (grand mean)
      # residual deviance
      # Akaike Information Criterion (AIC)   for between model comparisonss. based on the Deviance, but penalizes you for making the model more complicated. 
      
      
      #x_prob <- seq(0,1,0.01)
      #y_prob <- predict(model.i, list(Avg_Probability = x_prob), type = "response")
      #y_prob <- predict(model.i, list(Avg_Probability = x_prob, Agreement = x_prob), type = "response")
      
      # plot(x_prob, y_prob, pch = 16, xlab = "Avg_Probability_Agreement", ylab = prob, main = paste0(c, " Calibration Model Simulation (",n_data_pt,")"))
      
    }  else  if (method == "glmnet") {
      cal_labels.df.i[cal_labels.df.i$obs != c,]$obs <- "NonClass"
      cal_labels.df.i[cal_labels.df.i$obs == c,]$obs <- "Class"
      
      cal_labels.df.i$obs <-factor(cal_labels.df.i$obs, levels = c("NonClass","Class")) 
      
      cal_labels.df.i.mat <- cal_labels.df.i[,c("Avg_Probability"),drop=FALSE]
      cal_labels.df.i.mat$Constant <- rep(1,nrow(cal_labels.df.i.mat)) # glmnet requires multicolumn, add 1 Constant 
      cal_labels.df.i.class <- cal_labels.df.i$obs
      
      set.seed(849)
      ctrl.glmnet <- trainControl(method = "repeatedcv",
                                  repeats = 10,  # rep20 feature 30-2 CV10, takes about 40 minutes for rf
                                  number = 10,
                                  classProbs = TRUE, # save class probability for things like ROC curve
                                  #seeds = seeds, # see above
                                  savePredictions= TRUE,
                                  summaryFunction=twoClassSummary)  # for CaretEnsemble
      
      #train_model.glmnet <- train(cal_labels.df.i[,c("Avg_Probability","Agreement")],
      train_model.glmnet <- train(as.matrix(cal_labels.df.i.mat),
                                  cal_labels.df.i.class ,
                                  method = "glmnet",
                                  tuneLength = 10,# by default the function will tune through three values of each tuning parameter.
                                  trControl = ctrl.glmnet,
                                  metric="ROC" #The metric "Accuracy" was not in the result set. ROC will be used instead.
                                  #preProc = c("center", "scale")
      )
      
      Cal_models[[c]] <- train_model.glmnet
    } # glmnet
    
    
  } # for levels loop
  #dev.off()
  Cal_models$Data <- cal_labels.df
  
  #saveRDS(Cal_models, file=paste0(prefix,"_calibration_models.RDS"))
  return(Cal_models)
}



##### Calibrate with platt scaling #####
# v2 - added agreement as calibration parameter and fixes bug with calculating avg_cal_probability for ALL and voting scheme
# added glmnet
# Calibrate "Avg_Probability"
# Required column: "obs" for ground truth, "Class" for multi-class prediction made from classifier, "Sample" for Names for where the Avg_Probability originated, Avg_Probability for avg. probability predicted from the classifier.
#' @param cal_labels.df
#' @param out_path outpath location
#' @return Cal_models list and trained models saved as RDS
#' @example
#' library(ResourceSelection)
#' nano.calibrate(cal_labels.df)
nano.calibrate <- function(data, Cal_models, print_report = FALSE, method=c("glm","glmnet"), out_path=NULL){
  if(is.null(grep("Avg_Cal_Probability", colnames(data$test_results_agg)))){
    stop("[MSG] Data set has been calibrated already")
  }
  
  if(is.null(out_path)){
    out_path = paste(getwd(),data$run_info$run_id[1], sep = "/")
  }
  
  print(paste0("[MSG] Calibrate Avg_Probability..."))
  test_results_agg <- data$test_results_agg
  
  test_results_agg$Avg_Cal_Probability <- NA
  test_results_agg$Avg_Probability[is.na(test_results_agg$Avg_Probability)] <- 0
  test_results_agg <- test_results_agg[test_results_agg$Num_Features != "ALL",] # calculate ALL after calibration
  
  for (c in unique(test_results_agg$Class)){
    print(c)
    if (method == "glm"){
      print("[MSG] Using glm model")
    #Calibrate using Avg_Probability only
    test_results_agg[test_results_agg$Class == c,"Avg_Cal_Probability"] <- predict(Cal_models[[c]],
                                                                                   test_results_agg[test_results_agg$Class == c,"Avg_Probability",drop=FALSE],
                                                                                     type = "response")
    
    # Calibrate using Avg_Probability and Agreement
    # test_results_agg[test_results_agg$Class == c,"Avg_Cal_Probability"] <- predict(Cal_models[[c]], 
    #                                                                                test_results_agg[test_results_agg$Class == c,c("Avg_Probability","Agreement"),drop=FALSE],
    #                                                                               type = "response")
    } else if(method == "glmnet"){
      print("[MSG] Using glmnet model")
      test_results_agg$Constant <- rep(1,nrow(test_results_agg))
      calibrated_results <- predict(Cal_models[[c]], newdata = test_results_agg[test_results_agg$Class == c,c("Avg_Probability","Constant"),drop=FALSE], type = "prob")
      test_results_agg[test_results_agg$Class == c,"Avg_Cal_Probability"]  <- calibrated_results$Class
      test_results_agg <- subset(test_results_agg, select=-Constant)

    } # glmnet
     
  }
  ##### Calculate ALL #####
  test_results_agg.ALL.i.max_full <- data.frame()
  for (SAMPLE in unique(test_results_agg$Sample)){
    test_results_agg.ALL.i <- test_results_agg[test_results_agg$Sample == SAMPLE,,drop=FALSE]
    test_results_agg.ALL.i.data <- test_results_agg.ALL.i[,c("Class","Avg_Probability","Avg_Cal_Probability"),drop=FALSE]
    test_results_agg.ALL.i.data.agg <- aggregate(test_results_agg.ALL.i.data[,-1] , by=list(test_results_agg.ALL.i.data[,1]), FUN = mean)
    Count.agg <- aggregate(test_results_agg.ALL.i[,c("Count","N_models")], by=list(test_results_agg.ALL.i$Class), FUN= sum)
    test_results_agg.ALL.i.max <- merge(test_results_agg.ALL.i.data.agg,Count.agg, by="Group.1")
    test_results_agg.ALL.i.max$Agreement <- test_results_agg.ALL.i.max$Count/test_results_agg.ALL.i.max$N_models
    test_results_agg.ALL.i.max <- data.frame(SAMPLE,test_results_agg.ALL.i.max,"ALL")
    # rearrange #
    test_results_agg.ALL.i.max <- test_results_agg.ALL.i.max[,c(1,2,5,8,3,6,7,4)]
    colnames(test_results_agg.ALL.i.max) <- colnames(test_results_agg.ALL.i)
    test_results_agg.ALL.i.max_full <- rbind(test_results_agg.ALL.i.max_full,test_results_agg.ALL.i.max)
  }
  test_results_agg <- rbind(test_results_agg,test_results_agg.ALL.i.max_full)
  ## Format ##
  test_results_agg$Avg_Cal_Probability <- round(test_results_agg$Avg_Cal_Probability, digits = 7)
  test_results_agg$Avg_Cal_Probability <- format(test_results_agg$Avg_Cal_Probability, scientific = FALSE)
  test_results_agg$Avg_Cal_Probability <- as.numeric(test_results_agg$Avg_Cal_Probability)
  idx.raw <- grep("Avg_Probability",colnames(test_results_agg))
  idx.cal <- grep("Avg_Cal_Probability",colnames(test_results_agg))
  # rearrange
  test_results_agg <- test_results_agg[,c(1:idx.raw,idx.cal,(idx.raw+1):(idx.cal-1))]
  
  ##### update test_results ######
  test_results_agg.ALL <- test_results_agg[test_results_agg$Num_Features == "ALL",]
  test_results_agg.ALL.full <- data.frame()
  for (SAMPLE in unique(test_results_agg.ALL$Sample)){
    test_results_agg.ALL.i <- test_results_agg.ALL[test_results_agg.ALL$Sample == SAMPLE,]
    # Hard voting (Count priority) - Sort by count then Avg_Probability
    test_results_agg.ALL.i <- test_results_agg.ALL.i[with(test_results_agg.ALL.i, order(Count,Avg_Cal_Probability, decreasing = TRUE)),][1,]
    # Soft voting (Probability priority) - Sort by Avg_Probability then count
    #test_results_agg.ALL.i <- test_results_agg.ALL.i[with(test_results_agg.ALL.i, order(Avg_Cal_Probability, Count,decreasing = TRUE)),][1,]
    test_results_agg.ALL.full <- rbind(test_results_agg.ALL.full,test_results_agg.ALL.i)
  }
  
  ### Output ###
  
  if(!is.null(prefix)){prefix <- paste0(prefix,"_")}
  
  if(print_report == TRUE){
    write.table(test_results_agg, file=paste0(out_path,"/",prefix,"test_summary_aggregate.txt"), sep = "\t", col.names = NA)
    write.table(test_results_agg.ALL.full, file=paste0(out_path,"/",prefix,"test_summary.txt"), sep = "\t", col.names = NA)
  }
  data$test_results_agg <- test_results_agg
  data$test_results <- test_results_agg.ALL.full
  return(data)
  
} # end of nano.calibrate

### Open Directory Interactively ###
choose_directory = function(caption = 'Select data directory') {
  
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption) 
  } else {
    tk_choose.dir(caption = caption)
  }
}

#### Plot MDS ####
# v2 added function to plot test_results Group colour if present
# General-purpose data MDS plotting.
# Require nano.set.colour() to be ran first to determine colour_code
#' @param prefix
#' @param data
#' @param plot_type
#' @param data_name
#' @return plot
#' @example
#'  library(limma) #plotMDS
#'  library(plotly)
#'  library(reshape2)
#' nano.MDS(prefix = project_prefix, data=train, plot_type = "ggplot",data_name = "norm.t")


nano.MDS <- function(prefix, data, plot_type = c("boxplot","plot","ggplot","plotly"), data_name = c("norm.t","train.data.main","train.data.validate")){
  
  data.df <- data[[data_name]]
  
  if (data_name %in% c("train.data.main","train.data.validate")){
    if(is.null(data[[data_name]]$Group)){
      stop("[MSG] Data must have Group labels. Did you run nano.trainsplit()?")
    }
    group <- data.df[,"Group",drop=FALSE]
    data.df <- data.frame(t(subset(data.df, selec = -Group)))
  } else if (data_name %in% "norm.t"){
    if(!is.null(data[["test_results"]])){
      print(paste0("[MSG] Applying test_results to plot..."))
      group <- data.frame(matrix(nrow=nrow(data.df), ncol=1, data$test_results$Class))
    } else {
      group <- data.frame(matrix(nrow=nrow(data.df), ncol=1, rep(NA,nrow(data.df))))
    } 
    colnames(group) <- "Group"
    row.names(group) <- row.names(data.df)
    data.df <- data.frame(t(data.df))
  }
  
  # assign colour by colour_code or "black" if not present
  if(is.null(data$colour_code)){
    print(paste0("[MSG] colour code not detected, using default colours"))
    groupcol <- rep("black",nrow(data.df))
    col_code <- data.frame(Group=as.factor(c("Samples")),Group_Colour=c("black"), stringsAsFactors = FALSE)
  } else {
    col_code <- data$colour_code
    groupcol <- col_code[group$Group,"Group_Colour"]
    #groupcol[is.na(groupcol)] <- "nogroup"
    groupcol[is.na(groupcol)] <- "black"
  }

  
  # sample N check #
  if(ncol(data.df) <3){
    print("[MSG] Data must have minimum 3 samples to run nano.MDS().")
  } else {
      
      library(limma) #plotMDS
      library(plotly)
      library(reshape2)
      
      ##### Plots #####
      PlotTitle <- prefix

      # obtain MDS matrix #
      mds <- limma::plotMDS(data.df,pch=19, main=PlotTitle)
      group[is.na(group)] <- "black"
      #mds.anno <- merge(mds$cmdscale.out,group, by="row.names") # limma depreciated since 3.48.0 
      mds_2d_matrix <- data.frame(x=mds$x, y=mds$y)
      row.names(mds_2d_matrix) <- row.names(mds$distance.matrix.squared)
      mds.anno <- merge(mds_2d_matrix,group, by="row.names")
      
      colnames(mds.anno) <- c("Sample","X","Y","Group")
      mds.p <- ggplot(mds.anno, aes(x=X, y=Y, label=Sample, color=Group)) + 
        geom_point(size=2) +
        scale_color_manual(values = as.character(col_code$Group_Colour)) +
        #scale_color_manual(values = as.character(groupcol)) + # not correct
        xlab(label = "MDS1")+
        ylab(label = "MDS2")+
        labs(title = paste0(PlotTitle," MDS Plot"))+
        theme_minimal() +
        theme(aspect.ratio=1) + #squre shape
        theme(text = element_text(size = 16)) +
        theme(axis.ticks = element_line(size = 0.5))+
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_rect(colour = "black", size=1))

      mds.p.plotly <- ggplotly(mds.p)
      ## box plot ##
      if (plot_type == "boxplot"){
        boxplot(data.df,range=0,ylab="log2 intensity") #intensity runs from 5 to 16 on log2 scale
      }
      ## MDS plot ##
      if (plot_type == "plot"){
        limma::plotMDS(data.df, pch=19, col = groupcol, main=PlotTitle) # point and color
        limma::plotMDS(data.df,labels=colnames(data.df), pch=19, cex=0.5) # labels and color
      }
      ## MDS ggplot ##
      if (plot_type == "ggplot"){
        #mds.p <- mds.p + geom_label(nudge_y=nudge_y_value)
        print(mds.p)
      } else if (plot_type == "ggplot_label"){
        mds.p <- mds.p +geom_text_repel(aes(label = Sample))
        print(mds.p)
      }
      
      
      ## MDS plotly ##
      if (plot_type == "plotly"){
        print(mds.p.plotly)
      }
  }
}


##### Set Colour Code #####
# Assign colour_code to data. Can assign Group_names or auto search by providing 'NULL' in Group_names which searches for the 'Group' column in the specified dataframe
nano.set.colour <- function(data, Group_names = NULL, data_name = c("train.data.main","train.data.validate")){
  Group_lables <- vector()
  if (is.null(Group_names)) {
      if(data_name == "train.data.main" & is.null(data$train.data.main)){
        stop("[MSG] Data must have Group labels. Did you run nano.trainsplit()?")
      }
      if(data_name == "train.data.validate" & is.null(data$train.data.validate)){
        stop("[MSG] Data must have Group labels. Did you run nano.trainsplit()?")
      }
      data.df <- data[[data_name]]
      if (is.null(data.df$Group)){
        stop("[MSG] Data frame don't have Group labels")
      }else {Group_lables <- unique(data.df$Group)}
      
  } else {
    Group_lables <- Group_names
  }
  if (all(Group_lables %in% c("Group1","Group2"))){
    print(paste0("[MSG] Torchia et. al., 2015 ATRT Subgroups detected:"))
    print(Group_lables)
    col_code <- data.frame(Group=as.factor(c("Group1","Group2")),Group_Colour=c("red","blue"), stringsAsFactors = FALSE)
  } else if (all(Group_lables %in% c("Group1","Group2A","Group2B"))){
    print(paste0("[MSG] Torchia et. al., 2016 ATRT Subgroups detected:"))
    print(Group_lables)
    col_code <- data.frame(Group=as.factor(c("Group1","Group2A","Group2B")),Group_Colour=c("red","blue","green"), stringsAsFactors = FALSE)
  } else if (all(Group_lables %in% c("SHH","TYR","MYC"))){
    print(paste0("[MSG] Ho et. al., 2019 ATRT Subgroups detected:"))
    print(Group_lables)
    col_code <- data.frame(Group=as.factor(c("SHH","TYR","MYC")),Group_Colour=c("4074E5","DD1D06","23AE2E"), stringsAsFactors = FALSE)
  } else {
    print(paste0("[MSG] Using custom lables:"))
    print(Group_lables)
    col_code <- data.frame(Group=as.factor(Group_lables),Group_Colour=as.numeric(Group_lables))
  }
  data$colour_code <- col_code
  return(data)
}



##### Evaluate results #####
# working in progress. problem with confusionmatrix require more than 2 classes to work properly
# creates confusion nmatrix and stats in excel
# requires nano.plot to work properly
# 'Class' reserved for prediction results
# 'Subgroup' reserved for ground truth
#' @param prefix
#' @param use_class
#' @param Prob_range
#' @param prob
#' @param anno_table
#' @param GeoMean_thres
#' @param out_path output location
#' @param in_path input location for *_test_summary.txt. Default to 

nano.eval.test <- function(prefix, use_class, Prob_range, prob = "Avg_Probability", anno_table, GeoMean_thres=0, out_path=NULL, input_path=getwd()){
  
  if(is.null(out_path)){
    out_path = paste(getwd())
  }
  
  anno <- anno_table # subgroup list
  summary_file.df <- summary_file.agg.df <- summary_file.full.df <- data.frame()
  Test_Summary_Overall <- Test_Summary_Stats <- Test_Summary <- conf_matrix_list <-list()
  WD <- prefix
  ## Load Recursively ##
  # test summary
  for (summary_file in list.files(path=in_path, pattern = paste0(prefix,"*_test_summary.txt"), recursive = FALSE)){
       if(is.null(summary_file)){
         stop("[MSG] Can't find *test_summary.txt. Did you run nano.plot()?")
       }      
    summary_file.df.i <- read.table(paste0(summary_file), sep="\t", header = TRUE, row.names = 1,stringsAsFactors = FALSE)
    summary_file.df <- rbind(summary_file.df,summary_file.df.i)
  }
  # test_summary_aggregate
  for (summary_file in list.files(path=in_path, pattern = paste0(prefix,"*_test_summary_aggregate.txt"), recursive = FALSE)){
    if(is.null(summary_file)){
      stop("[MSG] Can't find *_test_summary_aggregate.txt. Did you run nano.plot()?")
    } 
    summary_file.agg.df.i <- read.table(paste0(summary_file), sep="\t", header = TRUE, row.names = 1,stringsAsFactors = FALSE)
    summary_file.agg.df <- rbind(summary_file.agg.df,summary_file.agg.df.i)
  }
  #testing_summary_full
  for (summary_file in list.files(path=in_path, pattern = paste0(prefix,"*_testing_summary_full.txt"), recursive = FALSE)){
    if(is.null(summary_file)){
      stop("[MSG] Can't find *_testing_summary_full.txt. Did you run nano.plot()?")
    } 
    summary_file.full.df.i <- read.table(paste0(summary_file), sep="\t", header = TRUE, row.names = 1,stringsAsFactors = FALSE)
    summary_file.full.df <- rbind(summary_file.full.df,summary_file.full.df.i)
  }
  
  ## Annotate ##
  #summary_file.df.anno.full <- merge(summary_file.df, anno, by.x="Sample",by.y="Nanostring_ID_fixed")
  summary_file.df.anno.full <- merge(summary_file.df, anno, by.x="Sample",by.y="nano_filename")
  row.names(summary_file.df.anno.full) <- summary_file.df.anno.full[,1]
  
  # change factor order #

  ##  summary_file.df.anno.full$Class <- factor(summary_file.df.anno.full$Class, levels = c("Group1","Group2A","Group2B"))
  #summary_file.df.anno.full$Class <- factor(summary_file.df.anno.full$Class, levels = c("Group1","Group2","NA"))
  summary_file.df.anno.full$Class <- factor(summary_file.df.anno.full$Class, levels = use_class)
  
  ## summary_file.df.anno.full$Subgroup <- factor(summary_file.df.anno.full$Subgroup, levels = c("Group1","Group2A","Group2B"))
  #summary_file.df.anno.full$Subgroup <- factor(summary_file.df.anno.full$Subgroup, levels = c("Group1","Group2","NA"))
  summary_file.df.anno.full$Subgroup <- factor(summary_file.df.anno.full$Subgroup, levels = use_class)
  
  summary_file.df.anno.full$Matching_class <- ifelse(summary_file.df.anno.full$Subgroup == summary_file.df.anno.full$Class, 1, 0)
  
  #summary_file.agg.df.anno <- merge(summary_file.agg.df, anno, by.x="Sample", by.y="Nanostring_ID_fixed")
  summary_file.agg.df.anno <- merge(summary_file.agg.df, anno, by.x="Sample", by.y="nano_filename")  
  #summary_file.full.df.anno <- merge(summary_file.full.df, anno, by.x="Sample", by.y="Nanostring_ID_fixed")
  summary_file.full.df.anno <- merge(summary_file.full.df, anno, by.x="Sample", by.y="nano_filename")
  ## Remove Failed Samples ##
  
  # skip for now
  summary_file.df.anno <- summary_file.df.anno.full
  #summary_file.df.anno <- summary_file.df.anno.full[summary_file.df.anno.full$GeoMean >= GeoMean_thres,]
  
  
  # frozen / ffpe filter #
  #select_matrials <- c("Frozen","extracted_RNA")
  # select_matrials <- "FFPE"
  # select_matrials <- "Frozen"
  # select_matrials <- "extracted_RNA"
  #summary_file.df.anno <-summary_file.df.anno[summary_file.df.anno$RNA.Material.Used.FFPE.Frozen %in%select_matrials, ]
  
  N_Total <- nrow(summary_file.df.anno.full)
  N_QC_Pass <- nrow(summary_file.df.anno)
  
  ##### Calcuate TP/TF #####
  
  Accuracy_Table <- Subgroup_Accuracy_Table <- data.frame()
  for (Prob in Prob_range){
    summary_file.df.anno.i <- summary_file.df.anno[summary_file.df.anno[,prob] >= Prob,]
    
    N_Pass_Prob <- nrow(summary_file.df.anno.i)
    confmat.i <- confusionMatrix(summary_file.df.anno.i$Class, summary_file.df.anno.i$Subgroup, positive = NULL)
    #conf_matrix_list[[eval(paste0(alg,"_",i))]] <- confmat.i
    conf_matrix_list[[as.character(Prob)]] <- confmat.i
    
    # overall level #
    Accuracy_Table.i <- data.frame(Probability = Prob, N=N_Pass_Prob, t(data.frame(confmat.i$overall)))
    Accuracy_Table <- rbind(Accuracy_Table,Accuracy_Table.i)
    
    # subgroup level #
    Subgroup_Accuracy_Table.i <- data.frame(Probability = Prob, Class = gsub("Class: ","",row.names(confmat.i$byClass)), confmat.i$byClass)
    GRP_Count <- data.frame()
    
    for(GRP in sort(unique(Subgroup_Accuracy_Table.i$Class))){
      Summary_file.df.anno.i.GRP.True <- summary_file.df.anno.i[summary_file.df.anno.i$Subgroup == GRP,,drop=FALSE]
      N_Pass_Prob.GRP.True <- nrow(Summary_file.df.anno.i.GRP.True)
      TP <- nrow(Summary_file.df.anno.i.GRP.True[Summary_file.df.anno.i.GRP.True$Class == GRP,])
      FN <- nrow(Summary_file.df.anno.i.GRP.True[Summary_file.df.anno.i.GRP.True$Class != GRP,])
      
      Summary_file.df.anno.i.GRP.False <- summary_file.df.anno.i[summary_file.df.anno.i$Subgroup != GRP,,drop=FALSE]
      N_Pass_Prob.GRP.False <- nrow(Summary_file.df.anno.i.GRP.False)
      TN <- nrow(Summary_file.df.anno.i.GRP.False[Summary_file.df.anno.i.GRP.False$Class != GRP,])
      FP <- nrow(Summary_file.df.anno.i.GRP.False[Summary_file.df.anno.i.GRP.False$Class == GRP,])
      
      GRP_Count.i <- data.frame(Class=GRP,N_Class=N_Pass_Prob.GRP.True, TP=TP,TN=TN,FP=FP,FN=FN,FPR=FP/(FP+TN),TRP=TP/(TP+FN))
      GRP_Count <- rbind(GRP_Count, GRP_Count.i)
    }
    
    Subgroup_Accuracy_Table.i <- merge(Subgroup_Accuracy_Table.i, GRP_Count, by="Class")
    Subgroup_Accuracy_Table <- rbind(Subgroup_Accuracy_Table,Subgroup_Accuracy_Table.i)
    
    
    
  } # prob
  
  ##### Output to list #####
  
  Accuracy_Table$N_Total <- N_Total
  Accuracy_Table$N_QC_Pass <- N_QC_Pass
  Accuracy_Table$QC_Pass_Classified <- Accuracy_Table$N/Accuracy_Table$N_QC_Pass
  # save to list. Limit number of char in name
  Test_Summary_Stats[[substr(WD,start = 1, stop = 27)]] <- Accuracy_Table
  
  
  #Test_Summary_Stats[[paste0(WD,"_G")]] <- Subgroup_Accuracy_Table # save entire table
  # split table by subgroup #
  for(GRP in sort(unique(Subgroup_Accuracy_Table$Class))){
    GRP_name <- gsub("Group","",GRP)
    Test_Summary_Stats[[paste0(substr(WD,start = 1, stop = 27),"_G",GRP_name)]] <- Subgroup_Accuracy_Table[Subgroup_Accuracy_Table$Class == GRP,,drop=FALSE]
  }
  
  Test_Summary[[paste0(substr(WD,start = 1, stop = 27))]] <-summary_file.df.anno.full
  Test_Summary[[paste0(substr(WD,start = 1, stop = 27),"_agg")]] <-summary_file.agg.df.anno
  Test_Summary[[paste0(substr(WD,start = 1, stop = 27),"_full")]] <-summary_file.full.df.anno
  
  #} # end of folder loop
  
  length(Test_Summary_Stats)
  # export as excel #
  #Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")
  openxlsx::write.xlsx(Test_Summary_Stats, file = paste0(out_path,"/",prefix,"_Test_Summary_Stats.xlsx"))
  openxlsx::write.xlsx(Test_Summary, file = paste0(out_path,"/",prefix,"_Test_Summary.xlsx"))
  saveRDS(conf_matrix_list, file = paste0(out_path,"/",prefix,"_conf_matrix_list.RDS"))
  
  Test_Summary_Overall[["overall_accuracy"]] <- sum(Test_Summary[[1]]$Matching_class) / length(Test_Summary[[1]]$Matching_class)
  Test_Summary_Overall[["confusion_matrix"]] <- conf_matrix_list[[1]]$table

  return(Test_Summary_Overall)
}


#### Plot MDS ####
# v2 extracting group membership from test results. Supply "memberships" to override results from training and testing
# Plot both training (train.data.main) and testing (norm.t) data by merging main.data and test results.
# Require nano.set.colour() to be ran first to determine colour_code
# Require nano.test() to be able to obtain testing results from test.df
# membership and gene_list are optional
#' @param prefix
#' @param data_train - normalized data used for training
#' @param data_test - normalized data used for testing
#' @param memberships - dataframe of sample name (row.names) and subgroup. NULL (default)
#' @param gene_list - select genes for plotting. NULL (default) for no filtering
#' @param plot_type - different plot types
#' @param train_ellipse - add ellipse around training samples when set to TRUE. FALSE (default)
#' @return plot
#' @example
#'  library(limma) #plotMDS
#'  library(plotly)
#'  library(reshape2)
#' nano.MDS(prefix = prefix, data=train, plot_type = "ggplot",data_name = "norm.t")

nano.MDS.train.test <- function(prefix, train.data , test.data , colour_code, plot_type = c("plot","ggplot","ggplot_label", "ggplot_label_test", "plotly"), train_ellipse=FALSE, memberships=NULL, gene_list=NULL, omit_sample=NULL){

  data_train <- train.data$train.data.main
  data_test <- test.data$norm.t

  if(is.null(memberships)){
    
      # get memberships from training data
      data_train_memberships <- data_train[,"Group",drop=FALSE]
      
      # get memberships from testing data
      if (is.null(test.data$test_results)){
          stop("[MSG] Missing testing results. Did you run nano.test()?")
      } else {
          data_test_memberships <- data.frame(test.data$test_results[,"Class",drop=FALSE])
          colnames(data_test_memberships) <- "Group"
          row.names(data_test_memberships) <- test.data$test_results[,"Sample"]
      }
      # merge memberships
      group_memberships <- rbind(data_train_memberships, data_test_memberships)
      
  } else { # when membership is provided
      group_memberships <- memberships
  }
  ## strip away Group ##
  if(!is.null(data_train$Group)){data_train$Group = NULL}
  if(!is.null(data_test$Group)){data_test$Group = NULL}
  data_train$Type <- "Train"
  data_test$Type <- "Test"
  
  # check gene names in column and combine data #

  if(all(colnames(data_train) %in% colnames(data_test))){
    data_test <- data_test[,match(colnames(data_test),colnames(data_train))] # sort according to training
    data.df <- rbind(data_train,data_test)
  }else {
    stop("[MSG] Mismatch gene set in data_train and data_test!")
  }

  # anno data.df #
  data.df <- merge(data.df,group_memberships, by="row.names")
  row.names(data.df) <- data.df[,1]
  data.df <- data.df[,-1]
  data.df$Type <- factor(data.df$Type, levels = c("Train","Test"))

  # select genes #
  if (is.null(gene_list)){
      gene_list <- colnames(subset(data.df, select=-c(Type, Group)))
      print(paste0("[MSG] Using default gene list m=",length(gene_list)))
  } else {
      print(paste0("[MSG] Using gene list m=",length(gene_list)))
  }

  ##### Plots #####
  print(paste0(gene_list))
  data.df <- data.df[,c(gene_list,"Group","Type")]
    
  col_code <- colour_code
  groupcol <- col_code[data.df$Group,"Group_Colour"]
  group <- data.df[,c("Group","Type")]

  library(limma) #plotMDS
  library(plotly)
  library(reshape2)
  
  PlotTitle <- prefix
  if(!is.null(omit_sample)){
    print(paste0("[MSG] Omitting samples from plot:"))
    print(omit_sample)
    data.df <- data.df[!(row.names(data.df) %in% omit_sample),] 
  }
  
  pdf(file = NULL) # prevent writing file
  mds <- limma::plotMDS(t(subset(data.df, select = -c(Group,Type))),pch=19, main=PlotTitle)
  dev.off()
  #mds.anno <- merge(mds$cmdscale.out,group, by="row.names")# limma depreciated since 3.48.0 
  mds_2d_matrix <- data.frame(x=mds$x, y=mds$y)
  row.names(mds_2d_matrix) <- row.names(mds$distance.matrix.squared)
  mds.anno <- merge(mds_2d_matrix,group, by="row.names")
  colnames(mds.anno) <- c("Sample","X","Y","Group","Type")
  
  mds.p <- ggplot(mds.anno, aes(x=X, y=Y, label=Sample, color=Group,group=Type)) + 
    geom_point(aes(shape=Type), size=3) +
    scale_color_manual(values = as.character(col_code$Group_Colour)) +
    scale_shape_manual(values=c(1, 19)) +
    xlab(label = "MDS1")+
    ylab(label = "MDS2")+
    labs(title = paste0(PlotTitle," Combined Train and Test MDS Plot"))+
    theme_minimal() +
    theme(aspect.ratio=1) + #squre shape
    theme(axis.ticks = element_line(size = 0.5))+
    theme(text = element_text(size = 16)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1))
  
  ### Output ###
  
  #add Ellipse#
  if(train_ellipse==TRUE){
    # The t-distribution is used as an alternative to the normal distribution when sample sizes are small in order to estimate confidence or determine critical values that an observation is a given distance from the mean
    mds.p <- mds.p + stat_ellipse(data=mds.p$data[mds.p$data$Type == "Train",], aes(group=Group),type = "t", show.legend = FALSE, linetype=4) 
  }

  ## MDS plotly ##
  if (plot_type == "plotly"){
    mds.p.plotly <- ggplotly(mds.p)
    print(mds.p.plotly)
  } else if (plot_type == "ggplot"){
    print(mds.p)
  } else if (plot_type == "ggplot_label"){
    mds.p <- mds.p +geom_text_repel(aes(label = Sample), show.legend = FALSE)
  } else if (plot_type == "ggplot_label_test"){
    mds.p <- mds.p +geom_text_repel(data=mds.p$data[mds.p$data$Type == "Test",], aes(label = Sample), show.legend = FALSE)
  }

} 



###### feature selection ######
# assume have Group label in main 
nano.feat.select <- function(nanostring_data){
    library(Boruta)
  
    Boruta.obj <- list()
    
    # maxRuns   increase runs to resolve tentative features
    # doTrace   get report of the progress
    
    # features selection
    boruta.nanostring_data <- Boruta(Group~., data = nanostring_data$train.data.main, doTrace = 2)
    print(boruta.nanostring_data)
    
    # take a call on tentative features
    #  simple comparison of the median feature Z-score with the median Z-score of the most important shadow feature
    boruta.nanostring_data.fix <- TentativeRoughFix(boruta.nanostring_data)
    print(boruta.nanostring_data.fix)
    
    # plots
    plot(boruta.nanostring_data.fix, xlab = "", xaxt = "n")
    lz<-lapply(1:ncol(boruta.nanostring_data.fix$ImpHistory),function(i)
      boruta.nanostring_data.fix$ImpHistory[is.finite(boruta.nanostring_data.fix$ImpHistory[,i]),i])
    names(lz) <- colnames(boruta.nanostring_data.fix$ImpHistory)
    Labels <- sort(sapply(lz,median))
    axis(side = 1,las=2,labels = names(Labels),
         at = 1:ncol(boruta.nanostring_data.fix$ImpHistory), cex.axis = 0.7)
    
    ### get list of important attributes ###
    important_genes <- getSelectedAttributes(boruta.nanostring_data.fix, withTentative = F)

    ### Extract attribute statistics ###
    boruta.nanostring_data.fix_df <- attStats(boruta.nanostring_data.fix)
    print(boruta.nanostring_data.fix_df)
    
    Boruta.obj[["Important_Genes"]] <- important_genes
    Boruta.obj[["Boruta_obj_rough_fix"]] <- boruta.nanostring_data.fix
    Boruta.obj[["Stats"]] <- boruta.nanostring_data.fix_df
    return(Boruta.obj)
}




##### extract samples from nano.obj #####
# extract sample from nano.obj and update run_info accordingly
# data  nano.obj
# keep_samples_path csv file expect no header, col one that match sample name and col 2 Subgroup
nano.extract <- function(data, keep_samples_path=NULL){
  if(is.null(keep_samples_path)){
    stop(paste0("[MSG] keep_samples_path must be provided with. CSV file expect no header, col one that match sample name and col 2 Subgroup."))
  }
  
  keep_samples <- read.table(keep_samples_path, header = FALSE, sep=",")
  keep_samples <- keep_samples[,1]
  
  data$run_info
  data$raw <- raw.obj$raw[,c(colnames(raw.obj$raw)[1:3],keep_samples)]
  data$prenorm_qc <- raw.obj$prenorm_qc[keep_samples,]
  data$norm <- raw.obj$norm[,keep_samples]
  data$norm.t <- raw.obj$norm.t[keep_samples,]
  ori_n <- nrow(raw.obj$norm.t)
  data$run_info$samples_loaded <- paste0(nrow(raw.obj$norm.t), " samples loaded.")
  
  print(paste0("[MSG] ",nrow(raw.obj$norm.t)," samples extracted from ",ori_n, " samples."))
  return(data)
  
}

# convert training sample validate data frame for test
# data = train
convert2test <- function(data){
  test <- list()
  test$norm.t <- subset(data$train.data.validate, select =-Group)
  test$prenorm_qc <- data$prenorm_qc
  test$run_info <- data$run_info
  return(test)
}
