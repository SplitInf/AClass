##### AClass #####
# TODO:
# - Fix TYPE vs object naming
# - Add subgroup-specific thresholds
# - Save probes_rank paths automatically

##### initialize.prj #####
#' @title Initialize environment and install required packages
#' @description Installs and loads required packages for AClass project. This includes both CRAN and Bioconductor packages,
#' and ensures all dependencies are available. Note: archived version of NanoStringNorm is utilized
#'
#' @param install_missing Logical. Whether to install missing packages automatically. Default is TRUE.
#' @return Prints messages about package loading status.
#' @export
initialize.prj <- function(install_missing = TRUE) {

  message(paste0("[MSG] Checking for libraries..."))

  cran_packages <- c(
    "tcltk", "rlang", "vctrs", "pillar", "scales", "ggplot2", "gower",
    "ipred", "reshape2", "ggrepel", "grid", "gridExtra", "XML", "cowplot",
    "viridis", "limma", "plotly", "colorspace", "lazyeval"
  )

  bioc_packages <- c(
    "data.table", "caret", "randomForest", "glmnet", "pamr", "klaR",
    "vsn", "ggpubr", "ResourceSelection", "Boruta"
  )

  # Install BiocManager if needed
  if (!requireNamespace("BiocManager", quietly = TRUE) && install_missing) {
    install.packages("BiocManager")
  }

  # Load CRAN packages
  for (cran_pkg in cran_packages) {
    if(!requireNamespace(cran_pkg, quietly =TRUE)){
      if (install_missing) {
        install.packages(cran_pkg, dependencies = TRUE)
      } else {
        warning(paste0("[MSG] CRAN package ",cran_pkg," is missing"))
      }
    }
    library(cran_pkg, character.only = TRUE)
  }

  # Load Bioconductor packages
  for (bioc_pkg in bioc_packages) {
    if(!requireNamespace(bioc_pkg, quietly =TRUE)){
      if (install_missing) {
        BiocManager::install(bioc_pkg, ask = FALSE, update = FALSE)
      } else {
        warning(paste0("[MSG] Bioconductor package ",bioc_pkg," is missing"))
      }
    }
    library(bioc_pkg, character.only = TRUE)
  }

  # Special case for NanoStringNorm (archived version)
  if (!requireNamespace("XML", quietly = TRUE)) {
    if (install_missing) {
      if (.Platform$OS.type == "windows") {
        install.packages("XML", type = "binary") #	Usually works for Windows. Sometimes for MacOS if can't compile
      } else {
        install.packages("XML") #linux compatibility
      }
    } else {
      warning("[MSG] XML package is missing. If using macOS, try: install.packages('XML', type = 'binary')")
    }
  }
  if (requireNamespace("XML", quietly = TRUE)) {
    library(XML)
  }


  if (!requireNamespace("NanoStringNorm", quietly = TRUE)) {
    if (install_missing) {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
      }
      devtools::install_url(
        "https://cran.r-project.org/src/contrib/Archive/NanoStringNorm/NanoStringNorm_1.2.1.tar.gz"
      )
    }
  }

  if (requireNamespace("NanoStringNorm", quietly = TRUE)) {
    library(NanoStringNorm)
  }


  message("[MSG] All libraries loaded.")
  invisible(NULL)
}

##### process.raw #####
#' @title Function to load, normalize and pre-process NanoString data
#' @description Wrapper function to process NanoString data by running `nano.load()`, `nano.prenorm.qc()`, `nano.prep()` and `nano.norm()`
#'
#' @param work_path Directory where the output folder will be created. Defaults to current working directory.
#' @param raw_path Path to NanoString data file or folder (default *NormalizedData.csv) file. Required
#' @param keep_file_path Optional. File listing samples to keep (one sample ID per row, no header).
#' @param omit_file_path Optional. File listing samples to omit (one sample ID per row, no header).
#' @param prefix Optional string to prefix output directories and files.
#' @param SampleContent Normalization method passed to `nano.norm()`. Default is `"housekeeping.geo.mean"`.
#' @param recursive_read Logical. Whether to recursively scan directories for input files. Default is `FALSE`.
#' @return A normalized and pre-processed AClass object.
#' @export
process.raw <- function(work_path=getwd(), raw_path=NULL, keep_file_path=NULL, omit_file_path=NULL, prefix=NULL, SampleContent = "housekeeping.geo.mean", recursive_read=FALSE){

  # check paths #
  if(is.null(raw_path)){stop("[MSG] raw_path missing")}

  if(is.null(omit_file_path)) {omit_file_path=""}
  if(is.null(keep_file_path)) {keep_file_path=""}

  message("[MSG] Processing raw data with parameters:")
  message(paste("work_path:", work_path))
  message(paste("raw_path:", raw_path))
  message(paste("keep_file_path:", keep_file_path))
  message(paste("omit_file_path:", omit_file_path))
  message(paste("prefix:", prefix))

  # create prefix folder #
  if(!is.null(prefix)){
    prj_prefix <- paste(prefix,format(Sys.time(), "%Y%m%d-%H%M"), sep = "_")
  } else {
    prj_prefix <- format(Sys.time(), "%Y%m%d-%H%M")
  }

  out_path <- file.path(work_path, prj_prefix)
  dir.create(out_path, showWarnings = FALSE)

  # Step 1: Load raw data
  raw.obj <- list()

  message("=nano.load=")
  raw.obj <- nano.load(
    raw_path = raw_path,
    keep_file_path=keep_file_path,
    omit_file_path=omit_file_path,
    recursive_read = recursive_read
  )

  ## check if samples loaded ##
  raw_n <- ncol(raw.obj$raw)-3

  raw.obj$run_info$run_id <- prj_prefix

  # Step 2: Pre-normalization QC
  if(raw_n > 0){
    message("=nano.prenorm.qc=")
    raw.obj <- nano.prenorm.qc(data=raw.obj, prefix = prefix, out_path = out_path)
  }

  # Step 3: Normalization
  if(raw_n > 1){
    message("=nano.norm=")
    raw.obj <- nano.norm(data=raw.obj, SampleContent = SampleContent)
  }

  # Step 4: Data preparation
  if(raw_n > 0){
    message("=nano.prep=")
    raw.obj <- nano.prep(data=raw.obj)# transpose and check for zero's
  }

  return(raw.obj)
}


##### batch.process.raw #####
#' @title Batch Processing Function of NanoString Data
#' @description A batch wrapper for `process.raw()` that allows multiple NanoString datasets
#' to be processed either individually (batch mode) or together (combined mode)
#' All parameters matche with `process.raw()` except raw_dir_path and mode
#'
#' @param work_path Directory where the output folder will be created. Defaults to current working directory.
#' @param raw_path Path to NanoString data file or folder (default *NormalizedData.csv) file. Required
#' @param keep_file_path Optional. File listing samples to keep (one sample ID per row, no header).
#' @param omit_file_path Optional. File listing samples to omit (one sample ID per row, no header).
#' @param prefix Optional string to prefix output directories and files.
#' @param SampleContent Normalization method passed to `nano.norm()`. Default is `"housekeeping.geo.mean"`.
#' @param recursive_read Logical. Whether to recursively scan directories for input files. Default is `FALSE`.
#' @param raw_dir_path  expects multiple raw directories stored within this path.
#' @param mode Processing mode that controls how data are loaded. `"batch"` option reads each directory in raw_dir_path and normalize separately before merging. `"combined"` option reads all data within raw_dir_path recursively and processes as one. Default is `"combined"`.
#' @return A normalized and pre-processed AClass object containing merged results across multiple batches or combined dataset.
#' @export
batch.process.raw <- function(work_path=getwd(), raw_dir_path, keep_file_path=NULL, omit_file_path=NULL, prefix=NULL, SampleContent = "housekeeping.geo.mean", mode="combined"){

  if(!(mode %in% c("batch","combined"))){
    stop("[MSG] mode must be batch or combined")
  }

  if (is.null(raw_dir_path)) stop("[MSG] raw_dir_path is required.")

  if(mode=="batch"){
    # Get path to csv's
    raw_dir <- as.data.frame(list.files(path = raw_dir_path, pattern = ".*NormalizedData.*.csv", recursive = TRUE, full.names=TRUE))
    raw_dir <- apply(raw_dir,1,function(x) unlist(strsplit(gsub(paste0(raw_dir_path,"/"),"",x),split = "/"))[1] )
    raw_dir <- unique(raw_dir)

    raw.obj <- list()

    for (i in raw_dir) {
      message(paste0("[MSG] Processing dir: ", i))
      raw_path_i <- file.path(raw_dir_path,i)
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
        message(a)
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
  } else if (mode=="combined"){
    raw.obj <- process.raw(work_path=work_path, raw_path=raw_dir_path, keep_file_path=keep_file_path, omit_file_path=omit_file_path, prefix=prefix, SampleContent = SampleContent, recursive_read = TRUE)
  }

  return(raw.obj)
}


##### classify.data #####
#' @title Perform ensemble classification using pre-trained model
#' @description Wrapper function to classify pre-normalized expression data (NanoString or compatible) using a pre-trained AClass model. Wraps `nano.test`, `get.nano.test.results`, `nano.set.colour`, and `nano.plot`.
#' @param work_path Directory where the output folder will be created. Default is NULL
#' @param data Pre-processed AClass object
#' @param prefix Optional string to prefix output directories and files.
#' @param training_model_obj Pre-trained model
#' @param alg_list Vector of algorithms in the training_model_obj to be used for classification. Default is `c("rf","glmnet","pam", "nb", "knn")`
#' @param keep_file_path Optional. File listing samples to keep (one sample ID per row, no header).
#' @param omit_file_path Optional. File listing samples to omit (one sample ID per row, no header).
#' @param out_path Alternative output location instead of run_id under work_path
#' @param thres_geomean threshold for report to be considered passing. Default is `100`. `NULL` indicate skipping threshold check.
#' @param report_type report format options. `Summary` or `Detailed` description. Refer to `nano.plot()`
#' @param remap_to_atrt_consensus Logical. If TRUE, will remap Torchia et al. 2016 subgroup names to ATRT consensus names from Ho et al. 2019.
#' @return A classified AClass object with test results, probability scores, assigned groups, and create visualization outputs.
#' @examples
#' data.obj <- classify.data(data = data.obj,prefix = "demo",training_model_obj = training_models)
#' @export
classify.data <- function(work_path=NULL, data, prefix, training_model_obj, alg_list = c("rf","glmnet","pam", "nb", "knn"), keep_file_path = NULL, omit_file_path = NULL, out_path=NULL, thres_geomean=100, report_type="Summary", remap_to_atrt_consensus=FALSE){

  if(is.null(prefix)) {stop("[MSG] prefix is required but missing.")}
  if(is.null(work_path)){stop("[MSG] work_path missing")}
  if(is.null(training_model_obj$training_model_list)) {
    stop("[MSG] training_model_obj is missing expected structure (no training_model_list).")
  }

  test.obj <- data

  # default to use run_id as out_path unless otherwise provided
  # Ensure test.obj$run_info$run_id exists e.g. custom dataset that didn't begin with process.raw()
  if (is.null(test.obj$run_info$run_id)) {
    test.obj$run_info$run_id <- paste0("run_", format(Sys.Date(), "%Y%m%d"))
  }
  # Set out_path if not provided
  if (is.null(out_path)) {
    out_path <- file.path(work_path, test.obj$run_info$run_id)
  }
  # create out_path if it doesn't exist
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
  }

  message(paste0("[MSG] out_path: ", out_path))

  if (!is.null(keep_file_path)) {
    test.keep <- read.table(file = keep_file_path, stringsAsFactors = FALSE, sep = "\t")
    test.keep <- apply(test.keep, 1, make.names)
    message(paste0("[MSG] Keeping ",length(test.keep), " samples"))
    test.obj$norm.t <- test.obj$norm.t[row.names(test.obj$norm.t) %in% test.keep,,drop=FALSE]
  }

  if (!is.null(omit_file_path)) {
    test.omit <- read.table(file = omit_file_path, stringsAsFactors = FALSE, sep = "\t") # files omitted
    test.omit <- apply(test.omit, 1, make.names)
    message(paste0("[MSG] Omitting ",length(test.omit), " samples"))
    test.obj$norm.t <- test.obj$norm.t[!row.names(test.obj$norm.t) %in% test.omit,,drop=FALSE]
  }

  # Ensure same gene names are used in test.obj and model
  test.obj <- check_and_prepare_test_data(test.obj, model = training_model_obj)

  # Step 6: Test
  # choose algorithms
  test.obj <- nano.test(prefix = prefix, training_model_obj = training_model_obj, data = test.obj, alg_list = alg_list, out_path = out_path) # output text file to out_path

  # Step 7: Consolidate results
  # choose min max range based on model accuracy
  test.obj <- get.nano.test.results(prefix,test.obj, out_path = out_path)
  #saveRDS(test.obj, file = file.path(out_path,paste0(prefix,"_test.data.tested.RDS")))

  # Step 8: Set Colour Code based on pre-trained models and remap_to_atrt_consensus flag
  group <- unique(training_model_obj$train.data.main$Group)
  test.obj <- nano.set.colour(test.obj, group, remap_to_atrt_consensus=remap_to_atrt_consensus)

  # Step 9: Generate report
  test.obj <- nano.plot(prefix = prefix, data = test.obj, prob= "Avg_Probability", report_type=report_type, print_report = TRUE, thres_avg_prob=0, thres_geomean = thres_geomean, out_path=out_path)
  saveRDS(test.obj, file = file.path(out_path,paste0(prefix,"_test.data.tested.RDS")))
  return(test.obj)
}


##### Load Nanostring Data #####
#' @title Load NanoString NormalizedData.csv files and return structured raw object
#'
#' @description Loads all files with "NormalizedData" in the filename and ".csv" extension, starting at row 16.
#' @param raw_path path of the directory where raw nanostring data is kept. Read recursively.
#' @param keep_file_path Optional. File listing samples to keep (one sample ID per row, no header).
#' @param omit_file_path Optional. File listing samples to omit (one sample ID per row, no header).
#' @param recursive_read Logical. Whether to recursively search subdirectories. May introduce run-to-run variation. Default is FALSE.
#' @return AClass object, which is a list of dataframe with read data and number of samples loaded.
#' @examples
#' raw.obj <- nano.load(
#'   raw_path = "data/",
#'   keep_file_path = "samples_keep.txt",
#'   omit_file_path = "samples_omit.txt"
#' )
#' @export
nano.load <- function(raw_path = getwd(), keep_file_path="", omit_file_path="", recursive_read=FALSE) {
  raw.summary <- list()
  raw.merge <- data.frame()
  for (nanofile in list.files(path=raw_path, pattern = ".*NormalizedData.*.csv", recursive = recursive_read)){

    #check#
    filepath <- file.path(raw_path, nanofile)
    if (!file.exists(filepath)) {
      stop(paste0("[MSG] File not found: ", filepath))
    }

    raw <- read.table(file.path(raw_path,nanofile), sep = ",", skip = 15, stringsAsFactors = FALSE) # skip header
    header <- c("Code.Class", "Name", "Accession")
    info <- read.table(file.path(raw_path,nanofile), sep = ",", skip = 2, stringsAsFactors = FALSE) #add column info in col1-4
    Sample_names <- info[1,4:ncol(info)] # get sample names
    colnames(raw) <- c(header,Sample_names)

    # run details
    raw.summary[["csv"]][[nanofile]][["details"]] <- list()
    raw.summary[["csv"]][[nanofile]][["details"]]$found <- paste0(dim(raw)[1]," features ", dim(raw)[2]-3, " samples.")
    raw.summary[["csv"]][[nanofile]][["details"]]$samples <- colnames(raw[,-c(1:3)])

    message(paste0("[MSG] Loading ",raw.summary[["csv"]][[nanofile]][["details"]]$found))

    if (ncol(raw.merge)==0) {
      raw.merge <- raw
    } else {
      raw.merge <- merge(raw.merge, raw, by = header)
    }
  }

  n_sample <- ifelse(ncol(raw.merge)-3 <0,0,ncol(raw.merge)-3)
  raw.summary[["samples_found"]] <- paste0(n_sample, " samples found.")

  message(paste0("[MSG] ",raw.summary[["samples_found"]]))

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
  message(paste0("[MSG] ",raw.summary[["samples_loaded"]])) #7
  raw.loaded <- list()

  raw.loaded$run_info <- raw.summary
  raw.loaded$raw <- raw.merge
  if (n_sample_loaded == 0){
    stop("[MSG] No samples were loaded. Check file omit list.")
  }
  return(raw.loaded)
}



##### nano.prenorm.qc #####
#' @title Perform data QC
#' @description This step is performed before normalization to evaluate data quality. Produces QC report and is required to calculate HK_geomean
#' @param data AClass object output from nano.load.
#' @param code_class Default "Housekeeping".
#' @param prefix project prefix
#' @return Report of the raw data and GeoMean table.
#' @examples
#' data.obj <- nano.prenorm.qc(data.obj)
#' @export
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
  message(paste0("[MSG] max_val:", max_val))
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
  message(paste0("[MSG] Exporting QC report"))

  if(!is.null(prefix)){prefix <- paste0(prefix,"_")}

  pdf(file = file.path(out_path,paste0(prefix,"Prenorm_QC","_",data_name,".pdf")), width= 8,  height = 10.5)
  grid.arrange(p.genes, p.qc, p.cv, p.mean_cv, ncol=2,
               top=textGrob("Prenorm Housekeeping Genes", gp=gpar(fontsize=15,font=8))
  )
  dev.off()

  write.table(hk.prenorm_qc, file = file.path(out_path, paste0(prefix,"Prenorm_QC","_",data_name,".txt")), sep = "\t", col.names = NA, quote = FALSE)

  data$prenorm_qc <- hk.prenorm_qc
  return(data)
}


##### nano.prep #####
#' @title Pretrain QC
#' @description Transform and check dataframe for missing values
#'
#' @param data AClass object with normalized expression matrix in \code{$norm}.
#' @return AClass object with added \code{$norm.t} (transposed normalized data). Issues a warning if near-zero variance features are detected.
#' @examples
#' data.obj <- nano.prep(data.obj)
#' @export
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
    message(norm.t.nzv[norm.t.nzv$nzv == "TRUE",])
    #stopifnot(nrow(norm.t.nzv[norm.t.nzv$nzv == "TRUE",]) == 0)
  }

  data$norm.t <- norm.t
  message(paste0("[MSG] Dataset read for testing/training"))
  return(data)

}


##### nano.trainsplit #####
#' @title Split data into training and testing dataset
#' @description Split normalized training data (NanoString or transcriptomic data) and assign 'Group' label to output
#'
#' @param data AClass object with \code{$norm.t}.
#' @param training_memberships_path Tab-delimited file with sample IDs (column 1) and group labels (column 2)
#' @param N.train.per Proportion of data to use for training (e.g., 0.8 for 80%).
#' @param seed Set to ensure reproducibility. Note it is needed to wrap function in loop
#' @return AClass object with added train.data.main and train.data.validate components
#' @examples
#' data.obj <- nano.trainsplit(data.obj, training_memberships_path, 0.8)
#' @export
nano.trainsplit <- function(data,training_memberships_path,N.train.per, seed=NULL){

  if(missing(training_memberships_path)){
    stop("[MSG] training_memberships_path missing")
  }

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

  if(is.null(seed)){
    set.seed(as.numeric(Sys.time()))
  } else if(!is.null(seed)){
    # adding random seed doesn't solve loop bug issue #
    # rnd <- sample(1:.Machine$integer.max,1)
    rnd <- seed
    set.seed(rnd)
    message(paste0("[MSG] custom seed: ",rnd))
  }

  duplicate_flag = TRUE
  while (duplicate_flag==TRUE){
    train.idx <- createDataPartition(y = train.data$Group, ## the outcome data are needed for random sampling
                                     p = N.train.per,     ## The percentage of data in the training set
                                     list = FALSE)
    train.data.main <- train.data[train.idx,,drop = FALSE]     # main training
    train.data.validate <- train.data[-train.idx,,drop = FALSE]     # training validation

    message(paste0("[MSG] Training: N",nrow(train.data.main),
                   "     Training_validatation: N",nrow(train.data.validate)))

    if(!is.null(data$train.data.main)){
      if(identical(train.data.main,data$train.data.main) & identical(train.data.validate,data$train.data.validate)){
        duplicate_flag = TRUE
        message("[MSG] Identical results obtained. Repeating...")
      } else {
        duplicate_flag = FALSE
      }
    }else{
      duplicate_flag = FALSE
    }
  }
  data$train.data.main <- train.data.main
  data$train.data.validate <- train.data.validate
  return(data)
}

##### nano.train #####
#' @title Ensemble model training
#' @description Prenorm Nanostring Data and produce report
#'
#' @param data AClass object with \code{$train.norm}. Expects column `"Group"`.
#' @param alg_list list of alg
#' @param work_path project work path
#' @param probes_rank_path path to probe/gene ranking list, one gene per line. E.g. gene list sorted by p value (most significant to least). Looks for probes_list.txt in the working directory if probes_rank_path is missing.
#' @param min_test_features Minimum number of features tested.
#' @param max_test_features Maximum number of features tested.
#' @param out_path output path. When not provided out_path will be extracted from run_info (default) and work_path
#' @return Updated AClass object containing trained models (\code{$training_model_list}), training metadata, and the training dataset. Also saves model performance matrices and confusion matrices to disk.
#'
#' @examples
#' training_models <- nano.train(
#'   prefix = "demo",
#'   data = train.obj,  # AClass object created using nano.trainsplit()
#'   work_path = "path/to/output",
#'   probes_rank_path = "your_sorted_list.txt"
#' )
#' @export
nano.train <- function(prefix, data, work_path, alg_list = c("rf","glmnet","pam", "nb", "knn"), probes_rank_path=NULL, min_test_features=20, max_test_features=30, c.method = "repeatedcv", c.repeats = 5, c.number = 10, out_path=NULL) {

  if(is.null(data$train.data.main) || is.null(data$train.data.validate)){
    stop("[MSG] train.data.main / train.data.validate missing. Did you run nano.trainsplit()?")
  }

  if (!"Group" %in% colnames(data$train.data.main)) {
    stop("[MSG] Group column not found in train.data.main. Check if nano.trainsplit() has ben ran or add manually.")
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
    if(missing(work_path)){
      stop(paste0("[MSG] out_path not provided and work_path is not defined"))
    }
    out_path = file.path(work_path,data$run_info$run_id[1])
  }

  train.settings <- list(alg_list=alg_list, probes_rank=probes.list.full, min_test_features=min_test_features,max_test_features=max_test_features,c.method =c.method , c.repeats = c.repeats, c.number = c.number )

  training_model_obj <- list()
  ###### Training Settings #####

  ctrl <- trainControl(method = train.settings$c.method,
                       repeats = train.settings$c.repeats,  # rep20 feature 30-2 CV10, takes about 40 minutes for rf
                       number = train.settings$c.number,  # x fold CV
                       classProbs = TRUE, # save class probability for things like ROC curve
                       #seeds = seeds, # implement seed assignment system
                       savePredictions= TRUE)  # for CaretEnsemble

  ### reset dataframe ##

  training_model_mat_list <- training_model_list <- c()  # stores trained models and results
  training_model_mat_list.full <- data.frame()

  for (alg in alg_list) {

    message(paste0("[MSG] training and optimizing algorithm: ",alg))

    ## different number of features for loop ##
    for (probe.idx in max_test_features:min_test_features){  ## select features based on list given
      message(paste0("[MSG] ----- ",paste0(alg,"_",probe.idx)," -----"))
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
      #saveRDS(training_model_list, file=file.path(out_path, paste0(prefix,"_Training_Models_List.tmp.RDS"))) # replaces original RDS with newer one every loop rather than in the end in case crashing

      ## show variables that were used in the final model
      message(predictors(train_model))

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
        } else if (length(bestTuneModel) == 3){
          #alg = "nb"
          train.mat$grid_index_col <- train.mat[,ncol_total - 3]  # order swapped for nb
          train.mat$grid_index_col2 <- train.mat[,ncol_total - 4] # order swapped for nb
          train.mat$grid_index_col3 <- train.mat[,ncol_total - 2]
          train.mat.best <- train.mat[train.mat[,"grid_index_col"] == bestTuneModel[[1]] & train.mat[,"grid_index_col2"] == bestTuneModel[[2]] & train.mat[,"grid_index_col3"] == bestTuneModel[[3]],]
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

      } else if(alg == "glmnet") {
        # if using incompatible alg
        if (length(bestTuneModel) == 2){
          train.mat$grid_index_col <- train.mat[,ncol_total - 3]
          train.mat$grid_index_col2 <- train.mat[,ncol_total - 2]
          train.mat.best <- train.mat[train.mat[,"grid_index_col"] == bestTuneModel[[1]] & round(train.mat[,"grid_index_col2"],5) == round(bestTuneModel[[2]],5),]
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
      }
      # confusion matrix #
      train_model.conmat <- confusionMatrix(train_model)
      #training_model_mat_list[[paste0(alg,"_",probe.idx,"_N_match")]] <- train.mat.best.simple
      training_model_mat_list[[paste0(alg,"_",probe.idx,"_matrix")]] <- train.mat.best.cast
      training_model_mat_list[[paste0(alg,"_",probe.idx,"_confmatrix")]] <- train_model.conmat
      training_model_mat_list[["Training_Result_Matrix"]] <- training_model_mat_list.full
      saveRDS(training_model_mat_list, file=file.path(out_path,paste0(prefix,"_Training_Models_Mat_List.RDS"))) # replaces original RDS with newer one every loop rather than in the end in case crashing

    } # probe.idx loop

  } # alg loop

  message(paste0("[MSG] ",length(training_model_list)," models created."))

  ##### Output #####
  training_model_obj[["training_model_list"]] <- training_model_list
  training_model_obj[["train.data.main"]] <- train.data.training_main
  training_model_obj[["train.settings"]] <- train.settings
  training_model_obj[["run_info"]] <- data$run_info

  saveRDS(training_model_obj, file=file.path(out_path,paste0(prefix,"_Training_Models_List.RDS")))

  message(paste0("[MSG] Training complete. Models saved to ", out_path))

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
  # message(gg)
  # dev.off()

}

######### nano.train.report ##########
#' @title Model training performance summary
#' @description
#' Generates summary statistics and visualizations of classification accuracy across algorithms and feature subset sizes from `nano.train()` output.
#' Produces line plots, heatmaps, and optional PDF reports to assist in evaluating the effect of probe number and model choice on training performance.
#'
#' Designed to facilitate selection of optimal classifiers and feature sizes by visualizing average accuracy trends and algorithm-specific behaviors.
#' @param prefix - file prefix that goes to algorithms performance report
#' @param training_model_obj - training model object where training performance is pulled out from
#' @param feature_min - controls the number of min plotting range (x-axis). Should lie within trained range.
#' @param feature_max - controls the number of max plotting range (x-axis). Should lie within trained range.
#' @param print_report - binary option to print as pdf or not
#' @param out_path output path. When not provided out_path will be extracted from run_info (default)
#' @param feature_box_range - controls plotting box over selected feature range in overall combined plot. Expects 4 numeric values: xmin, xmax, ymin, ymax in the format: c(x1,x2,y1,y2). Default. NULL.
#' @param annotate_alg - whether to annotate algorithms when in overall combined plot. Default NULL.
#' @param adj_y_range - controls y-axis in overall combined plot. Format: c(y_min,y_max).
#' @param add_legend - binary option to add legend to overall combined plot. Default FALSE.
#' @return A ggplot object showing overall accuracy trends by number of features and algorithm.
#' @examples
#' nano.train.report(prefix = "demo", training_model_obj = training_models, feature_min = 20, feature_max = 30)
#' @export
nano.train.report <- function(prefix, training_model_obj, feature_min, feature_max, print_report=TRUE, out_path=NULL, feature_box_range=NULL, annotate_alg=FALSE, adj_y_range=NULL, add_legend=FALSE){

  library(ggrepel)

  if(is.null(out_path)){
    out_path = paste(getwd(),training_model_obj$run_info$run_id[1], sep = "/")
  }

  train_list <- training_model_obj[["training_model_list"]]
  full_internal_performance <- internal_performance <- data.frame()

  for (model_idx in 1:length(train_list)){
    train_model <- train_list[[model_idx]]
    alg <- as.character(train_model$method)
    probe.idx <- as.numeric(gsub(paste0(alg,"_"),"",names(train_list[model_idx])))

    ## print model name ##
    #message(paste0("[MSG] Parsing ",names(train_list[model_idx]))) # comment out for brevity

    ########## Evaluate Results ##########
    #message(train_model$results) # Training Results
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
    # message(predictors(train_model)) # comment out for brevity

  } # for loop


  ##### Prepare output #####
  probe_max <- max(internal_performance$Num_Features) # max in trained data
  probe_min <- min(internal_performance$Num_Features) # min in trained data
  stopifnot(feature_max >= probe_min & feature_max <= probe_max)
  stopifnot(feature_min >= probe_min & feature_min <= probe_max)

  if(!is.null(prefix)){prefix <- paste0(prefix,"_")}

  # calculate avg #

  internal_performance$Alg <- as.factor(internal_performance$Alg)
  stats <- stats.i <- data.frame()
  for(i in probe_min:probe_max){
    for(j in probe_min:probe_max){
      #message(paste0(i," vs ",j))
      # it is ok for i==j
      if(i<=j){
        internal_performance.i <- internal_performance[internal_performance$Num_Features >= i & internal_performance$Num_Features <= j,]
      }else if(i>j){
        internal_performance.i <- internal_performance[internal_performance$Num_Features <= i & internal_performance$Num_Features >= j,]
      }
      #internal_performance.i[internal_performance.i$Num_Features > feature_min & internal_performance.i$Num_Features < feature_max ,]
      stats.i <- aggregate.data.frame(x=internal_performance.i$Accuracy, by=list(internal_performance.i$Alg),FUN = mean)
      # stats.i$Num_Features <- paste0(i,"_",j)
      stats.i$Num_Features.i <- i
      stats.i$Num_Features.j <- j
      colnames(stats.i) <- c("Alg","Avg_accuracy","Num_Features.i","Num_Features.j")
      stats <- rbind(stats,stats.i)
    } # j
  } # i


  ## plot ##
  g_conf_mat.facet <- ggplot(data =stats, aes(Num_Features.i, Num_Features.j,Avg_accuracy))+
    geom_tile(aes(fill = Avg_accuracy),colour = "white")   +
    scale_x_continuous(breaks = seq(from=probe_min, to =probe_max, by=2))+
    scale_y_continuous(breaks = seq(from=probe_min, to =probe_max, by=2))+
    scale_fill_gradientn(colours = c("cyan", "black", "red"))+
    theme_minimal() +  coord_equal(ratio = 1) +
    #facet_grid(Alg~. )
    facet_wrap(~Alg, ncol=2)


  g_conf_mat <- ggplot(data =stats, aes(x=Num_Features.i, y=Num_Features.j,z=Avg_accuracy))+
    geom_tile(aes(fill = Avg_accuracy),colour = "white")   +
    scale_x_continuous(breaks = seq(from=probe_min, to =probe_max, by=2))+
    scale_y_continuous(breaks = seq(from=probe_min, to =probe_max, by=2))+
    #scale_fill_gradient2(low = "blue",  high = "red") +
    #scale_fill_gradient2() +
    scale_fill_gradientn(colours = c("cyan", "black", "red"))+
    theme_minimal() + coord_equal(ratio = 1) +  ggtitle("Overall Average") +
    theme(axis.text.x=element_text(angle = 90, hjust = 0))


  # g_conf_mat + geom_density_2d(stats, aes(x = Num_Features.i, y = Num_Features.j, z = Avg_accuracy))
  # g_conf_mat + geom_density_2d()
  # g_conf_mat + stat_contour(stats, aes(x = Num_Features.i, y = Num_Features.j, z = Avg_accuracy))

  ## line graph ##
  internal_performance.select <- internal_performance[internal_performance$Num_Features >= feature_min & internal_performance$Num_Features <= feature_max,,drop=FALSE]
  internal_performance.select.num_features.agg <- aggregate(internal_performance.select[,"Accuracy",drop=FALSE], by=list(internal_performance.select$Num_Features), FUN=mean)
  colnames(internal_performance.select.num_features.agg) <- c("Num_Features","Avg_accuracy")

  ## calculate alg performance avg accuracy ##
  internal_performance.select.alg.agg <- aggregate(internal_performance.select[,"Accuracy",drop=FALSE], by=list(internal_performance.select$Alg), FUN=mean)
  colnames(internal_performance.select.alg.agg) <- c("Alg","Avg_accuracy")
  # add final alg performance info #
  alg.final_performance <- data.frame()
  for(i in unique(internal_performance.select$Alg)){
    alg.final_performance.i <- internal_performance.select[internal_performance.select$Alg == i,]
    alg.final_performance.i <- alg.final_performance.i[alg.final_performance.i$Num_Features == max(alg.final_performance.i$Num_Features),]
    alg.final_performance <- rbind(alg.final_performance, alg.final_performance.i)
  }
  internal_performance.select.alg.agg <- merge(internal_performance.select.alg.agg,alg.final_performance , by="Alg")

  ## plots ##
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


  gg_line.combined <- ggplot(data=internal_performance, aes(x=Num_Features, y=Accuracy, colour=Alg)) +
    geom_point(show.legend = FALSE) +
    #facet_grid(Alg~.)+
    scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.02)) +
    scale_x_continuous(breaks = seq(from = feature_min, to = feature_max ,by = 2)) +
    coord_cartesian(xlim = c(feature_min, feature_max)) +
    theme_minimal()+ geom_path(aes(colour = Alg),show.legend = FALSE)  +
    geom_hline(data = internal_performance.select.alg.agg, aes(yintercept=Avg_accuracy, group=Alg, colour=Alg), linetype = "dashed", show.legend = FALSE)+
    geom_text(data = internal_performance.select.alg.agg, aes(x=0,y=Avg_accuracy,colour="black", group=Alg, label=paste(Alg,"avg_acc",round(Avg_accuracy,3))), nudge_x=probe_min*2, cex=3,show.legend = FALSE)


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

  ## overall plot ##

  # base plot #
  gg_line.combined.overall <- ggplot(data=internal_performance, aes(x=Num_Features, y=Accuracy, colour=Alg)) +
    geom_point(show.legend = FALSE) +
    geom_path(aes(colour = Alg)) +
    #scale_color_brewer(palette="Set1") +
    # modified Set1
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3" ,"#FF7F00" ,"yellow2" ,"#A65628" ,"#F781BF", "#999999"))
  #scale_color_viridis(discrete=TRUE)

  # format #
  gg_line.combined.overall <- gg_line.combined.overall  +
    xlab(label = "Number of genes in model")+
    ylab(label = "Training accuracy")+
    labs(title = paste0("Training Accuracy Plot"))+
    scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.02)) +
    #coord_cartesian(xlim = c(feature_min, feature_max),default = TRUE,expand = FALSE) +
    #coord_cartesian(ylim = c(min(internal_performance[internal_performance$Num_Features >= feature_min & internal_performance$Num_Features <= feature_max,]$Accuracy),
    #                         max(internal_performance[internal_performance$Num_Features >= feature_min & internal_performance$Num_Features <= feature_max,]$Accuracy)))+
    #ylim(0,1)+
    theme_minimal() +
    theme(axis.ticks = element_line(size = 0.5))+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(text = element_text(size = 16))

  #add border
  gg_line.combined.overall <- gg_line.combined.overall +
    theme(panel.background = element_rect(colour = "black", size=1))
  #add alg annotation
  if(annotate_alg==TRUE){
    # point.size controls gap between segment and data
    gg_line.combined.overall <- gg_line.combined.overall + ggrepel::geom_text_repel(data = internal_performance.select.alg.agg, aes(x=Num_Features,y=Accuracy,colour=Alg,group=Alg,
                                                                                                                                    point.size = 7,
                                                                                                                                    label=paste(Alg,"avg.",round(Avg_accuracy,3))),
                                                                                    nudge_x=probe_min*2, cex=3,show.legend = FALSE, segment.size  = 0.2, segment.color = "black",segment.linetype = 1, arrow = arrow(length = unit(0.005, "npc"), type = "closed"))

    # min expand by 1, max epand by 10
    gg_line.combined.overall <- gg_line.combined.overall + scale_x_continuous(expand = expansion(mult = c(0, 0),
                                                                                                 add = c(1, 5)),
                                                                              breaks = seq(from = feature_min, to = feature_max ,by = 2))

  } else if (annotate_alg==FALSE){
    gg_line.combined.overall <- gg_line.combined.overall + scale_x_continuous(breaks = seq(from = feature_min, to = feature_max ,by = 2))

  }

  if (!is.null(adj_y_range)){
    gg_line.combined.overall <- gg_line.combined.overall+ coord_cartesian(xlim = c(feature_min, feature_max), ylim = c(adj_y_range[1],adj_y_range[2]),default = TRUE, expand = TRUE)
    #ylim(0,1)+
  } else if (is.null(adj_y_range)){
    gg_line.combined.overall <- gg_line.combined.overall+ coord_cartesian(xlim = c(feature_min, feature_max), default = TRUE, expand = TRUE)
  }
  # add box #
  if(!is.null(feature_box_range)){
    gg_line.combined.overall <- gg_line.combined.overall + geom_rect(aes(xmin = feature_box_range[1], xmax = feature_box_range[2], ymin = feature_box_range[3], ymax =feature_box_range[4]),
                                                                     fill = "transparent", color = "red", size = 0.5)
  }

  if(add_legend==TRUE){
    gg_line.combined.overall <- gg_line.combined.overall + theme(legend.position = "bottom")
  } else if (add_legend==FALSE){
    gg_line.combined.overall <- gg_line.combined.overall + theme(legend.position = "none")
  }
  ## prepare output ##
  overview_internal_performance <- get.training.stats(train_list)
  overview_internal_performance <- merge(overview_internal_performance,internal_performance.select.alg.agg, by="Alg")

  message(overview_internal_performance)

  ##### Output #####

  if (print_report == TRUE){

    message(paste0("[MSG] Check directory for detailed reports."))

    write.table(internal_performance, file = file.path(out_path,paste0(prefix,"Optimal_Training_Attributes.txt")), col.names = NA, sep = "\t")
    write.table(full_internal_performance, file = file.path(out_path, paste0(prefix,"Full_Training_Attributes.txt")), col.names = NA, sep = "\t")
    write.table(overview_internal_performance,file = file.path(out_path, paste0(prefix,"Overview_Training_Attributes.txt")), col.names = NA, sep = "\t")

    pdf(file = file.path(out_path,paste0(prefix,"Accuracy_by_Alg.pdf")), width = 10.5, height = 8)
    message(g_conf_mat.facet)
    message(g_conf_mat)
    message(gg_line.facet)
    message(gg_line.combined)
    message(gg_line)
    message(gg_line.combined.overall)
    dev.off()
  }

  return(gg_line.combined.overall)
}



######### get.training.stats ##########
#' @title Summarize training configuration
#' @description
#' Extracts and summarizes key model training settings from a `nano.train()` output AClass object, including algorithm type, cross-validation method, number of folds, number of repeats, number of trained models, and feature range tested. Used internally in `nano.train.report()`.
#'
#' @param training_model_obj - training model object where training performance is pulled out from
#' @export
get.training.stats <- function(training_model_obj){

  models = training_model_obj$training_model_list

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

######### nano.norm #########
#' @title Normalize NanoString data
#' @description
#' Applies normalization to NanoString raw data using the `NanoStringNorm` package. This function wraps common parameters for expression normalization and handles edge cases such as single-sample input. For details, refer to the help documentation of the `NanoStringNorm` package.
#' @param data AClass object with raw NanoString data stored in $raw.
#' @param SampleContent Normalization method. Default is "housekeeping.geo.mean".
#' @param round.values Logical. Whether to round normalized values. Default is FALSE.
#' @param take.log Logical. Whether to log-transform the data. Default is TRUE.
#' @param return.matrix.of.endogenous.probes Logical. Whether to return a matrix of endogenous probes. Default is TRUE.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#' @return AClass object with normalized expression matrix stored in $norm.
#' @export
nano.norm <- function(data, SampleContent = "housekeeping.geo.mean", round.values = FALSE, take.log = TRUE, return.matrix.of.endogenous.probes = TRUE, verbose = TRUE){

  library(NanoStringNorm)

  # Check for data$raw existence and format
  if (is.null(data$raw) || !is.data.frame(data$raw)) {
    stop("[MSG] 'data$raw' is missing or not a dataframe. Did you run nano.load()?")
  }

  norm <- NanoStringNorm(x = data$raw, SampleContent = SampleContent,  round.values = round.values,  take.log = take.log,   return.matrix.of.endogenous.probes = return.matrix.of.endogenous.probes, verbose = verbose)
  norm <- as.data.frame(norm)
  if(ncol(norm) == 1) {
    sample_id <- colnames(data$raw)[4] # sample name gets dropped from NanoStringNorm when there's only one sample
    colnames(norm) <- sample_id
  }

  data$norm <- as.data.frame(norm)
  return(data)

}

######### nano.test #########
#' @title Test classification model on NanoString or other transcriptomic expression data
#'
#' @description
#' Applies trained models from `nano.train()` to test data for classification. If known class labels (`Group`) are present, a confusion matrix is computed to evaluate performance. Prediction results are saved to disk and appended to the AClass object for downstream analysis.
#'
#' @param training_model_obj Object returned from `nano.train()`, must include `$training_model_list`.
#' @param alg_list Optional. Subset of algorithms to test. If not provided, all algorithms in the training object are used.
#' @param data AClass object containing test data in `norm.t`
#' @param min_test_features Minimum number of features to test. If not specified, inferred from training model.
#' @param max_test_features Maximum number of features to test. If not specified, inferred from training model.
#' @param prefix Output file prefix (e.g., `paste(project, format(Sys.time(), "%Y-%m-%d_%H%M"), sep = "_")`).
#' @param out_path Output path. If not provided, it will be inferred from `run_info`.
#'
#' @return AClass object with `$test_results_full` and updated `run_info`. If `Group` is present, also outputs confusion matrix results to disk.
#' @examples
#' test.obj <- nano.test(prefix = "demo", training_model_obj = trained_models, data = test.obj)
#' @export
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
  #training_model_list.mat <- data.frame(matrix(unlist(strsplit(names(training_model_list),"_")), ncol=2,byrow = TRUE, dimnames = list(NULL,c("Alg","Model"))), stringsAsFactors = FALSE)
  #split by last _
  training_model_list.mat <- data.frame(matrix(unlist(strsplit(names(training_model_list),"_(?=[^_]+$)", perl=TRUE)), ncol=2,byrow = TRUE, dimnames = list(NULL,c("Alg","Model"))), stringsAsFactors = FALSE)
  if(is.null(min_test_features)&is.null(max_test_features)) {
    min_test_features <- min(as.numeric(training_model_list.mat$Model))
    max_test_features <- max(as.numeric(training_model_list.mat$Model))
    message(paste0("[MSG] Using min and max number of features from training model - min: ",min_test_features," and max: ",max_test_features))
  }

  if(is.null(alg_list)){
    alg_list <- unique(training_model_list.mat$Alg)
    message(paste0("[MSG] Using algorithm list from training model:"))
    message(paste0(alg_list))
  } else {
    alg_list <- unique(training_model_list.mat$Alg)[unique(training_model_list.mat$Alg) %in% alg_list]
    message(paste0("[MSG] Using algorithm:"))
    message(paste0(alg_list))
  }
  library(caret)

  message(paste0("[MSG] Testing ",nrow(test.df), " samples:"))
  message(row.names(test.df))

  # Initialize Report #
  conf_matrix_list <- list()
  testing_results <- data.frame()
  conf_matrix_results.full <- data.frame() # from conf_matrix
  testing_results_summary_groups_score <- data.frame() # testing class and probability


  for (alg in alg_list) {

    for (i in max_test_features:min_test_features) {

      # Select models #
      testing_model <- training_model_list[[eval(paste0(alg,"_",i))]]

      # Predict #
      #message(paste0("Testing... ",alg,"_",i))

      if (("Group" %in% colnames(test.df))){
        # If Group information is present (ie testing known samples) #
        ## message(paste0("[MSG] Have labels. Contruct Confusion Matrix."))

        #omit#
        #if (alg == "svmLinear"){
        #  preProcValues <- preProcess(subset(train.data.training_main.selected_probes, select = -Group), method= c("center", "scale") )
        #}

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
        ##message("[MSG] Group not detected in dataframe. Skip Confusion Matrix...")

        #omit#
        #if (alg == "svmLinear"){
        #  preProcValues <- preProcess(train.data.training_main.selected_probes, method= c("center", "scale") )
        #}

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
  message("[MSG] Printing results. Also check conf_matrix for matrix")

  if(!is.null(prefix)){prefix <- paste0(prefix,"_")}

  write.table(testing_results_summary_groups_score, file=file.path(out_path,paste0(prefix,"testing_summary_full.txt")), sep = "\t", col.names = NA)

  data$test_results_full <- testing_results_summary_groups_score

  data$run_info$test_settings <- c("models" = names(training_model_list), "alg_list" = alg_list,"min_test_features" = min_test_features, "max_test_features" = max_test_features)

  # if you have conf_matrix # (need if statement)
  if (("Group" %in% colnames(test.df))){
    write.table(conf_matrix_results.full, file=file.path(out_path,paste0(prefix,"_conf_matrix_results_full_Test.txt")), sep = "\t", col.names = NA)
    saveRDS(conf_matrix_list, file=file.path(out_path,paste0(prefix,"_Alg_Conf_Matrix_Test.RDS"))) # more details
  }

  return(data)
}


##### get.nano.test.results #####
#' @title Aggregate and summarize test results from AClass predictions
#' @description Aggregates ensemble classification results, computes model agreement, and exports summary tables and top predictions.
#' @param prefix file prefix (optional)
#' @param data A processed AClass object containing `test_results_full`, `prenorm_qc`, and `run_info`.
#' @param print_report print or not. default is FALSE
#' @param out_path output path. When not provided out_path will be extracted from run_info (default)
#' @return Modified AClass object with two new elements: `test_results_agg` (summary of all predictions) and `test_results` (top prediction per sample). If `print_report = TRUE`, tab-delimited files `*_test_summary_aggregate.txt` and `*_test_summary.txt` are written to `out_path`.
#' @examples
#' test.obj <- get.nano.test.results(prefix = "demo", data = test.obj, print_report = TRUE)
#' @export
get.nano.test.results <- function(prefix, data, print_report=FALSE, out_path=NULL) {

  if(is.null(out_path)){
    out_path = file.path(getwd(),data$run_info$run_id[1])
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
    write.table(testing_results_summary_groups_score_i_summary_full, file=file.path(out_path,paste0(prefix,"test_summary_aggregate.txt")), sep = "\t", col.names = NA)
    write.table(testing_results_summary_groups_score_i_max_full, file=file.path(out_path,paste0(prefix,"test_summary.txt")), sep = "\t", col.names = NA)
  }

  data$test_results_agg <- testing_results_summary_groups_score_i_summary_full
  data$test_results <- testing_results_summary_groups_score_i_max_full

  return(data)
}


##### nano.plot #####
#' @title Visualize classification results and QC metrics
#' @description Generates PDF plots summarizing prediction confidence, model agreement, and sample quality control. Includes both detailed and summary report options with support for calibrated or raw probabilities.
#' @param prefix file prefix (optional)
#' @param prob  "Avg_Probability" or "Avg_Cal_Probability"
#' @param thres_avg_prob Probability threshold below which predictions are considered low confidence. Related to "prob". Default to 0.
#' @param thres_geomean geomean threshold cut off below which will be considered fail. Default to NULL which indicates this test is skipped.
#' @param report_type report format options. "Summary" or "Detailed"
#' @param print_report  binary option to print txt summary "test_summary_aggregate.txt" and "test_summary.txt"
#' @param out_path output path. When not provided out_path will be extracted from run_info (default)
#' @examples
#' test.obj <- nano.plot(prefix = "demo", data = test.obj, report_type = "Summary", print_report = TRUE)
#' @export
nano.plot <- function(prefix, data, prob="Avg_Probability", thres_avg_prob=0, thres_geomean = NULL, report_type=c("Summary","Detailed"), print_report=FALSE, out_path=NULL){

  report_type <- match.arg(report_type)

  # Check #
  if(prob == "Avg_Cal_Probability" & is.null(data$test_results_agg$Avg_Cal_Probability)){
    stop("[MSG] Avg_Cal_Probability missing from data. Run nano.calibrate() first to plot results from calibrated probability or proceed with Avg_Probability.")
  }

  if(!is.numeric(thres_avg_prob)){stop("[MSG] thres_avg_prob must be numeric")}

  if(is.null(data$colour_code)){
    stop("[MSG] colour_code missing from data. Run nano.set.colour() first.")
  }

  if(is.null(out_path)){
    out_path = file.path(getwd(),data$run_info$run_id[1])
  }

  ## prepare all data frames ##
  prenorm_qc <- data$prenorm_qc # QC stats
  col_code <- data$colour_code

  test_results_full <- data$test_results_full # full results: all models across probes and results for all classes.
  test_results_agg <- data$test_results_agg # aggregate results: ensemble models across probes and results for all classes. Contain "ALL" which combines all analysis.
  test_results <- data$test_results # final test output

  test_results_agg$Class <- factor(test_results_agg$Class, levels = col_code$Group) # Class required to be factor for plotting

  # prepare color chart of ggpubr #
  color_chart <- col_code$Group_Colour
  names(color_chart) <- col_code$Group

  # library(gridExtra)
  # library(ggpubr)

  test_results$prob <- test_results[,prob]
  test_results_agg$prob <- test_results_agg[,prob]

  if(prob == "Avg_Probability"){
    plot_title <-"Average probability"
  } else if(prob == "Avg_Cal_Probability"){
    plot_title <-"Average calibrated probability"
  }

  ### functions ###

  get_qc_report <- function(test_results,prenorm_qc, sample, prob, thres_geomean, thres_avg_prob){
    result.i <- c(test_results[test_results$Sample == sample,], prenorm_qc[prenorm_qc$Sample == sample,])

    # get sample QC stats #
    if (is.null(thres_geomean) ){
      qc <- "PASS"
      remarks1 <- "QC - SKIPPED"
    } else if (result.i$GeoMean >= thres_geomean){
      qc <- "PASS"
      remarks1 <- "QC - PASS"
    } else if (result.i$GeoMean < thres_geomean) {
      qc <- "FAIL"
      remarks1 <- "QC - FAIL"
    }

    if (result.i[prob] >= thres_avg_prob && qc == "PASS"){
      remarks2 <- "High Confidence"
    } else if (result.i[prob] > thres_avg_prob && qc == "FAIL") {
      remarks2 <- "Caution"
    } else if (result.i[prob] < thres_avg_prob) {
      remarks2 <- "Low Confidence"
    }

    # QC Summary Table #

    # gather results #
    s <- result.i$Sample
    c <- result.i$Class
    cp <- paste0(round(as.numeric(result.i[prob]), 4) * 100, "%")
    mt <- result.i$N_models
    ma <- paste0(round(result.i$Agreement * 100, 4),"%")

    #catch exceptions for non-nanostring data
    if (is.null(result.i$GeoMean)) {
      gm <- NA
    } else {
      gm <- result.i$GeoMean
    }
    if (is.null(result.i$CV)){
      cv <- NA
    } else {
      cv <- result.i$CV
    }

    # combine results #
    t.result <- t(data.frame("Sample" = s, "Class" = c, "Class Prob" = cp, "Prob Type"=prob, "Models Tested" = mt, "Models Agreement" =  ma, "GeoMean" = gm, "CV" = cv, "Remarks1" = remarks1, "Remarks2" = remarks2))
    colnames(t.result) <- "Summary Table"

    return(t.result)
  }

  get_agg_line_p <- function(test_results_agg, sample, color_chart, plot_title){
    # Num_Features != "ALL" to remove combined results
    agg.i <- test_results_agg[test_results_agg$Sample == sample & test_results_agg$Num_Features !="ALL",]
    agg.m <- reshape2::melt(agg.i[c("Num_Features",prob,"Class")], variable_name = "Class", id.vars=c("Num_Features",prob,"Class"))
    agg.m[is.na(agg.m)] <- 0 # for cases with no predictions probability

    ## ggline ##
    agg.p <- ggpubr::ggline(agg.m, x = "Num_Features", y = prob,
                            linetype = "Class",
                            color = "Class",
                            palette = color_chart,

                            #facet.by = "Class",
                            nrow = 1,
                            size = 0.5,
                            main = paste0(sample," - ",plot_title),
                            xlab = "Number of Genes",
                            ylab = prob
    ) +
      scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.2)) +
      #scale_x_continuous(breaks = seq(from = 0, to = 30 ,by = 1)) +
      coord_cartesian(ylim = c(0, 1))
    return(agg.p)
  }

  ## ggdotchart ##
  get_agg_dot_p <- function(test_results_agg, prob,color_chart, sample, col_code){

    agg.ALL <- test_results_agg[test_results_agg$Sample == sample & test_results_agg$Num_Features == "ALL" ,]

    # merge with Group to ensure all groups are present
    Groups <- data.frame(Class = col_code$Group)
    agg.ALL <- merge(Groups, agg.ALL,by = "Class", all.x = TRUE)

    # fill NA's with 0 for cases with no predictions. Col 1 and 2 are Class and Sample (factors) #
    agg.ALL.tmp <- agg.ALL[,-c(1,2)]
    agg.ALL.tmp[is.na(agg.ALL.tmp)] <- 0
    agg.ALL <- cbind(agg.ALL[,c(1,2)],agg.ALL.tmp)

    agg.ALL$Agreement <- as.numeric(round(agg.ALL$Agreement, digits = 4) * 100)

    agg.ALL.dot.p <- ggpubr::ggdotchart(agg.ALL, x= "Class",y = prob,
                                        linetype = "Class",
                                        color = "Class",
                                        palette = color_chart,
                                        dot.size = 5,
                                        nrow = 1,
                                        size = 0.5,
                                        legend="",
                                        xlab = "",
                                        group = "Class"
    ) +
      scale_y_continuous(breaks = seq(from = 0, to = 1 ,by = 0.2)) +
      coord_cartesian(ylim = c(0, 1))+
      #theme(plot.title = element_text(hjust = 0.0)) +
      theme_minimal()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      #theme(axis.text = element_text(size=10))+
      geom_hline(yintercept=seq(0,1,0.2), linetype="dashed", colour = "grey70")
    return(agg.ALL.dot.p)
  }

  get_agg_stacked_bar_p <- function(test_results_agg,sample,color_chart){

    if(!is.factor(test_results_agg$Class)){
      stop("[MSG] Class must be factor")
    }

    agg.ALL <- test_results_agg[test_results_agg$Sample == sample & test_results_agg$Num_Features == "ALL" ,]
    # reverse order to order from bottom up#
    agg.ALL$Class <- factor(agg.ALL$Class, levels = rev(levels(agg.ALL$Class)))
    agg.ALL.bar <- ggplot2::ggplot(data = agg.ALL, aes(x= Sample,y=Agreement, fill= Class, label = paste0(signif(Agreement,digits=3)*100,"%"))) +
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

    return(agg.ALL.bar)
  }

  get_full_horizontal_bar_p <- function(test_results_full,sample,color_chart, col_code){
    test_results_full.i <- test_results_full[test_results_full$Sample == sample,]
    all.horbar <- test_results_full.i

    #sort bars#
    #use original order
    #all.horbar$Order <- row.names(all.horbar)
    #use Group factor from col_code and sort by Probability
    all.horbar$Class <- factor(all.horbar$Class, levels = col_code$Group)
    all.horbar <- all.horbar[with(all.horbar,order(Class,-Probability)),]
    all.horbar$Order <- 1:nrow(all.horbar)

    #plot
    all.horbar.p <- ggpubr::ggbarplot(all.horbar, x="Order",y="Probability",
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
    return(all.horbar.p)
  }

  ## raw line plot by algorithms ##
  get_full_line_p <- function(test_results_full, sample,color_chart,plot_title, col_code){
    raw.i <- test_results_full[test_results_full$Sample == sample,]
    raw.m <- reshape2::melt(raw.i[c("Num_Features","Probability","Class","Alg")], variable_name = "Class", id.vars=c("Num_Features","Probability","Class","Alg"))
    #apply order from col_code
    raw.m$Class <- factor(raw.m$Class, levels = col_code$Group)

    raw.p <- ggpubr::ggline(raw.m, x = "Num_Features", y = "Probability",
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

    return(raw.p)
  }
  ### Output ###

  t.result.list <- list()
  if(!is.null(prefix)){prefix <- paste0(prefix,"_")}

  pdf(file = file.path(out_path,paste0(prefix,"test_result_",report_type,"_plots.pdf")), width = 10.5, height = 8)
  for (SAMPLE in unique(test_results_full$Sample)) {
    message(SAMPLE)

    #### Result stats and plots ####
    t.result <- get_qc_report(test_results=test_results,prenorm_qc=prenorm_qc, sample=SAMPLE, prob = prob, thres_geomean, thres_avg_prob)
    t.result.list[[SAMPLE]] <- t.result

    ## algorithms line plot ##
    agg.line.p <- get_agg_line_p(test_results_agg=test_results_agg, sample=SAMPLE, plot_title=plot_title, color_chart=color_chart)
    ## cleveland plot ##
    agg.dot.p <- get_agg_dot_p(test_results_agg=test_results_agg, sample=SAMPLE, prob=prob, color_chart=color_chart, col_code=col_code )
    ## get_agg_stacked_bar_p ##
    agg.stacked_bar.p <- get_agg_stacked_bar_p(test_results_agg=test_results_agg, sample=SAMPLE, color_chart=color_chart)
    # horiz bar plot #
    full.horizontal_bar.p <- get_full_horizontal_bar_p(test_results_full=test_results_full, sample=SAMPLE, color_chart=color_chart, col_code=col_code)
    # algorithms full line plot #
    full.line.p <- get_full_line_p(test_results_full=test_results_full, sample=SAMPLE, color_chart=color_chart, plot_title=plot_title, col_code=col_code)

    ##### Format Report ####

    if (report_type == "Detailed"){
      grid.arrange(
        top = paste0("Test Result Details - ",SAMPLE),
        tableGrob(t.result),
        #agg.line.p,
        agg.dot.p,
        agg.stacked_bar.p,
        full.line.p,
        ncol = 4,
        #widths = c(1,2)
        widths = c(1,1,1,1),
        bottom=paste0("ATRT Classifier - For Research Only")

      )
    } else if(report_type == "Summary"){
      #plot_grid(tableGrob(t.result), agg.dot.p, agg.stacked_bar.p = c('A', 'B','C'))
      grid.arrange(
        top = paste0("Test Result Summary - ",SAMPLE),
        tableGrob(t.result),
        arrangeGrob(arrangeGrob(agg.stacked_bar.p,full.horizontal_bar.p, ncol=2, widths = c(0.8,1),top=textGrob(paste0("Fig.1 - Models agreement"), gp = gpar(fontface=2)), bottom= textGrob("Probability distribution", gp = gpar(fontsize=11))),
                    arrangeGrob(agg.dot.p, widths = 1, top=textGrob(paste0("Fig.2 - ",plot_title," by class"), gp = gpar(fontface=2))),
                    ncol=1),
        ncol = 2,
        #widths = c(2,1,0.5),
        bottom=paste0("ATRT Classifier - For Research Only")
        #layout_matrix=matrix(c(1,2,1,3), byrow=TRUE,ncol = 2)

      )
    }
  } # SAMPLE for loop

  dev.off()

  if(print_report == TRUE){
    write.table(test_results_agg, file=file.path(out_path, paste0(prefix,"test_summary_aggregate.txt")), sep = "\t", col.names = NA)
    write.table(test_results, file=file.path(out_path,paste0(prefix,"test_summary.txt")), sep = "\t", col.names = NA)
  }

  data$test_summary <- t.result.list
  return(data)
}

##### nano.cal_model #####
#' @title Train class-wise recalibration models for prediction confidence
#' @description
#' Fits regression models to recalibrate classifier scores for each class using true labels and predicted probabilities, inspired by Capper et al. (2018) doi: 10.1038/nature26000. Supports both `glm` and `glmnet` methods.
#' It is recommended to use a held-out validation set or cross-validation for calibration, to avoid data leakage and overfitting.
#'
#' @param cal_labels.df Data frame. Must contain columns: "obs" for ground truth, "Class" for multi-class prediction made from classifier, "Sample" for Names for where the Avg_Probability originated, and "Avg_Probability" for average probability predicted from the classifier.
#' @param method Algorithm to use for score recalibration. Options are `glm` or `glmnet`
#' @return Named list of trained calibration models, one per class. Includes original calibration data under `$Data`.
#' @examples
#' cal_models <- nano.cal_model(cal_labels.df, method = "glm")
#' @export
nano.cal_model <- function(cal_labels.df, method=c("glm","glmnet")) {

  method <- match.arg(method)

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

      #set.seed(849)
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


##### nano.calibrate #####
#' @title Apply recalibration models for prediction confidence to calibrate "Avg_Probability". Refer to nano.cal_model().
#' @description
#' Fits regression models to recalibrate classifier scores for each class using true labels and predicted probabilities, inspired by Capper et al. (2018) doi: 10.1038/nature26000. Supports both `glm` and `glmnet` methods.
#' @param data AClass object with `test_results_agg` and `run_info`.
#' @param Cal_models A list of trained models returned from `nano.cal_model()`.
#' @param print_report Logical, whether to write test summary outputs.
#' @param method Character, either "glm" or "glmnet".
#' @param out_path Character, path to write outputs. If NULL, inferred from `data$run_info`.
#' @param prefix Optional character prefix for output files.
#' @return AClass data object with updated `test_results_agg` and `test_results` after calibration.
#' @examples
#' data <- nano.calibrate(data, Cal_models)
#' @export
nano.calibrate <- function(prefix=NULL, data, Cal_models, print_report = FALSE, method=c("glm","glmnet"), out_path=NULL){

  method <- match.arg(method)

  if("Avg_Cal_Probability" %in% colnames(data$test_results_agg)){
    stop("[MSG] Data set has been calibrated already")
  }

  if(is.null(out_path)){
    out_path = file.path(getwd(),data$run_info$run_id[1])
  }

  message(paste0("[MSG] Calibrate Avg_Probability..."))
  test_results_agg <- data$test_results_agg

  test_results_agg$Avg_Cal_Probability <- NA
  test_results_agg$Avg_Probability[is.na(test_results_agg$Avg_Probability)] <- 0
  test_results_agg <- test_results_agg[test_results_agg$Num_Features != "ALL",] # calculate ALL after calibration

  for (c in unique(test_results_agg$Class)){
    message(c)
    if (method == "glm"){
      message("[MSG] Using glm model")
      #Calibrate using Avg_Probability only
      test_results_agg[test_results_agg$Class == c,"Avg_Cal_Probability"] <- predict(Cal_models[[c]],
                                                                                     test_results_agg[test_results_agg$Class == c,"Avg_Probability",drop=FALSE],
                                                                                     type = "response")

      # Calibrate using Avg_Probability and Agreement
      # test_results_agg[test_results_agg$Class == c,"Avg_Cal_Probability"] <- predict(Cal_models[[c]],
      #                                                                                test_results_agg[test_results_agg$Class == c,c("Avg_Probability","Agreement"),drop=FALSE],
      #                                                                               type = "response")
    } else if(method == "glmnet"){
      message("[MSG] Using glmnet model")
      test_results_agg$Constant <- rep(1,nrow(test_results_agg))
      calibrated_results <- predict(Cal_models[[c]], newdata = test_results_agg[test_results_agg$Class == c,c("Avg_Probability","Constant"),drop=FALSE], type = "prob")
      test_results_agg[test_results_agg$Class == c,"Avg_Cal_Probability"]  <- calibrated_results[, "Class", drop = TRUE]
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
    write.table(test_results_agg, file=file.path(out_path,(paste0(prefix,"test_summary_aggregate.txt"))), sep = "\t", col.names = NA)
    write.table(test_results_agg.ALL.full, file=file.path(out_path,paste0(prefix,"test_summary.txt")), sep = "\t", col.names = NA)
  }
  data$test_results_agg <- test_results_agg
  data$test_results <- test_results_agg.ALL.full
  return(data)

}


### Open Directory Interactively ###
choose_directory = function(caption = 'Select data directory') {

  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption)
  } else {
    tk_choose.dir(caption = caption)
  }
}

##### Plot MDS #####
#' @title General-purpose data MDS plotting.
#' @description
#' Requires nano.set.colour() to be run first to determine colour_code
#' @param prefix Character string used in plot titles.
#' @param data AClass object with normalized expression data.
#' @param plot_type Type of plot to generate. Options: "boxplot", "plot", "ggplot", "ggplot_label", "ggplot_label_batch", "plotly", "plotly_batch"
#' @param data_name Name of the data slot to use. Options: "norm.t", "train.data.main", or "train.data.validate".
#' @return No object is returned; a plot (ggplot2 or plotly) is printed to the active device based on the selected type.
#' @examples
#'  library(limma) #plotMDS
#'  library(plotly)
#'  library(reshape2)
#' nano.MDS(prefix = project_prefix, data=train, plot_type = "ggplot",data_name = "norm.t")
#' @export
nano.MDS <- function(prefix, data, plot_type = c("boxplot","plot","ggplot","ggplot_label","ggplot_label_batch","plotly","plotly_batch"), data_name = c("norm.t","train.data.main","train.data.validate")){
  library(ggrepel)

  #apply and check default
  plot_type <- match.arg(plot_type)
  data_name <- match.arg(data_name)


  data_df <- data[[data_name]]
  data_csv <- data[["run_info"]]$csv

  if (data_name %in% c("train.data.main","train.data.validate")){
    if(is.null(data[[data_name]]$Group)){
      stop("[MSG] Data must have Group labels. Did you run nano.trainsplit()?")
    }
    group <- data_df[,"Group",drop=FALSE]
    data_df <- data.frame(t(subset(data_df, select = -Group)))
  } else if (data_name %in% "norm.t"){
    if(!is.null(data[["test_results"]])){
      message(paste0("[MSG] Applying test_results to plot..."))
      group <- data.frame(matrix(nrow=nrow(data_df), ncol=1, data$test_results$Class))
    } else {
      group <- data.frame(matrix(nrow=nrow(data_df), ncol=1, rep(NA,nrow(data_df))))
    }
    colnames(group) <- "Group"
    row.names(group) <- row.names(data_df)
    data_df <- data.frame(t(data_df))
  }

  # assign colour by colour_code or "black" if not present
  if(is.null(data$colour_code)){
    message(paste0("[MSG] colour code not detected, using default colours"))
    groupcol <- rep("black",nrow(data_df))
    col_code <- data.frame(Group=as.factor(c("Samples")),Group_Colour=c("black"), stringsAsFactors = FALSE)
  } else {
    col_code <- data$colour_code
    groupcol <- col_code[group$Group,"Group_Colour"]
    #groupcol[is.na(groupcol)] <- "nogroup"
    groupcol[is.na(groupcol)] <- "black"
  }

  #batch info#
  batch_details <- data.frame()
  for(i in 1:length(data_csv)){
    batch_details_i <- data.frame(Sample=data_csv[[i]]$details$samples,Batch=names(data_csv[i]))
    if(is.null(nrow(batch_details))){
      batch_details <- batch_details_i
    } else {
      batch_details <- rbind(batch_details,batch_details_i)
    }
  }

  # sample N check #
  if(ncol(data_df) <3){
    message("[MSG] Data must have minimum 3 samples to run nano.MDS().")
    return(invisible(NULL))
  } else {

    library(limma) #plotMDS
    library(plotly)
    library(reshape2)

    ##### Plots #####
    PlotTitle <- prefix

    # obtain MDS matrix #
    pdf(file = NULL) # prevent writing file
    mds <- limma::plotMDS(data_df,pch=19, main=PlotTitle, plot=FALSE)
    dev.off()
    group[is.na(group)] <- "black"
    #mds.anno <- merge(mds$cmdscale.out,group, by="row.names") # limma depreciated since 3.48.0
    mds_2d_matrix <- data.frame(x=mds$x, y=mds$y)
    row.names(mds_2d_matrix) <- row.names(mds$distance.matrix.squared)
    mds.anno <- merge(mds_2d_matrix,group, by="row.names")
    colnames(mds.anno) <- c("Sample","X","Y","Group")

    mds.anno <- merge(mds.anno, batch_details, by="Sample")

    dot_size <- 4

    mds.p <- ggplot(mds.anno, aes(x=X, y=Y, label=Sample, color=Group)) +
      geom_point(size=dot_size) +
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
            panel.background = element_rect(colour = "black", size=1))+
      theme(legend.position="bottom")

    mds.p.plotly <- plotly::ggplotly(mds.p)
    ## box plot ##
    if (plot_type == "boxplot"){
      boxplot(data_df,range=0,ylab="log2 intensity") #intensity runs from 5 to 16 on log2 scale
    }
    ## MDS plot ##
    if (plot_type == "plot"){
      limma::plotMDS(data_df, pch=19, col = groupcol, main=PlotTitle) # point and color
      limma::plotMDS(data_df,labels=colnames(data_df), pch=19, cex=0.5) # labels and color
    }
    ## MDS ggplot ##
    if (plot_type == "ggplot"){
      #mds.p <- mds.p + geom_label(nudge_y=nudge_y_value)
      print(mds.p)
    } else if (plot_type == "ggplot_label"){
      mds.p <- mds.p +geom_text_repel(aes(label = Sample))
      print(mds.p)
    } else if (plot_type == "ggplot_label_batch"){
      mds.p <- mds.p + geom_point(data = mds.anno, aes(color=Batch, size = dot_size)) + scale_color_manual(values = as.factor(c(unique(mds.anno$Batch),"black")))
      print(mds.p)
    } else if (plot_type == "plotly_batch"){
      mds.p <- mds.p + geom_point(data = mds.anno, aes(color=Batch, size = dot_size)) + scale_color_manual(values = as.factor(c(unique(mds.anno$Batch),"black")))
      print(plotly::ggplotly(mds.p)) # doesn't work when save to anther variable

    }


    ## MDS plotly ##
    if (plot_type == "plotly"){
      print(mds.p.plotly)
    }
  }
}

##### nano.set.colour #####
#' @title Set colour code
#' @description
#' Assign colour_code to data. Can assign Group_names or auto search by providing 'NULL' in Group_names which searches for the 'Group' column in the specified dataframe
#' @param data AClass object with a training data slot.
#' @param Group_names Optional character vector of group names. If `NULL`, will use labels in the data.
#' @param data_name Character. Name of the data slot to pull group labels from. Options: "train.data.main" or "train.data.validate".
#' @param remap_to_atrt_consensus Logical. If TRUE, will remap Torchia et al. 2016 subgroup names to ATRT consensus names from Ho et al. 2019.
#' @return Modified AClass object with a new `colour_code` data frame for plotting.
#' @examples
#' # Case 1: Auto-detect groups from train.data.main
#' data <- nano.set.colour(data)
#'
#' # Case 2: Auto-detect groups from train.data.validate
#' data <- nano.set.colour(data, data_name = "train.data.validate")
#'
#' # Case 3: Set group names manually (ATRT subgroups)
#' data <- nano.set.colour(data, Group_names = c("SHH", "TYR", "MYC"))
#'
#' # Case 4: Custom groups not matching preset schemes
#' data <- nano.set.colour(data, Group_names = c("A", "B", "C"))  # default numeric colors applied
#'
#' @export
nano.set.colour <- function(data, Group_names = NULL, data_name = c("train.data.main","train.data.validate"), remap_to_atrt_consensus = FALSE){

  #check and assign default
  data_name <- match.arg(data_name)

  Group_labels <- vector()
  if (is.null(Group_names)) {
    if(data_name == "train.data.main" && is.null(data$train.data.main)){
      stop("[MSG] Data must have Group labels. Did you run nano.trainsplit()?")
    }
    if(data_name == "train.data.validate" & is.null(data$train.data.validate)){
      stop("[MSG] Data must have Group labels. Did you run nano.trainsplit()?")
    }
    data.df <- data[[data_name]]
    if (is.null(data.df[["Group"]])){
      stop("[MSG] Data frame don't have Group labels")
    }else {Group_labels <- unique(data.df$Group)}

  } else {
    Group_labels <- Group_names
  }

  # Remap Torchia 2016 to Ho 2019, only if valid
  group_remap <- c("Group1" = "SHH", "Group2A" = "TYR", "Group2B" = "MYC")
  if (remap_to_atrt_consensus) {
    if (setequal(sort(Group_labels), sort(names(group_remap)))) {
      message("[MSG] Remapping Torchia 2016 subgroups to Ho et al. 2019 labels.")
      if (!is.null(data[[data_name]]$Group)) {
        data[[data_name]]$Group <- as.character(group_remap[as.character(data[[data_name]]$Group)])
      }
      Group_labels <- as.character(group_remap[Group_labels])
    } else {
      warning("[MSG] Remapping skipped. Group labels do not match expected Torchia et al., 2016 format.")
    }
  }

  #detect color schemes from subgroups#
  if (all(Group_labels %in% c("Group1","Group2"))){
    message(paste0("[MSG] Torchia et al., 2015 ATRT Subgroups detected:"))
    message(Group_labels)
    col_code <- data.frame(Group=as.factor(c("Group1","Group2")),Group_Colour=c("red","blue"), stringsAsFactors = FALSE)
    col_code$Group <- factor(col_code$Group, levels = c("Group1","Group2"))
  } else if (all(Group_labels %in% c("Group1","Group2A","Group2B"))){
    message(paste0("[MSG] Torchia et al., 2016 ATRT Subgroups detected:"))
    message(Group_labels)
    col_code <- data.frame(Group=as.factor(c("Group1","Group2A","Group2B")),Group_Colour=c("red","blue","green"), stringsAsFactors = FALSE)
    col_code$Group <- factor(col_code$Group, levels = c("Group1","Group2A","Group2B"))
  } else if (all(Group_labels %in% c("SHH","TYR","MYC"))){
    message(paste0("[MSG] Ho et al., 2019 ATRT Subgroups detected:"))
    message(Group_labels)
    col_code <- data.frame(Group=as.factor(c("SHH","TYR","MYC")),Group_Colour=c("#4074E5","#DD1D06","#23AE2E"), stringsAsFactors = FALSE)
    col_code$Group <- factor(col_code$Group, levels = c("SHH","TYR","MYC"))
  } else {
    message(paste0("[MSG] Using custom labels:"))
    message(Group_labels)
    col_code <- data.frame(Group=as.factor(Group_labels),Group_Colour=as.numeric(Group_labels))
  }
  data$colour_code <- col_code
  return(data)
}


##### nano.eval.test #####
#' @title Evaluate test results across probability thresholds
#' @description This function aggregates and evaluates model predictions from *_test_summary.txt files, generating performance metrics across specified probability thresholds and producing confusion matrices (if `Group` column exist in dataframe) and accuracy summaries. Note: `confusionMatrix()` from the `caret` package requires more than two classes to compute multiclass metrics reliably.
#' @param prefix Prefix string to identify result files.
#' @param use_class Character vector specifying class order for factor alignment.
#' @param Prob_range Numeric vector of thresholds to evaluate (e.g., seq(0,1,0.01)).
#' @param prob Column name in *_test_summary.txt indicating prediction probability.
#' @param training_memberships_path Path to a tab-delimited file with two columns: sample ID (column 1) and ground truth class (column2).
#' @param GeoMean_thres Optional numeric threshold to filter on GeoMean values.
#' @param out_path Directory to write results. Defaults to working directory.
#' @param in_path Directory to read *_test_summary.txt files. Defaults to working directory.
#' @param recursive_read Logical, whether to search in_path recursively.
#' @return A list with two entries: overall accuracy and confusion matrix.
#' @examples
#' # Run evaluation (requires *_test_summary.txt and mapping file)
#' # nano.eval.test(
#' #   prefix = "demo",
#' #   use_class = c("SHH", "TYR", "MYC"),
#' #   training_memberships_path = "sample_to_class.txt"
#' # )
#' @export
nano.eval.test <- function(prefix, use_class=NULL, Prob_range=seq(from=0,to=1,by=0.01), prob = "Avg_Probability", anno_table=NULL, training_memberships_path=NULL, GeoMean_thres=NULL, out_path=NULL, in_path=getwd(), recursive_read=FALSE){

  ### check ###
  if(is.null(out_path)){
    #out_path = paste(getwd(),data$run_info$run_id[1], sep = "/") # don't need data
    out_path = getwd()
  }

  if(!is.null(anno_table)){
    stop("[MSG] anno_table is no longer supported, use training_memberships_path instead")
  }

  if(is.null(use_class)){
    stop("[MSG] use_class is required. Hint: test_obj$colour_code$Group")
  }

  if(is.null(training_memberships_path)){
    stop("[MSG] training_memberships_path is required")
  } else if (!is.null(training_memberships_path)){
    anno <- read.table(training_memberships_path, header = FALSE, sep = "\t")
    if(ncol(anno) != 2){
      stop("[MSG] Expects 2 columns and header with nano_filename in column 1 and and Class in column 2. ")
    }
    colnames(anno) <- c("Sample","Ground_Truth")
  }


  summary_file.df <- summary_file.agg.df <- summary_file.full.df <- data.frame()
  Test_Summary_Overall <- Test_Summary_Stats <- Test_Summary <- conf_matrix_list <-list()

  if(nchar(prefix)>=20){
    message(paste0("[MSG] Trucating output names to 15 chars: ",substr(prefix,0,20)))
    WD <- substr(prefix,0,20)
  } else {
    WD <- prefix
  }
  ## Load Recursively ##
  # test summary
  for (summary_file in list.files(path=in_path, pattern = paste0(prefix,".*_test_summary.txt"), recursive = recursive_read, full.names = TRUE)){
    if(is.null(summary_file)){
      stop("[MSG] Can't find *test_summary.txt. Did you run nano.plot()?")
    }
    summary_file.df.i <- read.table(paste0(summary_file), sep="\t", header = TRUE, row.names = 1,stringsAsFactors = FALSE)
    summary_file.df <- rbind(summary_file.df,summary_file.df.i)
  }
  # test_summary_aggregate
  for (summary_file in list.files(path=in_path, pattern = paste0(prefix,".*_test_summary_aggregate.txt"), recursive = recursive_read, full.names = TRUE)){
    if(is.null(summary_file)){
      stop("[MSG] Can't find *_test_summary_aggregate.txt. Did you run nano.plot()?")
    }
    summary_file.agg.df.i <- read.table(paste0(summary_file), sep="\t", header = TRUE, row.names = 1,stringsAsFactors = FALSE)
    summary_file.agg.df <- rbind(summary_file.agg.df,summary_file.agg.df.i)
  }
  #testing_summary_full
  for (summary_file in list.files(path=in_path, pattern = paste0(prefix,".*_testing_summary_full.txt"), recursive = recursive_read, full.names = TRUE)){
    if(is.null(summary_file)){
      stop("[MSG] Can't find *_testing_summary_full.txt. Did you run nano.plot()?")
    }
    summary_file.full.df.i <- read.table(paste0(summary_file), sep="\t", header = TRUE, row.names = 1,stringsAsFactors = FALSE)
    summary_file.full.df <- rbind(summary_file.full.df,summary_file.full.df.i)
  }

  ## Annotate ##

  # annotate all ensemble models
  summary_file.df.anno.full <- merge(summary_file.df, anno, by.x="Sample",by.y="Sample")
  #row.names(summary_file.df.anno.full) <- summary_file.df.anno.full[,1] # not necessary and commenting out allow duplicate row names

  # change factor order #

  summary_file.df.anno.full$Class <- factor(summary_file.df.anno.full$Class, levels = use_class)

  summary_file.df.anno.full$Ground_Truth <- factor(summary_file.df.anno.full$Ground_Truth, levels = use_class)

  summary_file.df.anno.full$Matching_class <- ifelse(summary_file.df.anno.full$Ground_Truth == summary_file.df.anno.full$Class, 1, 0)

  summary_file.agg.df.anno <- merge(summary_file.agg.df, anno, by.x="Sample", by.y="Sample")
  summary_file.full.df.anno <- merge(summary_file.full.df, anno, by.x="Sample", by.y="Sample")

  ## Remove Failed Samples ##

  # skip for now
  summary_file.df.anno <- summary_file.df.anno.full

  if(!is.null(GeoMean_thres)){
    summary_file.df.anno <- summary_file.df.anno.full[summary_file.df.anno.full$GeoMean >= GeoMean_thres,]
  }

  # frozen / ffpe filter #
  #select_matrials <- c("Frozen","extracted_RNA")
  # select_matrials <- "FFPE"
  # select_matrials <- "Frozen"
  # select_matrials <- "extracted_RNA"
  #summary_file.df.anno <-summary_file.df.anno[summary_file.df.anno$RNA.Material.Used.FFPE.Frozen %in%select_matrials, ]

  N_Total <- nrow(summary_file.df.anno.full)
  N_QC_Pass <- nrow(summary_file.df.anno)

  ##### Calcuate TP/TF #####

  # Retrieving stats from confusion matrix generated from caret::confusionMatrix()
  confusionMatrix_stats <- function(conf_mat){
    cm <- conf_mat$table
    all_grps <- colnames(cm)

    conf_mat.stats <- data.frame()
    for (grp in all_grps){
      N_test <- sum(cm[,grp])
      TP <- cm[grp,grp]
      FN <- N_test - TP
      TN <- sum(cm[-grep(grp,colnames(cm)),-grep(grp,row.names(cm))])
      FP <- sum(cm[grp,])-TP
      conf_mat.stats.i <- data.frame(grp=grp,N_test,TP,FN,FP,TN)

      if(nrow(conf_mat.stats)==0){
        conf_mat.stats <- conf_mat.stats.i
      } else {
        conf_mat.stats <- rbind(conf_mat.stats,conf_mat.stats.i)
      }
    }
    return(conf_mat.stats)
  }

  Accuracy_Table <- Subgroup_Accuracy_Table <- data.frame()
  for (Prob in Prob_range){
    summary_file.df.anno.i <- summary_file.df.anno[summary_file.df.anno[,prob] >= Prob,]
    summary_file.df.anno.filtered.i <- summary_file.df.anno[summary_file.df.anno[,prob] < Prob,]

    N_Pass_Prob <- nrow(summary_file.df.anno.i)
    N_Filtered_Prob <- nrow(summary_file.df.anno.filtered.i)

    confmat.i <- caret::confusionMatrix(summary_file.df.anno.i$Class, summary_file.df.anno.i$Ground_Truth, positive = NULL)
    confmat.filtered.i <- caret::confusionMatrix(summary_file.df.anno.filtered.i$Class, summary_file.df.anno.filtered.i$Ground_Truth, positive = NULL)

    # stats <- confusionMatrix_stats(conf_mat = confmat.i)
    # stats_filtered <- confusionMatrix_stats(conf_mat = confmat.filtered.i)

    #overall stats taking into account of threshold#

    # TP <- stats[,"TP"]
    # FN <- stats[,"FN"] + stats_filtered[,c("TP","FN")]
    # TN <- stats[,"TN"] + stats_filtered[,c("TN","FP")]
    # FP <- stats[,"FP"]

    conf_matrix_list[[as.character(Prob)]] <- confmat.i

    # overall level #
    Accuracy_Table.i <- data.frame(Probability = Prob, N=N_Pass_Prob, N_Thres_Filtered=N_Filtered_Prob, t(data.frame(confmat.i$overall)))
    #Accuracy_Table.i <- data.frame(Probability = Prob, N=N_Pass_Prob, N_Thres_Filtered=N_Filtered_Prob, t(data.frame(confmat.i$overall)), TP=TP,TN=TN,FP=FP,FN=FN,FPR=FP/(FP+TN),TRP=TP/(TP+FN), Thres_Accuracy=(TP+TN)/(TP+TN+FP+FN), Thres_Specificity=TN/(TN+FP), Thres_Precision=TP/(TP+FP), Thres_Sensitivity=TP/(TP+FN), Youden_Index=(TN/(TN+FP))+(TP/(TP+FN))-1)
    Accuracy_Table <- rbind(Accuracy_Table,Accuracy_Table.i)

    # subgroup level #
    Subgroup_Accuracy_Table.i <- data.frame(Probability = Prob, Class = gsub("Class: ","",row.names(confmat.i$byClass)), confmat.i$byClass)
    GRP_Count <- data.frame()

    # calculate TP FN TN FN using group vs all
    for(GRP in sort(unique(Subgroup_Accuracy_Table.i$Class))){
      Summary_file.df.anno.i.GRP.True <- summary_file.df.anno.i[summary_file.df.anno.i$Ground_Truth == GRP,,drop=FALSE]
      N_Pass_Prob.GRP.True <- nrow(Summary_file.df.anno.i.GRP.True)
      summary_file.df.anno.filtered.i.GRP.True <- summary_file.df.anno.filtered.i[summary_file.df.anno.filtered.i$Ground_Truth == GRP,,drop=FALSE]
      N_Filtered_Prob.GRP.True <- nrow(summary_file.df.anno.filtered.i.GRP.True)

      Summary_file.df.anno.i.GRP.False <- summary_file.df.anno.i[summary_file.df.anno.i$Ground_Truth != GRP,,drop=FALSE]
      N_Pass_Prob.GRP.False <- nrow(Summary_file.df.anno.i.GRP.False)
      summary_file.df.anno.filtered.i.GRP.False <- summary_file.df.anno.filtered.i[summary_file.df.anno.filtered.i$Ground_Truth != GRP,,drop=FALSE]
      N_Filtered_Prob.GRP.False <- nrow(summary_file.df.anno.filtered.i.GRP.False)


      TP <- nrow(Summary_file.df.anno.i.GRP.True[Summary_file.df.anno.i.GRP.True$Class == GRP,])
      #FN <- nrow(Summary_file.df.anno.i.GRP.True[Summary_file.df.anno.i.GRP.True$Class != GRP,]) + nrow(summary_file.df.anno.filtered.i.GRP.True[summary_file.df.anno.filtered.i.GRP.True$Class == GRP,]) + nrow(summary_file.df.anno.filtered.i.GRP.False[summary_file.df.anno.filtered.i.GRP.False$Class == GRP,])
      FN <- nrow(Summary_file.df.anno.i.GRP.True[Summary_file.df.anno.i.GRP.True$Class != GRP,]) + nrow(summary_file.df.anno.filtered.i.GRP.True) # fixed

      #TN <- nrow(Summary_file.df.anno.i.GRP.False[Summary_file.df.anno.i.GRP.False$Class != GRP,]) + nrow(summary_file.df.anno.filtered.i.GRP.False[summary_file.df.anno.filtered.i.GRP.False$Class != GRP,]) + nrow(summary_file.df.anno.filtered.i.GRP.True[summary_file.df.anno.filtered.i.GRP.True$Class != GRP,])
      TN <- nrow(Summary_file.df.anno.i.GRP.False[Summary_file.df.anno.i.GRP.False$Class != GRP,]) + nrow(summary_file.df.anno.filtered.i.GRP.False)  #fixed
      FP <- nrow(Summary_file.df.anno.i.GRP.False[Summary_file.df.anno.i.GRP.False$Class == GRP,])

      GRP_Count.i <- data.frame(Class=GRP,N_Class=N_Pass_Prob.GRP.True, N_Other_Class=N_Pass_Prob.GRP.False, N_Thres_Filtered =nrow(summary_file.df.anno.filtered.i), TP=TP,TN=TN,FP=FP,FN=FN,FPR=FP/(FP+TN),TRP=TP/(TP+FN), Thres_Accuracy=(TP+TN)/(TP+TN+FP+FN), Thres_Specificity=TN/(TN+FP), Thres_Precision=TP/(TP+FP), Thres_Sensitivity=TP/(TP+FN), Youden_Index=(TN/(TN+FP))+(TP/(TP+FN))-1)
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

  Test_Summary_Overall[["overall_accuracy"]] <- sum(Test_Summary[[1]]$Matching_class) / length(Test_Summary[[1]]$Matching_class)
  Test_Summary_Overall[["confusion_matrix"]] <- conf_matrix_list[[1]]$table

  ### Output ###

  if(!is.null(prefix)){prefix <- paste0(prefix,"_")}

  # export as excel #
  #Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")
  openxlsx::write.xlsx(Test_Summary_Stats, file = file.path(out_path,paste0(prefix,"Test_Summary_Stats.xlsx")), overwrite = TRUE)
  openxlsx::write.xlsx(Test_Summary, file = file.path(out_path,paste0(prefix,"Test_Summary.xlsx")), overwrite = TRUE)
  saveRDS(conf_matrix_list, file = file.path(out_path,paste0(prefix,"conf_matrix_list.RDS")))
  saveRDS(Test_Summary_Overall, file = file.path(out_path,paste0(prefix,"conf_matrix_overall.RDS")))

  return(Test_Summary_Overall)
}

##### batch.nano.eval.test #####
#' @title Batch evaluation of classification results
#' @description
#' Batch processing `nano.eval.test()` and create an overall summary for all runs within dir
#' @param prefix Prefix string to identify result files.
#' @param use_class custom class output order. Required field.
#' @param Prob_range vector of probability intervals to be used in analysis. Default 0 to 1 by 0.01.
#' @param prob column name for probability present in *_test_summary.txt file. Default "Avg_Probability".
#' @param training_memberships_path no header, expects nano_filename in column 1 and and Class in column 2.
#' @param GeoMean_thres Housekeeping gene geometric mean threshold to be considered in analysis. Default NULL for no filtering and is same as using 0.
#' @param run_dir_path  expects multiple test results folders stored within this path.
#' @return No object returned. Evaluation results are written to Excel and RDS files in each subfolder and as a combined summary.
#' @export
batch.nano.eval.test <- function(prefix, use_class=NULL, Prob_range=seq(from=0,to=1,by=0.01), prob="Avg_Probability", training_memberships_path, GeoMean_thres=NULL, run_dir_path){


  for (path_i in list.dirs(run_dir_path,full.names = TRUE)){
    if(path_i == run_dir_path){next}
    in_path_i <- path_i
    out_path_i <- path_i
    nano.eval.test(prefix=prefix, use_class=use_class, Prob_range=Prob_range, prob = prob, anno_table=NULL, training_memberships_path=training_memberships_path, GeoMean_thres=GeoMean_thres, out_path=out_path_i, in_path=in_path_i)
  }

  in_path <- run_dir_path
  out_path <- run_dir_path
  nano.eval.test(prefix=prefix, use_class=use_class, Prob_range=Prob_range, prob = prob, anno_table=NULL, training_memberships_path=training_memberships_path, GeoMean_thres=GeoMean_thres, out_path=out_path, in_path =in_path, recursive_read = TRUE)

}


#### nano.MDS.train.test ####
#' @title Plot both training (train.data.main) and testing (norm.t) data by merging training and test results.
#' @description
#' Requires nano.set.colour() to be run first to determine colour_code. Requires nano.test() to obtain testing results.
#' Memberships and gene_list are optional.
#'
#' @param prefix Character string used in plot titles.
#' @param train.data AClass object containing normalized training data (must include train.data.main).
#' @param test.data AClass object containing normalized test data (must include norm.t and test_results).
#' @param colour_code Data frame specifying group-to-colour mappings. Typically set by nano.set.colour().
#' @param plot_type Type of plot to generate. Options: "plot", "ggplot", "ggplot_label", "ggplot_label_test", "plotly".
#' @param train_ellipse Logical. If TRUE, adds confidence ellipse around training samples. Default is FALSE.
#' @param memberships Optional data frame of sample names (as rownames) and subgroup labels. Overrides default assignments from training and test data.
#' @param gene_list Optional character vector of genes to include in the plot. Default is NULL (use all).
#' @param omit_sample Optional character vector of sample names to omit from the plot. Default is NULL.
#' @param prob Optional column name in test_results to use for coloring test samples by prediction score. Default is NULL.
#'
#' @return A plot is printed to the active device.
#' @examples
#' nano.MDS.train.test(prefix = "demo", train.data = train, test.data = test, colour_code = train$colour_code)
#' @export
nano.MDS.train.test <- function(prefix, train.data , test.data , colour_code, plot_type = c("plot","ggplot","ggplot_label", "ggplot_label_test", "plotly"), train_ellipse=FALSE, memberships=NULL, gene_list=NULL, omit_sample=NULL, prob=NULL){

  plot_type <- match.arg(plot_type)

  data_train <- train.data$train.data.main
  data_test <- test.data$norm.t
  data_test_results <- test.data$test_results

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
    message(paste0("[MSG] Using default gene list m=",length(gene_list)))
  } else {
    message(paste0("[MSG] Using gene list m=",length(gene_list)))
  }

  ##### Plots #####
  message("[MSG] First 50 genes...")
  message(head(paste0(gene_list),n=50))
  data.df <- data.df[,c(gene_list,"Group","Type")]

  col_code <- colour_code
  groupcol <- col_code[col_code$Group %in% data.df$Group,"Group_Colour"]
  group <- data.df[,c("Group","Type")]
  group$Group <- factor(group$Group, levels = colour_code$Group)

  library(limma) #plotMDS
  library(plotly)
  library(reshape2)

  PlotTitle <- prefix
  if(!is.null(omit_sample)){
    message(paste0("[MSG] Omitting samples from plot:"))
    message(omit_sample)
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

  #add test results #
  if(!is.null(prob)){
    data_test_results$pred_score <- data_test_results[,prob]
    mds.anno <- merge(mds.anno, data_test_results, by="Sample", all.x = TRUE)
  }

  # base plot #
  mds.p <- ggplot(mds.anno, aes(x=X,y=Y, label=Sample))

  # add colors #
  if(is.null(prob)){

    mds.p <- mds.p + geom_point(aes(color=Group,shape=Type), size=4)  + scale_color_manual(values = as.character(col_code$Group_Colour))
  } else if(!is.null(prob)){

    #qn = quantile(mds.anno$prob, c(0.01, 0.99), na.rm = TRUE)
    #qn01 <- rescale(c(qn, range(mds.anno$pred_score)))
    #fill.colors <- colorRampPalette(c("darkblue", "white", "darkred"))(20)
    #mds.p <-
    #mds.p + geom_point(aes(color=pred_score), size=3) +
    #scale_colour_gradientn(colours = fill.colors, breaks=seq(0,1,0.1), values = c(0,seq(qn01[1], qn01[2], length.out = 18),1), na.value = "whitesmoke", limits=c(0,1))

    mds.p<-  mds.p + geom_point(aes(color=pred_score,shape=Type), size=4) +
      scale_colour_gradientn(colours = c("red", "yellow", "darkgreen"), breaks=seq(0,1,0.1), values =c(0,0.7,1), na.value = "grey", limits=c(0,1))
    #scale_colour_gradient2(low="red",mid="yellow",high="green",  midpoint = 0.7, breaks=c(0,0.7,1))
    #scale_colour_gradient(low="red",high="green", midpoint = 0.7)
  }

  # add image format #
  mds.p <-mds.p +
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
          panel.background = element_rect(colour = "black", size=1)) +
    theme(legend.key.height = unit(0.75,"inches"))

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
    print(mds.p)
  } else if (plot_type == "ggplot_label_test"){
    mds.p <- mds.p +geom_text_repel(data=mds.p$data[mds.p$data$Type == "Test",], aes(label = Sample), show.legend = FALSE)
    print(mds.p)
  }

}


###### nano.feat.select ######
#' @title Feature selection using Boruta
#' @description
#' Performs wrapper-based feature selection using the Boruta algorithm on the training data.
#' Assumes the input object contains a `train.data.main` data frame with a `Group` column for classification.
#'
#' @param data An AClass object with a `train.data.main` slot. Must include a `Group` column for supervised classification.
#'
#' @return A list with selected important genes, the Boruta object after tentative fix, and selection statistics.
#' @export
nano.feat.select <- function(data){
  library(Boruta)

  if(!"train.data.main" %in% names(data)){
    stop("[MSG] Input file not supported. Check to ensure it contains train.data.main.")
  }

  boruta_data <- data$train.data.main
  boruta_data$Group <- as.factor(boruta_data$Group)
  Boruta.obj <- list()

  # maxRuns   increase runs to resolve tentative features
  # doTrace   get report of the progress

  # features selection
  boruta.data <- Boruta(Group~., data = boruta_data, doTrace = 2)
  print(boruta.data)

  # take a call on tentative features
  #  simple comparison of the median feature Z-score with the median Z-score of the most important shadow feature
  boruta.data.fix <- TentativeRoughFix(boruta.data)
  print(boruta.data.fix)

  # plots
  plot(boruta.data.fix, xlab = "", xaxt = "n")
  lz<-lapply(1:ncol(boruta.data.fix$ImpHistory),function(i)
    boruta.data.fix$ImpHistory[is.finite(boruta.data.fix$ImpHistory[,i]),i])
  names(lz) <- colnames(boruta.data.fix$ImpHistory)
  Labels <- sort(sapply(lz,median))
  axis(side = 1,las=2,labels = names(Labels),
       at = 1:ncol(boruta.data.fix$ImpHistory), cex.axis = 0.7)

  ### get list of important attributes ###
  important_genes <- getSelectedAttributes(boruta.data.fix, withTentative = F)

  ### Extract attribute statistics ###
  boruta.data.fix_df <- attStats(boruta.data.fix)
  print(boruta.data.fix_df)

  Boruta.obj[["Important_Genes"]] <- important_genes
  Boruta.obj[["Boruta_obj_rough_fix"]] <- boruta.data.fix
  Boruta.obj[["Stats"]] <- boruta.data.fix_df
  return(Boruta.obj)
}


##### nano.extract #####
#' @title Extract selected samples from an AClass object
#' @description
#' Subsets an AClass object to include only the specified samples.
#' Also updates the `run_info` metadata to reflect the new sample count.
#'
#' @param data AClass object containing slots: `raw`, `prenorm_qc`, `norm`, `norm.t`, and `run_info`.
#' @param keep_samples_path Path to a CSV file (no header) with two columns: sample name (column 1) and subgroup label (column 2). Default is NULL.
#'
#' @return A subsetted AClass-compatible list containing the selected samples.
#' @export
nano.extract <- function(data, keep_samples_path = NULL) {
  if (is.null(keep_samples_path)) {
    stop("[MSG] keep_samples_path must be provided. CSV file expects no header, with sample names in column 1 and subgroup labels in column 2.")
  }

  keep_samples <- read.table(keep_samples_path, header = FALSE, sep = ",")[, 1]

  data.obj <- list()
  data.obj$raw <- data$raw[, c(colnames(data$raw)[1:3], keep_samples)]
  data.obj$prenorm_qc <- data$prenorm_qc[keep_samples, , drop = FALSE]
  data.obj$norm <- data$norm[, keep_samples, drop = FALSE]
  data.obj$norm.t <- data$norm.t[keep_samples, , drop = FALSE]
  data.obj$run_info <- data$run_info
  data.obj$run_info$samples_loaded <- paste0(nrow(data.obj$norm.t), " samples loaded.")

  message("[MSG] ", nrow(data.obj$norm.t), " samples extracted from ", nrow(data$norm.t), " samples.")
  return(data.obj)
}

##### convert2test #####
#' @title Convert validation set to test input format
#' @description
#' Converts the validation portion of an AClass training object into a format compatible with test-mode functions.
#' It removes the `Group` column from `train.data.validate` if present, and carries over `prenorm_qc` and `run_info`.
#'
#' @param data An AClass object containing `train.data.validate`, `prenorm_qc`, and `run_info` slots.
#' @return A list with `norm.t`, `prenorm_qc`, and `run_info` slots, suitable for use with testing functions.
#' @export
convert2test <- function(data){

  if (!"train.data.validate" %in% names(data)) {
    stop("[MSG] train.data.validate not found in input data.")
  }

  df <- data$train.data.validate
  test <- list()
  # Handle labeled and unlabeled data
  if ("Group" %in% colnames(df)) {
    test$norm.t <- subset(df, select =-Group)
  } else {
    message("[MSG] Group column not found in train.data.validate. Proceeding with unlabeled test data.")
    test$norm.t <- df
  }
  test$prenorm_qc <- data$prenorm_qc
  test$run_info <- data$run_info
  return(test)
}

##### df2nano #####
#' @title Convert annotated data frame to AClass object format
#' @description
#' Converts a data frame into an AClass-style list object with support for
#' training/validation splits and group color annotations. This function is useful for importing custom transcriptomic data into the AClass workflow. By default, the input data frame is assigned to `train.data.main`.
#' @param df A data frame with genes/features as columns and samples as rows. If `Group` column is present it will be passed on to `norm.t`, which is required for `nano.trainsplit()`
#' @param colour_code Colour code data frame. Must contain `Group` and `Group_Colour` columns.
#' @param add_to Either a character ("train.data.main" or "train.data.validate") or a data frame with two columns.
#' In the data frame version, the first column should contain sample names (matching rownames in `df`), and the second column should contain split labels: "train" or "validate".
#' Samples labeled "train" will be assigned to `$train.data.main`, and those labeled "validate" to `$train.data.validate`. Default is "train.data.main".
#' @return A list representing a formatted transcriptomic object compatible with AClass tools.
#' @export
df2nano <- function(df, colour_code=NULL, add_to=c("train.data.main","train.data.validate")){

  # check
  if(!(is.data.frame(add_to)) & !(is.vector(add_to) & length(add_to) ==1 )){
    stop("[MSG] add_to must be either train.data.main or train.data.validate, or data frame with 2 columns e.g. sample_name1 test/validate")
  }

  if( is.vector(add_to) & length(add_to) ==1 ){
    if(add_to != "train.data.main" & add_to !="train.data.validate"){
      stop("[MSG] add_to must be either train.data.main or train.data.validate, or data frame with 2 columns e.g. sample_name1 test/validate")
    }
  }

  # check if there is Group column
  has_group <- "Group" %in% colnames(df)

  # If assigning to train.data.main, Group column is required
  if (!has_group && identical(add_to, "train.data.main")) {
    stop("[MSG] Group column is required when assigning to train.data.main.")
  }

  #convert dataframe to nano.obj #
  nano.obj <- list()
  nano.obj$norm.t <- if (has_group) subset(df, select = -Group) else df  #adaptive assignment depending on presence of Group

  nano.obj$run_info <- list()

  details <- list()
  n_features <- if (has_group) ncol(df) - 1 else ncol(df)
  details$found <- paste0(n_features, " features ", nrow(df), " samples.")  #adaptive assignment depending on presence of Group

  details$samples <- row.names(df)

  nano.obj$run_info$csv <- list()
  nano.obj$run_info$csv[["dummy_data"]] <- list()
  nano.obj$run_info$csv[["dummy_data"]][["details"]] <- details

  if(!is.null(colour_code)){
    #add colour_code
    nano.obj$colour_code <- colour_code
  }
  #add train.data.main

  #add to main or validate based on keywords
  if(is.data.frame(add_to)  ){
    stopifnot(ncol(add_to) ==2)
    nano.obj[["train.data.main"]]  <- df[row.names(df)%in% add_to[add_to[,2] == "train",1],]
    nano.obj[["train.data.validate"]]  <- df[row.names(df)%in% add_to[add_to[,2] == "validate",1],]
  } else if (is.vector(add_to) & length(add_to) ==1){
    nano.obj[[add_to]] <- df
  } else {
    stop("[MSG] add_to must be either train.data.main or train.data.validate")
  }

  return(nano.obj)
}

##### check_and_prepare_test_data #####
#' @title Check if test data object uses the same genes as pre-trained model
#' @description
#' Ensures that the test data (`norm.t`) contains all genes required by the pre-trained model.
#' Optionally renames gene symbols to match training model feature names (e.g., HGNC-approved to aliases).
#' A default mapping (e.g., CTSV → CTSL2, MIR9-1HG → C1orf61) is applied unless a custom `rename_map` is provided.
#'
#' @param test_data An AClass-style object containing `norm.t`.
#' @param model A pre-trained AClass model object with `training_model_list` and preprocessing info.
#' @param rename_map Optional named list of gene symbol remapping (e.g., list("CTSV" = "CTSL2")).
#' @return Modified test_data object with matched and reordered features.
#' @export
check_and_prepare_test_data <- function(test_data, model, rename_map = NULL) {
  # Genes expected by the trained model
  expected_features <- names(model$training_model_list[[1]]$preProcess$mean)
  test_features <- colnames(test_data$norm.t)

  # Use default rename_map if none supplied
  if (is.null(rename_map)) {
    rename_map <- list("CTSV" = "CTSL2", "MIR9-1HG" = "C1orf61")
  }

  # Apply renaming to test_data$norm.t
  for (old in names(rename_map)) {
    new <- rename_map[[old]]
    if (old %in% test_features && !(new %in% test_features)) {
      colnames(test_data$norm.t)[colnames(test_data$norm.t) == old] <- new
    }
  }

  # Re-collect feature names after renaming
  test_features <- colnames(test_data$norm.t)

  # Check for missing genes
  missing <- setdiff(expected_features, test_features)
  if (length(missing) > 0) {
    stop("[MSG] Missing genes in test data after renaming: ", paste(missing, collapse = ", "))
  }

  # Reorder columns in test data to match training model
  test_data$norm.t <- test_data$norm.t[, expected_features]

  return(test_data)
}


