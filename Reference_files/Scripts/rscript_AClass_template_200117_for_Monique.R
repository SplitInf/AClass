### ATRT Classifier Demo ###

# load ATRT Classifier package (AClass.R from Reference_files) #
source(file.choose()) # load package AClass_200117.R

# enter paths
work_path <- tcltk::tk_choose.dir() # choose "Projects" folder
raw_path <- tcltk::tk_choose.dir() # choose "Data" folder
ref_path <- tcltk::tk_choose.dir() # choose "Reference_files" folder

training_memberships_link <- file.choose() # labels used for training *optional* only required when needing to retrain models.

keep_file_path <- file.choose() # *optional* files to include from process.raw() and classify.data() If not used all samples in raw_path will be loaded. Works in opposite way as omit_file_path
omit_file_path <- file.choose() # *optional*
probes_rank_path <- file.choose() # *optional* probe ranking list used for training. By default uses "probes_list.txt" 

# name of project #
project_name <- "demo"

###### Code ######
# {1} Loading all required libraries
# initialize.prj() # install and load all required libraries

# Training. If you have trained models already you can skip all step 2
# {2} load training data
# By default, samples_omit.txt or samples_keep.txt can be supplied to control which samples to exclude or include respectively.
#train <- process.raw(work_path, project = project_name, raw_path = raw_path, ref_path = ref_path, keep_file_path = keep_file_path)

# {3.1} nano.trainsplit()
# splits up data into training and the remaining for testing
# N.train.per is the percentage for training 1=100% training, 0.6=60% training ..etc
#train <- nano.trainsplit(data = train, training_memberships_link, N.train.per=1)

# {3.2} nano.train() & nano.train.report()
# alg_list controls which algrithm to use. By default uses alg_list = c("rf","glmnet","pam", "nb", "knn").
# default training settings: 
# min_test_features=20, 
# max_test_features=30, 
# alg_list = c("rf","glmnet","pam", "nb", "knn"),  
# c.method = "repeatedcv", c.repeats = 5, c.number = 10
#models <- nano.train(prefix=project_name, data = train , alg_list = c("rf","pam"), min_test_features=28, max_test_features=30, c.repeats = 2)

# {4} load trained models for testing. 
# If there is no need for training the process begins here load trained model RDS *_Training_Models_List.RDS
models <- readRDS(file.choose())

# {5} load testing data
# sample loading process for training and testing is the same (refer to {1})


#test <- process.raw(work_path = work_path, raw_path = raw_path, ref_path = ref_path, project = project_name, SampleContent = "none")  # SampleContent default to "housekeeping.geo.mean" when not included. Options are none, housekeeping.sum, housekeeping.geo.mean, total.sum, low.cv.geo.mean, top.mean and top.geo.mean

test <- process.raw(work_path = work_path, raw_path = raw_path, ref_path = ref_path, project = project_name)

# {6} classify.data()
# use pre-trained models for classification and generates report output
# geomean value of housekeeping gene threshold of 100 is used by default
test <- classify.data(data = test, prefix = project_name, training_model_obj = models, omit_file_path = omit_file_path)



# Other tools # 
#1) Visualization with MDS plot

# static plot # (minimum 3 samples)
nano.MDS(prefix = project_name, data = test, plot_type = "ggplot",data_name = "norm.t")

# interactive plot # (minimum 3 samples)
nano.MDS(prefix = project_name, data = test, plot_type = "plotly",data_name = "norm.t")

#2) Visualize data in reference to training data
# can use training data or using the training models if training data is not available.
nano.MDS.train.test(prefix= project_name, train.data = train, test.data=test , colour_code=test$colour_code, plot_type = "ggplot")
nano.MDS.train.test(prefix= project_name, train.data = models, test.data=test , colour_code=test$colour_code, plot_type = "plotly")

#3) Get training information from model ##
get.training.stats(training_model_obj = models) # obtain model information

#4) Generate training report.
# command to pull out information about the models.
nano.train.report(prefix=project_name, training_model_obj=models, feature_min=28, feature_max=30)
