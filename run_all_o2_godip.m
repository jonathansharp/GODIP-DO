%% Run all scripts to make GOBAI-O2 with WOD profile data for DOMIP-1
t_whole_script=tic; % time entire script

%% Set configuration parameters
start_year = 1965;
end_year = 2025;
% system-specific worker configuration
numWorkers_train = 24;
numWorkers_predict = 24;
% float snapshot configuration
snap_download = 1;
snap_date = 202601;
file_date = datestr(datenum(floor(snap_date/1e2),...
    mod(snap_date,1e2),1),'mmm-yyyy');
glodap_year = 2023;
data_modes = {'D' 'A'};
float_file_ext = '_D_A';
% cluster configuration
num_clusters = 15;
clust_vars = {'Temperature' 'Salinity' 'depth'};
glodap_only = false; % EDIT THIS TO 'true' TO TEST WITH GLODAP DATA ONLY
thresh = 0.05;
num_folds = 5;
% algorithm training configuration
variables = ... % variables for algorithms
    {'lat' 'lon_cos_1' 'lon_cos_2' 'depth' 'sigma' ...
    'Temperature' 'Salinity' 'day_sin' 'day_cos' 'year'};
% % random forest regression configuration
% numtrees = 100;
% minLeafSize = 10;
% shallow neural network configuration
train_ratio = 0.8;
val_ratio = 0.1;
test_ratio = 0.1;
% % gradient boosting configuration
% numstumps = 500;
% numbins = 50;
% data and parameter configuration
data_per = 1; % set data reduction
data_per_osse = 0.01; % set data reduction for osse
param = 'o2';
param_props = param_config(param);
% base grid
base_grid = 'EN4';
fpaths = path_config(base_grid,param);
% osse parameters
model_types = {'GFDL-ESM4' 'CanESM5' 'IPSL-CM6A-LR' 'ACCESS-ESM1-5' 'MPI-ESM1-2-LR'};
model_folders = {'GFDL-ESM4' 'CanESM5' 'IPSL-CM6A-LR' 'ACCESS-ESM1-5' 'MPI-ESM1-2-LR'};
realizations = {'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1'};
grid_labels = {'gr' 'gn' 'gn' 'gn' 'gn'};
grid_types = {'regridded' 'native_grid' 'native_grid' 'native_grid' 'native_grid'};

%% Import data
% process_WOD_profile_data(start_year,end_year);
% % plot data histogram by year
% plot_data_hist(start_year,end_year);

%% create time-varying clusters and assign data points to them
% % form clusters
% gmm_clustering(param_props,fpaths,base_grid,start_year,...
%     end_year,snap_date,float_file_ext,clust_vars,num_clusters,...
%     numWorkers_predict);
% % plot cluster animations
% % plot_cluster_animation(param_props,fpaths,base_grid,num_clusters,...
% %     start_year,snap_date,numWorkers_train);
% % cluster data
% assign_data_to_clusters(param_props,base_grid,snap_date,...
%     float_file_ext,clust_vars,num_clusters,start_year,end_year);
% % plot clustered data points
% % plot_data_by_cluster(param_props,base_grid,file_date,float_file_ext,...
% %     num_clusters,numWorkers_predict,start_year,end_year);
% % develop k-fold evaluation indices
% kfold_split_data(param_props,file_date,float_file_ext,...
%     start_year,end_year,glodap_only,num_clusters,num_folds,thresh);

%% k-fold train models for evaluation stats
% % feed-forward neural networks
% train_gobai('FFNN',param_props,file_date,float_file_ext,...
%     num_clusters,variables,thresh,numWorkers_train,snap_date,...
%     start_year,end_year,'reduce_data',data_per,'train_ratio',...
%     train_ratio,'val_ratio',val_ratio,'test_ratio',test_ratio,...
%     'num_folds',num_folds);

%% train models to create GOBAI product
% feed-forward neural networks
train_gobai('FFNN',param_props,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,...
    start_year,end_year,'reduce_data',data_per,'train_ratio',...
    train_ratio,'val_ratio',val_ratio,'test_ratio',test_ratio);

%% estimate parameter on grid to create GOBAI product
% % feed-forward neural networks
% predict_gobai('FFNN',param_props,fpaths,base_grid,file_date,float_file_ext,...
%     num_clusters,variables,thresh,numWorkers_predict,clust_vars,start_year,...
%     end_year,snap_date,'train_ratio',train_ratio,'val_ratio',val_ratio,...
%     'test_ratio',test_ratio);
% plot_gobai_animation(param_props,fpaths.param_path,base_grid,num_clusters,'FFNN',...
%     file_date,float_file_ext,numWorkers_predict,'train_ratio',train_ratio,...
%     'val_ratio',val_ratio,'test_ratio',test_ratio);

%% extract data from models
run_osse(fpaths,model_types,model_folders,realizations,...
    grid_labels,grid_types,param_props,file_date,snap_date,...
    float_file_ext,start_year,end_year,num_clusters,variables,...
    clust_vars,train_ratio,val_ratio,test_ratio,thresh,...
    data_per_osse,numWorkers_train,numWorkers_predict,1,1,1);

%% finish timing
toc(t_whole_script)
