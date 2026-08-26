%% Run all scripts to make GOBAI-O2 with subsampled data from Dhruv
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
clust_vars = {'thetao' 'sal_abs' 'depth'};
glodap_only = false; % EDIT THIS TO 'true' TO TEST WITH GLODAP DATA ONLY
thresh = 0.05;
num_folds = 5;
% algorithm training configuration
variables = ... % variables for algorithms
    {'lat' 'lon_cos_1' 'lon_cos_2' 'depth' 'sigma' ...
    'thetao' 'sal_abs' 'day_sin' 'day_cos' 'year'};
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
data_per = 0.01; % set data reduction
data_per_osse = 0.01; % set data reduction for osse
param = 'o2';
param_props = param_config(param);
% base grid
base_grid = 'IAP';
fpaths = path_config(base_grid,param);
fpaths.model_path = '/data4/model/';
fpaths.param_path = '/data4/o2/GODIP-DO/';
% define models to use
model_names = {'ACCESS-ESM1-5' 'GFDL-ESM4' 'NorESM2-MM' 'MPI-ESM1-2-HR' ...
    'CanESM5-CanOE' 'MRI-ESM2' 'MIROC-ES2L' 'UKESM1.0LL' 'GFDL-CM4'};
model_file_names = {'ACCESS_ESM1.5' 'GFDL_ESM4' 'NorESM2_MM' 'MPI_ESM1_2_HR' ...
    'CanESM5_CanOE' 'MRI_ESM2_0' 'MIROC_ES2L' 'UKESM1.0LL' 'GFDL_CM4'};
model_realizations = {'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1' ...
    'r1i1p2f1' 'r1i2p1f1' 'r1i1p1f2' 'r1i1p1f2' 'r1i1p1f1'};
% suppress figure warning
warning('off','MATLAB:hg:AutoSoftwareOpenGL');

%% process subsampled model data
% process_subsampled_data(fpaths,model_names,model_file_names,...
%     model_realizations,start_year,end_year);

% Class 1 (more agreement with observed decadal variability):
% GFDL-CM4
% NorESM2-LM
% MPIESM-LR
% 
% Class 2:
% CanESM5
% CNRM-ESM2_1
% IPSL-CM6A-LR

%% loop through each model
for m = 5%1:length(model_names)

%% form clusters
gmm_clustering(param_props,fpaths,start_year,end_year,...
    snap_date,float_file_ext,clust_vars,num_clusters,numWorkers_predict,...
    model_names{m});

%% cluster data
assign_data_to_clusters(param_props,snap_date,float_file_ext,clust_vars,...
    num_clusters,start_year,end_year,model_names{m});
plot_data_by_cluster(param_props,file_date,float_file_ext,num_clusters,...
    numWorkers_train,start_year,end_year,model_names{m});
% plot_data_over_clusters(param,model_types{m},file_date,float_file_ext,num_clusters,numWorkers_train);

%% train algorithms
% train feed forward neural networks
train_gobai('FFNN',param_props,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,...
    start_year,end_year,model_names{m},'reduce_data',data_per_osse,'train_ratio',...
    train_ratio,'val_ratio',val_ratio,'test_ratio',test_ratio);

%% predict on grid
% predict on grid using feed forward neural networks
predict_gobai('FFNN',param_props,fpaths,...
    model_names{m},file_date,float_file_ext,num_clusters,...
    variables,thresh,numWorkers_predict,clust_vars,start_year,end_year,...
    snap_date,model_file_names{m},'train_ratio',train_ratio,'val_ratio',val_ratio,...
    'test_ratio',test_ratio,'rlz',model_realizations{m});

%% print finished model
disp(['Finished creating GOBAI for ' model_names{m}]);

end

%% compare reconstructions



%% finish timing
toc(t_whole_script)
