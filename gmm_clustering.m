function gmm_clustering(param_props,fpaths,base_grid,start_year,end_year,snap_date,...
    float_file_ext,clust_vars,num_clusters,numWorkers_predict,vrs,varargin)

%% process optional input arguments
% pre-allocate
rlz = NaN;
% process inputs
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'rlz')
        rlz = varargin{i+1};
    end
end

%% process date
date_str = num2str(snap_date);
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');

%% load basin mask file


%% load data
load(['O2/Data/' vrs '_data_' num2str(start_year) '_' num2str(end_year) '.mat']);

%% fit GMM from data points themselves
tic
% transform to normalized arrays
idx = false(length(all_data.Temperature),1);
for v = 1:length(clust_vars)
    idx = idx | ~isnan(all_data.(clust_vars{v}));
end
predictor_matrix = [];
for v = 1:length(clust_vars)
    predictor_matrix = [predictor_matrix all_data.(clust_vars{v})];
end
[X_norm,C,S] = normalize(predictor_matrix);
% reduce inputs for model training to 100,000 random data points
idx_rand = randperm(length(X_norm),100000)';
X_norm = X_norm(idx_rand,:);
% fit GMM
options = statset('MaxIter',1000); % increase max iterations to ensure convergence
gmm = fitgmdist(X_norm,num_clusters,...
    'Options',options,...
    'CovarianceType','diagonal',...
    'SharedCovariance',true,'Replicates',10);
% save GMM model
if ~isfolder([param_props.dir_name '/Data/GMM_' vrs '_' num2str(num_clusters)])
    mkdir([param_props.dir_name '/Data/GMM_' vrs '_' num2str(num_clusters)]);
end
save([param_props.dir_name '/Data/GMM_' vrs '_' num2str(num_clusters) '/model_' date_str],...
    'gmm','num_clusters','C','S','-v7.3');
clear gmm C S
toc

%% assign grid cells and probabilities to clusters
folder_name = [fpaths.param_path 'GMM_' base_grid '_' num2str(num_clusters)];
% determine length of cluster file if it exists
if exist([folder_name '/clusters.nc'],'file') == 2
    inf = ncinfo([folder_name '/clusters.nc']);
    for n = 1:length(inf.Dimensions)
        if strcmp(inf.Dimensions(n).Name,'time')
            time_idx = n;
        end
    end
end
% determine expected length
% year = str2double(date_str(1:4));
% month = str2double(date_str(5:6));
% length_expt = (year-start_year)*12 + month;
length_expt = (end_year-start_year)*12 + 12;

% delete file if it exists
if exist([folder_name '/clusters.nc'],'file') == 2
    delete([folder_name '/clusters.nc']);
end

% start timing cluster assignment
tic

% %% create netCDF file that will be end result
TS = load_EN4_dim(fpaths.temp_path,start_year,end_year);
% create file
create_nc_files(TS,num_clusters,base_grid,TS.xdim,TS.ydim,TS.zdim,folder_name);

% load GMM model
load([param_props.dir_name '/Data/GMM_' vrs '_' num2str(num_clusters) '/model_' ...
    date_str],'gmm','C','S');

% set up parallel pool
tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;

%% assign clusters for each timestep
parfor t = 1:TS.tdim

    % load dimensions
    TS = load_EN4_dim(fpaths.temp_path,start_year,end_year);

    % get salinity and temperature
    try
        TS.Temperature = ncread([fpaths.temp_path 'EN.4.2.2.f.analysis.c14.' ...
            num2str(TS.years(t)) sprintf('%02d',TS.months(t)) '.nc'],...
            'temperature');
        TS.Salinity = ncread([fpaths.temp_path 'EN.4.2.2.f.analysis.c14.' ...
            num2str(TS.years(t)) sprintf('%02d',TS.months(t)) '.nc'],...
            'salinity');
    catch
        TS.Temperature = ncread([fpaths.temp_path 'EN.4.2.2.p.analysis.c14.' ...
            num2str(TS.years(t)) sprintf('%02d',TS.months(t)) '.nc'],...
            'temperature');
        TS.Salinity = ncread([fpaths.temp_path 'EN.4.2.2.p.analysis.c14.' ...
            num2str(TS.years(t)) sprintf('%02d',TS.months(t)) '.nc'],...
            'salinity');
    end
    TS.depth = double(repmat(permute(TS.Depth,[3 2 1]),length(TS.Longitude),...
        length(TS.Latitude)));

    % transform to normalized arrays
    idx = ~isnan(TS.Temperature) & ~isnan(TS.Salinity);
    predictor_matrix = [];
    for v = 1:length(clust_vars)
        predictor_matrix = [predictor_matrix TS.(clust_vars{v})(idx)];
    end
    X_norm = normalize(predictor_matrix,'Center',C,'Scale',S);

    % assign data points to clusters
    assign_to_gmm_clusters(TS,base_grid,gmm,num_clusters,idx,X_norm,t,folder_name);

end

% end parallel session
delete(gcp('nocreate'));

%% concatenate cluster information in main file
for t = 1:TS.tdim
    % define file names
    filename = [folder_name '/clusters.nc'];
    filename_temp = [folder_name '/clust_' num2str(t) '.nc'];
    % read information from temporary file and write it to main file
    time = ncread(filename_temp,'time'); % read
    ncwrite(filename,'time',time,t); % write
    GMM_clusters = ncread(filename_temp,'GMM_clusters'); % read
    ncwrite(filename,'clusters',GMM_clusters,[1 1 1 t]); % write
    for c = 1:num_clusters
        GMM_cluster_probs = ncread(filename_temp,...
            ['GMM_cluster_probs_' num2str(c)]); % read
        ncwrite(filename,['cluster_probs_c' num2str(c)],...
            GMM_cluster_probs,[1 1 1 t]); % write
    end
    % delete temporary file
    delete(filename_temp);
end

% display information
toc
disp([num2str(num_clusters) ' clusters formed using ' base_grid ' grid']);

end

%% embedded function to assign points to GMM clusters
function assign_to_gmm_clusters(TS,base_grid,gmm,num_clusters,idx,X_norm,t,folder_name)
    % assign to clusters and obtain probabilities
    [clusters,~,p] = cluster(gmm,X_norm);
    % fill 3D clusters (highest probability cluster)
    GMM_clusters = nan(TS.xdim,TS.ydim,TS.zdim);
    GMM_clusters(idx) = uint8(clusters);
    % save cluster properties in temporary files
    filename = [folder_name '/clust_' num2str(t) '.nc'];
    if isfile(filename); delete(filename); end
    nccreate(filename,'time','Dimensions',{'time' 1});
    ncwrite(filename,'time',TS.Time(t));
    nccreate(filename,'GMM_clusters','Dimensions',{'lon' size(GMM_clusters,1) ...
        'lat' size(GMM_clusters,2) 'pres' size(GMM_clusters,3)});
    ncwrite(filename,'GMM_clusters',GMM_clusters);
    if strcmp(base_grid,'RFROM')
        % do nothing
    else
        % fill 3D probabilities (for each cluster)
        for c = 1:num_clusters
            GMM_cluster_probs = nan(TS.xdim,TS.ydim,TS.zdim);
            GMM_cluster_probs(idx) = int16(p(:,c)*10000);
            nccreate(filename,['GMM_cluster_probs_' num2str(c)],...
                'Dimensions',{'lon' size(GMM_cluster_probs,1) 'lat' ...
                size(GMM_cluster_probs,2) 'pres' size(GMM_cluster_probs,3)});
            ncwrite(filename,['GMM_cluster_probs_' num2str(c)],...
                GMM_cluster_probs);
    end
    end
end

% for creating netCDF file
function create_nc_files(TS,num_clusters,base_grid,xdim,ydim,zdim,folder_name)

    % define file name and create file
    if ~isfolder(folder_name); mkdir(folder_name); end
    filename = [folder_name '/clusters.nc'];

    % longitude
    nccreate(filename,'lon','Dimensions',{'lon',xdim},...
        'DataType','single','FillValue',NaN,'Format','netcdf4');
    ncwrite(filename,'lon',TS.Longitude);
    ncwriteatt(filename,'lon','units','degrees_east');
    ncwriteatt(filename,'lon','axis','X');
    ncwriteatt(filename,'lon','long_name','longitude');
    ncwriteatt(filename,'lon','_CoordinateAxisType','Lon');
    
    % latitude
    nccreate(filename,'lat','Dimensions',{'lat',ydim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'lat',TS.Latitude);
    ncwriteatt(filename,'lat','units','degrees_north');
    ncwriteatt(filename,'lat','axis','Y');
    ncwriteatt(filename,'lat','long_name','latitude');
    ncwriteatt(filename,'lat','_CoordinateAxisType','Lat');
    
    % pressure (or depth)
    nccreate(filename,'depth','Dimensions',{'depth',zdim},...
    'DataType','single','FillValue',NaN);
    ncwrite(filename,'depth',TS.Depth);
    ncwriteatt(filename,'depth','units','meters');
    ncwriteatt(filename,'depth','axis','Z');
    ncwriteatt(filename,'depth','long_name','depth');
    ncwriteatt(filename,'depth','_CoordinateAxisType','Depth');
    
    % time
    nccreate(filename,'time','Dimensions',{'time',Inf},...
        'DataType','single','FillValue',NaN);
    ncwriteatt(filename,'time','units','days since 0000-01-01');
    ncwriteatt(filename,'time','axis','T');
    ncwriteatt(filename,'time','long_name','time');
    ncwriteatt(filename,'time','_CoordinateAxisType','Time');

    % clusters
    if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
        nccreate(filename,'clusters','Dimensions',{'lon',xdim,'lat',ydim,...
            'pres',zdim,'time',Inf},'DataType','uint8','FillValue',NaN);
        for c = 1:num_clusters
            nccreate(filename,['cluster_probs_c' num2str(c)],'Dimensions',{'lon',xdim,'lat',ydim,...
                'pres',zdim,'time',Inf},'DataType','int16','FillValue',NaN);
        end
    else
        nccreate(filename,'clusters','Dimensions',{'lon',xdim,'lat',ydim,...
            'depth',zdim,'time',Inf},'DataType','uint8','FillValue',NaN);
        for c = 1:num_clusters
            nccreate(filename,['cluster_probs_c' num2str(c)],'Dimensions',{'lon',xdim,'lat',ydim,...
                'depth',zdim,'time',Inf},'DataType','int16','FillValue',NaN);
        end
    end

end
