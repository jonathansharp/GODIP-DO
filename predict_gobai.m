% predict_gobai
%
% DESCRIPTION:
% This function uses trained models to predict
% on a grid.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 2/6/2025

function predict_gobai(alg_type,param_props,param_path,temp_path,sal_path,...
    base_grid,file_date,float_file_ext,num_clusters,variables,thresh,...
    numWorkers_predict,clust_vars,start_year,end_year,snap_date,varargin)

%% process optional input arguments
% pre-allocate
rlz = NaN;
% process inputs
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'rlz')
        rlz = varargin{i+1};
    end
end

%% process necessary input arguments for model parameters
% pre-allocate
train_ratio = NaN;
val_ratio = NaN;
test_ratio = NaN;
numtrees = NaN;
minLeafSize = NaN;
numstumps = NaN;
numbins = NaN;
% process based on algorithm type
if strcmp(alg_type,'FFNN')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'train_ratio')
            train_ratio = varargin{i+1};
        elseif strcmpi(varargin{i}, 'val_ratio')
            val_ratio = varargin{i+1};
        elseif strcmpi(varargin{i}, 'test_ratio')
            test_ratio = varargin{i+1};
        end
    end
elseif strcmp(alg_type,'RFR')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'numtrees')
            numtrees = varargin{i+1};
        elseif strcmpi(varargin{i}, 'minLeafSize')
            minLeafSize = varargin{i+1};
        end
    end
elseif strcmp(alg_type,'GBM')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'numstumps')
            numstumps = varargin{i+1};
        elseif strcmpi(varargin{i}, 'numbins')
            numbins = varargin{i+1};
        end
    end
else
    disp('"alg_type" must be "FFNN", "RFR", or "GBM"')
end

%% process date
date_str = num2str(snap_date);

%% directory base
if strcmp(alg_type,'FFNN')
    dir_base = create_dir_base(alg_type,{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
elseif strcmp(alg_type,'RFR')
    dir_base = create_dir_base(alg_type,{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
elseif strcmp(alg_type,'GBM')
    dir_base = create_dir_base(alg_type,{base_grid;num_clusters;file_date;...
        float_file_ext;numstumps;numbins});
end

%% create directory and file names
alg_dir = [param_props.dir_name '/Models/' dir_base];
alg_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    alg_fnames(c) = ...
        {[alg_type '_' param_props.file_name '_C' num2str(c)]};
end
if strcmp(alg_type,'FFNN')
    gobai_alg_dir = ...
        [param_path 'GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) ...
        '_' file_date float_file_ext '/train' num2str(100*train_ratio) ...
        '_val' num2str(100*val_ratio) '_test' num2str(100*test_ratio) '/'];
elseif strcmp(alg_type,'RFR')
    gobai_alg_dir = ...
        [param_path 'GOBAI/' base_grid '/RFR/c' num2str(num_clusters) ...
        '_' file_date float_file_ext '/tr' num2str(numtrees) '_lf' ...
        num2str(minLeafSize) '/'];
elseif strcmp(alg_type,'GBM')
    gobai_alg_dir = ...
        [param_path 'GOBAI/' base_grid '/GBM/c' num2str(num_clusters) ...
        '_' file_date float_file_ext '/tr' num2str(numstumps) ...
        '_bin' num2str(numbins) '/'];
end

%% define variables for predictions
variables_TS = cell(size(variables));
for v = 1:length(variables)
    variables_TS{v} = [variables{v} '_array'];
end

%% create netCDF file that will be end result
TS = load_EN4_dim(temp_path,start_year,end_year);
% create file
create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_alg_dir,param_props);

%% set up parallel pool
tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;

%% start timing predictions
tStart = tic;

%% compute and save estimates for each month
parfor m = 1:length(TS.months)
    % load dimensions
    TS = load_EN4_dim(temp_path,start_year,end_year);
    % index to above 2000m
    idx_depth = find(TS.Depth < 2000);
    TS.Depth = TS.Depth(idx_depth);
    TS.zdim = length(idx_depth);
    TS = replicate_dims(base_grid,TS,1);
    % get monthly T and S
    TS.Temperature = ncread([temp_path '/EN.4.2.2.f.analysis.c14.' ...
        num2str(TS.years(m)) sprintf('%02d',TS.months(m)) '.nc'],...
        'temperature',[1 1 1 1],[Inf Inf idx_depth(end) 1])-273.15;
    TS.Salinity = ncread([sal_path '/EN.4.2.2.f.analysis.c14.' ...
        num2str(TS.years(m)) sprintf('%02d',TS.months(m)) '.nc'],...
        'salinity',[1 1 1 1],[Inf Inf idx_depth(end) 1]);
    % covert RG T and S to conservative Temperature and absolute Salinity
    % lon_3d = repmat(TS.Longitude,1,TS.ydim,TS.zdim);
    % lat_3d = repmat(TS.Latitude',TS.xdim,1,TS.zdim);
    % TS.Salinity_abs = gsw_SA_from_SP(TS.Salinity,TS.pressure,lon_3d,lat_3d);
    % TS.Temperature_cns = gsw_CT_from_t(TS.Salinity_abs,TS.Temperature,TS.pressure);
    % get time variables for just this timestep
    date_temp = datevec(TS.Time(m));
    date_temp0 = date_temp;
    date_temp0(:,2:3) = 1; % Jan. 1 of each year
    TS.year = date_temp(:,1);
    TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
    % transform day
    TS.day_sin = sin((2.*pi.*TS.day)/365.25);
    TS.day_cos = cos((2.*pi.*TS.day)/365.25);
    % apply model
    apply_model(alg_type,TS,num_clusters,alg_dir,alg_fnames,...
        base_grid,m,1,m,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
        thresh,gobai_alg_dir,param_props,param_path,date_str,clust_vars);
end

% end parallel session
delete(gcp('nocreate'));

%% concatenate gobai in main file
files = dir([gobai_alg_dir '/*.nc']); % count files in folder
filename = [gobai_alg_dir 'gobai-' param_props.file_name '.nc'];
for cnt = 1:length(files)-1
    % define file name
    filename_temp = [gobai_alg_dir 'gobai-' param_props.file_name '-' num2str(cnt) '.nc'];
    % read information from temporary file and write it to main file
    time = ncread(filename_temp,'time'); % read
    ncwrite(filename,'time',time,cnt); % write
    gobai_3d = ncread(filename_temp,param_props.file_name); % read
    ncwrite(filename,param_props.file_name,gobai_3d,[1 1 1 cnt]); % write
    % delete temporary file
    delete(filename_temp);
end

%% create monthly,one degree file
% for x = 1:length(TS.Longitude)
%     for y = 1:length(TS.Latitude)
%         for z = 1:length(TS.Pressure)
%             var = squeeze(ncread(filename,param_props.file_name,[x y z 1],[1 1 1 Inf]));
%             if any(~isnan(var))
% 
%             end
%         end
%     end
% end

%% stop timing predictions
fprintf([alg_type ' Prediction (' num2str(start_year) ' to ' date_str(1:4) '): ']);

tElapsed = toc(tStart);
disp(['Elapsed time is ' num2str(tElapsed/60) ' minutes.'])

end

%% for processing 3D grids and applying trained models to them
function apply_model(alg_type,TS,num_clusters,alg_dir,alg_fnames,...
    base_grid,m,w,cnt,xdim,ydim,zdim,variables_TS,thresh,gobai_alg_dir,...
    param_props,param_path,date_str,clust_vars)
    
    % define folder name
    folder_name = [param_path 'GMM_' base_grid '_' num2str(num_clusters)];

    % convert to arrays
    TS_index = ~isnan(TS.Temperature);
    vars = fieldnames(TS);
    for v = 1:length(vars)
        if ndims(TS.(vars{v})) == 3
            TS.([vars{v} '_array']) = TS.(vars{v})(TS_index);
            TS = rmfield(TS,vars{v});
        end
    end

    % replicate time variables as arrays
    vars = fieldnames(TS);
    for v = 1:length(vars)
        if isscalar(TS.(vars{v}))
            TS.([vars{v} '_array']) = repmat(TS.(vars{v}),size(TS.Temperature_array));
            TS = rmfield(TS,vars{v});
        end
    end

    % calculate absolute Salinity, conservative Temperature, potential density
    TS.Salinity_abs_array = gsw_SA_from_SP(TS.Salinity_array,TS.pressure_array,TS.lon_array,TS.lat_array);
    TS.Temperature_cns_array = gsw_CT_from_pt(TS.Salinity_abs_array,TS.Temperature_array);
    TS.sigma_array = gsw_sigma0(TS.Salinity_abs_array,TS.Temperature_cns_array);

    % pre-allocate
    gobai_matrix = single(nan(length(TS.Temperature_array),num_clusters));
    probs_matrix = single(nan(length(TS.Temperature_array),num_clusters));

    % apply GMM model for RFROM basegrid
    if strcmp(base_grid,'RFROM')
        % load GMM model
        load([param_props.dir_name '/Data/GMM_' base_grid '_' ...
            num2str(num_clusters) '/model_' date_str],'gmm','C','S');
        % transform to normalized arrays
        predictor_matrix = [];
        for v = 1:length(clust_vars)
            predictor_matrix = [predictor_matrix TS.([clust_vars{v} '_array'])];
        end
        X_norm = normalize(predictor_matrix,'Center',C,'Scale',S);
        % assign to clusters and obtain probabilities
        [~,~,p] = cluster(gmm,X_norm);
    end
    clear predictor_matrix gmm X_norm

    % apply models for each cluster
    for c = 1:num_clusters

      % check for existence of model
      if isfile([alg_dir '/' alg_fnames{c} '.mat'])

        % load GMM cluster probabilities for this cluster and month, and convert to array
        % load GMM model
        load([param_props.dir_name '/Data/GMM_' base_grid '_' ...
            num2str(num_clusters) '/model_' date_str],'gmm','C','S');
        % transform predictors to normalized arrays
        predictor_matrix = [];
        for v = 1:length(clust_vars)
            predictor_matrix = [predictor_matrix TS.([clust_vars{v} '_array'])];
        end
        X_norm = normalize(predictor_matrix,'Center',C,'Scale',S);
        % assign predictors to clusters and obtain probabilities
        [~,~,p] = cluster(gmm,X_norm);
        GMM_probs = nan(size(TS_index));
        GMM_probs(TS_index) = int16(p(:,c)*10000);
        probabilities_array = GMM_probs(TS_index)./10000; % convert to array
        probabilities_array(probabilities_array < thresh) = NaN; % remove probabilities below thresh
        probs_matrix(:,c) = probabilities_array; % add to probability matrix

        % load model for this cluster
        alg_struct = load([alg_dir '/' alg_fnames{c}],alg_type);
    
        % predict data for each cluster
        if strcmp(alg_type,'FFNN')
            gobai_matrix(:,c) = ...
                run_FFNN(alg_struct.FFNN,TS,probabilities_array,...
                true(size(TS.Temperature_array)),variables_TS,thresh);
        elseif strcmp(alg_type,'RFR')
            gobai_matrix(:,c) = ...
                run_RFR(alg_struct.RFR,TS,probabilities_array,...
                true(size(TS.Temperature_array)),variables_TS,thresh);
        elseif strcmp(alg_type,'GBM')
            gobai_matrix(:,c) = ...
                run_GBM(alg_struct.GBM,TS,probabilities_array,...
                true(size(TS.Temperature_array)),variables_TS,thresh);
        end

      end

    end

    % calculate weighted average over each cluster using probabilities
    gobai_array = ...
        double(sum(gobai_matrix.*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
    
    % convert back to 3D grid
    gobai_3d = nan(xdim,ydim,zdim);
    gobai_3d(TS_index) = gobai_array;

    % Write output in temporary files
    filename = [gobai_alg_dir 'gobai-' param_props.file_name '-' num2str(cnt) '.nc'];
    if exist(filename,'file')==2; delete(filename); end
    nccreate(filename,'time','Dimensions',{'time' 1});
    if strcmp(base_grid,'RFROM')
        ncwrite(filename,'time',TS.Time(cnt));
    else
        ncwrite(filename,'time',datenum(TS.years(cnt),TS.months(cnt),15)-datenum(1950,0,0));
    end
    nccreate(filename,param_props.file_name,'Dimensions',{'lon' xdim 'lat' ydim 'pres' zdim});
    ncwrite(filename,param_props.file_name,gobai_3d);

    % display information
    fprintf([alg_type ' Prediction (Month ' num2str(m) ', Week ' num2str(w) ')\n']);

end

%% for creating main netCDF file
function create_nc_file(TS,base_grid,xdim,ydim,zdim,gobai_alg_dir,...
    param_props)

% define file name
filename = [gobai_alg_dir 'gobai-' param_props.file_name '.nc'];

% create folder and file
if ~isfolder([gobai_alg_dir]); mkdir(gobai_alg_dir); end
if isfile(filename); delete(filename); end % delete file if it exists
% bgc parameter
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    nccreate(filename,param_props.file_name,'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
        'DataType','single','FillValue',NaN);
else
    nccreate(filename,param_props.file_name,'Dimensions',{'lon',xdim,'lat',ydim,'depth',zdim,'time',Inf},...
        'DataType','single','FillValue',NaN);
end
ncwriteatt(filename,param_props.file_name,'units',param_props.units);
ncwriteatt(filename,param_props.file_name,'long_name',param_props.long_param_name);
% longitude
nccreate(filename,'lon','Dimensions',{'lon',xdim},...
    'DataType','single','FillValue',NaN);
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
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    nccreate(filename,'pres','Dimensions',{'pres',zdim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'pres',TS.Pressure);
    ncwriteatt(filename,'pres','units','decibars');
    ncwriteatt(filename,'pres','axis','Z');
    ncwriteatt(filename,'pres','long_name','pressure');
    ncwriteatt(filename,'pres','_CoordinateAxisType','Pres');
else
    nccreate(filename,'depth','Dimensions',{'depth',zdim},...
    'DataType','single','FillValue',NaN);
    ncwrite(filename,'depth',TS.Depth);
    ncwriteatt(filename,'depth','units','meters');
    ncwriteatt(filename,'depth','axis','Z');
    ncwriteatt(filename,'depth','long_name','depth');
    ncwriteatt(filename,'depth','_CoordinateAxisType','Depth');
end
% time
nccreate(filename,'time','Dimensions',{'time',Inf},...
    'DataType','single','FillValue',NaN);
ncwriteatt(filename,'time','units','days since 1950-0-0');
ncwriteatt(filename,'time','axis','T');
ncwriteatt(filename,'time','long_name','time');
ncwriteatt(filename,'time','_CoordinateAxisType','Time');

end
