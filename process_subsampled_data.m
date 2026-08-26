function process_subsampled_data(fpaths,model_names,model_file_names,...
    model_realizations,start_year,end_year)

% loop through models
for m = 1:length(model_names)

% download from online
url = 'https://tigress-web.princeton.edu';
dir_sub = '/GEOCLIM/db9274/Monthly_subsampled_fields_1965_2021/';
o2_ext = '/o2_1x1bin_';
so_ext = '/so_1x1bin_';
thetao_ext = '/thetao_1x1bin_';
if ~isfolder([fpaths.model_path model_names{m}])
    mkdir([fpaths.model_path model_names{m}]);
end
% oxygen
web_name_o2 = [url dir_sub 'o2' o2_ext model_file_names{m} '_' model_realizations{m} '.nc'];
local_name_o2 = [fpaths.model_path model_names{m} o2_ext model_file_names{m} '_' model_realizations{m} '.nc'];
websave(local_name_o2,web_name_o2);
% salinity
web_name_so = [url dir_sub 'so' so_ext model_file_names{m} '_' model_realizations{m} '.nc'];
local_name_so = [fpaths.model_path model_names{m} so_ext model_file_names{m} '_' model_realizations{m} '.nc'];
websave(local_name_so,web_name_so);
% temperature
web_name_thetao = [url dir_sub 'thetao' thetao_ext model_file_names{m} '_' model_realizations{m} '.nc'];
local_name_thetao = [fpaths.model_path model_names{m} thetao_ext model_file_names{m} '_' model_realizations{m} '.nc'];
websave(local_name_thetao,web_name_thetao);

% load dimesions
lon = ncread(local_name_o2,'x');
lat = ncread(local_name_o2,'y');
depth = ncread(local_name_o2,'lev');
time = ncread(local_name_o2,'time');

% load and process variables and dimensions
vars = {'o2' 'thetao' 'so'};

% loop through each year and extract data from 4d NetCDF
for t = 1:12:length(time)
    % read each variable
    for v = 1:length(vars)
        if strcmp(vars{v},'o2')
            model.([vars{v} '_temp']) = ...
                single(ncread(local_name_o2,vars{v},...
                [1 1 1 t],[Inf Inf Inf 12]));
        elseif strcmp(vars{v},'so')
            model.([vars{v} '_temp']) = ...
                single(ncread(local_name_so,vars{v},...
            [1 1 1 t],[Inf Inf Inf 12]));
        elseif strcmp(vars{v},'thetao')
            model.([vars{v} '_temp']) = ...
                single(ncread(local_name_thetao,vars{v},...
            [1 1 1 t],[Inf Inf Inf 12]));
        end
    end
    % read time
    model.time = single(ncread(local_name_o2,'time',t,12));
    % index where temp, sal, and o2 are available
    idx = ~isnan(model.o2_temp) & ~isnan(model.thetao_temp) & ~isnan(model.so_temp);
    % expand dimensions to 4d
    [lon4d,lat4d,depth4d,time4d] = ndgrid(lon,lat,depth,model.time);
    lon4d = single(lon4d);
    lat4d = single(lat4d);
    depth4d = single(depth4d);
    time4d = single(time4d);
    % extract variables in arrays
    for v = 1:length(vars)
        data.(vars{v}) = model.([vars{v} '_temp'])(idx);
    end
    % extract dimensions in arrays
    data.lon = lon4d(idx);
    data.lat = lat4d(idx);
    data.depth = depth4d(idx);
    data.time = time4d(idx);
    % display year completed
    date = datevec(datenum(1850,1,double(time(t+11))));
    disp([model_names{m} ' subsampled for ' num2str(date(1)) '.'])
    % make directory to save netcdf
    if ~exist('O2/Data','dir'); mkdir('O2/Data'); end
    % define variables to save
    nc_vars = fieldnames(data);
    % define filename and delete if it exists
    filename = ['O2/Data/' model_names{m} '_' num2str(date(1)) '_data.nc'];
    if exist(filename,'file'); delete(filename); end
    % save variables in netcdf
    for v = 1:length(nc_vars)
        nccreate(filename,nc_vars{v},'Dimensions',{'observations' length(data.(nc_vars{v}))});
        ncwrite(filename,nc_vars{v},data.(nc_vars{v}));
    end
    % clear variables
    clear data lon4d lat4d depth4d time4d model idx
end

% remove 4d netcdf files
delete(local_name_o2);
delete(local_name_so);
delete(local_name_thetao);

% download full model fields (for use later)
dir_full = '/GEOCLIM/db9274/Monthly_full_fields_1965_2021/';
% oxygen
web_name_o2 = [url dir_full 'o2' o2_ext model_file_names{m} '_' model_realizations{m} '.nc'];
local_name_o2 = [fpaths.model_path model_names{m} o2_ext model_file_names{m} '_' model_realizations{m} '.nc'];
websave(local_name_o2,web_name_o2);
% salinity
web_name_so = [url dir_full 'so' so_ext model_file_names{m} '_' model_realizations{m} '.nc'];
local_name_so = [fpaths.model_path model_names{m} so_ext model_file_names{m} '_' model_realizations{m} '.nc'];
websave(local_name_so,web_name_so);
% temperature
web_name_thetao = [url dir_full 'thetao' thetao_ext model_file_names{m} '_' model_realizations{m} '.nc'];
local_name_thetao = [fpaths.model_path model_names{m} thetao_ext model_file_names{m} '_' model_realizations{m} '.nc'];
websave(local_name_thetao,web_name_thetao);

% pre-allocate all data structure
for v = 1:length(nc_vars)
    all_data.(nc_vars{v}) = [];
end

% add annual arrays to all data structure
for t = 1:12:length(time)
    date = datevec(datenum(1850,1,double(time(t+11))));
    filename = ['O2/Data/' model_names{m} '_' num2str(date(1)) '_data.nc'];
    for v = 1:length(nc_vars)
        all_data.(nc_vars{v}) = [all_data.(nc_vars{v});...
            ncread(filename,nc_vars{v})];
    end
    delete(filename);
end

% calculate pressure
all_data.pressure = gsw_p_from_z(-all_data.depth,all_data.lat);

% calculate potential density, conservative temperature, absolute salinity
all_data.sal_abs = gsw_SA_from_SP(all_data.so,all_data.pressure,all_data.lon,all_data.lat);
all_data.tmp_cns = gsw_CT_from_pt(all_data.sal_abs,all_data.thetao);
all_data.Temperature = gsw_t_from_pt0(all_data.sal_abs,all_data.thetao,all_data.pressure);
all_data.sigma = gsw_sigma0(all_data.sal_abs,all_data.tmp_cns);

% calculate sin and cos of lon and day
all_data.lon_cos_1 = cosd(all_data.lon-20);
all_data.lon_cos_2 = cosd(all_data.lon-110);
date = datevec(datenum(1850,1,15+all_data.time));
date0 = date; date0(:,2:3) = 0;
all_data.day_of_year = datenum(date) - datenum(date0);
all_data.day_sin = sin((2.*pi.*all_data.day_of_year)/365.25);
all_data.day_cos = cos((2.*pi.*all_data.day_of_year)/365.25);
all_data.year = date(:,1);
clear date date0

% plot all subsampled data (binned)
lon_edges = 0:360; lon = 0.5:359.5;
lat_edges = -90:90; lat = -89.5:89.5;
[~,~,Xnum] = histcounts(all_data.lon,lon_edges);
[~,~,Ynum] = histcounts(all_data.lat,lat_edges);
% accumulate 3D grid of test data point errors
subs = [Xnum, Ynum];
idx_subs = any(subs==0,2);
sz = [length(lon),length(lat)];
all_data.num_obs = accumarray(subs(~idx_subs,:),1,sz);
clear subs sz
% make plot
figure('visible','on','Position',[100 100 800 400]); hold on
worldmap([-90 90],[20 380]);
set(gca,'ColorScale','log');
title('number of observations per 1x1 degee grid cell');
[lon_temp,z] = reformat_lon(lon,all_data.num_obs,20);
pcolorm(lat,[lon_temp lon_temp(end)+1],[z;z(end,:)]');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar; 
cmap = [1,1,1;parula]; colormap(cmap);
mlabel off; plabel off;
if ~exist(['O2/Figures/' model_names{m}],'dir'); mkdir(['O2/Figures/' model_names{m}]); end
export_fig(gcf,['O2/Figures/' model_names{m} '/data_density_fig.png']); hold off;
close;

% define variables to save
nc_vars = fieldnames(all_data);
% define filename and delete if it exists
filename = ['O2/Data/' model_names{m} '_data_' num2str(start_year) ...
    '_' num2str(end_year) '.nc'];
if exist(filename,'file'); delete(filename); end
% save all subsampled data in netcdf
for v = 1:length(nc_vars)
    if size(all_data.(nc_vars{v}),2) == 1
        nccreate(filename,nc_vars{v},'Dimensions',{'observations' ...
            length(all_data.(nc_vars{v}))});
        ncwrite(filename,nc_vars{v},all_data.(nc_vars{v}));
    end
end

disp(['Finished processing subsampled data for ' model_names{m}]);

end