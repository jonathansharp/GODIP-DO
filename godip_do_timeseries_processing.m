% this compares annual means of 

%% file properties
file_date = 'Jan2026';
fpath = '/fast4/o2/GODIP-DO/O2_Maps/';
fname = ['GODIP-DO_NCEI_JDS_' file_date '.nc'];
fname_seas = ['GODIP-DO_NCEI_SEASONAL_JDS_' file_date '.nc'];
products = {'ncei' 'iap' 'gt_oi' 'rb' 'sjtu_gr' 'gobai' 'gt_ml' 'jingwei' 'han_zhou'};
labels = {'NCEI' 'IAP' 'GT-OI' 'RB' 'SJTU-GR' 'GOBAI' 'GT-ML' 'Jingwei' 'HZ'};

%% import ncei dimensions
datasets.ncei.name = 'o_19652022_5yr_objetively_analyzed_NCEIdata_NCEImapping.nc';
datasets.ncei.lat = ncread([fpath datasets.ncei.name],'lat');
datasets.ncei.lon = ncread([fpath datasets.ncei.name],'lon');
datasets.ncei.time = double(ncread([fpath datasets.ncei.name],'time'));
datasets.ncei.depth = ncread([fpath datasets.ncei.name],'depth');
datasets.ncei.date = datevec(datenum(1900,1,1)+datasets.ncei.time);
datasets.ncei.param_name = 'o_anom_an';

%% import iap dimensions
datasets.iap.name = 'IAPmapping_applied_to_NCEIdata_5yr_mean_1967to2020_4D.nc';
datasets.iap.lat = ncread([fpath datasets.iap.name],'lat');
datasets.iap.lon = ncread([fpath datasets.iap.name],'lon');
datasets.iap.time = double(ncread([fpath datasets.iap.name],'time'));
datasets.iap.depth = ncread([fpath datasets.iap.name],'depth_std');
datasets.iap.date = datevec(datetime(datasets.iap.time,6,1));
datasets.iap.param_name = 'oxygen';

%% import gt-oi dimensions
datasets.gt_oi.name = 'Ito22mapping_NCEIdata_5yrmean_67C0.nc';
datasets.gt_oi.lat = ncread([fpath datasets.gt_oi.name],'lat');
datasets.gt_oi.lon = ncread([fpath datasets.gt_oi.name],'lon');
datasets.gt_oi.time = double(ncread([fpath datasets.gt_oi.name],'time'));
datasets.gt_oi.depth = ncread([fpath datasets.gt_oi.name],'depth');
datasets.gt_oi.date = datevec(datenum(1967,1,1)+datasets.gt_oi.time);
datasets.gt_oi.param_name = 'o2';

%% import rb dimensions
datasets.rb.name = 'o_annual_mean_depth_NCEI_Roach_Bindoff_1965_2022_v2.nc';
datasets.rb.lat = ncread([fpath datasets.rb.name],'lat');
datasets.rb.lon = ncread([fpath datasets.rb.name],'lon');
datasets.rb.time = double(ncread([fpath datasets.rb.name],'year'));
datasets.rb.depth = ncread([fpath datasets.rb.name],'depth');
datasets.rb.date = datevec(datetime(datasets.rb.time,6,1));
datasets.rb.param_name = 'oxy_conc';

%% import sjtu-gr dimensions
datasets.sjtu_gr.name = '5yrmeanSJTU/O2anom_1965to1969_mean_1x1_NCEIdata0_2000m_UK_SJTU_Zhou250717.nc';
datasets.sjtu_gr.lat = ncread([fpath datasets.sjtu_gr.name],'lat');
datasets.sjtu_gr.lon = ncread([fpath datasets.sjtu_gr.name],'lon');
datasets.sjtu_gr.depth = ncread([fpath datasets.sjtu_gr.name],'depth');
datasets.sjtu_gr.time = (1967:2020)';
datasets.sjtu_gr.date = datevec(datenum(datasets.sjtu_gr.time,6,1));
datasets.sjtu_gr.files = dir([fpath '5yrmeanSJTU']);
datasets.sjtu_gr.param_name = 'o2anom_estimate';

%% import gobai dimensions
datasets.gobai.name = 'gobai-o2-godip-do-v1.2-with-float-correction.nc';
datasets.gobai.lat = ncread([fpath datasets.gobai.name],'lat');
datasets.gobai.lon = ncread([fpath datasets.gobai.name],'lon');
datasets.gobai.time = ncread([fpath datasets.gobai.name],'time');
datasets.gobai.depth = ncread([fpath datasets.gobai.name],'depth');
datasets.gobai.date = datevec(datenum(1950,0,0)+datasets.gobai.time);
datasets.gobai.param_name = 'o2';

%% import gt-ml dimensions
datasets.gt_ml.name1 = 'GT_ML_NN_updated_01042026.nc';
datasets.gt_ml.name2 = 'GT_ML_RF_updated_01042026.nc';
datasets.gt_ml.lat = ncread([fpath datasets.gt_ml.name1],'lat');
datasets.gt_ml.lon = ncread([fpath datasets.gt_ml.name1],'lon');
datasets.gt_ml.time = double(ncread([fpath datasets.gt_ml.name1],'time'));
datasets.gt_ml.depth = ncread([fpath datasets.gt_ml.name1],'depth');
datasets.gt_ml.date = datevec(datenum(1965,0,15)+datasets.gt_ml.time);
datasets.gt_ml.param_name = 'o2';

%% import jingwei dimensions
datasets.jingwei.name = 'o2map_v1.2_5500m_sjtu_jingwei.nc';
datasets.jingwei.lat = ncread([fpath datasets.jingwei.name],'lat');
datasets.jingwei.lon = ncread([fpath datasets.jingwei.name],'lon');
datasets.jingwei.time = double(ncread([fpath datasets.jingwei.name],'time'));
datasets.jingwei.depth = ncread([fpath datasets.jingwei.name],'depth');
datasets.jingwei.date = datevec(datenum(1965,0,15)+datasets.jingwei.time);
datasets.jingwei.param_name = 'o2';

%% import han and zhou dimensions
datasets.han_zhou.name = 'DO_DOMIP_SJTU_Han&Zhou.nc';
datasets.han_zhou.lat = double(ncread([fpath datasets.han_zhou.name],'Latitude'));
datasets.han_zhou.lon = double(ncread([fpath datasets.han_zhou.name],'Longitude'));
datasets.han_zhou.time = double(ncread([fpath datasets.han_zhou.name],'Time'));
datasets.han_zhou.depth = double(ncread([fpath datasets.han_zhou.name],'Depth'));
datasets.han_zhou.date = datevec(datenum([repmat(1965,length(datasets.han_zhou.time),1) ...
    1+datasets.han_zhou.time repmat(15,length(datasets.han_zhou.time),1)]));
datasets.han_zhou.param_name = 'DO_Reconstruction';

%% define common grid
d = nan(105,4);
for p = 1:length(products)
    d(1:length(datasets.(products{p}).depth),p) = ...
        datasets.(products{p}).depth;
end
common.lon = (-179.5:1:179.5)';
common.lat = (-89.5:1:89.5)';
common.depth = datasets.gt_ml.depth;
common.time = datenum(1965:2025,6,1)';

%% create nc file for interpolated annual means and seasonal cycle
if exist([fpath fname],'file'); delete([fpath fname]); end
nccreate([fpath fname],'lon','Dimensions',{'lon' length(common.lon)},'Format','netcdf4');
ncwrite([fpath fname],'lon',common.lon);
nccreate([fpath fname],'lat','Dimensions',{'lat' length(common.lat)});
ncwrite([fpath fname],'lat',common.lat);
nccreate([fpath fname],'depth','Dimensions',{'depth' length(common.depth)});
ncwrite([fpath fname],'depth',common.depth);
nccreate([fpath fname],'time','Dimensions',{'time' length(common.time)});
ncwrite([fpath fname],'time',common.time);
nccreate([fpath fname],'months','Dimensions',{'months' 12});
ncwrite([fpath fname],'months',(1:12)');
nccreate([fpath fname],'products','Dimensions',{'products' length(products)},'Datatype','string');
ncwrite([fpath fname],'products',products);
disp('Common NetCDF file created.');

%% interpolate annual means to common grid
% pre-allocate
common.mask = true(length(common.lon),length(common.lat),...
    length(common.depth),length(common.time));
years = datevec(common.time); years = years(:,1);
% loop through and interpolate each product
for p = 1:length(products)
    % pre-allocate
    common.(products{p}) = nan(length(common.lon),length(common.lat),...
        length(common.depth),length(common.time));
    tic % start timer
    years_temp = datasets.(products{p}).date(:,1); % extract years
    for t = 1:length(common.time) % loop through years
        if strcmp(products{p},'gt_ml')
            common.(products{p})(:,:,:,t) = ...
                interp_to_common_grid({[fpath datasets.(products{p}).name1];...
                [fpath datasets.(products{p}).name2]},...
                datasets.(products{p}).param_name,years_temp,years(t),...
                datasets.(products{p}).lon,datasets.(products{p}).lat,...
                datasets.(products{p}).depth,common.lon,common.lat,common.depth,...
                products{p},datasets.(products{p}));
        else
            common.(products{p})(:,:,:,t) = ...
                interp_to_common_grid([fpath datasets.(products{p}).name],...
                datasets.(products{p}).param_name,years_temp,years(t),...
                datasets.(products{p}).lon,datasets.(products{p}).lat,...
                datasets.(products{p}).depth,common.lon,common.lat,common.depth,...
                products{p},datasets.(products{p}));
        end
        % if start or end longitude is all nan, replace with mean of
        % bounding longitudes
        if all(all(isnan(common.(products{p})(1,:,:,t))))
            common.(products{p})(1,:,:,t) = ...
                mean([common.(products{p})(end,:,:,t);...
                common.(products{p})(2,:,:,t)],1);
        elseif all(all(isnan(common.(products{p})(end,:,:,t))))
            common.(products{p})(end,:,:,t) = ...
                mean([common.(products{p})(end-1,:,:,t);...
                common.(products{p})(1,:,:,t)],1);
        end
        disp([products{p} ' interpolated for ' num2str(years(t))]);
    end
    % adjust mask only for each year where there is data
    yr_idx = squeeze(any(any(any(common.(products{p}),1),2),3));
    mask_temp = common.mask(:,:,:,yr_idx);
    mask_temp(isnan(common.(products{p})(:,:,:,yr_idx))) = false;
    common.mask(:,:,:,yr_idx) = mask_temp;
    clear mask_temp
    % plot mask
    figure;
    pcolor(common.lon,common.lat,double(common.mask(:,:,11,52))');
    shading flat; colorbar;
    title([labels{p} ' in ' datestr(common.time(52),'yyyy')]);
    % plot o2
    figure;
    pcolor(common.lon,common.lat,mean(common.(products{p})(:,:,11,52),4,'omitnan')');
    shading flat; colorbar;
    title([labels{p} ' in ' datestr(common.time(52),'yyyy')]);
    % save
    nccreate([fpath fname],products{p},'Dimensions',{'lon' 'lat' 'depth' 'time'});
    ncwrite([fpath fname],products{p},common.(products{p}));
    % clean up
    common = rmfield(common,products{p});
    toc
end
% save common mask
nccreate([fpath fname],'mask','Dimensions',{'lon' 'lat' 'depth' 'time'});
ncwrite([fpath fname],'mask',int8(common.mask));

%% interpolate to seasonal cycle for ML products
% pre-allocate
common.mask = true(length(common.lon),length(common.lat),...
    length(common.depth),length(common.time));
% loop through and interpolate each product
for p = 6:7
    % pre-allocate
    common.([products{p} '_baseline_seas']) = nan(length(common.lon),length(common.lat),...
        length(common.depth),length(common.time));
    common.([products{p} '_anom_seas']) = nan(length(common.lon),length(common.lat),...
        length(common.depth),length(common.time));
    tic % start timer
    years_temp = datasets.(products{p}).date(:,1); % extract years
    y1 = 1990; y2 = 2020; y_anom = 2025;
    if strcmp(products{p},'gt_ml')
    [common.([products{p} '_baseline_seas']),common.([products{p} '_anom_seas']),...
        common.([products{p} '_baseline_seas_var'])] = ...
        interp_seasonal_cycle_to_common_grid({[fpath datasets.(products{p}).name1];...
        [fpath datasets.(products{p}).name2]},...
        datasets.(products{p}).param_name,years_temp,y1,y2,y_anom,...
        datasets.(products{p}).lon,datasets.(products{p}).lat,...
        datasets.(products{p}).depth,common.lon,common.lat,common.depth,...
        products{p},datasets.(products{p}));
    else
    [common.([products{p} '_baseline_seas']),common.([products{p} '_anom_seas']),...
        common.([products{p} '_baseline_seas_var'])] = ...
        interp_seasonal_cycle_to_common_grid([fpath datasets.(products{p}).name],...
        datasets.(products{p}).param_name,years_temp,y1,y2,y_anom,...
        datasets.(products{p}).lon,datasets.(products{p}).lat,...
        datasets.(products{p}).depth,common.lon,common.lat,common.depth,...
        products{p},datasets.(products{p}));
    end
    % if start or end longitude is all nan, replace with mean of
    % bounding longitudes
    for m = 1:12
        if all(all(isnan(common.([products{p} '_baseline_seas'])(1,:,:,m))))
            common.([products{p} '_baseline_seas'])(1,:,:,m) = ...
                mean([common.([products{p} '_baseline_seas'])(end,:,:,m);...
                common.([products{p} '_baseline_seas'])(2,:,:,m)],1);
        elseif all(all(isnan(common.([products{p} '_baseline_seas'])(end,:,:,m))))
            common.([products{p} '_baseline_seas'])(end,:,:,m) = ...
                mean([common.([products{p} '_baseline_seas'])(end-1,:,:,m);...
                common.([products{p} '_baseline_seas'])(1,:,:,m)],m);
        end
        if all(all(isnan(common.([products{p} '_anom_seas'])(1,:,:,m))))
            common.([products{p} '_anom_seas'])(1,:,:,m) = ...
                mean([common.([products{p} '_anom_seas'])(end,:,:,m);...
                common.([products{p} '_anom_seas'])(2,:,:,m)],1);
        elseif all(all(isnan(common.([products{p} '_anom_seas'])(end,:,:,m))))
            common.([products{p} '_anom_seas'])(end,:,:,m) = ...
                mean([common.([products{p} '_anom_seas'])(end-1,:,:,m);...
                common.([products{p} '_anom_seas'])(1,:,:,m)],1);
        end
    end
    disp([products{p} ' seasonal cycle interpolated']);
    % save
    nccreate([fpath fname],[products{p} '_baseline_seas'],...
        'Dimensions',{'lon' 'lat' 'depth' 'months'});
    ncwrite([fpath fname],[products{p} '_baseline_seas'],common.([products{p} '_baseline_seas']));
    nccreate([fpath fname],[products{p} '_anom_seas'],...
        'Dimensions',{'lon' 'lat' 'depth' 'months'});
    ncwrite([fpath fname],[products{p} '_anom_seas'],common.([products{p} '_anom_seas']));
    nccreate([fpath fname],[products{p} '_baseline_seas_var'],...
        'Dimensions',{'lon' 'lat' 'depth' 'months'});
    ncwrite([fpath fname],[products{p} '_baseline_seas_var'],common.([products{p} '_baseline_seas_var']));
    % clean up
    common = rmfield(common,[products{p} '_baseline_seas']);
    common = rmfield(common,[products{p} '_anom_seas']);
    toc
end

%% embedded interpolation function
function o2_interp = interp_to_common_grid(file_directory,param_name,years_temp,...
    yr,lon,lat,depth,lon_interp,lat_interp,depth_interp,prod,struct)
    t_idx = find(years_temp == yr); % index to year of loop
    if ~isempty(t_idx)
        lon_temp = convert_lon(lon,'format','-180:180');
        [lon_temp,lon_idx] = sort(lon_temp); % sort longitude in order
        if strcmp(prod,'iap')
            o2_temp = ncread(file_directory,param_name,... % extract o2 for year
                [t_idx(1) 1 1 1],[t_idx(end)-t_idx(1)+1 Inf Inf Inf]);
            o2_temp = permute(o2_temp,[3 4 2 1]);
        elseif strcmp(prod,'rb')
            o2_temp = ncread(file_directory,param_name,... % extract o2 for year
                [1 1 t_idx(1) 1],[Inf Inf t_idx(end)-t_idx(1)+1 Inf]);
            o2_temp = permute(o2_temp,[1 2 4 3]);
        elseif strcmp(prod,'sjtu_gr')
            o2_temp = ncread([extractBefore(file_directory,'O2anom_1965to1969')...
                struct.files(t_idx+2).name],param_name); % extract o2 for year
        elseif strcmp(prod,'gt_ml')
            % average NN and RF
            o2_temp1 = ncread(file_directory{1},param_name,... % extract o2 for year
                [1 1 1 t_idx(1)],[Inf Inf Inf t_idx(end)-t_idx(1)+1]);
            o2_temp2 = ncread(file_directory{2},param_name,... % extract o2 for year
                [1 1 1 t_idx(1)],[Inf Inf Inf t_idx(end)-t_idx(1)+1]);
            o2_temp = mean(cat(5,o2_temp1,o2_temp2),5,'omitnan');
            clear o2_temp1 o2_temp2
        else
            o2_temp = ncread(file_directory,param_name,... % extract o2 for year
                [1 1 1 t_idx(1)],[Inf Inf Inf t_idx(end)-t_idx(1)+1]);
        end
        o2_temp = mean(o2_temp,4,'omitnan'); % average across year
        o2_temp = o2_temp(lon_idx,:,:); % sort o2 by longitude index
        o2_interp = interp3(lat,lon_temp,depth,... % interpolate to common grid
            o2_temp,lat_interp,lon_interp',depth_interp);
    else
        o2_interp = nan(length(lon_interp),length(lat_interp),...
            length(depth_interp));
    end
%     figure; pcolor(lon_temp,lat,o2_temp(:,:,6)'); shading flat; colorbar; title('Original');
%     figure; pcolor(lon_interp,lat_interp,o2_interp(:,:,11)'); shading flat; colorbar; title('Interpolated');
end

%% embedded function to interpolate to baseline and anomaly seasonal cycle
function [o2_interp_seas_baseline,o2_interp_seas_anom,o2_interp_seas_baseline_var] = ...
    interp_seasonal_cycle_to_common_grid(file_directory,param_name,years_temp,...
        yr1,yr2,yr_anom,lon,lat,depth,lon_interp,lat_interp,depth_interp,prod,struct)
    baseline_idx = find(years_temp >= yr1 & years_temp <= yr2);
    anom_idx = find(years_temp == yr_anom);
    lon_temp = convert_lon(lon,'format','-180:180');
    [lon_temp,lon_idx] = sort(lon_temp); % sort longitude in order
    if strcmp(prod,'gt_ml')
        % obtain baseline seasonal cycle for NN
        o2_temp_seas_baseline1 = ncread(file_directory{1},param_name,... % extract o2 for baseline years
            [1 1 1 baseline_idx(1)],[Inf Inf Inf baseline_idx(end)-baseline_idx(1)+1]);
        [o2_seas_baseline1,o2_seas_baseline1_var] = seasonal_cycle(o2_temp_seas_baseline1);
        clear o2_temp_seas_baseline1
        % obtain anomalous seasonal cycle for NN
        o2_seas_anom1 = ncread(file_directory{1},param_name,... % extract o2 for anomaly year
            [1 1 1 anom_idx(1)],[Inf Inf Inf anom_idx(end)-anom_idx(1)+1]);
        clear o2_temp_seas_anom1
        % obtain baseline seasonal cycle for RF
        o2_temp_seas_baseline2 = ncread(file_directory{2},param_name,... % extract o2 for baseline years
            [1 1 1 baseline_idx(1)],[Inf Inf Inf baseline_idx(end)-baseline_idx(1)+1]);
        [o2_seas_baseline2,o2_seas_baseline2_var] = seasonal_cycle(o2_temp_seas_baseline2);
        clear o2_temp_seas_baseline2
        % obtain anomalous seasonal cycle for RF
        o2_seas_anom2 = ncread(file_directory{2},param_name,... % extract o2 for anomaly year
            [1 1 1 anom_idx(1)],[Inf Inf Inf anom_idx(end)-anom_idx(1)+1]);
        % average NN and RF
        o2_seas_baseline = mean(cat(5,o2_seas_baseline1,o2_seas_baseline2),5,'omitnan');
        o2_seas_baseline_var = mean(cat(5,o2_seas_baseline1_var,o2_seas_baseline2_var),5,'omitnan');
        o2_seas_anom = mean(cat(5,o2_seas_anom1,o2_seas_anom2),5,'omitnan');
    else
        % obtain baseline seasonal cycle
        o2_temp_seas_baseline = ncread(file_directory,param_name,... % extract o2 for baseline years
            [1 1 1 baseline_idx(1)],[Inf Inf Inf baseline_idx(end)-baseline_idx(1)+1]);
        [o2_seas_baseline,o2_seas_baseline_var] = seasonal_cycle(o2_temp_seas_baseline);
        clear o2_temp_seas_baseline
        % obtain anomalous seasonal cycle
        o2_seas_anom = ncread(file_directory,param_name,... % extract o2 for anomaly year
            [1 1 1 anom_idx(1)],[Inf Inf Inf anom_idx(end)-anom_idx(1)+1]);
    end
    % interpolate each month to common grid
    o2_interp_seas_baseline = nan(length(lon_interp),length(lat_interp),...
        length(depth_interp),12);
    o2_interp_seas_anom = nan(length(lon_interp),length(lat_interp),...
        length(depth_interp),12);
    o2_interp_seas_baseline_var = nan(length(lon_interp),length(lat_interp),...
        length(depth_interp),12);
    for m = 1:12
        o2_interp_seas_baseline(:,:,:,m) = interp3(lat,lon_temp,depth,... % interpolate to common grid
            o2_seas_baseline(lon_idx,:,:,m),lat_interp,lon_interp',depth_interp);
        o2_interp_seas_baseline_var(:,:,:,m) = interp3(lat,lon_temp,depth,... % interpolate to common grid
            o2_seas_baseline_var(lon_idx,:,:,m),lat_interp,lon_interp',depth_interp);
        try
            o2_interp_seas_anom(:,:,:,m) = interp3(lat,lon_temp,depth,... % interpolate to common grid
                o2_seas_anom(lon_idx,:,:,m),lat_interp,lon_interp',depth_interp);
        catch
        end
    end
end

%% embedded function to calculate seasonal cycle
function [data_seas,data_seas_var] = seasonal_cycle(data)
    data = reshape(data,size(data,1),size(data,2),size(data,3),12,size(data,4)/12);
    data_seas = mean(data,5,'omitnan');
    data_seas_var = std(data,[],5,'omitnan');
end