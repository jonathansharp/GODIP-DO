% this compares annual means of 

%% file properties
file_date = char(datetime('now','Format','MMMuuuu')); % e.g., 'Nov2026'
fpath = '/fast4/o2/GODIP-DO/O2_Maps/';
fname = ['GODIP-DO_SotC-2025_' file_date '.nc'];
fname_dens = ['GODIP-DO_SotC-2025_' file_date '_density.nc'];
% fname_seas = ['GODIP-DO_NCEI_SEASONAL_JDS_' file_date '.nc'];
products = {'ncei' 'iap' 'gt_oi' 'rb' 'sjtu_gr' 'gobai' 'gobai_iap' ...
    'gt_ml' 'gt_nn_en4' 'gt_nn_iap' 'gt_rf_en4' 'gt_rf_iap' 'jingwei' 'han_zhou'};
labels = {'NCEI' 'IAP' 'GT-OI' 'RB' 'SJTU-GR' 'GOBAI' 'GOBAI-IAP' 'GT-ML' ...
    'GT-NN-EN4' 'GT-NN-IAP' 'GT-RF-EN4' 'GT-RF-IAP' 'Jingwei' 'HZ'};
products_dens = {'iap' 'en4'};

%% import ncei dimensions
datasets.ncei.name = 'w25_o_19652025_yearly_anomaly_NCEI.nc';
datasets.ncei.lat = ncread([fpath datasets.ncei.name],'lat');
datasets.ncei.lon = ncread([fpath datasets.ncei.name],'lon');
datasets.ncei.time = double(ncread([fpath datasets.ncei.name],'time'));
datasets.ncei.depth = ncread([fpath datasets.ncei.name],'depth');
datasets.ncei.date = datevec(datenum(1900,1,1)+datasets.ncei.time);
datasets.ncei.param_name = 'o_anom_an';

%% import iap dimensions
% datasets.iap.name = 'IAPmapping_applied_to_IAPdata_IAPmapping_34_1960to2025_4D.nc';
datasets.iap.name = 'IAP_oxygen_a_1970_2025.nc';
datasets.iap.lat = ncread([fpath datasets.iap.name],'lat');
datasets.iap.lon = ncread([fpath datasets.iap.name],'lon');
datasets.iap.time = double(ncread([fpath datasets.iap.name],'time'));
datasets.iap.depth = ncread([fpath datasets.iap.name],'depth_std');
datasets.iap.date = datevec(datetime(datasets.iap.time,7,1));
datasets.iap.param_name = 'oxygen';

%% import gt-oi dimensions
datasets.gt_oi.name = 'o2anom_g1x1z67_NCEI_Ito22_1965-2025.nc';
datasets.gt_oi.lat = ncread([fpath datasets.gt_oi.name],'lat');
datasets.gt_oi.lon = ncread([fpath datasets.gt_oi.name],'lon');
datasets.gt_oi.time = double(ncread([fpath datasets.gt_oi.name],'time'));
datasets.gt_oi.depth = ncread([fpath datasets.gt_oi.name],'depth');
datasets.gt_oi.date = datevec(datenum(1965,7,1)+datasets.gt_oi.time);
datasets.gt_oi.param_name = 'o2';

%% import rb dimensions
datasets.rb.name = 'o_annual_mean_depth_NCEI_Roach_Bindoff_1965_2025.nc';
datasets.rb.lat = ncread([fpath datasets.rb.name],'lat');
datasets.rb.lon = ncread([fpath datasets.rb.name],'lon');
datasets.rb.time = double(ncread([fpath datasets.rb.name],'year'));
datasets.rb.depth = ncread([fpath datasets.rb.name],'depth');
datasets.rb.date = datevec(datetime(datasets.rb.time,7,1));
datasets.rb.param_name = 'oxy_conc';

%% import sjtu-gr dimensions
datasets.sjtu_gr.name = '5yrmeanSJTU/O2anom_1965to1969_mean_1x1_NCEIdata0_2000m_UK_SJTU_Zhou250717.nc';
datasets.sjtu_gr.lat = ncread([fpath datasets.sjtu_gr.name],'lat');
datasets.sjtu_gr.lon = ncread([fpath datasets.sjtu_gr.name],'lon');
datasets.sjtu_gr.depth = ncread([fpath datasets.sjtu_gr.name],'depth');
datasets.sjtu_gr.time = (1967:2020)';
datasets.sjtu_gr.date = datevec(datenum(datasets.sjtu_gr.time,7,1));
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

%% import gobai-iap dimensions
datasets.gobai_iap.name = 'gobai-o2-IAP-godip-do-v1.2-with-float-correction.nc';
datasets.gobai_iap.lat = ncread([fpath datasets.gobai_iap.name],'lat');
datasets.gobai_iap.lon = ncread([fpath datasets.gobai_iap.name],'lon');
datasets.gobai_iap.time = ncread([fpath datasets.gobai_iap.name],'time');
datasets.gobai_iap.depth = ncread([fpath datasets.gobai_iap.name],'depth');
datasets.gobai_iap.date = datevec(datenum(1950,0,0)+datasets.gobai_iap.time);
datasets.gobai_iap.param_name = 'o2';

%% import gt-ml dimensions (OLD)
% datasets.gt_ml.name1 = 'GT_ML_NN_updated_01042026.nc';
% datasets.gt_ml.name2 = 'GT_ML_RF_updated_01042026.nc';
% datasets.gt_ml.lat = ncread([fpath datasets.gt_ml.name1],'lat');
% datasets.gt_ml.lon = ncread([fpath datasets.gt_ml.name1],'lon');
% datasets.gt_ml.time = double(ncread([fpath datasets.gt_ml.name1],'time'));
% datasets.gt_ml.depth = ncread([fpath datasets.gt_ml.name1],'depth');
% datasets.gt_ml.date = datevec(datenum(1965,0,15)+datasets.gt_ml.time);
% datasets.gt_ml.param_name = 'o2';

%% import gt-ml dimensions
datasets.gt_ml.name = 'o2_g1x1z102_yearly_NCEI_EN4C14_I24_1965-2025.nc';
datasets.gt_ml.lat = ncread([fpath datasets.gt_ml.name],'lat');
datasets.gt_ml.lon = ncread([fpath datasets.gt_ml.name],'lon');
datasets.gt_ml.year = double(ncread([fpath datasets.gt_ml.name],'year'));
datasets.gt_ml.depth = ncread([fpath datasets.gt_ml.name],'depth');
datasets.gt_ml.date = datevec(datetime(datasets.gt_ml.year,7,1));
datasets.gt_ml.param_name = 'o2';

%% import gt-ml (NN, EN4) dimensions
datasets.gt_nn_en4.name = 'o2_g1x1z102_NCEI_EN4C14_I24NN_yearly_1965-2025.nc';
datasets.gt_nn_en4.lat = ncread([fpath datasets.gt_nn_en4.name],'lat');
datasets.gt_nn_en4.lon = ncread([fpath datasets.gt_nn_en4.name],'lon');
datasets.gt_nn_en4.year = double(ncread([fpath datasets.gt_nn_en4.name],'year'));
datasets.gt_nn_en4.depth = ncread([fpath datasets.gt_nn_en4.name],'depth');
datasets.gt_nn_en4.date = datevec(datetime(datasets.gt_nn_en4.year,7,1));
datasets.gt_nn_en4.param_name = 'o2';

%% import gt-ml (NN, IAP) dimensions
datasets.gt_nn_iap.name = 'o2_g1x1z102_NCEI_IAP_I24NN_yearly_1965-2025.nc';
datasets.gt_nn_iap.lat = ncread([fpath datasets.gt_nn_iap.name],'lat');
datasets.gt_nn_iap.lon = ncread([fpath datasets.gt_nn_iap.name],'lon');
datasets.gt_nn_iap.year = double(ncread([fpath datasets.gt_nn_iap.name],'year'));
datasets.gt_nn_iap.depth = ncread([fpath datasets.gt_nn_iap.name],'depth');
datasets.gt_nn_iap.date = datevec(datetime(datasets.gt_nn_iap.year,7,1));
datasets.gt_nn_iap.param_name = 'o2';

%% import gt-ml (RF, EN4) dimensions
datasets.gt_rf_en4.name = 'o2_g1x1z102_NCEI_EN4C14_I24RF_yearly_1965-2025.nc';
datasets.gt_rf_en4.lat = ncread([fpath datasets.gt_rf_en4.name],'lat');
datasets.gt_rf_en4.lon = ncread([fpath datasets.gt_rf_en4.name],'lon');
datasets.gt_rf_en4.year = double(ncread([fpath datasets.gt_rf_en4.name],'year'));
datasets.gt_rf_en4.depth = ncread([fpath datasets.gt_rf_en4.name],'depth');
datasets.gt_rf_en4.date = datevec(datetime(datasets.gt_rf_en4.year,7,1));
datasets.gt_rf_en4.param_name = 'o2';

%% import gt-ml (RF, IAP) dimensions
datasets.gt_rf_iap.name = 'o2_g1x1z102_NCEI_IAP_I24RF_yearly_1965-2025.nc';
datasets.gt_rf_iap.lat = ncread([fpath datasets.gt_rf_iap.name],'lat');
datasets.gt_rf_iap.lon = ncread([fpath datasets.gt_rf_iap.name],'lon');
datasets.gt_rf_iap.year = double(ncread([fpath datasets.gt_rf_iap.name],'year'));
datasets.gt_rf_iap.depth = ncread([fpath datasets.gt_rf_iap.name],'depth');
datasets.gt_rf_iap.date = datevec(datetime(datasets.gt_rf_iap.year,7,1));
datasets.gt_rf_iap.param_name = 'o2';

%% import jingwei dimensions
datasets.jingwei.name = 'o2map_v1.2_WOD25_ShipArgo_with_correction_5500m_sjtu_jingwei.nc';
datasets.jingwei.lat = ncread([fpath datasets.jingwei.name],'lat');
datasets.jingwei.lon = ncread([fpath datasets.jingwei.name],'lon');
datasets.jingwei.time = double(ncread([fpath datasets.jingwei.name],'time'));
datasets.jingwei.depth = ncread([fpath datasets.jingwei.name],'depth');
datasets.jingwei.date = datevec(datenum(1965,0,15)+datasets.jingwei.time);
datasets.jingwei.param_name = 'o2';

%% import han and zhou dimensions
datasets.han_zhou.name = 'DO_DOMIP_SJTU_Han_Zhou_correctedArgo_19652025.nc';
datasets.han_zhou.lat = double(ncread([fpath datasets.han_zhou.name],'Latitude'));
datasets.han_zhou.lon = double(ncread([fpath datasets.han_zhou.name],'Longitude'));
datasets.han_zhou.time = double(ncread([fpath datasets.han_zhou.name],'Time'));
datasets.han_zhou.depth = double(ncread([fpath datasets.han_zhou.name],'Depth'));
datasets.han_zhou.date = datevec(datenum([repmat(1965,length(datasets.han_zhou.time),1) ...
    1+datasets.han_zhou.time repmat(15,length(datasets.han_zhou.time),1)]));
datasets.han_zhou.param_name = 'DO_Reconstruction';

%% load NCEI mapped climatology
ncei_mean.name = 'DOMIP1_NCEI_19652022_LTM_annual_apr2025.nc';
ncei_mean.lat = ncread([fpath ncei_mean.name],'lat');
ncei_mean.lon = ncread([fpath ncei_mean.name],'lon');
% sort longitude
ncei_mean.lon = convert_lon(ncei_mean.lon,'format','-180:180');
[ncei_mean.lon,lon_idx] = sort(ncei_mean.lon);
ncei_mean.depth = ncread([fpath ncei_mean.name],'depth');
ncei_mean.o2 = mean(ncread([fpath ncei_mean.name],'o2'),4,'omitnan');
% ncei_mean.o2 = ncei_mean.o2(:,:,d1_idx:d2_idx);
ncei_mean.o2 = ncei_mean.o2(lon_idx,:,:);

%% define common grid
d = nan(105,4);
for p = 1:length(products)
    d(1:length(datasets.(products{p}).depth),p) = ...
        datasets.(products{p}).depth;
end
common.lon = (-179.5:1:179.5)';
common.lat = (-89.5:1:89.5)';
common.depth = datasets.gt_nn_en4.depth;
common.time = datenum(1965:2025,7,1)';

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
% nccreate([fpath fname],'months','Dimensions',{'months' 12});
% ncwrite([fpath fname],'months',(1:12)');
nccreate([fpath fname],'products','Dimensions',{'products' length(products)},'Datatype','string');
ncwrite([fpath fname],'products',products);
disp('Common NetCDF file created.');

%% create nc file for interpolated density
if exist([fpath fname_dens],'file'); delete([fpath fname_dens]); end
nccreate([fpath fname_dens],'lon','Dimensions',{'lon' length(common.lon)},'Format','netcdf4');
ncwrite([fpath fname_dens],'lon',common.lon);
nccreate([fpath fname_dens],'lat','Dimensions',{'lat' length(common.lat)});
ncwrite([fpath fname_dens],'lat',common.lat);
nccreate([fpath fname_dens],'depth','Dimensions',{'depth' length(common.depth)});
ncwrite([fpath fname_dens],'depth',common.depth);
nccreate([fpath fname_dens],'time','Dimensions',{'time' length(common.time)});
ncwrite([fpath fname_dens],'time',common.time);
% nccreate([fpath fname_dens],'months','Dimensions',{'months' 12});
% ncwrite([fpath fname_dens],'months',(1:12)');
% nccreate([fpath fname_dens],'products','Dimensions',{'products' length(products_dens)},'Datatype','string');
% ncwrite([fpath fname_dens],'products',products_dens);
for v = {'temperature_en4','salinity_en4','density_en4', ...
        'temperature_iap','salinity_iap','density_iap'}
    nccreate([fpath fname_dens],v{1}, ...
        'Dimensions',{'lon' 'lat' 'depth' 'time'}, ...
        'Datatype','single');
end
disp('Common NetCDF file for density created.');

%% pre-allocate mask and years
common.mask = true(length(common.lon),length(common.lat),...
    length(common.depth),length(common.time));
years = datevec(common.time); years = years(:,1);

%% interpolate annual O2 means to common grid
% loop through and interpolate each O2 product
for p = 1:length(products)
    % pre-allocate
    common.(products{p}) = nan(length(common.lon),length(common.lat),...
        length(common.depth),length(common.time));
    tic % start timer
    years_temp = datasets.(products{p}).date(:,1); % extract years
    for t = 1:length(common.time) % loop through years
        if strcmp(products{p},'gt_ml_old')
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
    % add woa mean to o2 anomalies
    if p <= 5; common.(products{p}) = common.(products{p}) + ncei_mean.o2; end
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

%% interpolate annual density means to common grid
en4_paths = path_config('EN4','o2');
en4_ref = fullfile(en4_paths.temp_path,'/EN.4.2.2.f.analysis.c14.200001.nc');
en4_lon = ncread(en4_ref,'lon');
en4_lat = ncread(en4_ref,'lat');
en4_depth = ncread(en4_ref,'depth');

% Sort EN4 longitudes and retain the ordering for the data arrays.
en4_lon = convert_lon(en4_lon,'format','-180:180');
[en4_lon,lon_idx] = sort(en4_lon);

% Native EN4 grid and pressure for the GSW calculations.
nlon_en4 = length(en4_lon);
nlat_en4 = length(en4_lat);
nz_en4 = length(en4_depth);
[lon_en4_2d,lat_en4_2d] = ndgrid(en4_lon,en4_lat);
lon_en4_3d = repmat(lon_en4_2d,1,1,nz_en4);
lat_en4_3d = repmat(lat_en4_2d,1,1,nz_en4);
depth_en4_3d = repmat(reshape(en4_depth,1,1,[]), ...
    nlon_en4,nlat_en4,1);
pressure_en4 = gsw_p_from_z(-depth_en4_3d,lat_en4_3d);

% Process each year.
for t = 1:length(common.time)
    y = years(t);
    monthly_files = cell(12,1);
    % create list of monthly files for the year
    for m = 1:12
        monthly_files{m} = fullfile(en4_paths.temp_path,...
        sprintf('EN.4.2.2.f.analysis.c14.%04d%02d.nc',y,m));
    end
    % add warning if one or more monthly files are missing
    if ~all(cellfun(@isfile,monthly_files))
        warning('Skipping EN4 %d: one or more monthly files are missing.',y);
        continue
    end
    % pre-allocate sums and counts for averaging
    temp_sum = zeros(nlon_en4,nlat_en4,nz_en4,'single');
    sal_sum = zeros(nlon_en4,nlat_en4,nz_en4,'single');
    rho_sum = zeros(nlon_en4,nlat_en4,nz_en4,'single');
    count = zeros(nlon_en4,nlat_en4,nz_en4,'single');
    % average monthly fields on the native EN4 grid
    for m = 1:12
        % read temperature and salinity from the monthly file
        temp_src = single(ncread(monthly_files{m},'temperature')) - 273.15; % convert from kelvin to celsius
        sal_src = single(ncread(monthly_files{m},'salinity'));
        temp_src = temp_src(lon_idx,:,:);
        sal_src = sal_src(lon_idx,:,:);
        % calculate absolute salinity, conservative temperature, and density using GSW
        sa = gsw_SA_from_SP(sal_src,pressure_en4, ...
            lon_en4_3d,lat_en4_3d);
        ct = gsw_CT_from_t(sa,temp_src,pressure_en4);
        rho_month = single(gsw_rho(sa,ct,pressure_en4));
        % accumulate sum for calculation ofannual mean, only where all three fields are finite
        valid = isfinite(temp_src) & isfinite(sal_src) & isfinite(rho_month);
        temp_sum(valid) = temp_sum(valid) + temp_src(valid);
        sal_sum(valid) = sal_sum(valid) + sal_src(valid);
        rho_sum(valid) = rho_sum(valid) + rho_month(valid);
        count(valid) = count(valid) + 1;
    end
    % pre-allocate annual mean arrays
    annual_temp_native = nan(nlon_en4,nlat_en4,nz_en4,'single');
    annual_sal_native = nan(nlon_en4,nlat_en4,nz_en4,'single');
    annual_rho_native = nan(nlon_en4,nlat_en4,nz_en4,'single');
    valid = count > 0;
    % calculate annual mean on native grid only where there is at least one valid month
    annual_temp_native(valid) = temp_sum(valid)./count(valid);
    annual_sal_native(valid) = sal_sum(valid)./count(valid);
    annual_rho_native(valid) = rho_sum(valid)./count(valid);
    % interpolate each annual mean to the common grid
    annual_temp = single(interp3(en4_lat,en4_lon,en4_depth, ...
        annual_temp_native,common.lat,common.lon',common.depth));
    annual_sal = single(interp3(en4_lat,en4_lon,en4_depth, ...
        annual_sal_native,common.lat,common.lon',common.depth));
    annual_rho = single(interp3(en4_lat,en4_lon,en4_depth, ...
        annual_rho_native,common.lat,common.lon',common.depth));
    % write common grid annual means to the density NetCDF file
    ncwrite([fpath fname_dens],'temperature_en4',annual_temp,[1 1 1 t]);
    ncwrite([fpath fname_dens],'salinity_en4',annual_sal,[1 1 1 t]);
    ncwrite([fpath fname_dens],'density_en4',annual_rho,[1 1 1 t]);
    disp(['EN4 interpolated for ' num2str(y)]);
end

%% interpolate annual IAP temperature, salinity, and density means
iap_paths = path_config('IAP','o2');
iap_temp_ref = fullfile(iap_paths.temp_path, ...
    'IAPv4_Temp_monthly_1_6000m_year_2025_month_09.nc');
iap_sal_ref = fullfile(iap_paths.sal_path, ...
    'IAPv2_Salinity_monthly_1_6000m_year_2025_month_09.nc');

iap_lon = ncread(iap_temp_ref,'lon');
iap_lat = ncread(iap_temp_ref,'lat');
iap_depth = ncread(iap_temp_ref,'depth_std');
iap_sal_lon = ncread(iap_sal_ref,'lon');
if ~isequal(iap_lat,ncread(iap_sal_ref,'lat')) || ...
        ~isequal(iap_depth,ncread(iap_sal_ref,'depth_std'))
    error('IAP temperature and salinity latitude/depth grids do not match.');
end

% Sort IAP longitudes and retain the ordering for each data array.
iap_lon = convert_lon(iap_lon,'format','-180:180');
[iap_lon,lon_idx_iap] = sort(iap_lon);
iap_sal_lon = convert_lon(iap_sal_lon,'format','-180:180');
[iap_sal_lon,lon_idx_sal] = sort(iap_sal_lon);
iap_sal_lon_ext = [iap_sal_lon(end)-360; iap_sal_lon; iap_sal_lon(1)+360];

% Native IAP grid and pressure for the GSW calculations.
nlon_iap = length(iap_lon);
nlat_iap = length(iap_lat);
nz_iap = length(iap_depth);
[lon_iap_2d,lat_iap_2d] = ndgrid(iap_lon,iap_lat);
lon_iap_3d = repmat(lon_iap_2d,1,1,nz_iap);
lat_iap_3d = repmat(lat_iap_2d,1,1,nz_iap);
depth_iap_3d = repmat(reshape(iap_depth,1,1,[]), ...
    nlon_iap,nlat_iap,1);
pressure_iap = gsw_p_from_z(-depth_iap_3d,lat_iap_3d);

% Process each year.
for t = 1:length(common.time)
    y = years(t);
    monthly_temp_files = cell(12,1);
    monthly_sal_files = cell(12,1);
    for m = 1:12
        monthly_temp_files{m} = fullfile(iap_paths.temp_path,sprintf( ...
            'IAPv4_Temp_monthly_1_6000m_year_%04d_month_%02d.nc',y,m));
        monthly_sal_files{m} = fullfile(iap_paths.sal_path,sprintf( ...
            'IAPv2_Salinity_monthly_1_6000m_year_%04d_month_%02d.nc',y,m));
    end
    if ~all(cellfun(@isfile,monthly_temp_files)) || ...
            ~all(cellfun(@isfile,monthly_sal_files))
        warning('Skipping IAP %d: one or more monthly files are missing.',y);
        continue
    end

    temp_sum = zeros(nlon_iap,nlat_iap,nz_iap,'single');
    sal_sum = zeros(nlon_iap,nlat_iap,nz_iap,'single');
    rho_sum = zeros(nlon_iap,nlat_iap,nz_iap,'single');
    count = zeros(nlon_iap,nlat_iap,nz_iap,'single');

    % Average monthly fields on the native IAP grid.
    for m = 1:12
        % IAP fields are stored depth x lon x lat; convert to lon x lat x depth.
        temp_src = permute(single(ncread(monthly_temp_files{m},'temp')),[2 3 1]);
        fill_warning = warning('off', ...
            'MATLAB:imagesci:netcdf:fillValueTypeMismatch');
        sal_src = permute(single(ncread(monthly_sal_files{m},'salinity')),[2 3 1]);
        warning(fill_warning);
        temp_src(temp_src == 999) = NaN;
        sal_src(sal_src == 999) = NaN;
        temp_src = temp_src(lon_idx_iap,:,:);
        sal_src = sal_src(lon_idx_sal,:,:);
        % Salinity longitudes are offset by 0.5 degrees; remap to the
        % temperature grid, adding cyclic points across the longitude seam.
        sal_src_ext = cat(1,sal_src(end,:,:),sal_src,sal_src(1,:,:));
        sal_src = interp3(iap_lat,iap_sal_lon_ext,iap_depth,sal_src_ext, ...
            iap_lat,iap_lon',iap_depth);

        % IAP salinity is already Absolute Salinity (g/kg).
        ct = gsw_CT_from_t(sal_src,temp_src,pressure_iap);
        rho_month = single(gsw_rho(sal_src,ct,pressure_iap));

        valid = isfinite(temp_src) & isfinite(sal_src) & isfinite(rho_month);
        temp_sum(valid) = temp_sum(valid) + temp_src(valid);
        sal_sum(valid) = sal_sum(valid) + sal_src(valid);
        rho_sum(valid) = rho_sum(valid) + rho_month(valid);
        count(valid) = count(valid) + 1;
    end

    annual_temp_native = nan(nlon_iap,nlat_iap,nz_iap,'single');
    annual_sal_native = nan(nlon_iap,nlat_iap,nz_iap,'single');
    annual_rho_native = nan(nlon_iap,nlat_iap,nz_iap,'single');
    valid = count > 0;
    annual_temp_native(valid) = temp_sum(valid)./count(valid);
    annual_sal_native(valid) = sal_sum(valid)./count(valid);
    annual_rho_native(valid) = rho_sum(valid)./count(valid);

    % Interpolate each annual mean to the common grid.
    annual_temp = single(interp3(iap_lat,iap_lon,iap_depth, ...
        annual_temp_native,common.lat,common.lon',common.depth));
    annual_sal = single(interp3(iap_lat,iap_lon,iap_depth, ...
        annual_sal_native,common.lat,common.lon',common.depth));
    annual_rho = single(interp3(iap_lat,iap_lon,iap_depth, ...
        annual_rho_native,common.lat,common.lon',common.depth));

    ncwrite([fpath fname_dens],'temperature_iap',annual_temp,[1 1 1 t]);
    ncwrite([fpath fname_dens],'salinity_iap',annual_sal,[1 1 1 t]);
    ncwrite([fpath fname_dens],'density_iap',annual_rho,[1 1 1 t]);
    disp(['IAP interpolated for ' num2str(y)]);
end

%% test plotting of density fields
ncfile = [fpath fname_dens];
plot_year = 1994;
target_depth = 500;  % metres; plots the nearest common-grid level
lon = ncread(ncfile,'lon');
lat = ncread(ncfile,'lat');
depth = ncread(ncfile,'depth');
time_parts = datevec(ncread(ncfile,'time'));
t_idx = find(time_parts(:,1) == plot_year,1);
if isempty(t_idx)
    error('Year %d is not in the output time coordinate.',plot_year);
end
[~,z_idx] = min(abs(depth - target_depth));
fprintf('Plotting %d m (nearest to requested %d m)\n', ...
    depth(z_idx),target_depth);
vars = {'temperature_en4','salinity_en4','density_en4';...
    'temperature_iap','salinity_iap','density_iap'};
titles = {'EN4, Temperature (deg C)','EN4, Salinity','EN4, Density (kg m^{-3})';...
    'IAP, Temperature (deg C)','IAP, Salinity','IAP, Density (kg m^{-3})'};
zlim = {[0 30], [33 40], [1025 1031];...
    [0 30], [33 40], [1025 1031]};
for p = 1:size(vars,1)
    figure('Name',sprintf('Density Fields for %d m in %d',depth(z_idx),plot_year), ...
        'NumberTitle','off','Position',[100 100 1600 400]);
    for v = 1:size(vars,2)
        field = squeeze(ncread(ncfile,vars{p,v}, ...
            [1 1 z_idx t_idx],[Inf Inf 1 1]));
        subplot(1,3,v);
        pcolor(lon,lat,field');
        shading flat;
        axis tight;
        xlabel('Longitude');
        ylabel('Latitude');
        title(sprintf('%s, %d, %g m',titles{p,v},plot_year,depth(z_idx)));
        colorbar;
        caxis(zlim{p,v});
        vals = field(isfinite(field));
        if isempty(vals)
            warning('%s has no finite values for this year and depth.',vars{p,v});
        else
            fprintf('%s range: %.3f to %.3f\n',vars{p,v},min(vals),max(vals));
        end
    end
    exportgraphics(gcf,sprintf('Figures/%s/density_plot_%s_%d_%dm.png',file_date,vars{p,v},plot_year,depth(z_idx)));
    close;
end

%% interpolate to seasonal cycle for ML products
% % pre-allocate
% common.mask = true(length(common.lon),length(common.lat),...
%     length(common.depth),length(common.time));
% % loop through and interpolate each product
% for p = [6 8]
%     % pre-allocate
%     common.([products{p} '_baseline_seas']) = nan(length(common.lon),length(common.lat),...
%         length(common.depth),length(common.time));
%     common.([products{p} '_anom_seas']) = nan(length(common.lon),length(common.lat),...
%         length(common.depth),length(common.time));
%     tic % start timer
%     years_temp = datasets.(products{p}).date(:,1); % extract years
%     y1 = 1990; y2 = 2020; y_anom = 2025;
%     if strcmp(products{p},'gt_ml_old')
%     [common.([products{p} '_baseline_seas']),common.([products{p} '_anom_seas']),...
%         common.([products{p} '_baseline_seas_var'])] = ...
%         interp_seasonal_cycle_to_common_grid({[fpath datasets.(products{p}).name1];...
%         [fpath datasets.(products{p}).name2]},...
%         datasets.(products{p}).param_name,years_temp,y1,y2,y_anom,...
%         datasets.(products{p}).lon,datasets.(products{p}).lat,...
%         datasets.(products{p}).depth,common.lon,common.lat,common.depth,...
%         products{p},datasets.(products{p}));
%     else
%     [common.([products{p} '_baseline_seas']),common.([products{p} '_anom_seas']),...
%         common.([products{p} '_baseline_seas_var'])] = ...
%         interp_seasonal_cycle_to_common_grid([fpath datasets.(products{p}).name],...
%         datasets.(products{p}).param_name,years_temp,y1,y2,y_anom,...
%         datasets.(products{p}).lon,datasets.(products{p}).lat,...
%         datasets.(products{p}).depth,common.lon,common.lat,common.depth,...
%         products{p},datasets.(products{p}));
%     end
%     % if start or end longitude is all nan, replace with mean of
%     % bounding longitudes
%     for m = 1:12
%         if all(all(isnan(common.([products{p} '_baseline_seas'])(1,:,:,m))))
%             common.([products{p} '_baseline_seas'])(1,:,:,m) = ...
%                 mean([common.([products{p} '_baseline_seas'])(end,:,:,m);...
%                 common.([products{p} '_baseline_seas'])(2,:,:,m)],1);
%         elseif all(all(isnan(common.([products{p} '_baseline_seas'])(end,:,:,m))))
%             common.([products{p} '_baseline_seas'])(end,:,:,m) = ...
%                 mean([common.([products{p} '_baseline_seas'])(end-1,:,:,m);...
%                 common.([products{p} '_baseline_seas'])(1,:,:,m)],m);
%         end
%         if all(all(isnan(common.([products{p} '_anom_seas'])(1,:,:,m))))
%             common.([products{p} '_anom_seas'])(1,:,:,m) = ...
%                 mean([common.([products{p} '_anom_seas'])(end,:,:,m);...
%                 common.([products{p} '_anom_seas'])(2,:,:,m)],1);
%         elseif all(all(isnan(common.([products{p} '_anom_seas'])(end,:,:,m))))
%             common.([products{p} '_anom_seas'])(end,:,:,m) = ...
%                 mean([common.([products{p} '_anom_seas'])(end-1,:,:,m);...
%                 common.([products{p} '_anom_seas'])(1,:,:,m)],1);
%         end
%     end
%     disp([products{p} ' seasonal cycle interpolated']);
%     % save
%     nccreate([fpath fname],[products{p} '_baseline_seas'],...
%         'Dimensions',{'lon' 'lat' 'depth' 'months'});
%     ncwrite([fpath fname],[products{p} '_baseline_seas'],common.([products{p} '_baseline_seas']));
%     nccreate([fpath fname],[products{p} '_anom_seas'],...
%         'Dimensions',{'lon' 'lat' 'depth' 'months'});
%     ncwrite([fpath fname],[products{p} '_anom_seas'],common.([products{p} '_anom_seas']));
%     nccreate([fpath fname],[products{p} '_baseline_seas_var'],...
%         'Dimensions',{'lon' 'lat' 'depth' 'months'});
%     ncwrite([fpath fname],[products{p} '_baseline_seas_var'],common.([products{p} '_baseline_seas_var']));
%     % clean up
%     common = rmfield(common,[products{p} '_baseline_seas']);
%     common = rmfield(common,[products{p} '_anom_seas']);
%     toc
% end

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
        elseif strcmp(prod,'gt_ml_old')
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
        for d = 5:-1:1 % upwards in depth from 20
            if sum(sum(double(~isnan(o2_interp(:,:,d))))) == 0 % if all are NaN
                o2_interp(:,:,d) = o2_interp(:,:,d+1); % fill with lower layer
            end
        end
    else
        o2_interp = nan(length(lon_interp),length(lat_interp),...
            length(depth_interp));
    end
    % figure; pcolor(lon_temp,lat,o2_temp(:,:,6)'); shading flat; colorbar; title('Original');
    % figure; pcolor(lon_interp,lat_interp,o2_interp(:,:,11)'); shading flat; colorbar; title('Interpolated');
end

%% embedded function to interpolate to baseline and anomaly seasonal cycle
function [o2_interp_seas_baseline,o2_interp_seas_anom,o2_interp_seas_baseline_var] = ...
    interp_seasonal_cycle_to_common_grid(file_directory,param_name,years_temp,...
        yr1,yr2,yr_anom,lon,lat,depth,lon_interp,lat_interp,depth_interp,prod,struct)
    baseline_idx = find(years_temp >= yr1 & years_temp <= yr2);
    anom_idx = find(years_temp == yr_anom);
    lon_temp = convert_lon(lon,'format','-180:180');
    [lon_temp,lon_idx] = sort(lon_temp); % sort longitude in order
    if strcmp(prod,'gt_ml_old')
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