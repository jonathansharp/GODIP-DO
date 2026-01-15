% this compares annual means of 

% file path
fpath = '/raid/sharp/matlab/O2_Intercomparison/O2_Maps/';

% taka's compilation
compilation.name = 'GODIP-DO_NCEI_67C0_Oct2025.nc';
products = ncread([fpath compilation.name],'prod');
lat = ncread([fpath compilation.name],'lat');
lon = ncread([fpath compilation.name],'lon');
time = ncread([fpath compilation.name],'time');
depth = ncread([fpath compilation.name],'depth');

% calculate long-term mean inventory
d1 = 0; d2 = 1000; % depth limits
[~,d1_idx] = min(abs(depth-d1));
[~,d2_idx] = min(abs(depth-d2));
o2_inv_mean_0_1000 = nan(length(lon),length(lat),length(products));
for p = 1:length(products)
    o2_temp = squeeze(ncread([fpath compilation.name],'o2anom',...
        [1 1 d1_idx 1 p],[Inf Inf d2_idx Inf 1]));
    o2_temp_mean = mean(o2_temp,4,'omitnan');
    dv = ncread([fpath compilation.name],'dv',...
        [1 1 d1_idx],[Inf Inf d2_idx]);
    o2_inv_mean_0_1000(:,:,p) = ...
        (sum(o2_temp_mean,3,'omitnan').*sum(dv,3,'omitnan'))./sum(dv,3,'omitnan');
    clear o2_temp o2_temp_mean dv
end

% plot mean inventories 
% for p = 1%:length(products)
%     figure; worldmap([-90 90],[-180 180]);
%     pcolorm(lat,lon,o2_inv_mean_0_1000(:,:,p)');
%     colorbar; clim([-10 10])
% end

% plot standard deviation among OI products
o2_inv_std_OI_0_1000 = std(o2_inv_mean_0_1000(:,:,1:5),[],3,'omitnan');
figure; worldmap([-90 90],[-180 180]);
pcolorm(lat,lon,o2_inv_std_OI_0_1000');
colorbar;

% plot standard deviation among ML products
o2_inv_std_ML_0_1000 = std(o2_inv_mean_0_1000(:,:,6:9),[],3,'omitnan');
figure; worldmap([-90 90],[-180 180]);
pcolorm(lat,lon,o2_inv_std_ML_0_1000');
colorbar;