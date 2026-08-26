% this compares annual means of 
close all;

%% file and plot properties
file_date = 'Feb2026';
% fpath = '/fast4/o2/GODIP-DO/';
% fpath_EN4 = '/med2/';
fpath = '/raid/sharp/matlab/GODIP-DO/';
fpath_EN4 = '/raid/sharp/matlab/GODIP-DO/';
compilation.name = ['GODIP-DO_NCEI_JDS_' file_date '.nc'];
% products = {'ncei' 'iap' 'gt_oi' 'rb' 'sjtu_gr' 'gobai' 'gt_ml' 'jingwei' 'han_zhou'};
products = ncread([fpath 'O2_Maps/' compilation.name],'products');
labels = {'NCEI' 'IAP' 'GT-OI' 'UTAS' 'SJTU-GR' 'GOBAI' 'GT-ML' 'Jingwei' 'SJTU-HZ'};
lat_lims = [-75 75]; lon_lims = [20 380];
o2_lims = [0 400]; o2_levels = o2_lims(1):12.5:o2_lims(end);
anom_lims = [-40 40]; anom_levels = anom_lims(1):5:anom_lims(end);
trend_lims = [-6 6]; trend_levels = trend_lims(1):0.5:trend_lims(end);
temp_anom_lims = [-3.5 3.5]; temp_anom_levels = temp_anom_lims(1):0.25:temp_anom_lims(end);
% inventories
layers = {'50-1000' '0-1800' '0-50' '100-600'};
lyr_top = [50 0 0 100]; lyr_bot = [1000 1800 50 600];
inv_anom_lims = [-50 50]; inv_anom_levels = inv_anom_lims(1):10:inv_anom_lims(end);
o2_inv_lims = [0,300;0,600;10,20;0,200];
o2_inv_levels = {o2_inv_lims(1,1):20:o2_inv_lims(1,end);...
    o2_inv_lims(2,1):40:o2_inv_lims(2,end);...
    o2_inv_lims(3,1):0.5:o2_inv_lims(3,end);...
    o2_inv_lims(4,1):10:o2_inv_lims(4,end)};
inv_trend_lims = [-3,3;-5,5;-0.2,0.2;-3,3];
inv_trend_levels = {inv_trend_lims(1,1):0.25:inv_trend_lims(1,end);...
    inv_trend_lims(2,1):0.5:inv_trend_lims(2,end);...
    inv_trend_lims(3,1):0.02:inv_trend_lims(3,end);...
    inv_trend_lims(4,1):0.25:inv_trend_lims(4,end)};
inv_anom_lims = [-15,15;-20,20;-1,1;-10,10];
inv_anom_levels = {inv_anom_lims(1,1):1.5:inv_anom_lims(1,end);...
    inv_anom_lims(2,1):4:inv_anom_lims(2,end);...
    inv_anom_lims(3,1):0.125:inv_anom_lims(3,end);...
    inv_anom_lims(4,1):2:inv_anom_lims(4,end)};

%% import compilation dimensions
lat = ncread([fpath 'O2_Maps/' compilation.name],'lat');
lon = ncread([fpath 'O2_Maps/' compilation.name],'lon');
time = ncread([fpath 'O2_Maps/' compilation.name],'time');
depth = ncread([fpath 'O2_Maps/' compilation.name],'depth');
mask = logical(ncread([fpath 'O2_Maps/' compilation.name],'mask'));

%% depth limits for inventory calculations
d1 = 0; d2 = 1800;
[~,d1_idx] = min(abs(depth-d1));
[~,d2_idx] = min(abs(depth-d2));
depth = depth(d1_idx:d2_idx);
depths_to_plot = [10 250 1000];
depth_indices = find(any(depth == depths_to_plot,2));
date = datevec(time); year = date(:,1);
% time indices
anom_idx_2025 = find(year == 2025);
anom_idx_2024 = find(year == 2024);
baseline_idx = find(year >= 1990 & year <= 2020);

%% load NCEI mapped climatology
ncei_mean.name = '/DOMIP1_NCEI_19652022_LTM_annual_apr2025.nc';
ncei_mean.lat = ncread([fpath 'O2_Maps/' ncei_mean.name],'lat');
ncei_mean.lon = ncread([fpath 'O2_Maps/' ncei_mean.name],'lon');
% sort longitude
ncei_mean.lon = convert_lon(ncei_mean.lon,'format','-180:180');
[ncei_mean.lon,lon_idx] = sort(ncei_mean.lon);
ncei_mean.depth = ncread([fpath 'O2_Maps/' ncei_mean.name],'depth');
ncei_mean.o2 = mean(ncread([fpath 'O2_Maps/' ncei_mean.name],'o2'),4,'omitnan');
ncei_mean.o2 = ncei_mean.o2(:,:,d1_idx:d2_idx);
ncei_mean.o2 = ncei_mean.o2(lon_idx,:,:);

%% calculate o2 inventories, trends, anomalies and other statistics
% define volume and height
[dv,da,dh] = weights3d(lon,lat,depth);
mask = mask(:,:,d1_idx:d2_idx,:);
mask = sum(mask,4) > 10; % set constant mask over time
% * weird randon NaNs in Bay of Bengal for one product? It's GT-OI
% * fixed (I think) by just setting "constant" mask to > 50 years

if ~exist([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_STATS_' file_date '.mat'],'file')

% calculate statistics
o2_mean = nan(length(lon),length(lat),d2_idx-d1_idx+1,length(products));
o2_trend = nan(length(lon),length(lat),d2_idx-d1_idx+1,length(products));
o2_inv_mean = nan(length(lon),length(lat),length(products),length(layers));
o2_avg_mean = nan(length(lon),length(lat),length(products),length(layers));
o2_inv_trend = nan(length(lon),length(lat),length(products),length(layers));
o2_avg_trend = nan(length(lon),length(lat),length(products),length(layers));
o2_inv_ts = nan(length(time),length(products),length(layers));
o2_inv_depth_ts = nan(length(time),length(products),length(depth));
o2_inv_anom_2025 = nan(length(lon),length(lat),length(layers),length(products));
o2_inv_anom_2024 = nan(length(lon),length(lat),length(layers),length(products));
o2_inv_diff_2025_2024 = nan(length(lon),length(lat),length(layers),length(products));
o2_avg_anom_2025 = nan(length(lon),length(lat),length(layers),length(products));
o2_avg_anom_2024 = nan(length(lon),length(lat),length(layers),length(products));
o2_avg_diff_2025_2024 = nan(length(lon),length(lat),length(layers),length(products));

% loop through products
for p = 1:length(products)
    % o2 maps over time in umol/kg
    o2 = squeeze(ncread([fpath 'O2_Maps/' compilation.name],products{p},...
        [1 1 d1_idx 1],[Inf Inf d2_idx-d1_idx+1 Inf]));
    if p <= 5 % add woa mean to o2 anomalies
        o2 = o2 + ncei_mean.o2;
    end
    % zero trap
    o2(o2<0) = 0;
    % convert to mol/m3 (approximate seawater density)
    o2_mean(:,:,:,p) = mean(o2,4,'omitnan'); % mean for each grid cell
    o2_m3 = o2.*1026.*(1e-6); % umol/kg * kg/m3 * mol/umol = mol/m3
    % remove values outside common mask
    o2(~repmat(mask,1,1,1,length(time))) = NaN;
    o2_m3(~repmat(mask,1,1,1,length(time))) = NaN;
    dv(~mask) = NaN; dh(~mask) = NaN; da(~mask) = NaN;
    % calculate long-term mean inventory
    for l = 1:length(layers)
        idx_top = find(depth==lyr_top(l));
        idx_bot = find(depth==lyr_bot(l));
        % mask that covers entire depth interval
        idx_full_mask = sum(mask(:,:,idx_top:idx_bot),3) == idx_bot-idx_top+1;
        % calculate inventory
        o2_inv_temp = squeeze(trapz(depth(idx_top:idx_bot),...
            o2_m3(:,:,idx_top:idx_bot,:),3)); % mol/m2
        o2_inv_temp(~repmat(idx_full_mask,1,1,length(time))) = NaN;
        % o2_inv_temp2 = squeeze(sum(o2_m3(:,:,idx_top:idx_bot,:).*...
        %     dh(:,:,idx_top:idx_bot),3,'omitnan')); % umol/kg
        % o2_inv_temp2(~idx_full_mask) = NaN;
        % o2_inv_temp2(o2_inv_temp2==0) = NaN;
        % calculate average concentration
        o2_avg_temp = squeeze(sum(o2(:,:,idx_top:idx_bot,:).*...
            dh(:,:,idx_top:idx_bot,:),3,'omitnan')./...
            sum(dh(:,:,idx_top:idx_bot,:),3,'omitnan')); % umol/kg
        o2_avg_temp(~repmat(idx_full_mask,1,1,length(time))) = NaN;
        % calculate mean inventory
        o2_inv_mean(:,:,p,l) = mean(o2_inv_temp,3,'omitnan'); % mol/m2
        % calculate mean average concentration
        o2_avg_mean(:,:,p,l) = mean(o2_avg_temp,3,'omitnan'); % umol/kg
        % 2025 inv. anomaly
        o2_inv_anom_temp = o2_inv_temp(:,:,anom_idx_2025) - mean(o2_inv_temp(:,:,...
            baseline_idx(1):baseline_idx(end)),3,'omitnan');
        o2_inv_anom_temp(~idx_full_mask) = NaN;
        o2_inv_anom_2025(:,:,l,p) = o2_inv_anom_temp;
        % 2024 inv. anomaly
        o2_inv_anom_temp = o2_inv_temp(:,:,anom_idx_2024) - mean(o2_inv_temp(:,:,...
            baseline_idx(1):baseline_idx(end)),3,'omitnan');
        o2_inv_anom_temp(~idx_full_mask) = NaN;
        o2_inv_anom_2024(:,:,l,p) = o2_inv_anom_temp;
        % 2025 - 2024 inv. difference
        o2_inv_diff_temp = o2_inv_temp(:,:,anom_idx_2025) - o2_inv_temp(:,:,anom_idx_2024);
        o2_inv_diff_temp(~idx_full_mask) = NaN;
        o2_inv_diff_2025_2024(:,:,l,p) = o2_inv_diff_temp;
        % 2025 average content anomaly
        o2_avg_anom_temp = o2_avg_temp(:,:,anom_idx_2025) - mean(o2_avg_temp(:,:,...
            baseline_idx(1):baseline_idx(end)),3,'omitnan');
        o2_avg_anom_temp(~idx_full_mask) = NaN;
        o2_avg_anom_2025(:,:,l,p) = o2_avg_anom_temp;
        % 2024 average content anomaly
        o2_avg_anom_temp = o2_avg_temp(:,:,anom_idx_2024) - mean(o2_avg_temp(:,:,...
            baseline_idx(1):baseline_idx(end)),3,'omitnan');
        o2_avg_anom_temp(~idx_full_mask) = NaN;
        o2_avg_anom_2024(:,:,l,p) = o2_avg_anom_temp;
        % 2025 - 2024 average content difference
        o2_avg_diff_temp = o2_avg_temp(:,:,anom_idx_2025) - o2_avg_temp(:,:,anom_idx_2024);
        o2_avg_diff_temp(~idx_full_mask) = NaN;
        o2_avg_diff_2025_2024(:,:,l,p) = o2_avg_diff_temp;
        % decadal inventory trend
        for x = 1:length(lon)
            for y =1:length(lat)
                if sum(~isnan(o2_inv_temp(x,y,:))) > length(time)/2 % if more than half the grid cells are filled
                    idx = ~isnan(o2_inv_temp(x,y,:));
                    fit_params = polyfit(time(idx),squeeze(o2_inv_temp(x,y,idx)),1); % mol/m2/day
                    o2_inv_trend(x,y,p,l) = 10.*365.25.*fit_params(1); % convert to mol/m2 per decade
                    % fit_params = fitlm(time(idx),squeeze(o2_inv_temp(x,y,idx))); % mol/m2/day
                    % o2_inv_trend(x,y,p,l) = 10.*365.25.*fit_params.Coefficients.Estimate(2); % conver to mol/m2 per decade
                end
                % if sum(~isnan(o2_inv_temp2(x,y,:))) > length(time)/2 % if more than half the grid cells are filled
                %     idx = ~isnan(o2_inv_temp2(x,y,:));
                %     fit_params = polyfit(time(idx),squeeze(o2_inv_temp2(x,y,idx)),1); % mol/m2/day
                %     o2_inv_trend2(x,y,p,l) = 10.*365.25.*fit_params(1); % conver to mol/m2 per decade
                % end
                if sum(~isnan(o2_avg_temp(x,y,:))) > length(time)/2 % if more than half the grid cells are filled
                    idx = ~isnan(o2_avg_temp(x,y,:));
                    fit_params = polyfit(time(idx),squeeze(o2_avg_temp(x,y,idx)),1); % umol/kg/day
                    o2_avg_trend(x,y,p,l) = 10.*365.25.*fit_params(1); % conver to umol/kg per decade
                end
            end
        end
    end
    % calculate trends in umol/kg
    for x = 1:length(lon)
        for y =1:length(lat)
            for z = 1:length(depth)
                if sum(~isnan(o2(x,y,z,:))) > length(time)/2 % if more than half the grid cells are filled
                    idx = ~isnan(o2(x,y,z,:));
                    fit_params = polyfit(time(idx),squeeze(o2(x,y,z,idx)),1); % per day
                    o2_trend(x,y,z,p) = 10.*365.25.*fit_params(1); % conver to umol/kg per decade
                end
            end
        end
    end
    % % plot trends
    % figure('Position',[100 100 1200 600]);
    % worldmap(lat_lims,lon_lims);
    % pcolorm(lat,lon,o2_inv_trend(:,:,p)');
    % c = colorbar;
    % mlabel off; plabel off; gridm off;
    % drawnow;
    % close;
    % inventory time series
    for l = 1:length(layers)
        idx_top = find(depth==lyr_top(l));
        idx_bot = find(depth==lyr_bot(l));
        for t = 1:length(time)
            o2_temp = o2_m3(:,:,idx_top:idx_bot,t);
            o2_temp = o2_temp(mask(:,:,idx_top:idx_bot)); % mol/m3
            dv_temp = dv(:,:,idx_top:idx_bot);
            dv_temp = dv_temp(mask(:,:,idx_top:idx_bot)); % m3
            o2_inv_ts(t,p,l) = sum(o2_temp.*dv_temp,'omitnan')./(10^15); % Pmol
        end
    end
    % average concentration time series w/ depth
    for t = 1:length(time)
        o2_temp = o2(:,:,:,t); % umol/kg
        o2_temp(~mask) = NaN; % umol/kg
        dv_temp = dv; % m3
        dv_temp(~mask) = NaN; % m3
        o2_inv_depth_ts(t,p,:) = squeeze(sum(sum(o2_temp.*dv_temp,1,'omitnan'),2,'omitnan'))./...
            squeeze(sum(sum(dv_temp,1,'omitnan'),2,'omitnan')); % umol/kg
    end
    disp(['Inventories calculated for ' products{p}]);
end
o2_inv_ts(o2_inv_ts == 0) = NaN;
o2_inv_depth_ts(o2_inv_depth_ts == 0) = NaN;

% save stats for each product
save([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_STATS_' file_date],...
    'o2_mean','o2_trend','o2_inv_mean','o2_inv_trend','o2_inv_ts',...
    'o2_inv_depth_ts','o2_inv_anom_2025','o2_inv_anom_2024',...
    'o2_avg_mean','o2_avg_trend','o2_inv_diff_2025_2024');

else

% load stats
load([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_STATS_' file_date '.mat']);

end

%% display trend statistics
% each prodct
for p = 1:length(products)
    for l = 1:length(layers)
         idx = ~isnan(o2_inv_ts(:,p,l));
         stats = polyfit(time(idx),o2_inv_ts(idx,p,l),1);
         trend(p,l) = 100*(stats(1)/mean(o2_inv_ts(idx,p,l))).*10.*365.25;
         iav(p,l) = std(o2_inv_ts(idx,p,l)-time(idx)*stats(1)+stats(2));
         disp([products{p} ' inventory trend  (' layers{l} ...
             ') = ' num2str(round(trend(p,l),3)) ' %/dec.']);
    end
end

for l = 1:length(layers)
    % volume
    idx_top = find(depth == lyr_top(l));
    idx_bot = find(depth == lyr_bot(l));
    dv_temp = dv(:,:,idx_top:idx_bot);
    mask_temp = mask(:,:,idx_top:idx_bot);
    volume = sum(dv_temp(mask_temp)).*1e-9*1e-6; %k m^3
    disp(['Volume (' layers{l} ') = ' num2str(volume) ' km^3 * 10^6']);
    % ensemble
    disp(['Ensemble mean inventory (' layers{l} ...
        ') = ' num2str(round(mean(mean(o2_inv_ts(:,:,l),'omitnan'),'omitnan'),3)) ' +/- ' ...
        num2str(round(std(mean(o2_inv_ts(:,:,l),'omitnan'),'omitnan'),3)) ' Pmol']);
    disp(['Ensemble mean inventory, 1990-2020 (' layers{l} ...
        ') = ' num2str(round(mean(mean(o2_inv_ts(baseline_idx,:,l),'omitnan'),'omitnan'),3)) ' +/- ' ...
        num2str(round(std(mean(o2_inv_ts(baseline_idx,:,l),'omitnan'),'omitnan'),3)) ' Pmol']);
    disp(['Ensemble mean inventory trend (' layers{l} ...
        ') = ' num2str(round(mean(trend(:,l)),3)) ' +/- ' ...
        num2str(round(std(trend(:,l)),3)) ' %/dec.']);
    disp(['Ensemble mean inventory IAV (' layers{l} ...
        ') = ' num2str(round(mean(iav(:,l)),3)) ' +/- ' ...
        num2str(round(std(iav(:,l)),3)) ' %/dec.']);
    % statistical interpolation
    disp(['SI ensemble mean inventory (' layers{l} ...
        ') = ' num2str(round(mean(mean(o2_inv_ts(:,1:5,l),'omitnan'),'omitnan'),3)) ' +/- ' ...
        num2str(round(std(mean(o2_inv_ts(:,1:5,l),'omitnan'),'omitnan'),3)) ' Pmol']);
    disp(['SI ensemble mean inventory trend (' layers{l} ...
        ') = ' num2str(round(mean(trend(1:5,l)),3)) ' +/- ' ...
        num2str(round(std(trend(1:5,l)),3)) ' %/dec.']);
    disp(['SI ensemble mean inventory IAV (' layers{l} ...
        ') = ' num2str(round(mean(iav(1:5,l)),3)) ' +/- ' ...
        num2str(round(std(iav(1:5,l)),3)) ' %/dec.']);
    % machine learning
    disp(['ML ensemble mean inventory (' layers{l} ...
        ') = ' num2str(round(mean(mean(o2_inv_ts(:,6:9,l),'omitnan'),'omitnan'),3)) ' +/- ' ...
        num2str(round(std(mean(o2_inv_ts(:,6:9,l),'omitnan'),'omitnan'),3)) ' Pmol']);
    disp(['ML ensemble mean inventory trend (' layers{l} ...
        ') = ' num2str(round(mean(trend(6:9,l)),3)) ' +/- ' ...
        num2str(round(std(trend(6:9,l)),3)) ' %/dec.']);
    disp(['ML ensemble mean inventory IAV (' layers{l} ...
        ') = ' num2str(round(mean(iav(6:9,l)),3)) ' +/- ' ...
        num2str(round(std(iav(6:9,l)),3)) ' %/dec.']);
end



%% calculate temp/sal anomalies
if ~exist([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_TEMP_SAL_' file_date '.mat'],'file')

% pre-allocate
EN4.temp = nan(360,173,42,length(time));
EN4.sal = nan(360,173,42,length(time));
EN4.temp_common = nan(length(lon),length(lat),length(depth),length(time));
EN4.sal_common = nan(length(lon),length(lat),length(depth),length(time));
% load dimensions
EN4.lon = ncread([fpath_EN4 'EN4.2.2/' 'EN.4.2.2.f.analysis.c14.196501.nc'],'lon');
EN4.lat = ncread([fpath_EN4 'EN4.2.2/' 'EN.4.2.2.f.analysis.c14.196501.nc'],'lat');
EN4.depth = ncread([fpath_EN4 'EN4.2.2/' 'EN.4.2.2.f.analysis.c14.196501.nc'],'depth');
% assemble annual means
for t = 1:length(time)
    temp_temp = []; sal_temp = [];
    for m = 1:12
        try
            temp_temp = cat(4,temp_temp,ncread([fpath_EN4 'EN4.2.2/' 'EN.4.2.2.f.analysis.c14.' ...
                num2str(year(t)) sprintf('%02d',m) '.nc'],'temperature'));
            sal_temp = cat(4,sal_temp,ncread([fpath_EN4 'EN4.2.2/' 'EN.4.2.2.f.analysis.c14.' ...
                num2str(year(t)) sprintf('%02d',m) '.nc'],'salinity'));
        catch
            try
                temp_temp = cat(4,temp_temp,ncread([fpath_EN4 'EN4.2.2/' 'EN.4.2.2.p.analysis.c14.' ...
                    num2str(year(t)) sprintf('%02d',m) '.nc'],'temperature'));
                sal_temp = cat(4,sal_temp,ncread([fpath_EN4 'EN4.2.2/' 'EN.4.2.2.p.analysis.c14.' ...
                    num2str(year(t)) sprintf('%02d',m) '.nc'],'salinity'));
            catch
                % do not replace NaNs
            end
        end
    end
    EN4.temp(:,:,:,t) = mean(temp_temp,4,'omitnan');
    EN4.sal(:,:,:,t) = mean(sal_temp,4,'omitnan');
end
% convert to celcius
EN4.temp = EN4.temp -273.15;
% interpolate to comon grid
[temp_lon,lon_idx] = sort(convert_lon(EN4.lon,'format','-180-180'));
for t = 1:length(time)
    EN4.temp_common(:,:,:,t) = interp3(EN4.lat,temp_lon',EN4.depth,...
        EN4.temp(lon_idx,:,:,t),lat,lon',depth);
    EN4.sal_common(:,:,:,t) = interp3(EN4.lat,temp_lon',EN4.depth,...
        EN4.sal(lon_idx,:,:,t),lat,lon',depth);
    % if start or end longitude is all nan, replace with mean of bounding longitudes
    if all(all(isnan(EN4.temp_common(1,:,:,t))))
        EN4.temp_common(1,:,:,t) = mean([EN4.temp_common(end,:,:,t);...
            EN4.temp_common(2,:,:,t)],1);
    elseif all(all(isnan(EN4.temp_common(end,:,:,t))))
        EN4.temp_common(end,:,:,t) = mean([EN4.temp_common(end-1,:,:,t);...
            EN4.temp_common(1,:,:,t)],1);
    end
    if all(all(isnan(EN4.sal_common(1,:,:,t))))
        EN4.sal_common(1,:,:,t) = mean([EN4.sal_common(end,:,:,t);...
            EN4.sal_common(2,:,:,t)],1);
    elseif all(all(isnan(EN4.sal_common(end,:,:,t))))
        EN4.sal_common(end,:,:,t) = mean([EN4.sal_common(end-1,:,:,t);...
            EN4.sal_common(1,:,:,t)],1);
    end
end
% replace zeros with NaN
EN4.temp_common(EN4.temp_common == 0) = NaN;
EN4.sal_common(EN4.sal_common == 0) = NaN;

% save temperature and salinity anomalies
save([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_TEMP_SAL_' file_date '.mat'],...
    'EN4','-v7.3');

else

% load stats
load([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_TEMP_SAL_' file_date '.mat']);

end

% %% calculate absolute o2 for oi products
% o2_anom_mean = o2_mean(:,:,:,1:5);
% o2_mean(:,:,:,1:5) = o2_anom_mean+ncei_mean.o2;

%% plot mean on depth levels
for d = 1:length(depth_indices)
    % plot each product individually
    for p = 1:length(products)
        figure('Position',[100 100 1200 600]);
        worldmap(lat_lims,lon_lims);
        title([labels{p} ' Oxygen at ' num2str(depth(depth_indices(d))) ...
            ' meters (\mumol kg^{-1})'],'FontSize',20);
        pcolorm(lat,[lon;lon(end)+1],...
            [o2_mean(1,:,depth_indices(d),p);o2_mean(:,:,depth_indices(d),p)]');
        contourm(lat,[lon;lon(end)+1],...
            [o2_mean(1,:,depth_indices(d),p);o2_mean(:,:,depth_indices(d),p)]',...
            o2_levels,'k','LineWidth',1,'ShowText','off');
        c = colorbar; clim(o2_lims); c.FontSize = 16; c.TickLength = 0;
        colormap(customcolormap([0;1],[0.7 0 1; 1 1 0],length(o2_levels)-1));
        plot_land('map',[1 1 1]);
        mlabel off; plabel off; gridm off;
        figname = ['Figures/mean_o2_' labels{p} '_' ...
            num2str(depth(depth_indices(d)))];
        export_fig([figname '.png'],'-transparent');
        figname = ['Figures/vectors/mean_o2_' labels{p} '_' ...
            num2str(depth(depth_indices(d)))];
        exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
        close;
    end
    % plot ensemble mean
    figure('Position',[100 100 1200 600]);
    worldmap(lat_lims,lon_lims);
    title(['Oxygen at ' num2str(depth(depth_indices(d))) ...
        ' meters (\mumol kg^{-1})'],'FontSize',20);
    pcolorm(lat,[lon;lon(end)+1],mean([o2_mean(1,:,depth_indices(d),:);...
        o2_mean(:,:,depth_indices(d),:)],4,'omitnan')');
    contourm(lat,[lon;lon(end)+1],mean([o2_mean(1,:,depth_indices(d),:);...
        o2_mean(:,:,depth_indices(d),:)],4,'omitnan')',...
        o2_levels,'k','LineWidth',1,'ShowText','off');
    c = colorbar; clim(o2_lims); c.FontSize = 16; c.TickLength = 0;
    colormap(customcolormap([0;1],[0.7 0 1; 1 1 0],length(o2_levels)-1));
    plot_land('map',[1 1 1]);
    mlabel off; plabel off; gridm off;
    figname = ['Figures/ensemble_mean_o2_' ...
        num2str(depth(depth_indices(d)))];
    export_fig([figname '.png'],'-transparent');
    figname = ['Figures/vectors/ensemble_mean_o2_' ...
        num2str(depth(depth_indices(d)))];
    exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
    close;
end

%% plot mean inventory
for l = 1:length(layers)
    % plot each product individually
    for p = 1:length(products)
        figure('Position',[100 100 1200 600]);
        worldmap(lat_lims,lon_lims);
        title([labels{p} ' Oxygen Inventory from ' layers{l} ...
            ' meters (\mumol kg^{-1})'],'FontSize',20);
        pcolorm(lat,[lon;lon(end)+1],...
            [o2_inv_mean(:,:,p,l);o2_inv_mean(end,:,p,l)]');
        contourm(lat,[lon;lon(end)+1],...
            [o2_inv_mean(:,:,p,l);o2_inv_mean(end,:,p,l)]',...
            o2_inv_levels{l,:},'k','LineWidth',1,'ShowText','off');
        c = colorbar; clim(o2_inv_lims(l,:)); c.FontSize = 16; c.TickLength = 0;
        colormap(customcolormap([0;1],[0.7 0 1; 1 1 0],length(o2_levels)-1));
        plot_land('map',[1 1 1]);
        mlabel off; plabel off; gridm off;
        figname = ['Figures/mean_o2_inv_' labels{p} '_' layers{l}];
        export_fig([figname '.png'],'-transparent');
        close;
    end
    % plot ensemble mean
    figure('Position',[100 100 1200 600]);
    worldmap(lat_lims,lon_lims);
    title(['Oxygen Inventory from ' layers{l} ...
        ' meters (\mumol kg^{-1})'],'FontSize',20);
    pcolorm(lat,[lon;lon(end)+1],[mean(o2_inv_mean(:,:,:,l),3,'omitnan');...
        mean(o2_inv_mean(end,:,:,l),3,'omitnan')]');
    contourm(lat,[lon;lon(end)+1],[mean(o2_inv_mean(:,:,:,l),3,'omitnan');...
        mean(o2_inv_mean(end,:,:,l),3,'omitnan')]',...
        o2_inv_levels{l,:},'k','LineWidth',1,'ShowText','off');
    c = colorbar; clim(o2_inv_lims(l,:)); c.FontSize = 16; c.TickLength = 0;
    colormap(customcolormap([0;1],[0.7 0 1; 1 1 0],length(o2_levels)-1));
    plot_land('map',[1 1 1]);
    mlabel off; plabel off; gridm off;
    figname = ['Figures/ensemble_mean_o2_inv_' layers{l}];
    export_fig([figname '.png'],'-transparent');
    close;
end

%% plot 2025 anomaly on depth levels
anom_idx = find(year == 2025);
baseline_idx = find(year >= 1990 & year <= 2020);
o2_anom = nan(length(lon),length(lat),length(depth_indices),length(products));
temp_anom = nan(length(lon),length(lat),length(depth_indices));
sal_anom = nan(length(lon),length(lat),length(depth_indices));
for d = 1:length(depth_indices)
    %% plot each product individually
    for p = 1:length(products)
        figure('Position',[100 100 1200 600]);
        worldmap(lat_lims,lon_lims);
        title([labels{p} ' Oxygen Anomaly in ' num2str(year(anom_idx)) ...
            ' at ' num2str(depth(depth_indices(d))) ...
            ' meters (\mumol kg^{-1})'],'FontSize',20);
        % calculate anomaly
        o2 = squeeze(ncread([fpath 'O2_Maps/' compilation.name],products{p},...
            [1 1 d1_idx 1],[Inf Inf d2_idx-d1_idx+1 Inf]));
        if ~any(any(~isnan(o2(:,:,depth_indices(d),anom_idx))))
            close; continue; end % check for any values in 2025
        o2_anom_temp = o2(:,:,depth_indices(d),anom_idx) - mean(o2(:,:,depth_indices(d),...
            baseline_idx(1):baseline_idx(end)),4,'omitnan');
        o2_anom_temp(~mask(:,:,depth_indices(d))) = NaN;
        o2_anom(:,:,d,p) = o2_anom_temp;
        pcolorm(lat-.25,[lon;lon(end)+1]-.25,[o2_anom(:,:,d,p);o2_anom(end,:,d,p)]');
        contourm(lat-.25,[lon;lon(end)+1]-.25,[o2_anom(:,:,d,p);o2_anom(end,:,d,p)]',...
            anom_levels,'k','LineWidth',1,'ShowText','off');
        c = colorbar; clim(anom_lims); c.FontSize = 16; c.TickLength = 0;
        colormap(cmocean('balance',length(anom_levels)-1,'pivot',0));
        plot_land('map',[1 1 1]);
        mlabel off; plabel off; gridm off;
        figname = ['Figures/o2_anom_' num2str(year(anom_idx)) '_' labels{p} '_' ...
            num2str(depth(depth_indices(d)))];
        export_fig([figname '.png'],'-transparent');
        close;
    end
    %% plot ensemble anomaly
    figure('Position',[100 100 1200 600]);
    worldmap(lat_lims,lon_lims);
    title(['Oxygen Anomaly at ' num2str(depth(depth_indices(d))) ...
        ' meters (\mumol kg^{-1})'],'FontSize',20);
    pcolorm(lat-.25,[lon;lon(end)+1]-.25,mean([o2_anom(:,:,d,:);o2_anom(end,:,d,:)],4,'omitnan')');
    contourm(lat-.25,[lon;lon(end)+1]-.25,mean([o2_anom(:,:,d,:);o2_anom(end,:,d,:)],4,'omitnan')',...
        anom_levels,'k','LineWidth',1,'ShowText','off');
    c = colorbar; clim(anom_lims); c.FontSize = 16; c.TickLength = 0;
    colormap(cmocean('balance',length(anom_levels)-1,'pivot',0));
    plot_land('map',[1 1 1]);
    mlabel off; plabel off; gridm off;
    figname = ['Figures/o2_ensemble_anom_' num2str(year(anom_idx)) '_' ...
        num2str(depth(depth_indices(d)))];
    export_fig([figname '.png'],'-transparent');
    figname = ['Figures/vectors/o2_ensemble_anom_' num2str(year(anom_idx)) '_' ...
        num2str(depth(depth_indices(d)))];
    exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
    close;
    %% plot temperature anomaly
    figure('Position',[100 100 1200 600]);
    worldmap(lat_lims,lon_lims);
    title(['Temperature Anomaly at ' num2str(depth(depth_indices(d))) ...
        ' meters (deg. C)'],'FontSize',20);
    % calculate anomaly
    temp_anom_temp = EN4.temp_common(:,:,depth_indices(d),anom_idx) - ...
        mean(EN4.temp_common(:,:,depth_indices(d),...
        baseline_idx(1):baseline_idx(end)),4,'omitnan');
    temp_anom_temp(~mask(:,:,depth_indices(d))) = NaN;
    temp_anom(:,:,d) = temp_anom_temp;
    pcolorm(lat-.25,[lon;lon(end)+1]-.25,mean([temp_anom(:,:,d);temp_anom(end,:,d)],4,'omitnan')');
    contourm(lat-.25,[lon;lon(end)+1]-.25,mean([temp_anom(:,:,d);temp_anom(end,:,d)],4,'omitnan')',...
        temp_anom_levels,'k','LineWidth',1,'ShowText','off');
    c = colorbar; clim(temp_anom_lims); c.FontSize = 16; c.TickLength = 0;
    colormap(cmocean('balance',length(temp_anom_levels)-1,'pivot',0));
    plot_land('map',[1 1 1]);
    mlabel off; plabel off; gridm off;
    figname = ['Figures/temp_ensemble_anom_' num2str(year(anom_idx)) '_' ...
        num2str(depth(depth_indices(d)))];
    export_fig([figname '.png'],'-transparent');
    figname = ['Figures/vectors/temp_ensemble_anom_' num2str(year(anom_idx)) '_' ...
        num2str(depth(depth_indices(d)))];
    exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
    close;
    %% calculate correlation between o2 and temp
    o2_anom_vec = mean(o2_anom(:,:,d,:),4,'omitnan');
    o2_anom_vec = o2_anom_vec(:);
    temp_anom_vec = temp_anom(:,:,d,:);
    temp_anom_vec = temp_anom_vec(:);
    idx = ~isnan(o2_anom_vec) & ~isnan(temp_anom_vec);
    [r2,pval] = corr(o2_anom_vec(idx),temp_anom_vec(idx));
    disp(['[O2] vs T correlation at 10 meters = ' sprintf('%.3g',r2) ...
        ' (p = ' sprintf('%.3g',pval) ')'])
end

%% plot 2025 inventory anomaly
for l = 1:length(layers)
    %% plot each product individually
    for p = 1:length(products)
        % plot inventory anomaly in 2025
        if ~any(any(~isnan(o2_inv_anom_2025)))
            continue; end % check for any values in 2025
        figure('Position',[100 100 1200 600]);
        worldmap(lat_lims,lon_lims);
        title([labels{p} ' Oxygen Anomaly in ' num2str(year(anom_idx_2025)) ...
            ' from ' layers{l} ' meters (mol m^{-2})'],'FontSize',20);
        pcolorm(lat-.25,[lon;lon(end)+1]-.25,...
            [o2_inv_anom_2025(:,:,l,p);o2_inv_anom_2025(end,:,l,p)]');
        contourm(lat-.25,[lon;lon(end)+1]-.25,...
            [o2_inv_anom_2025(:,:,l,p);o2_inv_anom_2025(end,:,l,p)]',...
            inv_anom_levels{l,:},'k','LineWidth',1,'ShowText','off');
        c = colorbar; clim(inv_anom_lims(l,:)); c.FontSize = 16; c.TickLength = 0;
        colormap(cmocean('balance',length(inv_anom_levels{l,:})-1,'pivot',0));
        plot_land('map',[1 1 1]);
        mlabel off; plabel off; gridm off;
        figname = ['Figures/o2_inv_anom_2025_' labels{p} '_' layers{l}];
        export_fig([figname '.png'],'-transparent');
        close;
    end
    %% plot ensemble anomaly
    figure('Position',[100 100 1200 600]);
    worldmap(lat_lims,lon_lims);
    title(['Oxygen Anomaly in 2025 from ' layers{l} ...
        ' meters (mol m^{-2})'],'FontSize',20);
    pcolorm(lat-.25,[lon;lon(end)+1]-.25,...
        mean([o2_inv_anom_2025(:,:,l,:);...
        o2_inv_anom_2025(end,:,l,:)],4,'omitnan')');
    contourm(lat-.25,[lon;lon(end)+1]-.25,...
            mean([o2_inv_anom_2025(:,:,l,:);...
            o2_inv_anom_2025(end,:,l,:)],4,'omitnan')',...
            inv_anom_levels{l,:},'k','LineWidth',1,'ShowText','off');
    c = colorbar; clim(inv_anom_lims(l,:)); c.FontSize = 16; c.TickLength = 0;
    colormap(flipud(cmocean('balance',length(inv_anom_levels{l,:})-1,'pivot',0)));
    plot_land('map',[1 1 1]);
    mlabel off; plabel off; gridm off;
    figname = ['Figures/o2_ensemble_inv_anom_2025_' layers{l}];
    export_fig([figname '.png'],'-transparent');
    close;
end

%% plot trends on depth levels
for d = 1:length(depth_indices)
    % plot each product trend individually
    for p = 1:length(products)
        figure('Position',[100 100 1200 600]);
        worldmap(lat_lims,lon_lims);
        set(gca,'Children');
        title([labels{p} ' Oxygen Trend at ' num2str(depth(depth_indices(d))) ...
            ' meters (\mumol kg^{-1} dec.^{-1})'],'FontSize',20);
        temp_mean_trend = o2_trend(:,:,depth_indices(d),p);
        pcolorm(lat,[lon;lon(end)+1],[temp_mean_trend(1,:);temp_mean_trend]');
        contourm(lat,[lon;lon(end)+1],[temp_mean_trend(1,:);temp_mean_trend]',...
            trend_levels,'k','LineWidth',1,'ShowText','off');
        c = colorbar; clim(trend_lims); c.FontSize = 16; c.TickLength = 0;
        colormap(cmocean('balance',length(trend_levels)-1,'pivot',0));
        plot_land('map',[1 1 1]);
        mlabel off; plabel off; gridm off;
        figname = ['Figures/o2_trend_' labels{p} '_' ...
            num2str(depth(depth_indices(d)))];
        export_fig([figname '.png'],'-transparent');
        figname = ['Figures/vectors/o2_trend_' labels{p} '_' ...
            num2str(depth(depth_indices(d)))];
        exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
        close;
    end
    % plot ensemble mean trend
    figure('Position',[100 100 1200 600]);
    worldmap(lat_lims,lon_lims);
    set(gca,'Children');
    title(['Oxygen Trend at ' num2str(depth(depth_indices(d))) ...
        ' meters (\mumol kg^{-1} dec.^{-1})'],'FontSize',20);
    ens_mean_trend = mean(o2_trend(:,:,depth_indices(d),:),4,'omitnan');
    sum_pos = sum(o2_trend(:,:,depth_indices(d),:) > 0,4);
    sum_neg = sum(o2_trend(:,:,depth_indices(d),:) < 0,4);
    idx = sum_pos >= 7 | sum_neg >= 7; ens_mean_trend(~idx) = NaN;
    pcolorm(lat,[lon;lon(end)+1],[ens_mean_trend(1,:);ens_mean_trend]');
    contourm(lat,[lon;lon(end)+1],[ens_mean_trend(1,:);ens_mean_trend]',...
        trend_levels,'k','LineWidth',1,'ShowText','off');
    stipplem(lat,lon,~idx','density',300);
    c = colorbar; clim(trend_lims); c.FontSize = 16; c.TickLength = 0;
    colormap(cmocean('balance',length(trend_levels),'pivot',0));
    plot_land('map',[1 1 1]);
    mlabel off; plabel off; gridm off;
    figname = ['Figures/ensemble_mean_o2_trend_' ...
        num2str(depth(depth_indices(d)))];
    export_fig([figname '.png'],'-transparent');
    figname = ['Figures/vectors/ensemble_mean_o2_trend_ ' ...
        num2str(depth(depth_indices(d)))];
    exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
    close;
end

%% plot seasonal cycle
% for d = 1:length(depth_indices)
%     % plot each product cycle individually
%     for p = 6:length(products)
%         % calculate mean seasonal cycle
% 
%         % plot
%         figure('Position',[100 100 1200 600]);
%         worldmap(lat_lims,lon_lims);
%         set(gca,'Children');
%         title([labels{p} ' Oxygen Trend at ' num2str(depth(depth_indices(d))) ...
%             ' meters (\mumol kg^{-1}) dec.^{-1})'],'FontSize',20);
%         temp_mean_trend = o2_trend(:,:,depth_indices(d),p);
%         pcolorm(lat,[lon;lon(end)+1],[temp_mean_trend(1,:);temp_mean_trend]');
%         contourm(lat,[lon;lon(end)+1],[temp_mean_trend(1,:);temp_mean_trend]',...
%             trend_levels,'k','LineWidth',1,'ShowText','off');
%         c = colorbar; clim(trend_lims); c.FontSize = 16; c.TickLength = 0;
%         colormap(cmocean('balance',length(trend_levels)-1,'pivot',0));
%         plot_land('map',[1 1 1]);
%         mlabel off; plabel off; gridm off;
%         figname = ['Figures/o2_trend_ ' labels{p} '_' ...
%             num2str(depth(depth_indices(d)))];
%         export_fig([figname '.png'],'-transparent');
%         figname = ['Figures/vectors/o2_trend_ ' labels{p} '_' ...
%             num2str(depth(depth_indices(d)))];
%         exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
%         close;
%     end
%     % plot ensemble mean trend
%     figure('Position',[100 100 1200 600]);
%     worldmap(lat_lims,lon_lims);
%     set(gca,'Children');
%     title(['Oxygen Trend at ' num2str(depth(depth_indices(d))) ...
%         ' meters (\mumol kg^{-1}) dec.^{-1})'],'FontSize',20);
%     ens_mean_trend = mean(o2_trend(:,:,depth_indices(d),:),4,'omitnan');
%     sum_pos = sum(o2_trend(:,:,depth_indices(d),:) > 0,4);
%     sum_neg = sum(o2_trend(:,:,depth_indices(d),:) < 0,4);
%     idx = sum_pos >= 7 | sum_neg >= 7; ens_mean_trend(~idx) = NaN;
%     pcolorm(lat,[lon;lon(end)+1],[ens_mean_trend(1,:);ens_mean_trend]');
%     contourm(lat,[lon;lon(end)+1],[ens_mean_trend(1,:);ens_mean_trend]',...
%         trend_levels,'k','LineWidth',1,'ShowText','off');
%     stipplem(lat,lon,~idx','density',300);
%     c = colorbar; clim(trend_lims); c.FontSize = 16; c.TickLength = 0;
%     colormap(cmocean('balance',length(trend_levels),'pivot',0));
%     plot_land('map',[1 1 1]);
%     mlabel off; plabel off; gridm off;
%     figname = ['Figures/ensemble_mean_o2_trend_ ' ...
%         num2str(depth(depth_indices(d)))];
%     export_fig([figname '.png'],'-transparent');
%     figname = ['Figures/vectors/ensemble_mean_o2_trend_ ' ...
%         num2str(depth(depth_indices(d)))];
%     exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
%     close;
% end

%% plot timeseries
for l = 1:length(layers)
    for included_products = 1:3
        figure('Position',[100 100 800 300]); hold on; box on; grid on;
        set(gca,'TitleHorizontalAlignment','left');
        if included_products == 1
            title('a) All products','FontWeight','normal');
        elseif included_products == 2
            title('b) Statistical interpolation products','FontWeight','normal');
        elseif included_products == 3
            title('c) Machine learning products','FontWeight','normal');
        end
        plot(time,repmat(0,length(time),1),'k--');
        for p = 1:length(products)
            if p <= 5
                if included_products == 1 | included_products == 2
                    plot(time,o2_inv_ts(:,p,l)-mean(o2_inv_ts(:,p,l),...
                        'omitnan'),'linewidth',2,'LineStyle','-');
                else
                   plot(time,o2_inv_ts(:,p,l)-mean(o2_inv_ts(:,p,l),...
                    'omitnan'),'linewidth',2,'LineStyle','none');
                end
            else
                if included_products == 1 | included_products == 3
                    plot(time,o2_inv_ts(:,p,l)-mean(o2_inv_ts(:,p,l),...
                        'omitnan'),'linewidth',2,'LineStyle','--');
                else
                    plot(time,o2_inv_ts(:,p,l)-mean(o2_inv_ts(:,p,l),...
                    'omitnan'),'linewidth',2,'LineStyle','none');
                end
            end
        end
        ylabel({'O_{2} Inventory Anomaly';['(Pmol, ' layers{l} 'm)']});
        xlim([datenum(1964,12,31) datenum(2026,1,1)]);
        datetick('x','keeplimits');
        figname = ['Figures/mean_inv_o2_' layers{l}];
        if included_products == 1
            legend([{''} labels],'NumColumns',3,'Location','southwest');
            export_fig([figname '.png'],'-transparent');
        elseif included_products == 2
            legend([{''} labels(1:5) {'' '' '' ''}],'NumColumns',2,'Location','southwest');
            export_fig([figname 'si_only.png'],'-transparent');
        elseif included_products == 3
            legend([{'' '' '' '' '' ''} labels(6:9)],'NumColumns',2,'Location','southwest');
            export_fig([figname 'ml_only.png'],'-transparent');
        end
        close;
    end
end

%% plot timeseries
for l = 1:length(layers)
    figure('Position',[100 100 800 300]); hold on;
    title(['O_{2} Inventory Anomaly (' layers{l} ' meters)']);
    plot(time,repmat(0,length(time),1),'k--');
    for p = 6:7
        if p <= 5
            plot(time,o2_inv_ts(:,p,l)-mean(o2_inv_ts(:,p,l),...
                'omitnan'),'linewidth',2,'LineStyle','-');
        else
            plot(time,o2_inv_ts(:,p,l)-mean(o2_inv_ts(:,p,l),...
                'omitnan'),'linewidth',2,'LineStyle','--');
        end
    end
    ylabel('O_{2} Inventory Anomaly (Pmol)');
    xlim([datenum(1964,12,31) datenum(2026,1,1)]);
    datetick('x','keeplimits');
    legend([{''} labels{6:7}],'NumColumns',3,'Location','southwest');
    figname = ['Figures/mean_inv_o2_EN4_vs_IAP_' layers{l}];
    export_fig([figname '.png'],'-transparent');
    close;
end


%% plot standard deviation products
% OI
o2_inv_std_OI_0_1000 = std(o2_inv_mean_10_1000(:,:,1:5),[],3,'omitnan');
figure; worldmap(lat_lims,lon_lims);
pcolorm(lat,lon,o2_inv_std_OI_0_1000');
colorbar;
% ML
o2_inv_std_ML_0_1000 = std(o2_inv_mean_10_1000(:,:,6:9),[],3,'omitnan');
figure; worldmap(lat_lims,lon_lims);
pcolorm(lat,lon,o2_inv_std_ML_0_1000');
colorbar;



%% plot SoTC Figure 2
% create figure
stip = false; stip_dens = 200; sig_prods = 7;
plt_alpha = 1; plt_bkgr = 'w';
figure('Position',[100 100 1200 1000],'visible','on','Color','w');
SotC_fig = tiledlayout(3,2,'TileSpacing','tight','Padding','none');
ax1 = nexttile; % inventory from 10-50m
idx_l = find(strcmp(layers,'0-50'));
worldmap(lat_lims,lon_lims);
set(ax1,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax1,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['a) Oxygen Inventory from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})'],...
    'FontWeight','normal','FontSize',12);
pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[mean(o2_inv_mean(:,:,:,idx_l),3,'omitnan');...
    mean(o2_inv_mean(end,:,:,idx_l),3,'omitnan')]');
c = colorbar; clim(o2_inv_lims(idx_l,:)); c.TickLength = 0;
colormap(ax1,customcolormap([0;1],[0.7 0 1; 1 1 0],length(o2_inv_levels{idx_l,:})-1));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax2 = nexttile; % inventory from 100-1000
idx_l = find(strcmp(layers,'50-1000'));
worldmap(lat_lims,lon_lims);
set(ax2,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax2,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['b) Oxygen Inventory from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})'],...
    'FontWeight','normal','FontSize',12);
pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[mean(o2_inv_mean(:,:,:,idx_l),3);...
    mean(o2_inv_mean(end,:,:,idx_l),3,'omitnan')]');
c = colorbar; clim(o2_inv_lims(idx_l,:)); c.TickLength = 0;
colormap(ax2,customcolormap([0;1],[0.7 0 1; 1 1 0],length(o2_inv_levels{idx_l,:})-1));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax3 = nexttile; % inventory trend from 10-50
idx_l = find(strcmp(layers,'0-50'));
worldmap(lat_lims,lon_lims);
set(ax3,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax3,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['c) Oxygen Inventory Trend from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2} dec.^{-1})'],...
    'FontWeight','normal','FontSize',12);
ens_mean_trend = mean(o2_inv_trend(:,:,:,idx_l),3,'omitnan');
ens_mean_trend_sig = mean(o2_inv_trend(:,:,:,idx_l),3,'omitnan');
ens_mean_trend_non = mean(o2_inv_trend(:,:,:,idx_l),3,'omitnan');
sum_pos = sum(o2_inv_trend(:,:,:,idx_l) > 0,3);
sum_neg = sum(o2_inv_trend(:,:,:,idx_l) < 0,3);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_mean_trend_sig(~idx) = NaN; ens_mean_trend_non(idx) = NaN; 
h1=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_mean_trend_sig;ens_mean_trend_sig(end,:)]');
h2=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_mean_trend_non;ens_mean_trend_non(end,:)]');
warning('off','MATLAB:nearlySingularMatrix')
hatchfill2(h2,'cross');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_mean_trend)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_trend_lims(idx_l,:)); c.TickLength = 0;
colormap(ax3,flipud(cmocean('balance',length(inv_trend_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax4 = nexttile; % inventory trend from 100 to 1000
idx_l = find(strcmp(layers,'50-1000'));
worldmap(lat_lims,lon_lims);
set(ax4,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax4,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['d) Oxygen Inventory Trend from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2} dec.^{-1})'],...
    'FontWeight','normal','FontSize',12);
ens_mean_trend = mean(o2_inv_trend(:,:,:,idx_l),3,'omitnan');
ens_mean_trend_sig = mean(o2_inv_trend(:,:,:,idx_l),3,'omitnan');
ens_mean_trend_non = mean(o2_inv_trend(:,:,:,idx_l),3,'omitnan');
sum_pos = sum(o2_inv_trend(:,:,:,idx_l) > 0,3);
sum_neg = sum(o2_inv_trend(:,:,:,idx_l) < 0,3);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_mean_trend_sig(~idx) = NaN; ens_mean_trend_non(idx) = NaN; 
h1=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_mean_trend_sig;ens_mean_trend_sig(end,:)]');
h2=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_mean_trend_non;ens_mean_trend_non(end,:)]');
warning('off','MATLAB:nearlySingularMatrix');
hatchfill2(h2,'cross');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_mean_trend)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_trend_lims(idx_l,:)); c.TickLength = 0;
colormap(ax4,flipud(cmocean('balance',length(inv_trend_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax5 = nexttile([1 2]); % timeseries
hold on; box on; grid on;
set(ax5,'FontSize',12,'LineWidth',2,'GridLineWidth',1,...
    'TitleHorizontalAlignment','left');
title('e) Global Oxygen Inventory Anomalies from 0 to 1800 meters (Pmol), [GOBAI inventory on right axis (Pmol)]',...
    'FontWeight','normal','FontSize',12);
plot(time,repmat(0,length(time),1),'k:');
idx_l = find(strcmp(layers,'0-1800'));
for p = 1:length(products)
    if p <= 5
        yyaxis left; ax5.YAxis(1).Color = 'k';
        colororder('default');
        idx_temp = ~isnan(o2_inv_ts(:,p,idx_l));
        plot(time(idx_temp),o2_inv_ts(idx_temp,p,idx_l)-mean(o2_inv_ts(baseline_idx,p,idx_l),'omitnan'),...
            'linewidth',2,'LineStyle','-');
    else
        yyaxis left; ax5.YAxis(1).Color = 'k';
        idx_temp = ~isnan(o2_inv_ts(:,p,idx_l));
        plot(time(idx_temp),o2_inv_ts(idx_temp,p,idx_l)-mean(o2_inv_ts(baseline_idx,p,idx_l),'omitnan'),...
            'linewidth',2,'LineStyle','--');
    end
    if p == 6
        yyaxis right; ax5.YAxis(2).Color = 'k';
        plot(time(idx_temp),o2_inv_ts(idx_temp,p,idx_l),'LineStyle','none');
        ylim([mean(o2_inv_ts(baseline_idx,p,idx_l))-2 mean(o2_inv_ts(baseline_idx,p,idx_l))+2]);
    end
end
% ylabel('O_{2} Inventory Anomaly (Pmol)');
xlim([min(time) max(time)]); datetick('x','keeplimits'); 
legend([{''} labels],'NumColumns',3,'Location','southwest');
hold off
% save figure
if stip; figname = 'Figures/SotC_Fig2_stipple';
else; figname = 'Figures/SotC_Fig2'; end
exportgraphics(gcf,[figname '.png'],'Resolution',600);
if stip; figname = 'Figures/vectors/SotC_Fig2_stipple';
else; figname = 'Figures/vectors/SotC_Fig2'; end
exportgraphics(gcf,[figname '.eps'],'Resolution',600,...
    'ContentType','vector');
close;

%% plot SoTC Figure 3
% create figure
stip = false; stip_dens = 200; sig_prods = 5;
plt_alpha = 0.6; plt_bkgr = 'w';
figure('Position',[100 100 1400 1050],'visible','off','Color','w');
SotC_fig = tiledlayout(3,4,'TileSpacing','tight','Padding','none');
ax1 = nexttile([1 2]); % anomaly at 10
idx_l = find(strcmp(layers,'0-50'));
worldmap(lat_lims,lon_lims);
set(ax1,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax1,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['a) 2025 Inventory Anomaly from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})'],...
    'FontWeight','normal','FontSize',12);
ens_inv_anom = mean(o2_inv_anom_2025(:,:,idx_l,:),4,'omitnan');
ens_inv_anom_sig = mean(o2_inv_anom_2025(:,:,idx_l,:),4,'omitnan');
ens_inv_anom_non = mean(o2_inv_anom_2025(:,:,idx_l,:),4,'omitnan');
sum_pos = sum(o2_inv_anom_2025(:,:,idx_l,:) > 0,4);
sum_neg = sum(o2_inv_anom_2025(:,:,idx_l,:) < 0,4);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_inv_anom_sig(~idx) = NaN; ens_inv_anom_non(idx) = NaN; 
h1=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_inv_anom_sig;ens_inv_anom_sig(end,:)]');
h2=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_inv_anom_non;ens_inv_anom_non(end,:)]');
warning('off','MATLAB:nearlySingularMatrix')
hatchfill2(h2,'cross');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_anom)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)); c.TickLength = 0;
colormap(ax1,flipud(cmocean('balance',length(inv_trend_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax2 = nexttile([1 2]); % anomaly at 10
idx_l = find(strcmp(layers,'50-1000'));
worldmap(lat_lims,lon_lims);
set(ax2,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax2,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['b) 2025 Inventory Anomaly from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})'],...
    'FontWeight','normal','FontSize',12);
ens_inv_anom = mean(o2_inv_anom_2025(:,:,idx_l,:),4,'omitnan');
ens_inv_anom_sig = mean(o2_inv_anom_2025(:,:,idx_l,:),4,'omitnan');
ens_inv_anom_non = mean(o2_inv_anom_2025(:,:,idx_l,:),4,'omitnan');
sum_pos = sum(o2_inv_anom_2025(:,:,idx_l,:) > 0,4);
sum_neg = sum(o2_inv_anom_2025(:,:,idx_l,:) < 0,4);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_inv_anom_sig(~idx) = NaN; ens_inv_anom_non(idx) = NaN; 
h1=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_inv_anom_sig;ens_inv_anom_sig(end,:)]');
h2=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_inv_anom_non;ens_inv_anom_non(end,:)]');
warning('off','MATLAB:nearlySingularMatrix')
hatchfill2(h2,'cross');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_anom)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)); c.TickLength = 0;
colormap(ax2,flipud(cmocean('balance',length(inv_anom_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax3 = nexttile([1 2]); % anomaly at 10
idx_l = find(strcmp(layers,'0-50'));
worldmap(lat_lims,lon_lims);
set(ax3,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax3,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['c) 2025' char(8722) '2024 Inventory from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})'],...
    'FontWeight','normal','FontSize',12);
ens_inv_diff = mean(o2_inv_diff_2025_2024(:,:,idx_l,:),4,'omitnan');
ens_inv_diff_sig = mean(o2_inv_diff_2025_2024(:,:,idx_l,:),4,'omitnan');
ens_inv_diff_non = mean(o2_inv_diff_2025_2024(:,:,idx_l,:),4,'omitnan');
sum_pos = sum(o2_inv_diff_2025_2024(:,:,idx_l,:) > 0,4);
sum_neg = sum(o2_inv_diff_2025_2024(:,:,idx_l,:) < 0,4);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_inv_diff_sig(~idx) = NaN; ens_inv_diff_non(idx) = NaN; 
h1=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_inv_diff_sig;ens_inv_diff_sig(end,:)]');
h2=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_inv_diff_non;ens_inv_diff_non(end,:)]');
warning('off','MATLAB:nearlySingularMatrix')
hatchfill2(h2,'cross');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_diff)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)); c.TickLength = 0;
colormap(ax3,flipud(cmocean('balance',length(inv_anom_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax4 = nexttile([1 2]); % anomaly at 10
idx_l = find(strcmp(layers,'50-1000'));
worldmap(lat_lims,lon_lims);
set(ax4,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax4,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['d) 2025' char(8722) '2024 Inventory from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})'],...
    'FontWeight','normal','FontSize',12);
ens_inv_diff = mean(o2_inv_diff_2025_2024(:,:,idx_l,:),4,'omitnan');
ens_inv_diff_sig = mean(o2_inv_diff_2025_2024(:,:,idx_l,:),4,'omitnan');
ens_inv_diff_non = mean(o2_inv_diff_2025_2024(:,:,idx_l,:),4,'omitnan');
sum_pos = sum(o2_inv_diff_2025_2024(:,:,idx_l,:) > 0,4);
sum_neg = sum(o2_inv_diff_2025_2024(:,:,idx_l,:) < 0,4);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_inv_diff_sig(~idx) = NaN; ens_inv_diff_non(idx) = NaN; 
h1=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_inv_diff_sig;ens_inv_diff_sig(end,:)]');
h2=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_inv_diff_non;ens_inv_diff_non(end,:)]');
warning('off','MATLAB:nearlySingularMatrix')
hatchfill2(h2,'cross');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_diff)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)); c.TickLength = 0; colormap(ax4,flipud(cmocean('balance',length(inv_trend_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
%
ax5 = nexttile([1 3]); % timeseries
hold on; box on;
contourf(time,depth,squeeze(mean(o2_inv_depth_ts,2,'omitnan'))'-...
    mean(squeeze(mean(o2_inv_depth_ts,2,'omitnan')),1,'omitnan')',-6:0.5:6);
set(ax5,'YDir','reverse','FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
ylim(ax5,[1 1800]);
title('e) Global Oxygen Anomalies (\mumol kg^{-1})',...
    'FontWeight','normal','FontSize',12);
c = colorbar; clim([-6 6]); c.TickLength = 0;
colormap(ax5,flipud(cmocean('balance',24,'pivot',0)));
ylabel('Depth (m)');
xlim([min(time) max(time)]); datetick('x','keeplimits'); 
hold off;
ax6 = nexttile([1 1]); hold on; box on;
set(ax6,'YDir','reverse','XColor','r','linewidth',1);
ens_mean_o2_by_depth = squeeze(mean(mean(o2_inv_depth_ts,1,'omitnan'),2,'omitnan'));
ens_std_o2_by_depth = squeeze(std(mean(o2_inv_depth_ts,1,'omitnan'),[],2,'omitnan'));
fill([ens_mean_o2_by_depth-ens_std_o2_by_depth;...
    flipud(ens_mean_o2_by_depth+ens_std_o2_by_depth)],...
    [depth;flipud(depth)],'r','FaceAlpha',0.5,'LineStyle','none');
plot(ens_mean_o2_by_depth,depth,'linewidth',2,'Color','r');
xlim([130 270]); ylim([1 1800]); hold off;
ax7 = axes('Position', ax6.Position, 'Color', 'none','XAxisLocation', 'top');
hold on; grid on;
set(ax7,'YDir','reverse','XColor','k','linewidth',1);
da(any(isnan(o2_trend),4)) = NaN;
ens_mean_trend_by_depth = squeeze(mean(sum(sum(o2_trend.*da,1,'omitnan'),2,'omitnan')./...
    sum(sum(da,1,'omitnan'),2,'omitnan'),4));
ens_std_trend_by_depth = squeeze(std(sum(sum(o2_trend.*da,1,'omitnan'),2,'omitnan')./...
    sum(sum(da,1,'omitnan'),2,'omitnan'),[],4));
fill([ens_mean_trend_by_depth-ens_std_trend_by_depth;...
    flipud(ens_mean_trend_by_depth+ens_std_trend_by_depth)],...
    [depth;flipud(depth)],'k','FaceAlpha',0.5,'LineStyle','none');
plot(ens_mean_trend_by_depth,depth,'linewidth',2,'Color','k');
xlim([-1.7 0]); ylim([1 1800]); 
text(-0.65,1400,{'ensemble';'mean trend';'(\mumol kg^{-1} dec.^{-1})'},...
    'FontSize',8,'HorizontalAlignment','right');
text(-1.5,1000,{'ensemble';'mean';'(\mumol kg^{-1})'},...
    'FontSize',8,'Color','r');
text(-1.9,-30,'f)','FontSize',12);
% save figure
if stip; figname = 'Figures/SotC_Fig3_stipple';
else; figname = 'Figures/SotC_Fig3'; end
exportgraphics(gcf,[figname '.png'],'Resolution',600);
if stip; figname = 'Figures/vectors/SotC_Fig3_stipple';
else; figname = 'Figures/vectors/SotC_Fig3'; end
exportgraphics(gcf,[figname '.eps'],'Resolution',600,...
    'ContentType','vector');
close;

%% plot SoTC Figure 2 (with average values)
% create figure
stip = false; stip_dens = 200; sig_prods = 7;
plt_alpha = 1; plt_bkgr = 'w';
figure('Position',[100 100 1200 1000],'visible','on','Color','w');
SotC_fig = tiledlayout(3,2,'TileSpacing','tight','Padding','none');
ax1 = nexttile; % inventory from 10-50m
idx_l = find(strcmp(layers,'0-50'));
worldmap(lat_lims,lon_lims);
set(ax1,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax1,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['a) Average [O_{2}] from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})'],...
    'FontWeight','normal','FontSize',12);
pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[mean(o2_avg_mean(:,:,:,idx_l),3,'omitnan');...
    mean(o2_avg_mean(end,:,:,idx_l),3,'omitnan')]');
c = colorbar; clim(o2_lims); c.TickLength = 0;
colormap(ax1,customcolormap([0;1],[0.7 0 1; 1 1 0],length(o2_levels)-1));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax2 = nexttile; % inventory from 100-1000
idx_l = find(strcmp(layers,'50-1000'));
worldmap(lat_lims,lon_lims);
set(ax2,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax2,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['b) Average [O_{2}] from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})'],...
    'FontWeight','normal','FontSize',12);
pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[mean(o2_avg_mean(:,:,:,idx_l),3,'omitnan');...
    mean(o2_avg_mean(end,:,:,idx_l),3,'omitnan')]');
c = colorbar; clim(o2_lims); c.TickLength = 0;
colormap(ax2,customcolormap([0;1],[0.7 0 1; 1 1 0],length(o2_levels)-1));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax3 = nexttile; % inventory trend from 10-50
idx_l = find(strcmp(layers,'0-50'));
worldmap(lat_lims,lon_lims);
set(ax3,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax3,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['c) Average [O_{2}] Trend from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2} dec.^{-1})'],...
    'FontWeight','normal','FontSize',12);
ens_mean_trend = mean(o2_avg_trend(:,:,:,idx_l),3,'omitnan');
ens_mean_trend_sig = mean(o2_avg_trend(:,:,:,idx_l),3,'omitnan');
ens_mean_trend_non = mean(o2_avg_trend(:,:,:,idx_l),3,'omitnan');
sum_pos = sum(o2_avg_trend(:,:,:,idx_l) > 0,3);
sum_neg = sum(o2_avg_trend(:,:,:,idx_l) < 0,3);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_mean_trend_sig(~idx) = NaN; ens_mean_trend_non(idx) = NaN; 
h1=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_mean_trend_sig;ens_mean_trend_sig(end,:)]');
h2=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_mean_trend_non;ens_mean_trend_non(end,:)]');
warning('off','MATLAB:nearlySingularMatrix')
hatchfill2(h2,'cross');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_mean_trend)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(trend_lims); c.TickLength = 0;
colormap(ax3,flipud(cmocean('balance',length(trend_levels)-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax4 = nexttile; % inventory trend from 100 to 1000
idx_l = find(strcmp(layers,'50-1000'));
worldmap(lat_lims,lon_lims);
set(ax4,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax4,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['d) Average [O_{2}] Trend from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2} dec.^{-1})'],...
    'FontWeight','normal','FontSize',12);
ens_mean_trend = mean(o2_avg_trend(:,:,:,idx_l),3,'omitnan');
ens_mean_trend_sig = mean(o2_avg_trend(:,:,:,idx_l),3,'omitnan');
ens_mean_trend_non = mean(o2_avg_trend(:,:,:,idx_l),3,'omitnan');
sum_pos = sum(o2_avg_trend(:,:,:,idx_l) > 0,3);
sum_neg = sum(o2_avg_trend(:,:,:,idx_l) < 0,3);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_mean_trend_sig(~idx) = NaN; ens_mean_trend_non(idx) = NaN; 
h1=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_mean_trend_sig;ens_mean_trend_sig(end,:)]');
h2=pcolorm(lat-0.5,[lon;lon(end)+1]-0.5,[ens_mean_trend_non;ens_mean_trend_non(end,:)]');
warning('off','MATLAB:nearlySingularMatrix');
hatchfill2(h2,'cross');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_mean_trend)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(trend_lims); c.TickLength = 0;
colormap(ax4,flipud(cmocean('balance',length(trend_levels)-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax5 = nexttile([1 2]); % timeseries
hold on; box on; grid on;
set(ax5,'FontSize',12,'LineWidth',2,'GridLineWidth',1,...
    'TitleHorizontalAlignment','left');
title('e) Global Oxygen Inventory Anomalies from 0 to 1800 meters (Pmol), [GOBAI inventory on right axis (Pmol)]',...
    'FontWeight','normal','FontSize',12);
plot(time,repmat(0,length(time),1),'k:');
idx_l = find(strcmp(layers,'0-1800'));
for p = 1:length(products)
    if p <= 5
        yyaxis left; ax5.YAxis(1).Color = 'k';
        colororder('default');
        idx_temp = ~isnan(o2_inv_ts(:,p,idx_l));
        plot(time(idx_temp),o2_inv_ts(idx_temp,p,idx_l)-mean(o2_inv_ts(baseline_idx,p,idx_l),'omitnan'),...
            'linewidth',2,'LineStyle','-');
    else
        yyaxis left; ax5.YAxis(1).Color = 'k';
        idx_temp = ~isnan(o2_inv_ts(:,p,idx_l));
        plot(time(idx_temp),o2_inv_ts(idx_temp,p,idx_l)-mean(o2_inv_ts(baseline_idx,p,idx_l),'omitnan'),...
            'linewidth',2,'LineStyle','--');
    end
    if p == 6
        yyaxis right; ax5.YAxis(2).Color = 'k';
        plot(time(idx_temp),o2_inv_ts(idx_temp,p,idx_l),'LineStyle','none');
        ylim([mean(o2_inv_ts(baseline_idx,p,idx_l))-2 mean(o2_inv_ts(baseline_idx,p,idx_l))+2]);
    end
end
% ylabel('O_{2} Inventory Anomaly (Pmol)');
xlim([min(time) max(time)]); datetick('x','keeplimits'); 
legend([{''} labels],'NumColumns',3,'Location','southwest');
hold off
% save figure
if stip; figname = 'Figures/SotC_Fig2_stipple';
else; figname = 'Figures/SotC_Fig2'; end
exportgraphics(gcf,[figname '.png'],'Resolution',600);
if stip; figname = 'Figures/vectors/SotC_Fig2_stipple';
else; figname = 'Figures/vectors/SotC_Fig2'; end
exportgraphics(gcf,[figname '.eps'],'Resolution',600,...
    'ContentType','vector');
close;