% this compares annual means of 
close all;

%% file and plot properties
file_date = 'Jan2026';
fpath = '/raid/sharp/matlab/GODIP-DO/';
compilation.name = ['GODIP-DO_NCEI_JDS_' file_date '.nc'];
% products = {'ncei' 'iap' 'gt_oi' 'rb' 'sjtu_gr' 'gobai' 'gt_ml' 'jingwei' 'han_zhou'};
products = ncread([fpath 'O2_Maps/' compilation.name],'products');
labels = {'NCEI' 'IAP' 'GT-OI' 'UTAS' 'SJTU-GR' 'PMEL' 'GT-ML' 'SJTU-JW' 'SJTU-HZ'};
lat_lims = [-75 75]; lon_lims = [20 380];
o2_lims = [0 400]; o2_levels = o2_lims(1):12.5:o2_lims(end);
anom_lims = [-40 40]; anom_levels = anom_lims(1):2.5:anom_lims(end);
trend_lims = [-6 6]; trend_levels = trend_lims(1):0.5:trend_lims(end);
temp_anom_lims = [-3.5 3.5]; temp_anom_levels = temp_anom_lims(1):0.25:temp_anom_lims(end);

%% import compilation dimensions
lat = ncread([fpath 'O2_Maps/' compilation.name],'lat');
lon = ncread([fpath 'O2_Maps/' compilation.name],'lon');
time = ncread([fpath 'O2_Maps/' compilation.name],'time');
depth = ncread([fpath 'O2_Maps/' compilation.name],'depth');
mask = logical(ncread([fpath 'O2_Maps/' compilation.name],'mask'));


%% depth limits for inventory calculations
d1 = 10; d2 = 1000;
[~,d1_idx] = min(abs(depth-d1));
[~,d2_idx] = min(abs(depth-d2));
depth = depth(d1_idx:d2_idx);
depths_to_plot = [10 250];
depth_indices = find(any(depth == depths_to_plot,2));
date = datevec(time); year = date(:,1);

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

%% calculate o2 inventories and other statistics
% define volume and height
[dv,da,dh] = weights3d(lon,lat,depth);
dv = repmat(dv,1,1,1,length(time));
mask = mask(:,:,d1_idx:d2_idx,:);
mask = all(mask,4); % set constant mask over time

if ~exist([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_STATS_' file_date '.mat'],'file')

% calculate statistics
o2_mean = nan(length(lon),length(lat),d2_idx-d1_idx+1,length(products));
o2_trend = nan(length(lon),length(lat),d2_idx-d1_idx+1,length(products));
o2_inv_mean_10_1000 = nan(length(lon),length(lat),length(products));
o2_inv_mean_10_100 = nan(length(lon),length(lat),length(products));
o2_inv_mean_100_1000 = nan(length(lon),length(lat),length(products));
o2_inv_ts_10_1000 = nan(length(time),length(products));
for p = 1:length(products)
    % o2 maps over time in umol/kg
    o2 = squeeze(ncread([fpath 'O2_Maps/' compilation.name],products{p},...
        [1 1 d1_idx 1],[Inf Inf d2_idx-d1_idx+1 Inf]));
    if p <= 5 % add woa mean to o2 anomalies
        o2 = o2+ncei_mean.o2;
    end
    o2_mean(:,:,:,p) = mean(o2,4,'omitnan');
    o2_m3 = o2.*1027.*(1e-6); % umol/kg * kg/m3 * mol/umol = mol/m3
    o2_m3(~mask) = NaN; dv(~mask) = NaN; % remove values outside common mask
    % calculate long-term mean inventory
    o2_inv_temp_10_1000 = squeeze((sum(o2_m3,3,'omitnan').*sum(dh,3,'omitnan'))); % mol/m2
    o2_inv_mean_10_1000(:,:,p) = mean(o2_inv_temp_10_1000,3); % mol/m2
    % calculate trends
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
    % plot trends
    figure('Position',[100 100 1200 600]);
    worldmap(lat_lims,lon_lims);
    pcolorm(lat,lon,o2_trend(:,:,1,p)');
    c = colorbar;
    mlabel off; plabel off; gridm off;
    drawnow;
    close;
    % inventory time series
    for t = 1:length(time)
        o2_temp = o2_m3(:,:,:,t); o2_temp = o2_temp(mask); % mol/m3
        dv_temp = dv(:,:,:,t); dv_temp = dv_temp(mask); % m3
        o2_inv_ts_10_1000(t,p) = sum(o2_temp.*dv_temp,'omitnan')./(10^15); % Pmol
    end
    disp(['Inventories calculated for ' products{p}]);
end
o2_inv_mean_10_1000(o2_inv_mean_10_1000==0) = NaN; % replace zeros with NaN
o2_inv_ts_10_1000(o2_inv_ts_10_1000==0) = NaN; % replace zeros with NaN

% save stats for each product
save([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_STATS_' file_date],'o2_inv_mean_10_1000',...
    'o2_inv_ts_10_1000','o2_mean','o2_trend');

else

% load stats
load([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_STATS_' file_date '.mat']);

end

%% display trend statistics
for p = 1:length(products)
     idx = ~isnan(o2_inv_ts_10_1000(:,p));
     stats = polyfit(time(idx),o2_inv_ts_10_1000(idx,p),1);
     trend(p) = 100*(stats(1)/mean(o2_inv_ts_10_1000(idx,p))).*10.*365.25;
     iav(p) = std(o2_inv_ts_10_1000(idx,p)-time(idx)*stats(1)+stats(2));
     disp([products{p} ' trend = ' num2str(round(trend(p),3)) ' %/dec.']);
end

mean(mean(o2_inv_ts_10_1000,'omitnan'))
std(mean(o2_inv_ts_10_1000,'omitnan'))

mean(trend)
std(trend)

mean(mean(o2_inv_ts_10_1000(:,1:5),'omitnan'))
std(mean(o2_inv_ts_10_1000(:,1:5),'omitnan'))

mean(trend(1:5))
std(trend(1:5))
mean(iav(1:5))
std(iav(1:5))

mean(mean(o2_inv_ts_10_1000(:,6:9),'omitnan'))
std(mean(o2_inv_ts_10_1000(:,6:9),'omitnan'))

mean(trend(6:9))
std(trend(6:9))
mean(iav(6:9))
std(iav(6:9))

%% calculate temp/sal anomalies
if ~exist([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_TEMP_SAL_' file_date '.mat'],'file')

% pre-allocate
EN4.temp = nan(360,173,42,length(time));
EN4.sal = nan(360,173,42,length(time));
EN4.temp_common = nan(length(lon),length(lat),length(depth),length(time));
EN4.sal_common = nan(length(lon),length(lat),length(depth),length(time));
% load dimensions
EN4.lon = ncread([fpath 'EN4.2.2/' 'EN.4.2.2.f.analysis.c14.196501.nc'],'lon');
EN4.lat = ncread([fpath 'EN4.2.2/' 'EN.4.2.2.f.analysis.c14.196501.nc'],'lat');
EN4.depth = ncread([fpath 'EN4.2.2/' 'EN.4.2.2.f.analysis.c14.196501.nc'],'depth');
% assemble annual means
for t = 1:length(time)
    temp_temp = []; sal_temp = [];
    for m = 1:12
        try
            temp_temp = cat(4,temp_temp,ncread([fpath 'EN4.2.2/' 'EN.4.2.2.f.analysis.c14.' ...
                num2str(year(t)) sprintf('%02d',m) '.nc'],'temperature'));
            sal_temp = cat(4,sal_temp,ncread([fpath 'EN4.2.2/' 'EN.4.2.2.f.analysis.c14.' ...
                num2str(year(t)) sprintf('%02d',m) '.nc'],'salinity'));
        catch
            try
                temp_temp = cat(4,temp_temp,ncread([fpath 'EN4.2.2/' 'EN.4.2.2.p.analysis.c14.' ...
                    num2str(year(t)) sprintf('%02d',m) '.nc'],'temperature'));
                sal_temp = cat(4,sal_temp,ncread([fpath 'EN4.2.2/' 'EN.4.2.2.p.analysis.c14.' ...
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

%% plot mean inventory
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
        figname = ['Figures/mean_o2_ ' labels{p} '_' ...
            num2str(depth(depth_indices(d)))];
        export_fig([figname '.png'],'-transparent');
        figname = ['Figures/vectors/mean_o2_ ' labels{p} '_' ...
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
    figname = ['Figures/ensemble_mean_o2_ ' ...
        num2str(depth(depth_indices(d)))];
    export_fig([figname '.png'],'-transparent');
    figname = ['Figures/vectors/ensemble_mean_o2_ ' ...
        num2str(depth(depth_indices(d)))];
    exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
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
        figname = ['Figures/vectors/o2_anom_' num2str(year(anom_idx)) '_' labels{p} '_' ...
            num2str(depth(depth_indices(d)))];
        exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
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
    figname = ['Figures/o2_ensemble_anom_' num2str(year(anom_idx)) '_' ...
        num2str(depth(depth_indices(d)))];
    export_fig([figname '.png'],'-transparent');
    figname = ['Figures/vectors/o2_ensemble_anom_' num2str(year(anom_idx)) '_' ...
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
        figname = ['Figures/vectors/o2_anom_' num2str(year(anom_idx)) '_' labels{p} '_' ...
            num2str(depth(depth_indices(d)))];
        exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
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
    figname = ['Figures/o2_ensemble_anom_' num2str(year(anom_idx)) '_' ...
        num2str(depth(depth_indices(d)))];
    export_fig([figname '.png'],'-transparent');
    figname = ['Figures/vectors/o2_ensemble_anom_' num2str(year(anom_idx)) '_' ...
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

%% plot trends
for d = 1:length(depth_indices)
    % plot each product trend individually
    for p = 1:length(products)
        figure('Position',[100 100 1200 600]);
        worldmap(lat_lims,lon_lims);
        set(gca,'Children');
        title([labels{p} ' Oxygen Trend at ' num2str(depth(depth_indices(d))) ...
            ' meters (\mumol kg^{-1}) dec.^{-1})'],'FontSize',20);
        temp_mean_trend = o2_trend(:,:,depth_indices(d),p);
        pcolorm(lat,[lon;lon(end)+1],[temp_mean_trend(1,:);temp_mean_trend]');
        contourm(lat,[lon;lon(end)+1],[temp_mean_trend(1,:);temp_mean_trend]',...
            trend_levels,'k','LineWidth',1,'ShowText','off');
        c = colorbar; clim(trend_lims); c.FontSize = 16; c.TickLength = 0;
        colormap(cmocean('balance',length(trend_levels)-1,'pivot',0));
        plot_land('map',[1 1 1]);
        mlabel off; plabel off; gridm off;
        figname = ['Figures/o2_trend_ ' labels{p} '_' ...
            num2str(depth(depth_indices(d)))];
        export_fig([figname '.png'],'-transparent');
        figname = ['Figures/vectors/o2_trend_ ' labels{p} '_' ...
            num2str(depth(depth_indices(d)))];
        exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
        close;
    end
    % plot ensemble mean trend
    figure('Position',[100 100 1200 600]);
    worldmap(lat_lims,lon_lims);
    set(gca,'Children');
    title(['Oxygen Trend at ' num2str(depth(depth_indices(d))) ...
        ' meters (\mumol kg^{-1}) dec.^{-1})'],'FontSize',20);
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
    figname = ['Figures/ensemble_mean_o2_trend_ ' ...
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
% figure('Position',[100 100 800 300]); hold on;
% for p = 1:length(products)
%     if p <= 4
%     plot(time,o2_inv_ts_10_1000(:,p),'linewidth',2);
%     else
%     plot(time,o2_inv_ts_10_1000(:,p)-mean(o2_inv_ts_10_1000(:,p),'omitnan'),...
%         'linewidth',2);
%     end
% end
% ylabel('O_{2} Inventory Anomaly (Pmol)');
% datetick('x'); legend(labels,'NumColumns',2);
% figname = 'Figures/global_mean_o2';
% export_fig([figname '.png'],'-transparent');
% figname = 'Figures/vectors/global_mean_o2';
% exportgraphics(gcf,[figname '.eps'],'ContentType','vector');

%% plot SoTC Figure 2
figure('Position',[100 100 1200 1000]);
SotC_fig = tiledlayout(3,2,'TileSpacing','tight','Padding','none');
ax1 = nexttile; % distribution at 10
worldmap(lat_lims,lon_lims);
set(ax1,'FontSize',12,'LineWidth',2);
setm(ax1,'MapProjection','robinson'); tightmap;
title(['[O_{2}] at ' num2str(depth(depth_indices(1))) ...
    ' meters (\mumol kg^{-1})']);
pcolorm(lat,[lon;lon(end)+1],mean([o2_mean(:,:,depth_indices(1),:);...
    o2_mean(end,:,depth_indices(1),:)],4,'omitnan')');
% contourm(lat,[lon;lon(end)+1],mean([o2_mean(:,:,depth_indices(1),:);...
%     o2_mean(end,:,depth_indices(1),:)],4,'omitnan')',o2_levels,'k',...
%             'LineWidth',1,'ShowText','off');
c = colorbar; clim(o2_lims); c.TickLength = 0;
colormap(ax1,customcolormap([0;1],[0.7 0 1; 1 1 0],length(o2_levels)-1));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax2 = nexttile; % distribution at 250
worldmap(lat_lims,lon_lims);
set(ax2,'FontSize',12,'LineWidth',2);
setm(ax2,'MapProjection','robinson'); tightmap;
title(['[O_{2}] at ' num2str(depth(depth_indices(2))) ...
    ' meters (\mumol kg^{-1})']);
pcolorm(lat,[lon;lon(end)+1],mean([o2_mean(:,:,depth_indices(2),:);...
    o2_mean(end,:,depth_indices(2),:)],4,'omitnan')');
% contourm(lat,[lon;lon(end)+1],mean([o2_mean(:,:,depth_indices(2),:);...
%     o2_mean(end,:,depth_indices(2),:)],4,'omitnan')',o2_levels,'k',...
%             'LineWidth',1,'ShowText','off');
c = colorbar; clim(o2_lims); c.TickLength = 0;
colormap(ax2,customcolormap([0;1],[0.7 0 1; 1 1 0],length(o2_levels)-1));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax3 = nexttile; % trend at 10
worldmap(lat_lims,lon_lims);
set(ax3,'FontSize',12,'LineWidth',2);
setm(ax3,'MapProjection','robinson'); tightmap;
title(['[O_{2}] Trend at ' num2str(depth(depth_indices(1))) ...
    ' meters (\mumol kg^{-1} dec.^{-1})']);
ens_mean_trend = mean(o2_trend(:,:,depth_indices(1),:),4,'omitnan');
ens_mean_trend_sig = mean(o2_trend(:,:,depth_indices(1),:),4,'omitnan');
ens_mean_trend_non = mean(o2_trend(:,:,depth_indices(1),:),4,'omitnan');
sum_pos = sum(o2_trend(:,:,depth_indices(1),:) > 0,4);
sum_neg = sum(o2_trend(:,:,depth_indices(1),:) < 0,4);
idx = sum_pos >= 7 | sum_neg >= 7;
ens_mean_trend_sig(~idx) = NaN; ens_mean_trend_non(idx) = NaN;
h1=pcolorm(lat,[lon;lon(end)+1],[ens_mean_trend_sig;ens_mean_trend_sig(end,:)]');
%h2=pcolorm(lat,[lon;lon(end)+1],[ens_mean_trend_non;ens_mean_trend_non(end,:)]');
%h2.FaceAlpha = 0.5;
% contourm(lat,[lon;lon(end)+1],[ens_mean_trend;ens_mean_trend(end,:)]',...
%     trend_levels,'k','LineWidth',1,'ShowText','off');
stipplem(lat,lon,~idx' & ~isnan(ens_mean_trend)','density',300,'color',[.5 .5 .5]);
c = colorbar; clim(trend_lims); c.TickLength = 0;
colormap(ax3,flipud(cmocean('curl',length(trend_levels)-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax4 = nexttile; % trend at 250
worldmap(lat_lims,lon_lims);
set(ax4,'FontSize',12,'LineWidth',2);
setm(ax4,'MapProjection','robinson'); tightmap;
title(['[O_{2}] Trend at ' num2str(depth(depth_indices(2))) ...
    ' meters (\mumol kg^{-1} dec.^{-1})']);
ens_mean_trend = mean(o2_trend(:,:,depth_indices(2),:),4,'omitnan');
ens_mean_trend_sig = mean(o2_trend(:,:,depth_indices(2),:),4,'omitnan');
ens_mean_trend_non = mean(o2_trend(:,:,depth_indices(2),:),4,'omitnan');
sum_pos = sum(o2_trend(:,:,depth_indices(2),:) > 0,4);
sum_neg = sum(o2_trend(:,:,depth_indices(2),:) < 0,4);
idx = sum_pos >=7 | sum_neg >= 7;
ens_mean_trend_sig(~idx) = NaN; ens_mean_trend_non(idx) = NaN;
h1=pcolorm(lat,[lon;lon(end)+1],[ens_mean_trend_sig;ens_mean_trend_sig(end,:)]');
%h2=pcolorm(lat,[lon;lon(end)+1],[ens_mean_trend_non;ens_mean_trend_non(end,:)]');
%h2.FaceAlpha = 0.5;
% contourm(lat,[lon;lon(end)+1],[ens_mean_trend;ens_mean_trend(end,:)]',...
%     trend_levels,'k','LineWidth',1,'ShowText','off');
stipplem(lat,lon,~idx' & ~isnan(ens_mean_trend)','density',300,'color',[.5 .5 .5]);
c = colorbar; clim(trend_lims); c.TickLength = 0;
colormap(ax4,flipud(cmocean('curl',length(trend_levels)-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax5 = nexttile([1 2]); % timeseries
hold on; box on;
set(ax5,'FontSize',12,'LineWidth',2,'YAxisLocation','right');
plot(time,repmat(0,length(time),1),'k--');
for p = 1:length(products)
    if p <= 5
        plot(time,o2_inv_ts_10_1000(:,p)-mean(o2_inv_ts_10_1000(:,p),'omitnan'),...
            'linewidth',2,'LineStyle','-');
    else
        plot(time,o2_inv_ts_10_1000(:,p)-mean(o2_inv_ts_10_1000(:,p),'omitnan'),...
            'linewidth',2,'LineStyle','--');
    end
end
ylabel('O_{2} Inventory Anomaly (Pmol)');
xlim([min(time) max(time)]); datetick('x','keeplimits'); 
legend([{''} labels],'NumColumns',3,'Location','southwest');
hold off
% add annotations
annotation('textbox',[0.01 0.965 0 0],'string','(a)','FontSize',16);
annotation('textbox',[0.51 0.965 0 0],'string','(b)','FontSize',16);
annotation('textbox',[0.01 .645 0 0],'string','(c)','FontSize',16);
annotation('textbox',[0.51 .645 0 0],'string','(d)','FontSize',16);
annotation('textbox',[0.01 .36 0 0],'string','(e)','FontSize',16);
% save figure
figname = 'Figures/SotC_Fig2';
export_fig([figname '.png'],'-transparent');
figname = 'Figures/vectors/SotC_Fig2';
exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
close;

%% plot SoTC Figure 3
figure('Position',[100 100 1200 1000]);
SotC_fig = tiledlayout(3,2,'TileSpacing','tight','Padding','none');
ax1 = nexttile; % anomaly at 10
worldmap(lat_lims,lon_lims);
set(ax1,'FontSize',12,'LineWidth',2);
setm(ax1,'MapProjection','robinson'); tightmap;
title(['2025 [O_{2}] Anomaly at ' num2str(depth(depth_indices(1))) ...
    ' meters (\mumol kg^{-1})']);
ens_mean_anom = mean(o2_anom(:,:,1,:),4,'omitnan');
ens_mean_anom_sig = mean(o2_anom(:,:,1,:),4,'omitnan');
ens_mean_anom_non = mean(o2_anom(:,:,1,:),4,'omitnan');
sum_pos = sum(o2_anom(:,:,1,:) > 0,4);
sum_neg = sum(o2_anom(:,:,1,:) < 0,4);
idx = sum_pos >= 2 | sum_neg >= 2;
ens_mean_anom_sig(~idx) = NaN; ens_mean_anom_non(idx) = NaN;
h1=pcolorm(lat,[lon;lon(end)+1],[ens_mean_anom_sig;ens_mean_anom_sig(end,:)]');
%h2=pcolorm(lat,[lon;lon(end)+1],[ens_mean_anom_non;ens_mean_anom_non(end,:)]');
%h2.FaceAlpha = 0.5;
stipplem(lat,lon,~idx' & ~isnan(ens_mean_anom)','density',300,'color',[.5 .5 .5]);
c = colorbar; clim(anom_lims); c.TickLength = 0;
colormap(ax1,flipud(cmocean('curl',length(anom_levels)-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax2 = nexttile; % anomaly at 250
worldmap(lat_lims,lon_lims);
set(ax2,'FontSize',12,'LineWidth',2);
setm(ax2,'MapProjection','robinson'); tightmap;
title(['2025 [O_{2}] Anomaly at ' num2str(depth(depth_indices(2))) ...
    ' meters (\mumol kg^{-1})']);
ens_mean_anom = mean(o2_anom(:,:,2,:),4,'omitnan');
ens_mean_anom_sig = mean(o2_anom(:,:,2,:),4,'omitnan');
ens_mean_anom_non = mean(o2_anom(:,:,2,:),4,'omitnan');
sum_pos = sum(o2_anom(:,:,2,:) > 0,4);
sum_neg = sum(o2_anom(:,:,2,:) < 0,4);
idx = sum_pos >= 2 | sum_neg >= 2;
ens_mean_anom_sig(~idx) = NaN; ens_mean_anom_non(idx) = NaN;
h1=pcolorm(lat,[lon;lon(end)+1],[ens_mean_anom_sig;ens_mean_anom_sig(end,:)]');
%h2=pcolorm(lat,[lon;lon(end)+1],[ens_mean_anom_non;ens_mean_anom_non(end,:)]');
%h2.FaceAlpha = 0.5;
stipplem(lat,lon,~idx' & ~isnan(ens_mean_anom)','density',300,'color',[.5 .5 .5]);
c = colorbar; clim(anom_lims); c.TickLength = 0;
colormap(ax2,flipud(cmocean('curl',length(anom_levels)-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax3 = nexttile; % temp anomaly at 10
worldmap(lat_lims,lon_lims);
set(ax3,'FontSize',12,'LineWidth',2);
setm(ax3,'MapProjection','robinson'); tightmap;
title(['2025 EN4.2.2 T Anomaly at ' num2str(depth(depth_indices(1))) ...
    ' meters (' char(176) ' C)']);
pcolorm(lat-.25,[lon;lon(end)+1]-.25,[temp_anom(:,:,1);temp_anom(end,:,1)]');
% contourm(lat-.25,[lon;lon(end)+1]-.25,[temp_anom(:,:,1);temp_anom(end,:,1)]',...
%     temp_anom_levels(3:3:end),'k','LineWidth',1,'ShowText','off');
c = colorbar; clim(temp_anom_lims); c.TickLength = 0;
colormap(ax3,cmocean('balance',length(temp_anom_levels)-1,'pivot',0));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax4 = nexttile; % temp anomaly at 250
worldmap(lat_lims,lon_lims);
set(ax4,'FontSize',12,'LineWidth',2);
setm(ax4,'MapProjection','robinson'); tightmap;
title(['2025 EN4.2.2 T Anomaly at ' num2str(depth(depth_indices(2))) ...
    ' meters (' char(176) ' C)']);
pcolorm(lat-.25,[lon;lon(end)+1]-.25,[temp_anom(:,:,2);temp_anom(end,:,2)]');
% contourm(lat-.25,[lon;lon(end)+1]-.25,[temp_anom(:,:,2);temp_anom(end,:,2)]',...
%     temp_anom_levels(3:3:end),'k','LineWidth',1,'ShowText','off');
c = colorbar; clim(temp_anom_lims); c.TickLength = 0;
colormap(ax4,cmocean('balance',length(temp_anom_levels)-1,'pivot',0));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
% calculate seasonal cycles
mean_o2_seas_baseline_N = nan(length(products),length(depth_indices),12);
mean_o2_seas_anom_N = nan(length(products),length(depth_indices),12);
mean_o2_seas_baseline_var_N = nan(length(products),length(depth_indices),12);
for p = 6:7
    % 10m
    o2_seas_baseline = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_baseline_seas'],[1 1 d1_idx 1],[Inf Inf 1 Inf]));
    o2_seas_anom = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_anom_seas'],[1 1 d1_idx 1],[Inf Inf 1 Inf]));
    o2_seas_baseline_var = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_baseline_seas_var'],[1 1 d2_idx 1],[Inf Inf 1 Inf]));
    idx = all(~isnan(o2_seas_baseline),3) & all(~isnan(o2_seas_baseline),3);
    [~,area,~] = weights3d(lon,lat,depth); area = area(:,:,1);
    idx = all(~isnan(o2_seas_baseline),3) & all(~isnan(o2_seas_baseline),3);
    o2_seas_baseline(~idx) = nan; o2_seas_anom(~idx) = nan; area(~idx) = nan;
    idx_lat = lat >= 0;
    mean_o2_seas_baseline_N(p,1,:) = ...
        (sum(reshape(o2_seas_baseline(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
    mean_o2_seas_anom_N(p,1,:) = ...
        (sum(reshape(o2_seas_anom(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
    mean_o2_seas_baseline_var_N(p,1,:) = ...
        (sum(reshape(o2_seas_baseline_var(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
    % 250m
    o2_seas_baseline = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_baseline_seas'],[1 1 d2_idx 1],[Inf Inf 1 Inf]));
    o2_seas_anom = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_anom_seas'],[1 1 d2_idx 1],[Inf Inf 1 Inf]));
    o2_seas_baseline_var = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_baseline_seas_var'],[1 1 d2_idx 1],[Inf Inf 1 Inf]));
    [~,area,~] = weights3d(lon,lat,depth); area = area(:,:,1);
    idx = all(~isnan(o2_seas_baseline),3) & all(~isnan(o2_seas_baseline),3);
    o2_seas_baseline(~idx) = nan; o2_seas_anom(~idx) = nan; area(~idx) = nan;
    mean_o2_seas_baseline_N(p,2,:) = ...
        (sum(reshape(o2_seas_baseline(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
    mean_o2_seas_anom_N(p,2,:) = ...
        (sum(reshape(o2_seas_anom(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
    mean_o2_seas_baseline_var_N(p,2,:) = ...
        (sum(reshape(o2_seas_baseline_var(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
end
ax5 = nexttile; % seasonal cycle at 10m
hold on; box on; set(gca,'FontSize',12,'LineWidth',2,'YAxisLocation','right');
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,1,:),1,'omitnan')),'LineWidth',3);
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,1,:)+mean_o2_seas_baseline_var_N(:,1,:),1,'omitnan')),'k--','LineWidth',1);
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,1,:)-mean_o2_seas_baseline_var_N(:,1,:),1,'omitnan')),'k--','LineWidth',1);
plot(1:11,squeeze(mean(mean_o2_seas_anom_N(:,1,1:11),1,'omitnan')),'LineWidth',3);
title('N. Hem. Average [O_{2}] at 10 meters (\mumol kg^{-1})');
xlim([0 13]); xticks(1:12);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
legend({'1990-2020 Mean' '' '' '2025'},'NumColumns',1);
hold off
ax6 = nexttile; % seasonal cycle at 250m
hold on; box on; set(gca,'FontSize',12,'LineWidth',2,'YAxisLocation','right');
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,2,:),1,'omitnan')),'LineWidth',3);
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,2,:)+mean_o2_seas_baseline_var_N(:,2,:),1,'omitnan')),'k--','LineWidth',1);
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,2,:)-mean_o2_seas_baseline_var_N(:,2,:),1,'omitnan')),'k--','LineWidth',1);
plot(1:11,squeeze(mean(mean_o2_seas_anom_N(:,2,1:11),1,'omitnan')),'LineWidth',3);
title('N. Hem. Average [O_{2}] at 250 meters (\mumol kg^{-1})');
xlim([0 13]); xticks(1:12); ylim([90 110]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
legend({'1990-2020 Mean' '' '' '2025'},'NumColumns',1);
hold off
% add annotations
annotation('textbox',[0.01 0.965 0 0],'string','(a)','FontSize',16);
annotation('textbox',[0.51 0.965 0 0],'string','(b)','FontSize',16);
annotation('textbox',[0.01 .64 0 0],'string','(c)','FontSize',16);
annotation('textbox',[0.51 .64 0 0],'string','(d)','FontSize',16);
annotation('textbox',[0.01 .3 0 0],'string','(e)','FontSize',16);
annotation('textbox',[0.51 .3 0 0],'string','(f)','FontSize',16);
% save figure
figname = 'Figures/SotC_Fig3';
export_fig([figname '.png'],'-transparent');
figname = 'Figures/vectors/SotC_Fig3';
exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
close;

%% plot standard deviation among OI products
% o2_inv_std_OI_0_1000 = std(o2_inv_mean_10_1000(:,:,1:5),[],3,'omitnan');
% figure; worldmap(lat_lims,lon_lims);
% pcolorm(lat,lon,o2_inv_std_OI_0_1000');
% colorbar;

%% plot standard deviation among ML products
% o2_inv_std_ML_0_1000 = std(o2_inv_mean_10_1000(:,:,6:9),[],3,'omitnan');
% figure; worldmap(lat_lims,lon_lims);
% pcolorm(lat,lon,o2_inv_std_ML_0_1000');
% colorbar;

%% plot SoTC Figure 3
figure('Position',[100 100 1200 1000]);
set(gcf,'Color','w');
[lat3d,lon3d] = ndgrid(lat,lon);
SotC_fig = tiledlayout(3,2,'TileSpacing','tight','Padding','none');
ax1 = nexttile; % anomaly at 10
worldmap(lat_lims,lon_lims);
set(ax1,'FontSize',12,'LineWidth',2);
setm(ax1,'MapProjection','robinson'); tightmap;
title(['2025 [O_{2}] Anomaly at ' num2str(depth(depth_indices(1))) ...
    ' meters (\mumol kg^{-1})']);
black_idx = double(isnan(mean(o2_anom(:,:,1,:),4,'omitnan')));
black_idx(black_idx==0) = NaN;
pcolorm(lat,[lon;lon(end)+1],double([black_idx;black_idx(end,:)])');
c = colorbar; c.Visible = 'off';
colormap(ax1,[0 0 0]); %caxis(ax1,[1 1]);
ax11 = axes('Position',ax1.Position,'color','none','visible','off');
worldmap(lat_lims,lon_lims);
setm(ax11,'MapProjection','robinson'); tightmap;
ens_mean_anom_sig = mean(o2_anom(:,:,1,:),4,'omitnan');
ens_mean_anom_non = mean(o2_anom(:,:,1,:),4,'omitnan');
sum_pos = sum(o2_anom(:,:,1,:) > 0,4);
sum_neg = sum(o2_anom(:,:,1,:) < 0,4);
idx = sum_pos >= 2 | sum_neg >= 2;
ens_mean_anom_sig(~idx) = NaN; ens_mean_anom_non(idx) = NaN;
pcolorm(lat,[lon;lon(end)+1],[ens_mean_anom_sig;ens_mean_anom_sig(end,:)]');
%h2=pcolorm(lat,[lon;lon(end)+1],[ens_mean_anom_non;ens_mean_anom_non(end,:)]');
%h2.FaceAlpha = 0.5;
%stipplem(lat,lon,~idx,'density',300);
c = colorbar; clim(anom_lims); c.TickLength = 0;
colormap(ax11,cmocean('diff',length(anom_levels)-1,'pivot',0));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax2 = nexttile; % anomaly at 250
worldmap(lat_lims,lon_lims);
set(ax2,'FontSize',12,'LineWidth',2);
setm(ax2,'MapProjection','robinson'); tightmap;
title(['2025 [O_{2}] Anomaly at ' num2str(depth(depth_indices(2))) ...
    ' meters (\mumol kg^{-1})']);
ens_mean_anom_sig = mean(o2_anom(:,:,2,:),4,'omitnan');
ens_mean_anom_non = mean(o2_anom(:,:,2,:),4,'omitnan');
sum_pos = sum(o2_anom(:,:,2,:) > 0,4);
sum_neg = sum(o2_anom(:,:,2,:) < 0,4);
idx = sum_pos >= 2 | sum_neg >= 2;
ens_mean_anom_sig(~idx) = NaN; ens_mean_anom_non(idx) = NaN;
h1=pcolorm(lat,[lon;lon(end)+1],[ens_mean_anom_sig;ens_mean_anom_sig(end,:)]');
%h2=pcolorm(lat,[lon;lon(end)+1],[ens_mean_anom_non;ens_mean_anom_non(end,:)]');
%h2.FaceAlpha = 0.5;
c = colorbar; clim(anom_lims); c.TickLength = 0;
colormap(ax2,cmocean('diff',length(anom_levels)-1,'pivot',0));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax3 = nexttile; % temp anomaly at 10
worldmap(lat_lims,lon_lims);
set(ax3,'FontSize',12,'LineWidth',2);
setm(ax3,'MapProjection','robinson'); tightmap;
title(['2025 EN4.2.2 T Anomaly at ' num2str(depth(depth_indices(1))) ...
    ' meters (' char(176) ' C)']);
pcolorm(lat-.25,[lon;lon(end)+1]-.25,[temp_anom(:,:,1);temp_anom(end,:,1)]');
% contourm(lat-.25,[lon;lon(end)+1]-.25,[temp_anom(:,:,1);temp_anom(end,:,1)]',...
%     temp_anom_levels(3:3:end),'k','LineWidth',1,'ShowText','off');
c = colorbar; clim(temp_anom_lims); c.TickLength = 0;
colormap(ax3,cmocean('balance',length(temp_anom_levels)-1,'pivot',0));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax4 = nexttile; % temp anomaly at 250
worldmap(lat_lims,lon_lims);
set(ax4,'FontSize',12,'LineWidth',2);
setm(ax4,'MapProjection','robinson'); tightmap;
title(['2025 EN4.2.2 T Anomaly at ' num2str(depth(depth_indices(2))) ...
    ' meters (' char(176) ' C)']);
pcolorm(lat-.25,[lon;lon(end)+1]-.25,[temp_anom(:,:,2);temp_anom(end,:,2)]');
% contourm(lat-.25,[lon;lon(end)+1]-.25,[temp_anom(:,:,2);temp_anom(end,:,2)]',...
%     temp_anom_levels(3:3:end),'k','LineWidth',1,'ShowText','off');
c = colorbar; clim(temp_anom_lims); c.TickLength = 0;
colormap(ax4,cmocean('balance',length(temp_anom_levels)-1,'pivot',0));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
% calculate seasonal cycles
mean_o2_seas_baseline_N = nan(length(products),length(depth_indices),12);
mean_o2_seas_anom_N = nan(length(products),length(depth_indices),12);
mean_o2_seas_baseline_var_N = nan(length(products),length(depth_indices),12);
for p = 6:7
    % 10m
    o2_seas_baseline = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_baseline_seas'],[1 1 d1_idx 1],[Inf Inf 1 Inf]));
    o2_seas_anom = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_anom_seas'],[1 1 d1_idx 1],[Inf Inf 1 Inf]));
    o2_seas_baseline_var = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_baseline_seas_var'],[1 1 d2_idx 1],[Inf Inf 1 Inf]));
    idx = all(~isnan(o2_seas_baseline),3) & all(~isnan(o2_seas_baseline),3);
    [~,area,~] = weights3d(lon,lat,depth); area = area(:,:,1);
    idx = all(~isnan(o2_seas_baseline),3) & all(~isnan(o2_seas_baseline),3);
    o2_seas_baseline(~idx) = nan; o2_seas_anom(~idx) = nan; area(~idx) = nan;
    idx_lat = lat >= 0;
    mean_o2_seas_baseline_N(p,1,:) = ...
        (sum(reshape(o2_seas_baseline(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
    mean_o2_seas_anom_N(p,1,:) = ...
        (sum(reshape(o2_seas_anom(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
    mean_o2_seas_baseline_var_N(p,1,:) = ...
        (sum(reshape(o2_seas_baseline_var(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
    % 250m
    o2_seas_baseline = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_baseline_seas'],[1 1 d2_idx 1],[Inf Inf 1 Inf]));
    o2_seas_anom = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_anom_seas'],[1 1 d2_idx 1],[Inf Inf 1 Inf]));
    o2_seas_baseline_var = squeeze(ncread([fpath 'O2_Maps/' compilation.name],...
        [products{p} '_baseline_seas_var'],[1 1 d2_idx 1],[Inf Inf 1 Inf]));
    [~,area,~] = weights3d(lon,lat,depth); area = area(:,:,1);
    idx = all(~isnan(o2_seas_baseline),3) & all(~isnan(o2_seas_baseline),3);
    o2_seas_baseline(~idx) = nan; o2_seas_anom(~idx) = nan; area(~idx) = nan;
    mean_o2_seas_baseline_N(p,2,:) = ...
        (sum(reshape(o2_seas_baseline(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
    mean_o2_seas_anom_N(p,2,:) = ...
        (sum(reshape(o2_seas_anom(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
    mean_o2_seas_baseline_var_N(p,2,:) = ...
        (sum(reshape(o2_seas_baseline_var(:,idx_lat,:),[length(lon)*length(lat(idx_lat)) 12]).*...
        reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan')./...
        sum(reshape(area(:,idx_lat),[length(lon)*length(lat(idx_lat)) 1]),'omitnan'));
end
ax5 = nexttile; % seasonal cycle at 10m
hold on; box on; set(gca,'FontSize',12,'LineWidth',2,'YAxisLocation','right');
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,1,:),1,'omitnan')),'LineWidth',3);
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,1,:)+mean_o2_seas_baseline_var_N(:,1,:),1,'omitnan')),'k--','LineWidth',1);
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,1,:)-mean_o2_seas_baseline_var_N(:,1,:),1,'omitnan')),'k--','LineWidth',1);
plot(1:11,squeeze(mean(mean_o2_seas_anom_N(:,1,1:11),1,'omitnan')),'LineWidth',3);
title('N. Hem. Average [O_{2}] at 10 meters (\mumol kg^{-1})');
xlim([0 13]); xticks(1:12);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
legend({'1990-2020 Mean' '' '' '2025'},'NumColumns',1);
hold off
ax6 = nexttile; % seasonal cycle at 250m
hold on; box on; set(gca,'FontSize',12,'LineWidth',2,'YAxisLocation','right');
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,2,:),1,'omitnan')),'LineWidth',3);
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,2,:)+mean_o2_seas_baseline_var_N(:,2,:),1,'omitnan')),'k--','LineWidth',1);
plot(1:12,squeeze(mean(mean_o2_seas_baseline_N(:,2,:)-mean_o2_seas_baseline_var_N(:,2,:),1,'omitnan')),'k--','LineWidth',1);
plot(1:11,squeeze(mean(mean_o2_seas_anom_N(:,2,1:11),1,'omitnan')),'LineWidth',3);
title('N. Hem. Average [O_{2}] at 250 meters (\mumol kg^{-1})');
xlim([0 13]); xticks(1:12); ylim([90 110]);
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});
legend({'1990-2020 Mean' '' '' '2025'},'NumColumns',1);
hold off
% add annotations
annotation('textbox',[0.01 0.965 0 0],'string','(a)','FontSize',16);
annotation('textbox',[0.51 0.965 0 0],'string','(b)','FontSize',16);
annotation('textbox',[0.01 .64 0 0],'string','(c)','FontSize',16);
annotation('textbox',[0.51 .64 0 0],'string','(d)','FontSize',16);
annotation('textbox',[0.01 .3 0 0],'string','(e)','FontSize',16);
annotation('textbox',[0.51 .3 0 0],'string','(f)','FontSize',16);
% save figure
figname = 'Figures/SotC_Fig3';
export_fig([figname '.png'],'-transparent');
figname = 'Figures/vectors/SotC_Fig3';
exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
close;