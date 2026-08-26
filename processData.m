% bin-aggregated data: 1.9 million data points
% raw discrete data: 13.4 million data points
% vertically interpolated discrete data: 7.6 million data points

%% for GLODAP only
% glodap_idx = OXY_table.platform < 10^6;
% OXY_table = OXY_table(glodap_idx,:);
% clear glodap_idx

%% split data into training/test sets
% rng(5); % for reproducibility
% test_per_oxy = 0.2; % percentage to use as test data
% plat_ids_oxy = unique(OXY_table.platform);
% num_plats_oxy = length(plat_ids_oxy);
% split_idx_oxy = randperm(num_plats_oxy) <= test_per_oxy.*num_plats_oxy;
% test_plat_oxy = plat_ids_oxy(split_idx_oxy);
% train_plat_oxy = plat_ids_oxy(~split_idx_oxy);
% test_idx_oxy = ismember(OXY_table.platform,test_plat_oxy);
% train_idx_oxy = ismember(OXY_table.platform,train_plat_oxy);
% clear plat_ids_oxy num_plats_oxy split_idx_oxy test_plat_oxy
% clear test_per_oxy train_plat_oxy train_plat_oxy

%% add ancillary data to table
% calculate date
date = datevec(double(OXY_table.time));
% extract year
OXY_table.year = date(:,1);
% extract day of year
date0 = date; date0(:,2:3) = 1;
OXY_table.day = datenum(date) - datenum(date0);
clear date date0
% calculate sine and cosine of day
OXY_table.day_sin = sin((2.*pi.*OXY_table.day)./365.25);
OXY_table.day_cos = cos((2.*pi.*OXY_table.day)./365.25);
% calculate cosines of longitude
OXY_table.lon_cos1 = cosd(OXY_table.longitude-20);
OXY_table.lon_cos2 = cosd(OXY_table.longitude-110);
% calculate log10 of chlorophyll
%OXY_table.log10_chl = log10(OXY_table.chlorophyll);
% determine bottom depth
OXY_table.bottom_depth = bottom_depth(OXY_table.latitude,OXY_table.longitude);

%% assemble data into table
OXY.all = OXY_table;
clear OXY_table train_idx_oxy test_idx_oxy

% re-order table
OXY.all = [OXY.all(:,1:3) OXY.all(:,16:17) OXY.all(:,4:6) ...
           OXY.all(:,12:15) OXY.all(:,7:10) OXY.all(:,18) ...
           OXY.all(:,11)];

% artificially thin data
data_per = 0.01;
idx = randperm(height(OXY.all));
OXY.all(idx >= height(OXY.all)*data_per,:) = [];
clear idx

% Predictor table
% 1 - platform
% 2 - latitude
% 3 - longitude
% 4 - cos1(lon)
% 5 - cos2(lon)
% 6 - sigma
% 7 - pressure
% 8 - time
% 9 - year
% 10 - day
% 11 - sin(day)
% 12 - cos(day)
% 13 - temperature
% 14 - cons. temperature
% 15 - salinity
% 16 - abs. salinity
%%%%%%%%%%%% 17 - log(chl)
% 17 - bottom depth
% 18 - oxygen

% Select predictors for oxygen
preds_oxy_RF = [2 4 5 6 7 9 11 12 14 16 17];
all_oxy_RF = [preds_oxy_RF 18];
preds_oxy_NN = [2 4 5 6 7 9 11 12 14 16 17];
all_oxy_NN = [preds_oxy_RF 18];

%% Calculate gridding uncertainty
% % establish edges of bins
% x_edges = -180:180;
% x_bins = -179.5:179.5;
% y_edges = -65:80;
% y_bins = -64.5:79.5;
% z_edges = [0 5:10:175 190:20:450 475:50:1375 1450:100:1950 2000];
% z_bins = [2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975];
% t_edges = [0 31 60 91 121 152 182 213 244 274 305 335 366];
% % get histogram counts in each bin
% [~,~,Xnum] = histcounts(OXY.all.longitude,x_edges);
% [~,~,Ynum] = histcounts(OXY.all.latitude,y_edges);
% [~,~,Znum] = histcounts(OXY.all.pressure,z_edges);
% [~,~,Tnum] = histcounts(OXY.all.day,t_edges);
% % accumulate index of counts
% subs = [Xnum,Ynum,Znum,Tnum];
% idx = ~any(subs==0,2);
% % determine size of 4D grid
% sz = [length(x_edges)-1,length(y_edges)-1,length(z_edges)-1,length(t_edges)-1];
% % get variability within each grid cell
% oxygen_std = accumarray(subs(idx,:),OXY.all.oxygen(idx),sz,@nanstd,nan);
% oxygen_pres = accumarray(subs(idx,:),OXY.all.pressure(idx),sz,@nanmean,nan);
% oxygen_lat = accumarray(subs(idx,:),OXY.all.latitude(idx),sz,@nanmean,nan);
% oxygen_lon = accumarray(subs(idx,:),OXY.all.longitude(idx),sz,@nanmean,nan);
% oxygen_temp = accumarray(subs(idx,:),OXY.all.temperature(idx),sz,@nanmean,nan);
% oxygen_sal = accumarray(subs(idx,:),OXY.all.salinity(idx),sz,@nanmean,nan);
% oxygen_sig = accumarray(subs(idx,:),OXY.all.sigma(idx),sz,@nanmean,nan);
% % remove standard deviations calculated from nine or fewer measurements
% oxygen_count = accumarray(subs(idx,:),1);
% oxygen_count_idx = oxygen_count < 10;
% oxygen_std(oxygen_count_idx) = nan;
% oxygen_pres(oxygen_count_idx) = nan;
% oxygen_lat(oxygen_count_idx) = nan;
% oxygen_lon(oxygen_count_idx) = nan;
% oxygen_sig(oxygen_count_idx) = nan;
% % calculate distance from shore
% oxygen_dist = dist2coast(oxygen_lat,oxygen_lon);
% % scatter different parameters against std
% % figure; scatter(oxygen_lon(:),oxygen_std(:));
% % fit model of variability vs depth, sigma, and distance from shore
% [b,~,~,~,stats] = ...
%     regress(oxygen_std(:),[ones(size(oxygen_pres(:))) ...
%             oxygen_pres(:) oxygen_pres(:).^2 oxygen_sig(:) ...
%             oxygen_sig(:).^2 oxygen_dist(:) oxygen_dist(:).^2]);
% stats(1) % R^2 statistic (0.0367)
% save('std_coefs','b')
% % plot figure of gridding uncertainty versus n
% grid_vs_n_disc
% % clean up
% clear oxygen_count oxygen_dist oxygen_lat oxygen_lon
% clear oxygen_pres oxygen_sig oxygen_bot oxygen_std stats b
% 
% %% Plot all data by months of year represented (do this if not binning)
% figure; hold on;
% worldmap([-90 90],[20 380]);
% setm(gca,'mapprojection','robinson');
% set(gcf,'units','inches','position',[0 5 20 10]);
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',16);
% x_bins = -180:180; y_bins = -90:90;
% [longitude,latitude] = ndgrid(-179.5:179.5,-89.5:89.5);
% if exist('float_data','var')
%     profiles = [float_data.OXY_PROF_ID;glodap_data.OXY_ID];
% else
%     profiles = [glodap_data.OXY_ID];
% end
% profiles_test = profiles(idx);
% [~,~,Xnum] = histcounts(convert_lon(OXY.all.longitude(idx)),x_bins);
% [~,~,Ynum] = histcounts(OXY.all.latitude(idx),y_bins);
% subs = [Xnum,Ynum]; sz = size(longitude);
% prof_counts = accumarray(subs,profiles_test,sz,@(x)length(unique(x)),0);
% pcolorm(double(latitude),double(longitude),prof_counts);
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar('location','northoutside');
% c.Label.String = 'Profile Count';
% c.FontSize = 22;
% caxis([-0.5 20.5]);
% cm=flipud(bone(21)); cm(1,:)=1; colormap(cm);
% c.Ticks = [0 5 10 15 20];
% c.TickLabels = {'0' '5' '10' '15' '20+'};
% c.TickLength = 0;
% mlabel off; plabel off;
% if ~exist('Figures','dir'); mkdir('Figures'); end
% if ~exist('Figures/Validation','dir'); mkdir('Figures/Validation'); end
% exportgraphics(gcf,'Figures/Validation/Data_all_map_binned.jpg');
% clear cm idx latitude longitude prof_counts profiles profiles_test
% close

%% Create basin shapes
basins = {'Atl' 'Pac' 'Ind' 'Arc' 'Med' 'NSou' 'SSou'};
lat = -89.5:89.5;
lon = -179.5:180.5;
[lat,lon] = meshgrid(lat,lon);
basin_shape = zeros(size(lat));
land = shaperead('landareas', 'UseGeoCoords', true);

for b = 1:7
    if b == 2
        basin_shape(basin_shape ~= 0 &...
            (inpolygon(lon,lat,Poly.([basins{b} '1'])(:,1),...
            Poly.([basins{b} '1'])(:,2)) | ...
            inpolygon(lon,lat,Poly.([basins{b} '2'])(:,1),...
            Poly.([basins{b} '2'])(:,2)))) = 8.5;
        basin_shape(basin_shape == 0 &...
            (inpolygon(lon,lat,Poly.([basins{b} '1'])(:,1),...
            Poly.([basins{b} '1'])(:,2)) | ...
            inpolygon(lon,lat,Poly.([basins{b} '2'])(:,1),...
            Poly.([basins{b} '2'])(:,2)))) = b+0.5;
    else
        basin_shape(basin_shape ~= 0 &...
            inpolygon(lon,lat,Poly.(basins{b})(:,1),...
            Poly.(basins{b})(:,2))) = 8.5;
        basin_shape(basin_shape == 0 &...
            inpolygon(lon,lat,Poly.(basins{b})(:,1),...
            Poly.(basins{b})(:,2))) = b+0.5;
    end
end
for n = 1:length(land)
    basin_shape(inpolygon(lon,lat,land(n).Lon,land(n).Lat)) = 0.5;
end
clear n land b basins

%% Plot data distribution by type of platform
figure; worldmap([-90 90],[20 380]); hold on;
pcolorm(lat,lon,basin_shape);
colormap([1,1,1;cbrewer('qual','Pastel2',8)]);
caxis([0.5 9.5]);
c=colorbar;
set(c,'YDir','reverse','Limits',[0.5 8.5]);
c.TickLabels = {'Atl.','Pac.','Ind.','Arc.','Med.','N. Sou.','S. Sou','Ovrlp.'};
c.FontSize = 20;
c.TickLength = 0;
set(gcf,'units','inches','position',[0 5 20 10]);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',14);
clrs = cmocean('thermal',7);
p1=scatterm(wod.lat(wod.type==3),wod.lon(wod.type==3),5,clrs(2,:),'filled');
p2=scatterm(wod.lat(wod.type==1),wod.lon(wod.type==1),5,clrs(4,:),'filled');
p3=scatterm(wod.lat(wod.type==2),wod.lon(wod.type==2),5,clrs(6,:),'filled');
[~,icons]=legend([p1 p2 p3],{'Argo' 'OSD' 'CTD'},...
    'location','northoutside','Orientation','horizontal','fontsize',24);
icons(4).Children.MarkerSize = 10;
icons(5).Children.MarkerSize = 10;
icons(6).Children.MarkerSize = 10;
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',[.7 .7 .7]);
mlabel off; plabel off;
if ~exist('Figures','dir'); mkdir('Figures'); end
if ~exist('Figures/Data','dir'); mkdir('Figures/Data'); end
exportgraphics(gcf,'Figures/Data/data_locations.png');
close

%% Plot all data distribution (do this either way)
figure; worldmap([-90 90],[20 380]); hold on;
set(gcf,'units','inches','position',[0 5 20 10]);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',14);
clrs = cmocean('thermal',7);
p1=scatterm(wod.lat(wod.type==3),wod.lon(wod.type==3),5,clrs(2,:),'filled');
p2=scatterm(wod.lat(wod.type==1),wod.lon(wod.type==1),5,clrs(4,:),'filled');
p3=scatterm(wod.lat(wod.type==2),wod.lon(wod.type==2),5,clrs(6,:),'filled');
[~,icons]=legend([p1 p2 p3],{'Argo' 'OSD' 'CTD'},...
    'location','northoutside','Orientation','horizontal','fontsize',24);
icons(4).Children.MarkerSize = 10;
icons(5).Children.MarkerSize = 10;
icons(6).Children.MarkerSize = 10;
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',[.7 .7 .7]);
mlabel off; plabel off;
exportgraphics(gcf,'Figures/Data/data_locations_all.jpg');
close

%% Clean up and save data
% if ~isfolder('Data'); mkdir('Data'); end
% save(['Data/Data-' date],'-v7.3','OXY','all_oxy_NN','all_oxy_RF','preds_oxy_NN','preds_oxy_RF');
