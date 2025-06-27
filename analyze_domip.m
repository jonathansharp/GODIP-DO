% analyze DOMIP

%% define filenames
gobai_filename = '/fast4/o2/domip/GOBAI/EN4/FFNN/c15_Jan-2025_D/train80_val10_test10/gobai-o2.nc';
gt_filename = '/fast4/o2/domip/O2map_v2.2.G.4.6.4.nc';
woa_filename = '/fast4/o2/woa23_all_o00_01.nc';

%% import mean climatologies
% load gobai dimensions
gobai.lon = ncread(gobai_filename,'lon');
gobai.lon = convert_lon(gobai.lon,'-180-180');
gobai.lat = ncread(gobai_filename,'lat');
gobai.depth = ncread(gobai_filename,'depth');
gobai.time = ncread(gobai_filename,'time');

% load gt dimensions
gt.lon = ncread(gt_filename,'lon');
gt.lat = ncread(gt_filename,'lat');
gt.depth = ncread(gt_filename,'depth');
gt.time = ncread(gt_filename,'time');

% load woa dimensions
woa.lon = ncread(woa_filename,'lon');
woa.lat = ncread(woa_filename,'lat');
woa.depth = ncread(woa_filename,'depth');
woa.o2 = ncread(woa_filename,'o_an');

% calculate long-term gobai mean
gobai.o2 = ncread(gobai_filename,'o2');
gobai.mean_o2 = mean(gobai.o2,4,'omitnan');
gobai = rmfield(gobai,'o2');

% calculate long-term gt mean
gt.o2 = ncread(gt_filename,'o2');
gt.mean_o2 = mean(gt.o2,4,'omitnan');
gt = rmfield(gt,'o2');

% create 3d dimensions
[gobai.lon3d,gobai.lat3d,gobai.depth3d] = ...
    ndgrid(gobai.lon,gobai.lat,gobai.depth(1:30));
[gt.lon3d,gt.lat3d,gt.depth3d] = ...
    ndgrid(gt.lon,gt.lat,gt.depth);

% interpolate gobai mean to woa/gt depth levels
gobai.mean_o2_interp = griddata(gobai.lon3d,gobai.lat3d,gobai.depth3d,...
    gobai.mean_o2(:,:,1:30),gt.lon3d,gt.lat3d,gt.depth3d,'linear');

%% plot long-term means
xgap = 0.00;
ygap = 0.03;
h_marg = [0 0];
w_marg = [0 0.1];
f=figure('Position',[100 100 1200 900]);
% 10m
h=subtightplot(4,3,1,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,3)',[0 350],'WOA, 10m',parula);
h=subtightplot(4,3,2,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gt.mean_o2(:,:,3)',[0 350],'GT, 10m',parula);
h=subtightplot(4,3,3,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gobai.mean_o2_interp(:,:,3)',[0 350],'GOBAI, 10m',parula);
% 100m
h=subtightplot(4,3,4,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,21)',[0 350],'WOA, 100m',parula)
h=subtightplot(4,3,5,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gt.mean_o2(:,:,21)',[0 350],'GT, 100m',parula)
h=subtightplot(4,3,6,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gobai.mean_o2_interp(:,:,21)',[0 350],'GOBAI, 100m',parula)
% 700m
h=subtightplot(4,3,7,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,41)',[0 350],'WOA, 700m',parula)
h=subtightplot(4,3,8,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gt.mean_o2(:,:,41)',[0 350],'GT, 700m',parula)
h=subtightplot(4,3,9,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gobai.mean_o2_interp(:,:,41)',[0 350],'GOBAI, 700m',parula)
% 1700m
h=subtightplot(4,3,10,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,61)',[0 350],'WOA, 1700m',parula)
h=subtightplot(4,3,11,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gt.mean_o2(:,:,61)',[0 350],'GT, 1700m',parula)
h=subtightplot(4,3,12,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gobai.mean_o2_interp(:,:,61)',[0 350],'GOBAI, 1700m',parula)
% colorbar
h=axes(f,'visible','off');
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,parula);
clim(h,[0 350]);
% export
export_fig(gcf,'O2/Figures/mean_comparisons.png','-transparent');
close

%% plot differences
xgap = 0.00;
ygap = 0.03;
h_marg = [0 0];
w_marg = [0 0.1];
f=figure('Position',[100 100 1200 900]);
% 10m
h=subtightplot(4,3,1,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,3)'-gt.mean_o2(:,:,3)',[-50 50],'WOA-GT, 10m',cmocean('balance'));
h=subtightplot(4,3,2,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,3)'-gobai.mean_o2_interp(:,:,3)',[-50 50],'WOA-GOBAI, 10m',cmocean('balance'));
h=subtightplot(4,3,3,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gobai.mean_o2_interp(:,:,3)'-gt.mean_o2(:,:,3)',[-50 50],'GOBAI-GT, 10m',cmocean('balance'));
% 100m
h=subtightplot(4,3,4,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,21)'-gt.mean_o2(:,:,21)',[-50 50],'WOA-GT, 100m',cmocean('balance'));
h=subtightplot(4,3,5,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,21)'-gobai.mean_o2_interp(:,:,21)',[-50 50],'WOA-GOBAI, 100m',cmocean('balance'));
h=subtightplot(4,3,6,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gobai.mean_o2_interp(:,:,21)'-gt.mean_o2(:,:,21)',[-50 50],'GOBAI-GT, 100m',cmocean('balance'));
% 700m
h=subtightplot(4,3,7,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,41)'-gt.mean_o2(:,:,41)',[-50 50],'WOA-GT, 700m',cmocean('balance'));
h=subtightplot(4,3,8,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,41)'-gobai.mean_o2_interp(:,:,41)',[-50 50],'WOA-GOBAI, 700m',cmocean('balance'));
h=subtightplot(4,3,9,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gobai.mean_o2_interp(:,:,41)'-gt.mean_o2(:,:,41)',[-50 50],'GOBAI-GT, 700m',cmocean('balance'));
% 1700m
h=subtightplot(4,3,10,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,61)'-gt.mean_o2(:,:,61)',[-50 50],'WOA-GT, 1700m',cmocean('balance'));
h=subtightplot(4,3,11,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,woa.o2(:,:,61)'-gobai.mean_o2_interp(:,:,61)',[-50 50],'WOA-GOBAI, 1700m',cmocean('balance'));
h=subtightplot(4,3,12,[xgap,ygap],h_marg,w_marg);
mean_map_plot(h,gt.lat,gt.lon,gobai.mean_o2_interp(:,:,61)'-gt.mean_o2(:,:,61)',[-50 50],'GOBAI-GT, 1700m',cmocean('balance'));
% colorbar
h=axes(f,'visible','off');
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,cmocean('balance'));
clim(h,[-50 50]);
% export
export_fig(gcf,'O2/Figures/mean_diff_comparisons.png','-transparent');
close

%% plot profile
[lon3d,lat3d] = meshgrid(woa.lon,woa.lat);
ll = (((lon3d + 0.5) - (lon3d - 0.5)) .* 111.320 .* cosd(lat3d)); % longitude distance (km)
ww = (((lat3d + 0.5) - (lat3d - 0.5)) .* 110.574); % latitude distance (km)
area = ll'.*ww';
% calculate means
for z = 1:length(gt.depth)
    idx = ~isnan(woa.o2(:,:,z)) & ~isnan(gt.mean_o2(:,:,z)) & ~isnan(gobai.mean_o2_interp(:,:,z));
    woa_o2_temp = woa.o2(:,:,z);
    mean_prof_woa(z) = sum(woa_o2_temp(idx).*area(idx))/sum(area(idx));
    gt_o2_temp = gt.mean_o2(:,:,z);
    mean_prof_gt(z) = sum(gt_o2_temp(idx).*area(idx))/sum(area(idx));
    gobai_o2_temp = gobai.mean_o2_interp(:,:,z);
    mean_prof_gobai(z) = sum(gobai_o2_temp(idx).*area(idx))/sum(area(idx));
end
% plot global mean O2
figure; hold on;
set(gcf,'position',[100 100 400 600]);
set(gca,'ydir','reverse');
plot(mean_prof_woa,gt.depth,'LineWidth',1);
plot(mean_prof_gt,gt.depth,'LineWidth',1);
plot(mean_prof_gobai,gt.depth,'LineWidth',1);
legend({'WOA' 'GT' 'GOBAI'},'Location','southeast');
xlim([130 250]);
ylabel('Depth (m)');
xlabel('[O_{2}] (\mumol kg^{-1})');
export_fig(gcf,'O2/Figures/profile_comparisons.png','-transparent');
close
% plot global mean delta O2
figure; hold on;
set(gcf,'position',[100 100 400 600]);
set(gca,'ydir','reverse');
plot(mean_prof_gt-mean_prof_woa,gt.depth);
plot(mean_prof_gobai-mean_prof_woa,gt.depth);
legend({'GT-WOA' 'GOBAI-WOA'},'Location','southeast');
xlim([-4 4]);
ylabel('Depth (m)');
xlabel('\Delta[O_{2}] (\mumol kg^{-1})');
export_fig(gcf,'O2/Figures/profile_diff_comparisons.png','-transparent');
close

%% plot timeseries
% calculate volume
idx_depth = gt.depth <= 700;
[lon3d,lat3d] = meshgrid(woa.lon,woa.lat);
ll = (((lon3d + 0.5) - (lon3d - 0.5)) .* 111.320 .* cosd(lat3d)); % longitude distance (km)
ww = (((lat3d + 0.5) - (lat3d - 0.5)) .* 110.574); % latitude distance (km)
area = ll'.*ww';
area = repmat(area,1,1,length(gt.depth(idx_depth)));
hh = diff(gt.depth);
hh(end) = [];
hh = repmat(permute(hh(idx_depth),[3 2 1]),length(woa.lon),length(woa.lat),1);
vol = area.*hh;
% calculate gt o2 timeseries
gt.o2 = ncread(gt_filename,'o2',[1 1 1 1],[Inf Inf sum(idx_depth) Inf]);
gt_o2_ts = nan(length(gt.time),1);
for t = 1:length(gt.time)
    o2_temp = gt.o2(:,:,:,t);
    idx = ~isnan(o2_temp);
    gt_o2_ts(t) = sum(o2_temp(idx).*vol(idx))./sum(vol(idx));
end
gt = rmfield(gt,'o2');
% calculate gobai volume
idx_depth = gobai.depth <= 700;
[lon3d,lat3d] = meshgrid(gobai.lon,gobai.lat);
ll = (((lon3d + 0.5) - (lon3d - 0.5)) .* 111.320 .* cosd(lat3d)); % longitude distance (km)
ww = (((lat3d + 0.5) - (lat3d - 0.5)) .* 110.574); % latitude distance (km)
area = ll'.*ww';
area = repmat(area,1,1,length(gt.depth(idx_depth)));
hh = diff(gobai.depth);
hh(end) = [];
hh = repmat(permute(hh(idx_depth),[3 2 1]),length(gobai.lon),length(gobai.lat),1);
vol = area.*hh;
% calculate gobai o2 timeseries
gobai.o2 = ncread(gobai_filename,'o2',[1 1 1 1],[Inf Inf sum(idx_depth) Inf]);
gobai_o2_ts = nan(length(gobai.time),1);
for t = 1:length(gobai.time)
    o2_temp = gobai.o2(:,:,:,t);
    idx = ~isnan(o2_temp);
    gobai_o2_ts(t) = sum(o2_temp(idx).*vol(idx))./sum(vol(idx));
end
% plot figure
figure; hold on;
set(gcf,'position',[100 100 800 400]);
title('Average [O_{2}] (0-700m)');
plot(double(datenum(1965,1,1)+gt.time),gt_o2_ts,'LineWidth',2);
plot(double(datenum(1950,1,1)+gobai.time),gobai_o2_ts,'LineWidth',2);
datetick('x');
ylabel('[O_{2}] (\mumol kg^{-1})');
export_fig(gcf,'O2/Figures/timeseries_comparisons.png','-transparent');
close

%% plot climatology

