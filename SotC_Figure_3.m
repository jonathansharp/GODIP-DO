%% plot SoTC Figure 3
% load data
file_date = 'Jan2026';
fpath = '/fast4/o2/GODIP-DO/';
load([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_STATS_' file_date '.mat']);
% create figure
stip = false; stip_dens = 200; sig_prods = 5;
plt_alpha = 0.6; plt_bkgr = 'w';
figure('Position',[100 100 1400 1050],'visible','off','Color','w');
SotC_fig = tiledlayout(3,2,'TileSpacing','tight','Padding','none');
ax1 = nexttile; % anomaly at 10
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
pcolorm(lat,[lon;lon(end)+1],[ens_inv_anom_sig;ens_inv_anom_sig(end,:)]');
h2=pcolorm(lat,[lon;lon(end)+1],[ens_inv_anom_non;ens_inv_anom_non(end,:)]');
h2.FaceAlpha = plt_alpha;
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_anom)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)); c.TickLength = 0;
colormap(ax1,flipud(cmocean('balance',length(inv_trend_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax2 = nexttile; % anomaly at 10
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
pcolorm(lat,[lon;lon(end)+1],[ens_inv_anom_sig;ens_inv_anom_sig(end,:)]');
h2=pcolorm(lat,[lon;lon(end)+1],[ens_inv_anom_non;ens_inv_anom_non(end,:)]');
h2.FaceAlpha = plt_alpha;
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_anom)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)/2); c.TickLength = 0;
colormap(ax2,flipud(cmocean('balance',length(inv_anom_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax3 = nexttile; % anomaly at 10
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
pcolorm(lat,[lon;lon(end)+1],[ens_inv_diff_sig;ens_inv_diff_sig(end,:)]');
h2=pcolorm(lat,[lon;lon(end)+1],[ens_inv_diff_non;ens_inv_diff_non(end,:)]');
h2.FaceAlpha = plt_alpha;
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_diff)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)/2); c.TickLength = 0;
colormap(ax3,flipud(cmocean('balance',length(inv_anom_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax4 = nexttile; % anomaly at 10
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
pcolorm(lat,[lon;lon(end)+1],[ens_inv_diff_sig;ens_inv_diff_sig(end,:)]');
h2=pcolorm(lat,[lon;lon(end)+1],[ens_inv_diff_non;ens_inv_diff_non(end,:)]');
h2.FaceAlpha = plt_alpha;
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_diff)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)); c.TickLength = 0; colormap(ax4,flipud(cmocean('balance',length(inv_trend_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
%
ax5 = nexttile([1 2]); % timeseries
hold on; box on;
contourf(time,depth,squeeze(mean(o2_inv_depth_ts,2,'omitnan'))'-...
    mean(squeeze(mean(o2_inv_depth_ts,2,'omitnan')),1,'omitnan')',-6:0.5:6);
set(ax5,'YDir','reverse','FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
title('e) Global Oxygen Anomalies (\mumol kg^{-1})',...
    'FontWeight','normal','FontSize',12);
c = colorbar; clim([-6 6]); c.TickLength = 0;
colormap(ax5,flipud(cmocean('balance',24,'pivot',0)));
ylabel('Depth (m)');
xlim([min(time) max(time)]); datetick('x','keeplimits'); 
hold off;
% save figure
if stip; figname = 'Figures/SotC_Fig3_stipple';
else; figname = 'Figures/SotC_Fig3'; end
exportgraphics(gcf,[figname '.png']);
if stip; figname = 'Figures/vectors/SotC_Fig3_stipple';
else; figname = 'Figures/vectors/SotC_Fig3'; end
exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
close;
