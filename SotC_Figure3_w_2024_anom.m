%% plot SoTC Figure 3
% load data
file_date = 'Jan2026';
fpath = '/fast4/o2/GODIP-DO/';
load([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_STATS_' file_date '.mat']);
% create figure
stip = false; stip_dens = 200; sig_prods = 5;
figure('Position',[100 100 1400 1050],'visible','off');
set(gcf,'Color','w');
[lat3d,lon3d] = ndgrid(lat,lon);
SotC_fig = tiledlayout(3,2,'TileSpacing','tight','Padding','none');
ax1 = nexttile; % anomaly at 10
idx_l = find(strcmp(layers,'0-50'));
worldmap(lat_lims,lon_lims);
set(ax1,'FontSize',12,'LineWidth',2);
setm(ax1,'MapProjection','miller'); tightmap;
title(['2024 Inv. Anom. from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})']);
ens_inv_anom = mean(o2_inv_anom_2024(:,:,idx_l,:),4,'omitnan');
ens_inv_anom_sig = mean(o2_inv_anom_2024(:,:,idx_l,:),4,'omitnan');
sum_pos = sum(o2_inv_anom_2024(:,:,idx_l,:) > 0,4);
sum_neg = sum(o2_inv_anom_2024(:,:,idx_l,:) < 0,4);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_inv_anom_sig(~idx) = NaN;
pcolorm(lat,[lon;lon(end)+1],[ens_inv_anom_sig;ens_inv_anom_sig(end,:)]');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_anom)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)); c.TickLength = 0;
colormap(ax1,flipud(cmocean('balance',length(inv_trend_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax2 = nexttile; % anomaly at 10
idx_l = find(strcmp(layers,'50-1000'));
worldmap(lat_lims,lon_lims);
set(ax2,'FontSize',12,'LineWidth',2);
setm(ax2,'MapProjection','miller'); tightmap;
title(['2024 Inv. Anom. from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})']);
ens_inv_anom = mean(o2_inv_anom_2024(:,:,idx_l,:),4,'omitnan');
ens_inv_anom_sig = mean(o2_inv_anom_2024(:,:,idx_l,:),4,'omitnan');
sum_pos = sum(o2_inv_anom_2024(:,:,idx_l,:) > 0,4);
sum_neg = sum(o2_inv_anom_2024(:,:,idx_l,:) < 0,4);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_inv_anom_sig(~idx) = NaN;
pcolorm(lat,[lon;lon(end)+1],[ens_inv_anom_sig;ens_inv_anom_sig(end,:)]');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_anom)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)); c.TickLength = 0;
colormap(ax2,flipud(cmocean('balance',length(inv_trend_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax3 = nexttile; % anomaly at 10
idx_l = find(strcmp(layers,'0-50'));
worldmap(lat_lims,lon_lims);
set(ax3,'FontSize',12,'LineWidth',2);
setm(ax3,'MapProjection','miller'); tightmap;
title(['2025' char(8722) '2024 Inv. from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})']);
ens_inv_diff = mean(o2_inv_diff_2025_2024(:,:,idx_l,:),4,'omitnan');
ens_inv_diff_sig = mean(o2_inv_diff_2025_2024(:,:,idx_l,:),4,'omitnan');
sum_pos = sum(o2_inv_diff_2025_2024(:,:,idx_l,:) > 0,4);
sum_neg = sum(o2_inv_diff_2025_2024(:,:,idx_l,:) < 0,4);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_inv_diff_sig(~idx) = NaN;
pcolorm(lat,[lon;lon(end)+1],[ens_inv_diff_sig;ens_inv_diff_sig(end,:)]');
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_inv_diff)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_anom_lims(idx_l,:)); c.TickLength = 0;
colormap(ax3,flipud(cmocean('balance',length(inv_trend_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax4 = nexttile; % anomaly at 10
idx_l = find(strcmp(layers,'50-1000'));
worldmap(lat_lims,lon_lims);
set(ax4,'FontSize',12,'LineWidth',2);
setm(ax4,'MapProjection','miller'); tightmap;
title(['2025' char(8722) '2024 Inv. from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})']);
ens_inv_diff = mean(o2_inv_diff_2025_2024(:,:,idx_l,:),4,'omitnan');
ens_inv_diff_sig = mean(o2_inv_diff_2025_2024(:,:,idx_l,:),4,'omitnan');
sum_pos = sum(o2_inv_diff_2025_2024(:,:,idx_l,:) > 0,4);
sum_neg = sum(o2_inv_diff_2025_2024(:,:,idx_l,:) < 0,4);
idx = sum_pos >= sig_prods | sum_neg >= sig_prods;
ens_inv_diff_sig(~idx) = NaN;
pcolorm(lat,[lon;lon(end)+1],[ens_inv_diff_sig;ens_inv_diff_sig(end,:)]');
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
set(ax5,'YDir','reverse','FontSize',12,'LineWidth',2);
title('Global [O_{2}] Anomalies (\mumol kg^{-1})');
c = colorbar; clim([-6 6]); c.TickLength = 0;
colormap(ax5,flipud(cmocean('balance',24,'pivot',0)));
ylabel('Depth (m)');
xlim([min(time) max(time)]); datetick('x','keeplimits'); 
hold off;
% add annotations
annotation('textbox',[0.083 0.94 0 0],'string','a','FontSize',14,'FontWeight','bold');
annotation('textbox',[0.55 0.94 0 0],'string','b','FontSize',14,'FontWeight','bold');
annotation('textbox',[0.083 0.614 0 0],'string','c','FontSize',14,'FontWeight','bold');
annotation('textbox',[0.55 0.614 0 0],'string','d','FontSize',14,'FontWeight','bold');
annotation('textbox',[0.07 0.363 0 0],'string','e','FontSize',14,'FontWeight','bold');
% save figure
if stip; figname = 'Figures/SotC_Fig3_stipple_w2024';
else; figname = 'Figures/SotC_Fig3_w2024'; end
export_fig([figname '.png'],'-transparent');
if stip; figname = 'Figures/vectors/SotC_Fig3_stipple_w2024';
else; figname = 'Figures/vectors/SotC_Fig3_w2024'; end
exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
close;
