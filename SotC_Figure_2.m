%% plot SoTC Figure 2
% load data
file_date = 'Jan2026';
fpath = '/fast4/o2/GODIP-DO/';
load([fpath 'O2_Maps/' 'GODIP-DO_NCEI_JDS_STATS_' file_date '.mat']);
% create figure
stip = false; stip_dens = 200; sig_prods = 9;
plt_alpha = 0.6; plt_bkgr = 'w';
figure('Position',[100 100 1200 1000],'visible','off','Color','w');
SotC_fig = tiledlayout(3,2,'TileSpacing','tight','Padding','none');
ax1 = nexttile; % inventory from 10-50m
idx_l = find(strcmp(layers,'0-50'));
worldmap(lat_lims,lon_lims);
set(ax1,'FontSize',12,'LineWidth',2,'TitleHorizontalAlignment','left');
setm(ax1,'MapProjection','miller','FFaceColor',plt_bkgr); tightmap;
title(['a) Oxygen Inventory from ' num2str(lyr_top(idx_l)) ...
    ' to ' num2str(lyr_bot(idx_l)) ' m (mol m^{-2})'],...
    'FontWeight','normal','FontSize',12);
pcolorm(lat,[lon;lon(end)+1],[mean(o2_inv_mean(:,:,:,idx_l),3,'omitnan');...
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
pcolorm(lat,[lon;lon(end)+1],[mean(o2_inv_mean(:,:,:,idx_l),3,'omitnan');...
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
idx = sum_pos >= 7 | sum_neg >= 7;
ens_mean_trend_sig(~idx) = NaN; ens_mean_trend_non(idx) = NaN; 
pcolorm(lat,[lon;lon(end)+1],[ens_mean_trend_sig;ens_mean_trend_sig(end,:)]');
h2=pcolorm(lat,[lon;lon(end)+1],[ens_mean_trend_non;ens_mean_trend_non(end,:)]');
h2.FaceAlpha = plt_alpha;
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
idx = sum_pos >= 7 | sum_neg >= 7;
ens_mean_trend_sig(~idx) = NaN; ens_mean_trend_non(idx) = NaN; 
pcolorm(lat,[lon;lon(end)+1],[ens_mean_trend_sig;ens_mean_trend_sig(end,:)]');
h2=pcolorm(lat,[lon;lon(end)+1],[ens_mean_trend_non;ens_mean_trend_non(end,:)]');
h2.FaceAlpha = plt_alpha;
if stip; stipplem(lat,lon,~idx' & ~isnan(ens_mean_trend)',...
        'density',stip_dens,'color',[.5 .5 .5]); end
c = colorbar; clim(inv_trend_lims(idx_l,:)); c.TickLength = 0;
colormap(ax4,flipud(cmocean('balance',length(inv_trend_levels{idx_l,:})-1,'pivot',0)));
plot_land('map',[1 1 1]);
mlabel off; plabel off; gridm off;
ax5 = nexttile([1 2]); % timeseries
hold on; box on;
set(ax5,'FontSize',12,'LineWidth',2,'YAxisLocation','right','TitleHorizontalAlignment','left');
title('e) Global Oxygen Inventory Anomalies from 0 to 1800 meters (Pmol)',...
    'FontWeight','normal','FontSize',12)
plot(time,repmat(0,length(time),1),'k:');
idx_l = find(strcmp(layers,'0-1800'));
for p = 1:length(products)
    if p <= 5
        plot(time,o2_inv_ts(:,p,idx_l)-mean(o2_inv_ts(baseline_idx,p,idx_l),'omitnan'),...
            'linewidth',2,'LineStyle','-');
    else
        plot(time,o2_inv_ts(:,p,idx_l)-mean(o2_inv_ts(baseline_idx,p,idx_l),'omitnan'),...
            'linewidth',2,'LineStyle','--');
    end
end
% ylabel('O_{2} Inventory Anomaly (Pmol)');
xlim([min(time) max(time)]); datetick('x','keeplimits'); 
legend([{''} labels],'NumColumns',3,'Location','southwest');
hold off
% save figure
if stip; figname = 'Figures/SotC_Fig2_stipple';
else; figname = 'Figures/SotC_Fig2'; end
exportgraphics(gcf,[figname '.png']);
if stip; figname = 'Figures/vectors/SotC_Fig2_stipple';
else; figname = 'Figures/vectors/SotC_Fig2'; end
exportgraphics(gcf,[figname '.eps'],'ContentType','vector');
close;
