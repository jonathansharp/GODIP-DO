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