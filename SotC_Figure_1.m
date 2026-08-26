% plot SotC Figure 1


% parameters
close all;
fpath = '/fast4/o2/GODIP-DO/O2_Maps/';
vrs = 'WOD25';
y1 = 1965;
y2 = 2025;

% create figure
h=figure('Position',[100 100 1200 1200],'visible','on');
a=tiledlayout(h,3,1,'TileSpacing','tight','Padding','tight');

% % load data from ncei binned file
% datasets.ncei.name = 'w25_o_19652025_yearly_anomaly_NCEI.nc';
% datasets.ncei.lat = ncread([fpath datasets.ncei.name],'lat');
% datasets.ncei.lon = ncread([fpath datasets.ncei.name],'lon');
% datasets.ncei.time = double(ncread([fpath datasets.ncei.name],'time'));
% datasets.ncei.depth = ncread([fpath datasets.ncei.name],'depth');
% datasets.ncei.date = datevec(datenum(1900,1,1)+datasets.ncei.time);
% datasets.ncei.param_name = 'o_anom_an';

% load processes profile data
load(['O2/Data/' vrs '_data_' num2str(y1) '_' num2str(y2)],'all_data');

% plot map
ax1 = nexttile(a,1,[2 1]);
worldmap([-75 75],[20 380]); hold on;
setm(gca,'MapProjection','miller');
ax1.TitleHorizontalAlignment = 'left';
title('a) Spatial distribution of O_{2} profile data by platform (1965-2025)',...
    'FontWeight','normal','FontSize',14);

% index to unique profiles
[~,prof_idx] = unique(all_data.profile);
lon = all_data.lon(prof_idx);
lat = all_data.lat(prof_idx);
type = all_data.type(prof_idx);

% 
clrs = colororder;

% plot map
ax_pfl = scatterm(lat(type==3),lon(type==3),5,clrs(3,:),'.');
mkr_pfl = scatterm(nan,nan,1000,clrs(3,:),'o','filled');
ax_osd = scatterm(lat(type==1),lon(type==1),5,clrs(1,:),'.');
mkr_osd = scatterm(nan,nan,1000,clrs(1,:),'o','filled');
ax_ctd = scatterm(lat(type==2),lon(type==2),5,clrs(2,:),'.');
mkr_ctd = scatterm(nan,nan,1000,clrs(2,:),'o','filled');
mlabel off; plabel on; gridm off;
plot_land('map','w');

% add legend
legend([mkr_osd,mkr_ctd,mkr_pfl],...
    {'Ocean Station Data (OSD)' 'CTD Sensor (CTD)' ...
    'Profiling Floats (PFL)'},'Location','northwest',...
    'FontSize',14,'NumColumns',1);

% plot histogram
ax2=nexttile(a);

% index to unique profiles
[~,prof_idx] = unique(all_data.profile);
vars = fieldnames(all_data);
for v = 1:length(vars)
    all_data.(vars{v}) = all_data.(vars{v})(prof_idx);
end

% index to each dataset type
y_ctd = all_data.year(all_data.type==2);
y_osd = all_data.year(all_data.type==1);
y_flt_a = all_data.year(all_data.type==3 & all_data.mode==3);
y_flt_d = all_data.year(all_data.type==3 & all_data.mode==2);

% count profiles in each year
counts_ctd = histc(y_ctd,y1:y2);
counts_osd = histc(y_osd,y1:y2);
counts_flt_a = histc(y_flt_a,y1:y2);
counts_flt_d = histc(y_flt_d,y1:y2);

% plot histogram
bar_chart = bar(y1:y2,[counts_osd counts_ctd counts_flt_a counts_flt_d]/1000,'stacked');

% adjust colors
set(bar_chart,'FaceColor','Flat');
bar_chart(4).CData = bar_chart(3).CData + 0.15;

% adjust properties
ylim([0 29]);
ax2.TitleHorizontalAlignment = 'left';
title('b) Annual number of O_{2} profiles by platform (1965-2025)',...
    'FontWeight','normal','FontSize',14);
ylabel({'Number of Profiles';'(in thousands)'},'FontSize',14);

% add legend
legend({'OSD' 'CTD' 'PFL[A]' 'PFL[D]'},'Location','northwest',...
    'FontSize',14,'NumColumns',4);

% save figure
exportgraphics(gcf,['Figures/SotC_Fig1_' vrs '_dataset_histogram_' num2str(y1) '_' ...
    num2str(y2) '.png'],'Resolution',600);
exportgraphics(gcf,['Figures/vectors/SotC_Fig1_' vrs '_dataset_histogram_' num2str(y1) '_' ...
    num2str(y2) '.eps'],'Resolution',600,'ContentType','vector');
close;
