% plot data histogram
function plot_data_map(y1,y2,vrs)

% load data
load(['O2/Data/' vrs '_data_' num2str(y1) '_' num2str(y2)],'all_data');

% establish figure
figure('visible','on'); 
worldmap([-90 90],[20 380]); hold on;
setm(gca,'MapProjection','robinson');
set(gcf,'Position',[100 100 1200 800]);

% index to unique profiles
[~,prof_idx] = unique(all_data.profile);
lon = all_data.lon(prof_idx);
lat = all_data.lat(prof_idx);
type = all_data.type(prof_idx);

% plot histogram
ax_pfl = scatterm(lat(type==3),lon(type==3),1,'k','filled');
ax_ctd = scatterm(lat(type==1),lon(type==1),1,'r','filled');
ax_osd = scatterm(lat(type==2),lon(type==2),1,'b','filled');
set(gca,'FontSize',20);
mlabel off; plabel off;

% add legend
legend([ax_pfl,ax_ctd,ax_osd],{'PFL' 'CTD' 'OSD'},'Location','northwest','FontSize',20,'NumColumns',3);

% save figure
export_fig(['O2/Figures/dataset_map_' num2str(y1) '_' ...
    num2str(y2) '.png'],'-transparent','-silent');
