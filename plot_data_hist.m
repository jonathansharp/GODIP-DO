% plot data histogram
function plot_data_hist(y1,y2,vrs)

% load data
load(['O2/Data/' vrs '_data_' num2str(y1) '_' num2str(y2)],'all_data');

% establish figure
figure('visible','off');
set(gcf,'Position',[100 100 1200 400]);

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
set(gca,'FontSize',20);
ylabel({'Number of Profiles';'(in thousands)'});

% add legend
legend({'OSD' 'CTD' 'PFL[A]' 'PFL[D]'},'Location','northwest'...
    ,'FontSize',20,'NumColumns',4);

% save figure
export_fig(['O2/Figures/' vrs '_dataset_histogram_' num2str(y1) '_' ...
    num2str(y2) '.png'],'-transparent','-silent');
