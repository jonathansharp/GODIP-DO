%% file information
file_date = 'May2026';
fpath = '/fast4/o2/GODIP-DO/';
fpath_EN4 = '/med2/';
compilation.name = ['GODIP-DO_SotC-2025_' file_date '.nc'];

% product indices
p_indices = [1:5,6,8,13,14];
products = ncread([fpath 'O2_Maps/' compilation.name],'products');
labels = {'NCEI' 'IAP' 'GT-OI' 'UTAS' 'SJTU-GR' 'GOBAI' 'GOBAI-IAP' ...
    'GT-ML' 'GT-NN-EN4' 'GT-NN-IAP' 'GT-RF-EN4' 'GT-RF-IAP' 'Jingwei' 'SJTU-HZ'};

% download dimensions
lat = ncread([fpath 'O2_Maps/' compilation.name],'lat');
lon = ncread([fpath 'O2_Maps/' compilation.name],'lon');
time = ncread([fpath 'O2_Maps/' compilation.name],'time');
depth = ncread([fpath 'O2_Maps/' compilation.name],'depth');
mask = logical(ncread([fpath 'O2_Maps/' compilation.name],'mask'));
dates = datevec(time); years = dates(:,1);

% define depth limits
d1 = 0; d2 = 1800;
[~,d1_idx] = min(abs(depth-d1));
[~,d2_idx] = min(abs(depth-d2));
depth = depth(d1_idx:d2_idx);

for p = 1:length(p_indices)

%% figure information
f = gcf;
f.Position(3) = f.Position(3)*2;
f.Position(4) = f.Position(4)*2;
tiledlayout(f,2,1);

% download oxygen
oxy = squeeze(ncread([fpath 'O2_Maps/' compilation.name],products{p},...
        [1 1 d1_idx 1],[Inf Inf d2_idx-d1_idx+1 Inf]));
% calculate weights
vol = weights3d(lon,lat,depth);
vol(isnan(mean(oxy,4,'omitnan'))) = NaN;
dens = 1025; mass = vol .* dens;
% calculate global mean
% global_mean = sum(reshape(oxy,[length(lon)*...
%     length(lat)*length(depth) length(time)]).*...
%     vol(:),'omitnan')./(sum(vol(:),'omitnan'));
% [yf,yr,x] = leastsq2(double(time),global_mean,double(time(1)),0);
% calculate global inv
global_inv = sum(reshape(oxy,[length(lon)*...
    length(lat)*length(depth) length(time)]).*...
    mass(:),'omitnan')./(10^18);
global_inv(global_inv == 0) = NaN;
[yf,yr,x] = leastsq2(double(time),global_inv,double(time(1)),0);

%% plot timeseries with fit
nexttile;
plot(double(time),global_inv,'LineWidth',2); hold on;
plot(double(time),yf,'LineWidth',2,'LineStyle','--','Color','k');
datetick('x','yyyy');
% ylim([149 153]);
title(['Global Mean ' labels{p} ' O_{2}']);
ylabel('O_{2} Inventory (Pmol)');
legend({'O_{2} Inventory' 'Least-squares fit'})

%% obtain annual mean enso index
enso_table = readtable('meiv2.data.txt');
enso_index_y = table2array(enso_table(:,2:13))';
enso_year = table2array(enso_table(:,1))';
mei_yr_idx = ismember(enso_year,years);
enso_index_y = mean(enso_index_y(:,mei_yr_idx))';
enso_year = enso_year(mei_yr_idx)';
% enso_index = table2array(enso_table(:,2:13))';
% enso_index = enso_index(:,mei_yr_idx);
% enso_index = enso_index(:);
% enso_index_sm = movmean(enso_index,12);
% enso_index_m = enso_index_sm(7:12:end);

%% obtain APO data
apo_ljo = readtable('monthly_o2_ljo.csv');
apo_yr = table2array(apo_ljo(2:end,1));
apo_mn = table2array(apo_ljo(2:end,2));
o2_atm = table2array(apo_ljo(2:end,5));
apo_time = datenum(apo_yr,apo_mn,15);
[o2_atm_yf,o2_atm_yr] = leastsq2(double(apo_time),o2_atm,...
    double(apo_time(1)),2,[365.2425 265.2425/2]);

%% plot anomaly
common_year_idx = ismember(years,enso_year);
corr_coeff = corr(enso_index_y,yr(common_year_idx));
nexttile;
title(['Global Mean ' labels{p} ' O_{2} Anomaly, ' ...
    '\itr\rm = ' num2str(corr_coeff)],'FontWeight','bold');
yyaxis left;
plot(double(time),yr,'LineWidth',2);
datetick('x','yyyy');
%ylim([-0.3 0.3]);
ylabel('O_{2} Inventory Anomaly (Tmol)');
yyaxis right;
plot(datenum(enso_year,7,1),enso_index_y,'LineWidth',2);
ylabel('Multivariate ENSO Index');
text(datenum(enso_year(5),7,1),max(yr)-200,...
    ['\itr\rm = ' num2str(corr_coeff)],...
    'VerticalAlignment','bottom','HorizontalAlignment','right')

%%
% nexttile;
% plot(apo_time,o2_atm_yr,'LineWidth',2);

%% export figure
exportgraphics(gcf,['Figures/global_' products{p} '_anomaly.png']);
close

end