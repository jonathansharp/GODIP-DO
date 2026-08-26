%% Climatological EN4 average
fpath = '/raid/Data/EN4/c14/';
fname_clim = 'EN.4.2.2.f.analysis.c14.1965-2023.climatology.nc';
% create climatology file
if isfile([fpath fname_clim]); delete([fpath fname_clim]); end
nccreate([fpath fname_clim],'temperature','Dimensions',{'lon',360,'lat',173,'depth',42,'month',12});
nccreate([fpath fname_clim],'salinity','Dimensions',{'lon',360,'lat',173,'depth',42,'month',12});
nccreate([fpath fname_clim],'lon','Dimensions',{'lon'});
nccreate([fpath fname_clim],'lat','Dimensions',{'lat'});
nccreate([fpath fname_clim],'depth','Dimensions',{'depth'});
nccreate([fpath fname_clim],'month','Dimensions',{'month'});
% write dimensions to file
fname = 'EN.4.2.2.f.analysis.c14.196501.nc';
lon = ncread([fpath fname],'lon'); ncwrite([fpath fname_clim],'lon',lon);
lat = ncread([fpath fname],'lat'); ncwrite([fpath fname_clim],'lat',lat);
depth = ncread([fpath fname],'depth'); ncwrite([fpath fname_clim],'depth',depth);
ncwrite([fpath fname_clim],'month',(1:12)');
% get average for each month of year
for m = 1:12
    % establish varaibles
    temperature = nan(360,173,42,2024-1965);
    salinity = nan(360,173,42,2024-1965);
    % download values for each year
    for y = 1965:2023
        fname = ['EN.4.2.2.f.analysis.c14.' num2str(y) sprintf('%02d',m) '.nc'];
        ncinfo([fpath fname]);
        temperature(:,:,:,y-1964) = ncread([fpath fname],'temperature');
        salinity(:,:,:,y-1964) = ncread([fpath fname],'salinity');
    end
    % write to climatology file
    ncwrite([fpath 'EN.4.2.2.f.analysis.c14.1965-2023.climatology.nc'],...
        'temperature',mean(temperature,4,'omitnan'),[1 1 1 m]);
    ncwrite([fpath 'EN.4.2.2.f.analysis.c14.1965-2023.climatology.nc'],...
        'salinity',mean(salinity,4,'omitnan'),[1 1 1 m]);
end

%% read EN4 coordinates
EN4.latitude = ncread('/raid/Data/EN4/c14/EN.4.2.2.f.analysis.c14.196501.nc','lat');
EN4.longitude = ncread('/raid/Data/EN4/c14/EN.4.2.2.f.analysis.c14.196501.nc','lon');
EN4.depth = ncread('/raid/Data/EN4/c14/EN.4.2.2.f.analysis.c14.196501.nc','depth');
% expand latitude and longitude to 2-D
latitude_3d = repmat(EN4.latitude',length(EN4.longitude),1);
longitude_3d = repmat(EN4.longitude,1,length(EN4.latitude));
longitude_3d(longitude_3d>180) = longitude_3d(longitude_3d>180) - 360;

%% Create basin index
createShapes;
basins = {'Atl' 'Pac' 'Ind' 'Arc' 'Med' 'NSou' 'SSou'};
for b = 1:7
    basin_index.EN4.(basins{b}) = single(zeros(size(EN4.temp)));
    if strcmp(basins{b},'Pac')
        basin_index.EN4.(basins{b}) = ...
            inpolygon(longitude_3d,latitude_3d,...
            Poly.([basins{b} '1'])(:,1),Poly.([basins{b} '1'])(:,2)) | ...
            inpolygon(longitude_3d,latitude_3d,...
            Poly.([basins{b} '2'])(:,1),Poly.([basins{b} '2'])(:,2));
    else
        basin_index.EN4.(basins{b}) = ...
            inpolygon(longitude_3d,latitude_3d,...
            Poly.(basins{b})(:,1),Poly.(basins{b})(:,2));
    end
end
clear longitude_3d latitude_3d

%% Predict O2 on EN4 grid 
EN4.latitude = ncread('/raid/Data/EN4/c14/EN.4.2.2.f.analysis.c14.196501.nc','lat');
EN4.longitude = ncread('/raid/Data/EN4/c14/EN.4.2.2.f.analysis.c14.196501.nc','lon');
EN4.depth = ncread('/raid/Data/EN4/c14/EN.4.2.2.f.analysis.c14.196501.nc','depth');

% load
EN4.temp = nan(length(EN4.longitude),length(EN4.latitude),length(EN4.depth),12);
EN4.sal = nan(length(EN4.longitude),length(EN4.latitude),length(EN4.depth));
EN4.time = nan(12,1);
for m = 1:12
    EN4.temp(:,:,:,m) = ncread(['/raid/Data/EN4/c14/EN.4.2.2.f.analysis.c14.' ...
        num2str(y1-1+y) sprintf('%02d',m) '.nc'],'temperature');
    EN4.sal(:,:,:,m) = ncread(['/raid/Data/EN4/c14/EN.4.2.2.f.analysis.c14.' ...
        num2str(y1-1+y) sprintf('%02d',m) '.nc'],'salinity');
    EN4.time(m) = ncread(['/raid/Data/EN4/c14/EN.4.2.2.f.analysis.c14.' ...
        num2str(y1-1+y) sprintf('%02d',m) '.nc'],'time');
end

% convert from day since 1/1/1800 to datenum
EN4.time = datenum(1800,1,EN4.time);
date = datevec(EN4.time);
date0 = date;
date0(:,2:3) = 1;
EN4.day = (datenum(date) - datenum(date0));
EN4.year = date(:,1);
clear date date0

% convert longitude
EN4.longitude(EN4.longitude>180) = EN4.longitude(EN4.longitude>180)-360;

% convert to singles to save memory
EN4.temp = single(EN4.temp);
EN4.sal = single(EN4.sal);
EN4.day = single(EN4.day);
EN4.year = single(EN4.year);

% expand latitude, longitude, depth, and time to 4-D
EN4.latitude = repmat(EN4.latitude',length(EN4.longitude),1,...
    length(EN4.depth),length(EN4.time));
EN4.longitude = repmat(EN4.longitude,1,size(EN4.latitude,2),...
    length(EN4.depth),length(EN4.time));
EN4.depth = repmat(permute(EN4.depth,[3 2 1]),size(EN4.longitude,1),...
    size(EN4.latitude,2),1,length(EN4.time));
EN4.year = repmat(permute(EN4.year,[4 3 2 1]),size(EN4.longitude,1),...
    size(EN4.latitude,2),size(EN4.depth,3),1);
EN4.day = repmat(permute(EN4.day,[4 3 2 1]),size(EN4.longitude,1),...
    size(EN4.latitude,2),size(EN4.depth,3),1);

% calculate pressure from depth
EN4.pressure = gsw_p_from_z(-EN4.depth,EN4.latitude);

% replace pressure with NaN where no T/S values exist
EN4.pressure(isnan(EN4.temp)) = NaN;

% ***** replace data with constant values to remove certain controls *****
% EN4.year(~isnan(EN4.year)) = mean(EN4.year(:),'omitnan');
% EN4.day(~isnan(EN4.day)) = 182;
% EN4.sal(~isnan(EN4.sal)) = mean(EN4.sal(:),'omitnan');
% EN4.temp(~isnan(EN4.temp)) = mean(EN4.temp(:),'omitnan');
% EN4.latitude(~isnan(EN4.latitude)) = 0;
% EN4.longitude(~isnan(EN4.longitude)) = 0;

% Calculate absolute salinity, conservative temperature, potential density,
% and spice from climatologies (do this using a loop for each timestep, to
% avoid memory over-usage)
EN4.sal_abs = single(nan(size(EN4.temp)));
EN4.temp_cns = single(nan(size(EN4.temp)));
EN4.sigma0 = single(nan(size(EN4.temp)));
% EN4.spice = single(nan(size(EN4.temp)));
for t = 1:length(EN4.time)
    EN4.sal_abs(:,:,:,t) = gsw_SA_from_SP(EN4.sal(:,:,:,t),EN4.pressure(:,:,:,t),EN4.longitude(:,:,:,t),EN4.latitude(:,:,:,t));
    EN4.temp_cns(:,:,:,t) = gsw_CT_from_t(EN4.sal_abs(:,:,:,t),EN4.temp(:,:,:,t),EN4.pressure(:,:,:,t));
    EN4.sigma0(:,:,:,t) = gsw_sigma0(EN4.sal_abs(:,:,:,t),EN4.temp_cns(:,:,:,t));
    % EN4.spice(:,:,:,t) = gsw_spiciness0(EN4.sal_abs(:,:,:,t),EN4.temp_cns(:,:,:,t));
end
clear t

% Transform day by sine and cosine:
EN4.day_sin = sin((2.*pi.*EN4.day)./366);
EN4.day_cos = cos((2.*pi.*EN4.day)./366);
EN4 = rmfield(EN4,'day');

% Transform longitude by sine and cosine:
EN4.lon_cos1 = cosd(EN4.longitude-20);
EN4.lon_cos2 = cosd(EN4.longitude-110);

% calculate bottom depth
EN4.bottom_depth = ...
    single(bottom_depth(EN4.latitude(:,:,1,1),EN4.longitude(:,:,1,1)));
EN4.bottom_depth = repmat(EN4.bottom_depth,1,1,...
    size(EN4.latitude,3),size(EN4.latitude,4));

EN4.oxy_ENS = nan(size(EN4.temp));
EN4.oxy_RF = nan(size(EN4.temp));
EN4.oxy_NN = nan(size(EN4.temp));

%% Loop through and make predictions for each basin

for b = 1:7

idx = basin_index.EN4.(basins{b});
idx = repmat(idx,1,1,size(EN4.oxy_ENS,3),size(EN4.oxy_ENS,4));

% Index to where P,T,S are available and matches basin index
index = ~isnan(EN4.pressure(:)) & idx(:);
clear idx

% Pre-allocate weights
weights = zeros(size(index));

% Concatenate predictors 
EN4.all = table(EN4.latitude(:),EN4.lon_cos1(:),EN4.lon_cos2(:),...
    EN4.sigma0(:),EN4.pressure(:),EN4.year(:),EN4.day_sin(:),...
    EN4.day_cos(:),EN4.temp_cns(:),EN4.sal_abs(:),EN4.bottom_depth(:),...
    'VariableNames',{'latitude','lon_cos1','lon_cos2','sigma','pressure',...
    'year','day_sin','day_cos','temperature_cns','salinity_abs','bottom_depth'});

% For Atlantic basin
if strcmp(basins{b},'Atl')
    % initially set all basin weights equal to 1
    weights(index) = 1;
    % adjust overlaps with Arctic
    idx_Arc = index & EN4.latitude(:) > 35;
    weights(idx_Arc) = (45 - EN4.latitude(idx_Arc))./10;
    % adjust overlaps with Southern
    idx_NSou = index & EN4.latitude(:) < -25;
    weights(idx_NSou) = (EN4.latitude(idx_NSou) + 35)./10;
    % calculate Arctic weights
    weights_Arc = calculateWeights(weights,idx_Arc);
    % calculate N. Southern weights
    weights_NSou = calculateWeights(weights,idx_NSou);
    % run data through random forest models
    oxy_RF_Atl_temp = runRF('All','Atl','Atl',EN4.all,index,1:11);
    oxy_RF_Arc_temp = runRF('All','Atl','Arc',EN4.all,index,1:11,idx_Arc);
    oxy_RF_NSou_temp = runRF('All','Atl','NSou',EN4.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_Atl_temp .* weights(index) + ...
        oxy_RF_Arc_temp .* weights_Arc(index) + ...
        oxy_RF_NSou_temp .* weights_NSou(index);
    % run data through neural networks
    oxy_NN_Atl_temp = runNN('All','Atl','Atl',EN4.all,index,1:11);
    oxy_NN_Arc_temp = runNN('All','Atl','Arc',EN4.all,index,1:11,idx_Arc);
    oxy_NN_NSou_temp = runNN('All','Atl','NSou',EN4.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_NN.(basins{b}) = ...
        oxy_NN_Atl_temp .* weights(index) + ...
        oxy_NN_Arc_temp .* weights_Arc(index) + ...
        oxy_NN_NSou_temp .* weights_NSou(index);
clear oxy_RF_Atl_temp oxy_RF_Arc_temp oxy_RF_NSou_temp
clear oxy_NN_Atl_temp oxy_NN_Arc_temp oxy_NN_NSou_temp
clear idx_Arc idx_NSou weights weights_Arc weights_NSou

% For Pacific basin
elseif strcmp(basins{b},'Pac')
    % initially set all basin weights equal to 1
    weights(index) = 1;
    % adjust overlaps with Arctic
    idx_Arc = index & EN4.latitude(:) > 64;
    weights(idx_Arc) = (70 - EN4.latitude(idx_Arc))./6;
    % adjust overlaps with Southern
    idx_NSou = index & EN4.latitude(:) < -25;
    weights(idx_NSou) = (EN4.latitude(idx_NSou) + 35)./10;
    % calculate Arctic weights
    weights_Arc = calculateWeights(weights,idx_Arc);
    % calculate N. Southern weights
    weights_NSou = calculateWeights(weights,idx_NSou);
    % run data through random forest models
    oxy_RF_Pac_temp = runRF('All','Pac','Pac',EN4.all,index,1:11);
    oxy_RF_Arc_temp = runRF('All','Pac','Arc',EN4.all,index,1:11,idx_Arc);
    oxy_RF_NSou_temp = runRF('All','Pac','NSou',EN4.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_Pac_temp .* weights(index) + ...
        oxy_RF_Arc_temp .* weights_Arc(index) + ...
        oxy_RF_NSou_temp .* weights_NSou(index);
    % run data through neural networks
    oxy_NN_Pac_temp = runNN('All','Pac','Pac',EN4.all,index,1:11);
    oxy_NN_Arc_temp = runNN('All','Pac','Arc',EN4.all,index,1:11,idx_Arc);
    oxy_NN_NSou_temp = runNN('All','Pac','NSou',EN4.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_NN.(basins{b}) = ...
        oxy_NN_Pac_temp .* weights(index) + ...
        oxy_NN_Arc_temp .* weights_Arc(index) + ...
        oxy_NN_NSou_temp .* weights_NSou(index);
clear oxy_RF_Pac_temp oxy_RF_Arc_temp oxy_RF_NSou_temp
clear oxy_NN_Pac_temp oxy_NN_Arc_temp oxy_NN_NSou_temp
clear idx_Arc idx_NSou weights weights_Arc weights_NSou

% For Indian basin
elseif strcmp(basins{b},'Ind')
    % initially set all basin weights equal to 1
    weights(index) = 1;
    % adjust overlaps with Southern
    idx_NSou = index & EN4.latitude(:) < -35;
    weights(idx_NSou) = (EN4.latitude(idx_NSou) + 45)./10;
    % calculate N. Southern weights
    weights_NSou = calculateWeights(weights,idx_NSou);
    % run data through random forest models
    oxy_RF_Ind_temp = runRF('All','Ind','Ind',EN4.all,index,1:11);
    oxy_RF_NSou_temp = runRF('All','Ind','NSou',EN4.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_Ind_temp .* weights(index) + ...
        oxy_RF_NSou_temp .* weights_NSou(index);
    % run data through neural networks
    oxy_NN_Ind_temp = runNN('All','Ind','Ind',EN4.all,index,1:11);
    oxy_NN_NSou_temp = runNN('All','Ind','NSou',EN4.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_NN.(basins{b}) = ...
        oxy_NN_Ind_temp .* weights(index) + ...
        oxy_NN_NSou_temp .* weights_NSou(index);
clear oxy_NN_Ind_temp oxy_NN_NSou_temp
clear oxy_RF_Ind_temp oxy_RF_NSou_temp
clear idx_NSou weights weights_NSou

% For Arctic basin
elseif strcmp(basins{b},'Arc')
    % initially set all basin weights equal to 1
    weights(index) = 1;
    % adjust overlaps with Pacific
    idx_Pac = index & EN4.latitude(:) < 70 & EN4.longitude(:) < -120;
    weights(idx_Pac) = (EN4.latitude(idx_Pac) - 64)./6;
    % adjust overlaps with Atlantic
    idx_Atl = index & EN4.latitude(:) < 45 & EN4.longitude(:) > -120;
    weights(idx_Atl) = (EN4.latitude(idx_Atl) - 35)./10;
    % calculate Pacific weights
    weights_Pac = calculateWeights(weights,idx_Pac);
    % calculate Atlantic weights
    weights_Atl = calculateWeights(weights,idx_Atl);
    % run data through random forest models
    oxy_RF_Arc_temp = runRF('All','Arc','Arc',EN4.all,index,1:11);
    oxy_RF_Pac_temp = runRF('All','Arc','Pac',EN4.all,index,1:11,idx_Pac);
    oxy_RF_Atl_temp = runRF('All','Arc','Atl',EN4.all,index,1:11,idx_Atl);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_Arc_temp .* weights(index) + ...
        oxy_RF_Pac_temp .* weights_Pac(index) + ...
        oxy_RF_Atl_temp .* weights_Atl(index);
    % run data through neural networks
    oxy_NN_Arc_temp = runNN('All','Arc','Arc',EN4.all,index,1:11);
    oxy_NN_Pac_temp = runNN('All','Arc','Pac',EN4.all,index,1:11,idx_Pac);
    oxy_NN_Atl_temp = runNN('All','Arc','Atl',EN4.all,index,1:11,idx_Atl);
    % calculate weighted averages
    oxy_NN.(basins{b}) = ...
        oxy_NN_Arc_temp .* weights(index) + ...
        oxy_NN_Pac_temp .* weights_Pac(index) + ...
        oxy_NN_Atl_temp .* weights_Atl(index);
clear oxy_NN_Arc_temp oxy_NN_Pac_temp oxy_NN_Atl_temp
clear oxy_RF_Arc_temp oxy_RF_Pac_temp oxy_RF_Atl_temp
clear idx_Pac idx_Atl weights weights_Pac weights_Atl

% For Mediterranean basin
elseif strcmp(basins{b},'Med')
    % run data through random forest model
    oxy_RF.(basins{b}) = runRF('All','Med','Med',EN4.all,index,1:11);
    % run data through neural networks
    oxy_NN.(basins{b}) = runNN('All','Med','Med',EN4.all,index,1:11);

% For N. Southern basin
elseif strcmp(basins{b},'NSou')
    % initially set all basin weights equal to 1
    weights(index) = 1;
    % adjust overlaps with Pacific
    idx_Pac = index & EN4.latitude(:) > -35 & ...
        (EN4.longitude(:) < -70 | EN4.longitude(:) > 145);
    weights(idx_Pac) = (-25 - EN4.latitude(idx_Pac))./10;
    % adjust overlaps with Atlantic
    idx_Atl = index & EN4.latitude(:) > -35 & ...
        EN4.longitude(:) < 15 & EN4.longitude(:) > -70;
    weights(idx_Atl) = (-25 - EN4.latitude(idx_Atl))./10;
    % adjust overlaps with Indian
    idx_Ind = index & EN4.latitude(:) > -35 & ...
        EN4.longitude(:) < 145 & EN4.longitude(:) > 15;
    weights(idx_Ind) = (-25 - EN4.latitude(idx_Ind))./10;
    % adjust overlaps with S. Southern
    idx_SSou = index & EN4.latitude(:) < -50;
    weights(idx_SSou) = (EN4.latitude(idx_SSou) + 60)./10;
    % calculate Pacific weights
    weights_Pac = calculateWeights(weights,idx_Pac);
    % calculate Atlantic weights
    weights_Atl = calculateWeights(weights,idx_Atl);
    % calculate Indian weights
    weights_Ind = calculateWeights(weights,idx_Ind);
    % calculate S. Southern weights
    weights_SSou = calculateWeights(weights,idx_SSou);
    % run data through random forest models
    oxy_RF_NSou_temp = runRF('All','NSou','NSou',EN4.all,index,1:11);
    oxy_RF_Pac_temp = runRF('All','NSou','Pac',EN4.all,index,1:11,idx_Pac);
    oxy_RF_Atl_temp = runRF('All','NSou','Atl',EN4.all,index,1:11,idx_Atl);
    oxy_RF_Ind_temp = runRF('All','NSou','Ind',EN4.all,index,1:11,idx_Ind);
    oxy_RF_SSou_temp = runRF('All','NSou','SSou',EN4.all,index,1:11,idx_SSou);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_NSou_temp .* weights(index) + ...
        oxy_RF_Pac_temp .* weights_Pac(index) + ...
        oxy_RF_Atl_temp .* weights_Atl(index) + ...
        oxy_RF_Ind_temp .* weights_Ind(index) + ...
        oxy_RF_SSou_temp .* weights_SSou(index);
    % run data through models
    oxy_NN_NSou_temp = runNN('All','NSou','NSou',EN4.all,index,1:11);
    oxy_NN_Pac_temp = runNN('All','NSou','Pac',EN4.all,index,1:11,idx_Pac);
    oxy_NN_Atl_temp = runNN('All','NSou','Atl',EN4.all,index,1:11,idx_Atl);
    oxy_NN_Ind_temp = runNN('All','NSou','Ind',EN4.all,index,1:11,idx_Ind);
    oxy_NN_SSou_temp = runNN('All','NSou','SSou',EN4.all,index,1:11,idx_SSou);
    % calculate weighted averages
    oxy_NN.(basins{b}) = ...
        oxy_NN_NSou_temp .* weights(index) + ...
        oxy_NN_Pac_temp .* weights_Pac(index) + ...
        oxy_NN_Atl_temp .* weights_Atl(index) + ...
        oxy_NN_Ind_temp .* weights_Ind(index) + ...
        oxy_NN_SSou_temp .* weights_SSou(index);
clear oxy_NN_NSou_temp oxy_NN_Pac_temp
clear oxy_RF_NSou_temp oxy_RF_Pac_temp oxy_NN_SSou_temp
clear oxy_NN_Atl_temp oxy_NN_Ind_temp
clear oxy_RF_Atl_temp oxy_RF_Ind_temp oxy_RF_SSou_temp
clear idx_Pac idx_Atl idx_Ind idx_SSou
clear weights weights_Pac weights_Atl weights_Ind weights_SSou
clear oxy_RF_NSou_SSou

% For S. Southern basin
elseif strcmp(basins{b},'SSou')
    % initially set all basin weights equal to 1
    weights(index) = 1;
    % adjust overlaps with N. Southern
    idx_NSou = index & EN4.latitude(:) > -60;
    weights(idx_NSou) = (-50 - EN4.latitude(idx_NSou))./10;
    % calculate N. Southern weights
    weights_NSou = calculateWeights(weights,idx_NSou);
    % run data through random forest models
    oxy_RF_SSou_temp = runRF('All','SSou','SSou',EN4.all,index,1:11);
    oxy_RF_NSou_temp = runRF('All','SSou','NSou',EN4.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_SSou_temp .* weights(index) + ...
        oxy_RF_NSou_temp .* weights_NSou(index);
    % run data through models
    oxy_NN_SSou_temp = runNN('All','SSou','SSou',EN4.all,index,1:11);
    oxy_NN_NSou_temp = runNN('All','SSou','NSou',EN4.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_NN.(basins{b}) = ...
        oxy_NN_SSou_temp .* weights(index) + ...
        oxy_NN_NSou_temp .* weights_NSou(index);
clear oxy_NN_NSou_temp oxy_NN_SSou_temp
clear oxy_RF_NSou_temp oxy_RF_SSou_temp
clear idx_NSou weights weights_NSou

end

% average predictions (ENS)
oxy_ENS.(basins{b}) = mean([oxy_RF.(basins{b}),oxy_NN.(basins{b})],2);

% reshape results
EN4.oxy_RF(index) = oxy_RF.(basins{b});
EN4.oxy_NN(index) = oxy_NN.(basins{b});
EN4.oxy_ENS(index) = oxy_ENS.(basins{b});

end

% Zero tran for ensemble average
EN4.oxy_ENS(EN4.oxy_ENS<0) = 0;

% clean up
clear b index oxy_RF oxy_NN oxy_ENS
EN4 = rmfield(EN4,'all');

% calculate percent O2 saturation
%EN4.oxy_sat = o2satv2b(EN4.sal,EN4.temp);
%EN4.oxy_ENS_per_sat = (EN4.oxy_ENS./EN4.oxy_sat) .* 100;
keyboard
% save chunks of data
if ~isfolder('Data'); mkdir('Data'); end
% save(['Data/EN4_' num2str(n)],'EN4','-v7.3');
save('Data/EN4_clim','EN4','-v7.3');
clear EN4

% disp([num2str(2003+n) ' finished']);
disp(['climatology finished']);

clear n basins basin_index

%% save product
% if ~isfolder('Data'); mkdir('Data'); end
% save('Data/EN4','EN4','-v7.3');
% clear

%% Plot O2 at 20 dbar
% figure; worldmap(latlim,lonlim);
% title('Oxygen at 20 dbar','fontsize',16)
% %set(gcf,'Position',[617, 599, 820, 820])
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',12);
% pcolorm(double(EN4.latitude(:,:,1,1)),double(EN4.longitude(:,:,1,1)),...
%     mean(EN4.oxy_ENS(:,:,3,:),4));
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar; caxis([0 350]);
% colormap(cmocean('haline',15));
% c.Label.String = 'Dissolved Oxygen (\mumol kg^{-1})';
% c.FontSize = 12;
% mlabel off; plabel off;
% exportgraphics(gcf,'Figures/Surface Plots/oxy_20_dbar_NN.jpg');

%% Plot O2 at 300 dbar
% figure; worldmap(latlim,lonlim);
% title('Oxygen at 300 dbar','fontsize',16)
% %set(gcf,'Position',[617, 599, 820, 820])
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',12);
% pcolorm(double(EN4.latitude(:,:,1,1)),double(EN4.longitude(:,:,1,1)),...
%     mean(EN4.oxy_ENS(:,:,25,:),4));
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar; caxis([0 350]);
% colormap(cmocean('haline',15));
% c.Label.String = 'Dissolved Oxygen (\mumol kg^{-1})';
% c.FontSize = 12;
% mlabel off; plabel off;
% exportgraphics(gcf,'Figures/Surface Plots/oxy_300_dbar_NN.jpg');

%% Plot O2 at 1000 dbar
% figure; worldmap(latlim,lonlim);
% title('Oxygen at 1000 dbar','fontsize',16)
% %set(gcf,'Position',[617, 599, 820, 820])
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',12);
% pcolorm(double(EN4.latitude(:,:,1,1)),double(EN4.longitude(:,:,1,1)),...
%     mean(EN4.oxy_ENS(:,:,44,:),4));
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar;
% colormap(cmocean('haline'));
% caxis([0 350]);
% c.Label.String = 'Dissolved Oxygen (\mumol kg^{-1})';
% c.FontSize = 12;
% exportgraphics(gcf,'Figures/Surface Plots/oxy_1000_dbar_ENS.jpg');

%% Plot O2 variability at 20 dbar
% figure; worldmap(latlim,lonlim);
% title('Oxygen variability at 20 dbar','fontsize',16)
% %set(gcf,'Position',[617, 599, 820, 820])
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',12);
% pcolorm(double(EN4.latitude(:,:,1,1)),double(EN4.longitude(:,:,1,1)),...
%     std(EN4.oxy_ENS(:,:,3,:),[],4));
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar;
% colormap(cmocean('thermal'));
% %caxis([200 300]);
% c.Label.String = 'Dissolved Oxygen variability (\mumol kg^{-1})';
% c.FontSize = 12;
% exportgraphics(gcf,'Figures/Surface Plots/oxy_var_20_dbar_ENS.jpg');

%% Plot O2 variability at 300 dbar
% figure; worldmap([20 44],[-140 -110]);
% title('Jon''s O2 std. dev. at 300 dbar','fontsize',16)
% %set(gcf,'Position',[617, 599, 820, 820])
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',12);
% pcolorm(double(EN4.latitude(:,:,1,1)),double(EN4.longitude(:,:,1,1)),...
%     std(EN4.oxy_ENS(:,:,25,:),[],4));
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar;
% colormap(cmocean('balance'));
% caxis([0 40]);
% c.Label.String = 'Dissolved Oxygen variability (\mumol kg^{-1})';
% c.FontSize = 12;
% exportgraphics(gcf,'Figures/Surface Plots/oxy_var_300_dbar_ENS.jpg');

%% Plot deseasonalized O2 variability at 300 dbar
% figure; worldmap([20 44],[-140 -110]);
% title('Jon''s deseasonalized O2 std. dev. at 300 dbar','fontsize',16)
% %set(gcf,'Position',[617, 599, 820, 820])
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',12);
% oxy_ENS_deseas = permute(cell2mat(permute(EN4.oxy_ENS_resid_notrend(:,:,25),[3 2 1])),[3 2 1]);
% pcolorm(double(EN4.latitude(:,:,1,1)),double(EN4.longitude(:,:,1,1)),...
%     std(oxy_ENS_deseas,[],3));
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar;
% colormap(cmocean('balance'));
% caxis([0 40]);
% c.Label.String = 'Dissolved Oxygen variability (\mumol kg^{-1})';
% c.FontSize = 12;
% exportgraphics(gcf,'Figures/Surface Plots/oxy_deseas_var_300_dbar_ENS.jpg');

%% Plot O2 variability at 1000 dbar
% figure; worldmap(latlim,lonlim);
% title('Oxygen variability at 1000 dbar','fontsize',16)
% %set(gcf,'Position',[617, 599, 820, 820])
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',12);
% pcolorm(double(EN4.latitude(:,:,1,1)),double(EN4.longitude(:,:,1,1)),...
%     std(EN4.oxy_ENS(:,:,44,:),[],4));
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar;
% colormap(cmocean('thermal'));
% %caxis([0 50]);
% c.Label.String = 'Dissolved Oxygen variability (\mumol kg^{-1})';
% c.FontSize = 12;
% exportgraphics(gcf,'Figures/Surface Plots/oxy_var_1000_dbar_ENS.jpg');
