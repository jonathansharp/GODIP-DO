%% Import Roemmich-Gilson T/S product (and satellite chlorophyll)
% importRG;
% RG = rmfield(RG,{'temp_mean' 'temp_anom' 'sal_mean' 'sal_anom'});
% if ~isfolder('Data'); mkdir('Data'); end
% save('Data/RG','RG','-v7.3'); clear RG

%% Break RG product into yearly parts
% load('Data/RG');
% RG_all = RG; clear RG;
% for n = 1:size(RG_all.temp,4)/12
%     RG.latitude = RG_all.latitude;
%     RG.longitude = RG_all.longitude;
%     RG.pressure = RG_all.pressure;
%     RG.temp = RG_all.temp(:,:,:,(n-1)*12+1:(n-1)*12+12);
%     RG.sal = RG_all.sal(:,:,:,(n-1)*12+1:(n-1)*12+12);
%     RG.time = RG_all.time((n-1)*12+1:(n-1)*12+12);
%     if ~isfolder('Data'); mkdir('Data'); end
%     save(['Data/RG_' num2str(n)],'RG','-v7.3');
%     clear RG
% end
% clear n

%% Climatological RG average
load('Data/RG');
RG_all = RG; clear RG;
RG.latitude = RG_all.latitude;
RG.longitude = RG_all.longitude;
RG.pressure = RG_all.pressure;
RG.temp = nan([size(RG_all.temp,[1 2 3]),12]);
RG.sal = nan([size(RG_all.sal,[1 2 3]),12]);
RG.time = nan(12,1);
for m = 1:12
    RG.temp(:,:,:,m) = mean(RG_all.temp(:,:,:,m:12:end),4,'omitnan');
    RG.sal(:,:,:,m) = mean(RG_all.sal(:,:,:,m:12:end),4,'omitnan');
    RG.time(m) = mean(RG_all.time(m:12:end),'omitnan');
end
save('Data/RG_clim','RG','-v7.3');
clear RG

%% Create basin index

createShapes
basins = {'Atl' 'Pac' 'Ind' 'Arc' 'Med' 'NSou' 'SSou'};

% expand latitude and longitude to 2-D
RG_all.latitude = repmat(RG_all.latitude',length(RG_all.longitude),1);
RG_all.longitude = repmat(RG_all.longitude,1,size(RG_all.latitude,2));
RG_all.longitude(RG_all.longitude>180) = ...
    RG_all.longitude(RG_all.longitude>180) - 360;

for b = 1:7
    basin_index.RG.(basins{b}) = single(zeros(size(RG_all.temp)));
    if strcmp(basins{b},'Pac')
        basin_index.RG.(basins{b}) = ...
            inpolygon(RG_all.longitude,RG_all.latitude,...
            Poly.([basins{b} '1'])(:,1),Poly.([basins{b} '1'])(:,2)) | ...
            inpolygon(RG_all.longitude,RG_all.latitude,...
            Poly.([basins{b} '2'])(:,1),Poly.([basins{b} '2'])(:,2));
    else
        basin_index.RG.(basins{b}) = ...
            inpolygon(RG_all.longitude,RG_all.latitude,...
            Poly.(basins{b})(:,1),Poly.(basins{b})(:,2));
    end
end

clear RG_all b Poly

%% Predict O2 on RG grid 

% for n = 1:20 % (in yearly parts)


% load(['Data/RG_' num2str(n)]);
load('Data/RG_clim');

% convert from month since 1/1/2004 to datenum
RG.time = datenum([repmat(2004,size(RG.time)) double(RG.time)+0.5 repmat(15,size(RG.time))]);
date = datevec(RG.time);
date0 = date;
date0(:,2:3) = 1;
RG.day = (datenum(date) - datenum(date0));
RG.year = date(:,1);
clear date date0

% convert longitude
RG.longitude(RG.longitude>180) = RG.longitude(RG.longitude>180)-360;

% convert to singles to save memory
RG.temp = single(RG.temp);
RG.sal = single(RG.sal);
RG.day = single(RG.day);
RG.year = single(RG.year);

% expand latitude, longitude, pressure, and time to 4-D
RG.latitude = repmat(RG.latitude',length(RG.longitude),1,...
    length(RG.pressure),length(RG.time));
RG.longitude = repmat(RG.longitude,1,size(RG.latitude,2),...
    length(RG.pressure),length(RG.time));
RG.pressure = repmat(permute(RG.pressure,[3 2 1]),size(RG.longitude,1),...
    size(RG.latitude,2),1,length(RG.time));
RG.year = repmat(permute(RG.year,[4 3 2 1]),size(RG.longitude,1),...
    size(RG.latitude,2),size(RG.pressure,3),1);
RG.day = repmat(permute(RG.day,[4 3 2 1]),size(RG.longitude,1),...
    size(RG.latitude,2),size(RG.pressure,3),1);

% replace pressure with NaN where no T/S values exist
RG.pressure(isnan(RG.temp)) = NaN;

% ***** replace data with constant values to remove certain controls *****
% RG.year(~isnan(RG.year)) = mean(RG.year(:),'omitnan');
% RG.day(~isnan(RG.day)) = 182;
% RG.sal(~isnan(RG.sal)) = mean(RG.sal(:),'omitnan');
% RG.temp(~isnan(RG.temp)) = mean(RG.temp(:),'omitnan');
% RG.latitude(~isnan(RG.latitude)) = 0;
% RG.longitude(~isnan(RG.longitude)) = 0;

% Calculate absolute salinity, conservative temperature, potential density,
% and spice from climatologies (do this using a loop for each timestep, to
% avoid memory over-usage)
RG.sal_abs = single(nan(size(RG.temp)));
RG.temp_cns = single(nan(size(RG.temp)));
RG.sigma0 = single(nan(size(RG.temp)));
% RG.spice = single(nan(size(RG.temp)));
for t = 1:length(RG.time)
    RG.sal_abs(:,:,:,t) = gsw_SA_from_SP(RG.sal(:,:,:,t),RG.pressure(:,:,:,t),RG.longitude(:,:,:,t),RG.latitude(:,:,:,t));
    RG.temp_cns(:,:,:,t) = gsw_CT_from_t(RG.sal_abs(:,:,:,t),RG.temp(:,:,:,t),RG.pressure(:,:,:,t));
    RG.sigma0(:,:,:,t) = gsw_sigma0(RG.sal_abs(:,:,:,t),RG.temp_cns(:,:,:,t));
    % RG.spice(:,:,:,t) = gsw_spiciness0(RG.sal_abs(:,:,:,t),RG.temp_cns(:,:,:,t));
end
clear t

% Transform day by sine and cosine:
RG.day_sin = sin((2.*pi.*RG.day)./366);
RG.day_cos = cos((2.*pi.*RG.day)./366);
RG = rmfield(RG,'day');

% Transform longitude by sine and cosine:
RG.lon_cos1 = cosd(RG.longitude-20);
RG.lon_cos2 = cosd(RG.longitude-110);

% calculate bottom depth
RG.bottom_depth = ...
    single(bottom_depth(RG.latitude(:,:,1,1),RG.longitude(:,:,1,1)));
RG.bottom_depth = repmat(RG.bottom_depth,1,1,...
    size(RG.latitude,3),size(RG.latitude,4));

RG.oxy_ENS = nan(size(RG.temp));
RG.oxy_RF = nan(size(RG.temp));
RG.oxy_NN = nan(size(RG.temp));

%% Loop through and make predictions for each basin

for b = 1:7

idx = basin_index.RG.(basins{b});
idx = repmat(idx,1,1,size(RG.oxy_ENS,3),size(RG.oxy_ENS,4));

% Index to where P,T,S are available and matches basin index
index = ~isnan(RG.pressure(:)) & idx(:);
clear idx

% Pre-allocate weights
weights = zeros(size(index));

% Concatenate predictors 
RG.all = table(RG.latitude(:),RG.lon_cos1(:),RG.lon_cos2(:),...
    RG.sigma0(:),RG.pressure(:),RG.year(:),RG.day_sin(:),...
    RG.day_cos(:),RG.temp_cns(:),RG.sal_abs(:),RG.bottom_depth(:),...
    'VariableNames',{'latitude','lon_cos1','lon_cos2','sigma','pressure',...
    'year','day_sin','day_cos','temperature_cns','salinity_abs','bottom_depth'});

% For Atlantic basin
if strcmp(basins{b},'Atl')
    % initially set all basin weights equal to 1
    weights(index) = 1;
    % adjust overlaps with Arctic
    idx_Arc = index & RG.latitude(:) > 35;
    weights(idx_Arc) = (45 - RG.latitude(idx_Arc))./10;
    % adjust overlaps with Southern
    idx_NSou = index & RG.latitude(:) < -25;
    weights(idx_NSou) = (RG.latitude(idx_NSou) + 35)./10;
    % calculate Arctic weights
    weights_Arc = calculateWeights(weights,idx_Arc);
    % calculate N. Southern weights
    weights_NSou = calculateWeights(weights,idx_NSou);
    % run data through random forest models
    oxy_RF_Atl_temp = runRF('All','Atl','Atl',RG.all,index,1:11);
    oxy_RF_Arc_temp = runRF('All','Atl','Arc',RG.all,index,1:11,idx_Arc);
    oxy_RF_NSou_temp = runRF('All','Atl','NSou',RG.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_Atl_temp .* weights(index) + ...
        oxy_RF_Arc_temp .* weights_Arc(index) + ...
        oxy_RF_NSou_temp .* weights_NSou(index);
    % run data through neural networks
    oxy_NN_Atl_temp = runNN('All','Atl','Atl',RG.all,index,1:11);
    oxy_NN_Arc_temp = runNN('All','Atl','Arc',RG.all,index,1:11,idx_Arc);
    oxy_NN_NSou_temp = runNN('All','Atl','NSou',RG.all,index,1:11,idx_NSou);
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
    idx_Arc = index & RG.latitude(:) > 64;
    weights(idx_Arc) = (70 - RG.latitude(idx_Arc))./6;
    % adjust overlaps with Southern
    idx_NSou = index & RG.latitude(:) < -25;
    weights(idx_NSou) = (RG.latitude(idx_NSou) + 35)./10;
    % calculate Arctic weights
    weights_Arc = calculateWeights(weights,idx_Arc);
    % calculate N. Southern weights
    weights_NSou = calculateWeights(weights,idx_NSou);
    % run data through random forest models
    oxy_RF_Pac_temp = runRF('All','Pac','Pac',RG.all,index,1:11);
    oxy_RF_Arc_temp = runRF('All','Pac','Arc',RG.all,index,1:11,idx_Arc);
    oxy_RF_NSou_temp = runRF('All','Pac','NSou',RG.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_Pac_temp .* weights(index) + ...
        oxy_RF_Arc_temp .* weights_Arc(index) + ...
        oxy_RF_NSou_temp .* weights_NSou(index);
    % run data through neural networks
    oxy_NN_Pac_temp = runNN('All','Pac','Pac',RG.all,index,1:11);
    oxy_NN_Arc_temp = runNN('All','Pac','Arc',RG.all,index,1:11,idx_Arc);
    oxy_NN_NSou_temp = runNN('All','Pac','NSou',RG.all,index,1:11,idx_NSou);
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
    idx_NSou = index & RG.latitude(:) < -35;
    weights(idx_NSou) = (RG.latitude(idx_NSou) + 45)./10;
    % calculate N. Southern weights
    weights_NSou = calculateWeights(weights,idx_NSou);
    % run data through random forest models
    oxy_RF_Ind_temp = runRF('All','Ind','Ind',RG.all,index,1:11);
    oxy_RF_NSou_temp = runRF('All','Ind','NSou',RG.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_Ind_temp .* weights(index) + ...
        oxy_RF_NSou_temp .* weights_NSou(index);
    % run data through neural networks
    oxy_NN_Ind_temp = runNN('All','Ind','Ind',RG.all,index,1:11);
    oxy_NN_NSou_temp = runNN('All','Ind','NSou',RG.all,index,1:11,idx_NSou);
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
    idx_Pac = index & RG.latitude(:) < 70 & RG.longitude(:) < -120;
    weights(idx_Pac) = (RG.latitude(idx_Pac) - 64)./6;
    % adjust overlaps with Atlantic
    idx_Atl = index & RG.latitude(:) < 45 & RG.longitude(:) > -120;
    weights(idx_Atl) = (RG.latitude(idx_Atl) - 35)./10;
    % calculate Pacific weights
    weights_Pac = calculateWeights(weights,idx_Pac);
    % calculate Atlantic weights
    weights_Atl = calculateWeights(weights,idx_Atl);
    % run data through random forest models
    oxy_RF_Arc_temp = runRF('All','Arc','Arc',RG.all,index,1:11);
    oxy_RF_Pac_temp = runRF('All','Arc','Pac',RG.all,index,1:11,idx_Pac);
    oxy_RF_Atl_temp = runRF('All','Arc','Atl',RG.all,index,1:11,idx_Atl);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_Arc_temp .* weights(index) + ...
        oxy_RF_Pac_temp .* weights_Pac(index) + ...
        oxy_RF_Atl_temp .* weights_Atl(index);
    % run data through neural networks
    oxy_NN_Arc_temp = runNN('All','Arc','Arc',RG.all,index,1:11);
    oxy_NN_Pac_temp = runNN('All','Arc','Pac',RG.all,index,1:11,idx_Pac);
    oxy_NN_Atl_temp = runNN('All','Arc','Atl',RG.all,index,1:11,idx_Atl);
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
    oxy_RF.(basins{b}) = runRF('All','Med','Med',RG.all,index,1:11);
    % run data through neural networks
    oxy_NN.(basins{b}) = runNN('All','Med','Med',RG.all,index,1:11);

% For N. Southern basin
elseif strcmp(basins{b},'NSou')
    % initially set all basin weights equal to 1
    weights(index) = 1;
    % adjust overlaps with Pacific
    idx_Pac = index & RG.latitude(:) > -35 & ...
        (RG.longitude(:) < -70 | RG.longitude(:) > 145);
    weights(idx_Pac) = (-25 - RG.latitude(idx_Pac))./10;
    % adjust overlaps with Atlantic
    idx_Atl = index & RG.latitude(:) > -35 & ...
        RG.longitude(:) < 15 & RG.longitude(:) > -70;
    weights(idx_Atl) = (-25 - RG.latitude(idx_Atl))./10;
    % adjust overlaps with Indian
    idx_Ind = index & RG.latitude(:) > -35 & ...
        RG.longitude(:) < 145 & RG.longitude(:) > 15;
    weights(idx_Ind) = (-25 - RG.latitude(idx_Ind))./10;
    % adjust overlaps with S. Southern
    idx_SSou = index & RG.latitude(:) < -50;
    weights(idx_SSou) = (RG.latitude(idx_SSou) + 60)./10;
    % calculate Pacific weights
    weights_Pac = calculateWeights(weights,idx_Pac);
    % calculate Atlantic weights
    weights_Atl = calculateWeights(weights,idx_Atl);
    % calculate Indian weights
    weights_Ind = calculateWeights(weights,idx_Ind);
    % calculate S. Southern weights
    weights_SSou = calculateWeights(weights,idx_SSou);
    % run data through random forest models
    oxy_RF_NSou_temp = runRF('All','NSou','NSou',RG.all,index,1:11);
    oxy_RF_Pac_temp = runRF('All','NSou','Pac',RG.all,index,1:11,idx_Pac);
    oxy_RF_Atl_temp = runRF('All','NSou','Atl',RG.all,index,1:11,idx_Atl);
    oxy_RF_Ind_temp = runRF('All','NSou','Ind',RG.all,index,1:11,idx_Ind);
    oxy_RF_SSou_temp = runRF('All','NSou','SSou',RG.all,index,1:11,idx_SSou);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_NSou_temp .* weights(index) + ...
        oxy_RF_Pac_temp .* weights_Pac(index) + ...
        oxy_RF_Atl_temp .* weights_Atl(index) + ...
        oxy_RF_Ind_temp .* weights_Ind(index) + ...
        oxy_RF_SSou_temp .* weights_SSou(index);
    % run data through models
    oxy_NN_NSou_temp = runNN('All','NSou','NSou',RG.all,index,1:11);
    oxy_NN_Pac_temp = runNN('All','NSou','Pac',RG.all,index,1:11,idx_Pac);
    oxy_NN_Atl_temp = runNN('All','NSou','Atl',RG.all,index,1:11,idx_Atl);
    oxy_NN_Ind_temp = runNN('All','NSou','Ind',RG.all,index,1:11,idx_Ind);
    oxy_NN_SSou_temp = runNN('All','NSou','SSou',RG.all,index,1:11,idx_SSou);
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
    idx_NSou = index & RG.latitude(:) > -60;
    weights(idx_NSou) = (-50 - RG.latitude(idx_NSou))./10;
    % calculate N. Southern weights
    weights_NSou = calculateWeights(weights,idx_NSou);
    % run data through random forest models
    oxy_RF_SSou_temp = runRF('All','SSou','SSou',RG.all,index,1:11);
    oxy_RF_NSou_temp = runRF('All','SSou','NSou',RG.all,index,1:11,idx_NSou);
    % calculate weighted averages
    oxy_RF.(basins{b}) = ...
        oxy_RF_SSou_temp .* weights(index) + ...
        oxy_RF_NSou_temp .* weights_NSou(index);
    % run data through models
    oxy_NN_SSou_temp = runNN('All','SSou','SSou',RG.all,index,1:11);
    oxy_NN_NSou_temp = runNN('All','SSou','NSou',RG.all,index,1:11,idx_NSou);
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
RG.oxy_RF(index) = oxy_RF.(basins{b});
RG.oxy_NN(index) = oxy_NN.(basins{b});
RG.oxy_ENS(index) = oxy_ENS.(basins{b});

end

% Zero tran for ensemble average
RG.oxy_ENS(RG.oxy_ENS<0) = 0;

% clean up
clear b index oxy_RF oxy_NN oxy_ENS
RG = rmfield(RG,'all');

% calculate percent O2 saturation
%RG.oxy_sat = o2satv2b(RG.sal,RG.temp);
%RG.oxy_ENS_per_sat = (RG.oxy_ENS./RG.oxy_sat) .* 100;

% save chunks of data
if ~isfolder('Data'); mkdir('Data'); end
% save(['Data/RG_' num2str(n)],'RG','-v7.3');
save('Data/RG_clim','RG','-v7.3');
clear RG

% disp([num2str(2003+n) ' finished']);
disp(['climatology finished']);

% end

clear n basins basin_index

%% Concatenate product (uses up to 50GB memory)

% % load yearly chunks
% extra_vars = {'year' 'sal_abs' 'temp_cns' 'day_sin' ...
%     'day_cos' 'lon_cos1' 'lon_cos2' 'bottom_depth'};
% load('Data/RG_1'); RG_1=RG; clear RG; RG_1 = rmfield(RG_1,extra_vars);
% load('Data/RG_2'); RG_2=RG; clear RG; RG_2 = rmfield(RG_2,extra_vars);
% load('Data/RG_3'); RG_3=RG; clear RG; RG_3 = rmfield(RG_3,extra_vars);
% load('Data/RG_4'); RG_4=RG; clear RG; RG_4 = rmfield(RG_4,extra_vars);
% load('Data/RG_5'); RG_5=RG; clear RG; RG_5 = rmfield(RG_5,extra_vars);
% load('Data/RG_6'); RG_6=RG; clear RG; RG_6 = rmfield(RG_6,extra_vars);
% load('Data/RG_7'); RG_7=RG; clear RG; RG_7 = rmfield(RG_7,extra_vars);
% load('Data/RG_8'); RG_8=RG; clear RG; RG_8 = rmfield(RG_8,extra_vars);
% load('Data/RG_9'); RG_9=RG; clear RG; RG_9 = rmfield(RG_9,extra_vars);
% load('Data/RG_10'); RG_10=RG; clear RG; RG_10 = rmfield(RG_10,extra_vars);
% load('Data/RG_11'); RG_11=RG; clear RG; RG_11 = rmfield(RG_11,extra_vars);
% load('Data/RG_12'); RG_12=RG; clear RG; RG_12 = rmfield(RG_12,extra_vars);
% load('Data/RG_13'); RG_13=RG; clear RG; RG_13 = rmfield(RG_13,extra_vars);
% load('Data/RG_14'); RG_14=RG; clear RG; RG_14 = rmfield(RG_14,extra_vars);
% load('Data/RG_15'); RG_15=RG; clear RG; RG_15 = rmfield(RG_15,extra_vars);
% load('Data/RG_16'); RG_16=RG; clear RG; RG_16 = rmfield(RG_16,extra_vars);
% load('Data/RG_17'); RG_17=RG; clear RG; RG_17 = rmfield(RG_17,extra_vars);
% load('Data/RG_18'); RG_18=RG; clear RG; RG_18 = rmfield(RG_18,extra_vars);
% load('Data/RG_19'); RG_19=RG; clear RG; RG_19 = rmfield(RG_19,extra_vars);
% load('Data/RG_20'); RG_20=RG; clear RG; RG_20 = rmfield(RG_20,extra_vars);
% 
% % add latitude, longitude, and pressure
% RG.longitude = squeeze(RG_1.longitude(:,1,1,1));
% RG.latitude = squeeze(RG_1.latitude(1,:,1,1))';
% RG.pressure = squeeze(RG_1.pressure(1,1,:,1));
% 
% % concatenate
% vars = fieldnames(RG_1);
% for n = 4:length(vars)
%     if ~strcmp((vars{n}),'time')
%         RG.(vars{n}) = ...
%             cat(4,RG_1.(vars{n}),RG_2.(vars{n}),RG_3.(vars{n}),...
%                   RG_4.(vars{n}),RG_5.(vars{n}),RG_6.(vars{n}),...
%                   RG_7.(vars{n}),RG_8.(vars{n}),RG_9.(vars{n}),...
%                   RG_10.(vars{n}),RG_11.(vars{n}),RG_12.(vars{n}),...
%                   RG_13.(vars{n}),RG_14.(vars{n}),RG_15.(vars{n}),...
%                   RG_16.(vars{n}),RG_17.(vars{n}),RG_18.(vars{n}),...
%                   RG_19.(vars{n}),RG_20.(vars{n}));
%     else
%         RG.(vars{n}) = ...
%             cat(1,RG_1.(vars{n}),RG_2.(vars{n}),RG_3.(vars{n}),...
%                   RG_4.(vars{n}),RG_5.(vars{n}),RG_6.(vars{n}),...
%                   RG_7.(vars{n}),RG_8.(vars{n}),RG_9.(vars{n}),...
%                   RG_10.(vars{n}),RG_11.(vars{n}),RG_12.(vars{n}),...
%                   RG_13.(vars{n}),RG_14.(vars{n}),RG_15.(vars{n}),...
%                   RG_16.(vars{n}),RG_17.(vars{n}),RG_18.(vars{n}),...
%                   RG_19.(vars{n}),RG_20.(vars{n}));
%     end
%     RG_1 = rmfield(RG_1,(vars{n}));
%     RG_2 = rmfield(RG_2,(vars{n}));
%     RG_3 = rmfield(RG_3,(vars{n}));
%     RG_4 = rmfield(RG_4,(vars{n}));
%     RG_5 = rmfield(RG_5,(vars{n}));
%     RG_6 = rmfield(RG_6,(vars{n}));
%     RG_7 = rmfield(RG_7,(vars{n}));
%     RG_8 = rmfield(RG_8,(vars{n}));
%     RG_9 = rmfield(RG_9,(vars{n}));
%     RG_10 = rmfield(RG_10,(vars{n}));
%     RG_11 = rmfield(RG_11,(vars{n}));
%     RG_12 = rmfield(RG_12,(vars{n}));
%     RG_13 = rmfield(RG_13,(vars{n}));
%     RG_14 = rmfield(RG_14,(vars{n}));
%     RG_15 = rmfield(RG_15,(vars{n}));
%     RG_16 = rmfield(RG_16,(vars{n}));
%     RG_17 = rmfield(RG_17,(vars{n}));
%     RG_18 = rmfield(RG_18,(vars{n}));
%     RG_19 = rmfield(RG_19,(vars{n}));
%     RG_20 = rmfield(RG_20,(vars{n}));
% end
% clear RG_1 RG_2 RG_3 RG_4 RG_5 RG_6 RG_7 RG_8 RG_9 RG_10
% clear RG_11 RG_12 RG_13 RG_14 RG_15 RG_16 RG_17 RG_18 RG_19 RG_20

%% save product
% if ~isfolder('Data'); mkdir('Data'); end
% save('Data/RG','RG','-v7.3');
% clear

%% Plot O2 at 20 dbar
% figure; worldmap(latlim,lonlim);
% title('Oxygen at 20 dbar','fontsize',16)
% %set(gcf,'Position',[617, 599, 820, 820])
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',12);
% pcolorm(double(RG.latitude(:,:,1,1)),double(RG.longitude(:,:,1,1)),...
%     mean(RG.oxy_ENS(:,:,3,:),4));
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
% pcolorm(double(RG.latitude(:,:,1,1)),double(RG.longitude(:,:,1,1)),...
%     mean(RG.oxy_ENS(:,:,25,:),4));
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
% pcolorm(double(RG.latitude(:,:,1,1)),double(RG.longitude(:,:,1,1)),...
%     mean(RG.oxy_ENS(:,:,44,:),4));
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
% pcolorm(double(RG.latitude(:,:,1,1)),double(RG.longitude(:,:,1,1)),...
%     std(RG.oxy_ENS(:,:,3,:),[],4));
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
% pcolorm(double(RG.latitude(:,:,1,1)),double(RG.longitude(:,:,1,1)),...
%     std(RG.oxy_ENS(:,:,25,:),[],4));
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
% oxy_ENS_deseas = permute(cell2mat(permute(RG.oxy_ENS_resid_notrend(:,:,25),[3 2 1])),[3 2 1]);
% pcolorm(double(RG.latitude(:,:,1,1)),double(RG.longitude(:,:,1,1)),...
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
% pcolorm(double(RG.latitude(:,:,1,1)),double(RG.longitude(:,:,1,1)),...
%     std(RG.oxy_ENS(:,:,44,:),[],4));
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar;
% colormap(cmocean('thermal'));
% %caxis([0 50]);
% c.Label.String = 'Dissolved Oxygen variability (\mumol kg^{-1})';
% c.FontSize = 12;
% exportgraphics(gcf,'Figures/Surface Plots/oxy_var_1000_dbar_ENS.jpg');
