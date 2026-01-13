function all_data = process_WOD_profile_data(y1,y2)

% load WOD profile data
folder = 'WOD25_Profiles_Data';
types = {'CTD' 'OSD' 'PFL'};

% for oxygen
vars_o2 = {'cruise' 'profile' 'time' 'year' 'month' 'day' 'lat' 'lon' 'depth' 'Oxygen'};
vars_other = {'Temperature' 'Salinity'};
vars_both = [vars_o2 vars_other];
% pre-allocate
for v = 1:length(vars_both); all_data.(vars_both{v}) = []; end
all_data.type = [];
% load oxygen variables
for y = y1:y2
    for x = 1:length(types)
        file = [folder '/Oxygen_' types{x} '_WOD25/Oxygen_' types{x} '_' num2str(y) '.nc'];
        if exist(file,'file')
            schema = ncinfo(file);
            pdim = schema.Dimensions(1).Length;
            zdim = schema.Dimensions(2).Length;
            % oxygen
            for v = 1:length(vars_o2)
                all_data_temp.(vars_o2{v}) = ncread(file,vars_o2{v});
            end
            % temp
            file = [folder '/Temperature_' types{x} '_WOD25/Temperature_' types{x} '_' num2str(y) '.nc'];
            all_data_temp.Temperature = ncread(file,'Temperature');
            all_data_temp.temp_profile = ncread(file,'profile');
            % sal
            file = [folder '/Salinity_' types{x} '_WOD25/Salinity_' types{x} '_' num2str(y) '.nc'];
            all_data_temp.Salinity = ncread(file,'Salinity');
            all_data_temp.sal_profile = ncread(file,'profile');
            % indices
            idx_o2 = ismember(all_data_temp.profile,all_data_temp.temp_profile) & ismember(all_data_temp.profile,all_data_temp.sal_profile);
            idx_temp = ismember(all_data_temp.temp_profile,all_data_temp.profile) & ismember(all_data_temp.temp_profile,all_data_temp.sal_profile);
            idx_sal = ismember(all_data_temp.sal_profile,all_data_temp.profile) & ismember(all_data_temp.sal_profile,all_data_temp.temp_profile);
            % filter to only profiles with oxygen, temperature, and salinity
            for v = 1:length(vars_o2)
                if strcmp(vars_o2{v},'depth')
                    all_data_temp.(vars_o2{v}) = repmat(all_data_temp.(vars_o2{v}),1,sum(idx_o2));
                elseif strcmp(vars_o2{v},'Oxygen')
                    all_data_temp.(vars_o2{v}) = all_data_temp.(vars_o2{v})(1:zdim,idx_o2);
                else
                    all_data_temp.(vars_o2{v}) = repmat(all_data_temp.(vars_o2{v})(idx_o2)',zdim,1);
                end
            end
            all_data_temp.Temperature = all_data_temp.Temperature(1:zdim,idx_temp);
            all_data_temp.Salinity = all_data_temp.Salinity(1:zdim,idx_sal);
            % add values to vector
            idx = ~isnan(all_data_temp.Oxygen);
            for v = 1:length(vars_both)
                all_data.(vars_both{v}) = [all_data.(vars_both{v});all_data_temp.(vars_both{v})(idx)];
            end
            % type: CTD=1 OSD=2 PFL=3
            all_data.type = [all_data.type;repmat(x,sum(idx(:)),1)];
        end
    end
end

% calculate ancillary parameters
all_data.pressure = gsw_p_from_z(-all_data.depth,all_data.lat);
all_data.Salinity_Abs = gsw_SA_from_SP(all_data.Salinity,all_data.pressure,all_data.lon,all_data.lat);
all_data.Temperature_Cons = gsw_CT_from_pt(all_data.Salinity_Abs,all_data.Temperature);
all_data.sigma = gsw_sigma0(all_data.Salinity_Abs,all_data.Temperature_Cons);

% transform longitude and day
all_data.lon_cos_1 = cosd(all_data.lon-20);
all_data.lon_cos_2 = cosd(all_data.lon-110);
date = datevec(datenum(double(all_data.year),double(all_data.month),double(all_data.day)));
date0 = date;
date0(:,2:3) = 0;
all_data.day_of_year = datenum(date) - datenum(date0);
all_data.day_sin = sin((2.*pi.*all_data.day_of_year)/365.25);
all_data.day_cos = cos((2.*pi.*all_data.day_of_year)/365.25);

% calculate float oxygen saturation
Ts = log((298.15 - all_data.Temperature)./(273.15 + all_data.Temperature)); 
% The coefficents below are from the second column of Table 1 of Garcia and
% Gordon (1992)
a0 =  5.80871; a1 =  3.20291; a2 =  4.17887; a3 =  5.10006;
a4 = -9.86643e-2; a5 =  3.80369; b0 = -7.01577e-3; b1 = -7.70028e-3;
b2 = -1.13864e-2; b3 = -9.51519e-3; c0 = -2.75915e-7;
all_data.Oxygen_sat = exp(a0 + Ts.*(a1 + Ts.*(a2 + Ts.*(a3 + Ts.*(a4 + a5*Ts)))) + ...
      all_data.Salinity.*(b0 + Ts.*(b1 + Ts.*(b2 + b3*Ts)) + c0*all_data.Salinity));
all_data.Oxygen_sat_per = 100.*(all_data.Oxygen./all_data.Oxygen_sat);
all_data.Oxygen_uncorrected = all_data.Oxygen;

% index to float data
idx_float = all_data.type == 3;

% adjust float data according to GOBAI-O2 float correction (July)
load('/home/sockeye/sharp/GOBAI/O2/Data/float_corr_Jul-2025_D');
all_data.Oxygen(idx_float) = double(((all_data.Oxygen_sat_per(idx_float) - ...
    (slp .* all_data.Oxygen_sat(idx_float) + int))./100) .* ...
    all_data.Oxygen_sat(idx_float));

% remove float data for test with ship data only
% vars = fieldnames(all_data);
% for v = 1:length(vars)
%     all_data.(vars{v})(idx_float) = [];
% end

% save data
save(['O2/Data/wod_data_' num2str(y1) '_' num2str(y2)],'all_data','-v7.3');

