function all_data = process_IAP_profile_data(y1,y2,vrs)

% load WOD profile data
nc_label = 'IAP';
folder = godip_data_path('IAP_Profiles_Data');
types = {'OSD' 'CTD' 'PFL'};

% for oxygen
vars_download = {'DOprofile_anomaly_stddepth' 'DOprofile_stddepth' 'DOprofile_info_record'...
    'DOprofile_info_string' 'DOprofile_info_record_name' 'DOprofile_info_string_name' 'stddepth' 'DOprofile_num' ...
    'CODC_depth_temperature'  'CODC_temperature' 'CODC_depth_salinity' 'CODC_salinity' 'matched_location'};
vars_o2 = {'cruise' 'profile' 'time' 'year' 'month' 'day' 'lat' 'lon' 'depth' 'Oxygen'};
vars_other = {'Temperature' 'Salinity'};
vars_both = [vars_o2 vars_other];
% pre-allocate
for v = 1:length(vars_both); all_data.(vars_both{v}) = []; end
all_data.type = []; all_data.mode = []; idx_mode = [];
% load oxygen variables
for y = y1:y2
    for m = 1:12
        file = [folder '/CAS_DO_profiles_' num2str(y) '_' num2str(m) '.nc'];
        if exist(file,'file')
            schema = ncinfo(file);
            pdim = schema.Dimensions(1).Length;
            zdim = schema.Dimensions(2).Length;
            % download data
            for v = 1:length(vars_download)
                all_data_temp.(vars_download{v}) = ncread(file,vars_download{v});
            end
            % iterpolate T/S/O2 to same depth levels
            all_data_temp.CODC_temperature_interp = nan(size(all_data_temp.DOprofile_stddepth));
            all_data_temp.CODC_salinity_interp = nan(size(all_data_temp.DOprofile_stddepth));
            for p = 1:all_data_temp.DOprofile_num
                % interpolate temperature
                idx_temp = find(~isnan(all_data_temp.CODC_depth_temperature(p,:)));
                if sum(idx_temp) > 1
                    all_data_temp.CODC_temperature_interp(p,:) = interp1(all_data_temp.CODC_depth_temperature(p,idx_temp)',...
                        all_data_temp.CODC_temperature(p,idx_temp)',all_data_temp.stddepth);
                end
                % interpolate salinity
                idx_temp = find(~isnan(all_data_temp.CODC_depth_salinity(p,:)));
                if sum(idx_temp) > 1
                    all_data_temp.CODC_salinity_interp(p,:) = interp1(all_data_temp.CODC_depth_salinity(p,idx_temp)',...
                        all_data_temp.CODC_salinity(p,idx_temp)',all_data_temp.stddepth);
                end
            end

            % process latitude, longitude, and time
            yr_idx = find(strcmp(all_data_temp.DOprofile_info_record_name,'Year'));
            mn_idx = find(strcmp(all_data_temp.DOprofile_info_record_name,'Month'));
            dy_idx = find(strcmp(all_data_temp.DOprofile_info_record_name,'Day'));
            lon_idx = find(strcmp(all_data_temp.DOprofile_info_record_name,'Longitude'));
            lat_idx = find(strcmp(all_data_temp.DOprofile_info_record_name,'Latitude'));
            all_data_temp.time = repmat(datenum([all_data_temp.DOprofile_info_record(:,yr_idx) ...
                    all_data_temp.DOprofile_info_record(:,mn_idx) ...
                    all_data_temp.DOprofile_info_record(:,dy_idx)]),...
                    1,length(all_data_temp.stddepth));
            all_data_temp.year = repmat(all_data_temp.DOprofile_info_record(:,yr_idx),1,length(all_data_temp.stddepth));
            all_data_temp.month = repmat(all_data_temp.DOprofile_info_record(:,mn_idx),1,length(all_data_temp.stddepth));
            all_data_temp.day = repmat(all_data_temp.DOprofile_info_record(:,dy_idx),1,length(all_data_temp.stddepth));
            all_data_temp.lat = repmat(all_data_temp.DOprofile_info_record(:,lat_idx),1,length(all_data_temp.stddepth));
            all_data_temp.lon = repmat(all_data_temp.DOprofile_info_record(:,lon_idx),1,length(all_data_temp.stddepth));
            all_data_temp.depth = repmat(all_data_temp.stddepth',all_data_temp.DOprofile_num,1);

            % process type, cruise, and profile
            type_idx = find(strcmp(all_data_temp.DOprofile_info_string_name,'ObservationType'));
            all_data_temp.type = nan(size(all_data_temp.DOprofile_stddepth));
            all_data_temp.type(strcmp(all_data_temp.DOprofile_info_string(:,type_idx),"OSD"),:) = 1;
            all_data_temp.type(strcmp(all_data_temp.DOprofile_info_string(:,type_idx),"CTD"),:) = 2;
            all_data_temp.type(strcmp(all_data_temp.DOprofile_info_string(:,type_idx),"ARG"),:) = 3;

            % process cruise
            cruise_idx = find(strcmp(all_data_temp.DOprofile_info_record_name,'Cruise_number(OSD,CTD)/WMO(ARG)'));
            all_data_temp.cruise = repmat(all_data_temp.DOprofile_info_record(:,cruise_idx),1,length(all_data_temp.stddepth));
            
            % process profile number
            stat_idx = find(strcmp(all_data_temp.DOprofile_info_record_name,'Stat_num(OSD,CTD)/Cycle_num(ARG)'));
            id_idx = find(strcmp(all_data_temp.DOprofile_info_record_name,'NODC_ID(OSD,CTD)/-9(ARG)'));
            all_data_temp.profile = nan(all_data_temp.DOprofile_num,1);
            % for OSD and CTD, assign NODC_ID as unique profile number
            all_data_temp.profile(all_data_temp.type(:,1) == 1 | all_data_temp.type(:,1) == 2) = ...
                all_data_temp.DOprofile_info_record(find(all_data_temp.type(:,1) == 1 | all_data_temp.type(:,1) == 2),id_idx);
            % for ARG, assign WMO+Cycle Number as unique profile number
            all_data_temp.profile(all_data_temp.type(:,1) == 3) = ...
                1000.*all_data_temp.DOprofile_info_record(all_data_temp.type(:,1) == 3,cruise_idx) + ...
                all_data_temp.DOprofile_info_record(all_data_temp.type(:,1) == 3,stat_idx);
            all_data_temp.profile = repmat(all_data_temp.profile,1,length(all_data_temp.stddepth));
            % indices to ensure all data points have O2, temp, and sal and are above 2250m
            idx_o2 = ~isnan(all_data_temp.DOprofile_stddepth);
            idx_temp = ~isnan(all_data_temp.CODC_temperature_interp);
            idx_sal = ~isnan(all_data_temp.CODC_salinity_interp);
            idx_depth = true(size(all_data_temp.depth)); % all_data_temp.depth <= 2250;
            idx_all = idx_o2 & idx_temp & idx_sal & idx_depth;

            % log data into large arrays
            all_data.cruise = [all_data.cruise;all_data_temp.cruise(idx_all)];
            all_data.profile = [all_data.profile;all_data_temp.profile(idx_all)];
            all_data.time = [all_data.time;all_data_temp.time(idx_all)];
            all_data.year = [all_data.year;all_data_temp.year(idx_all)];
            all_data.month = [all_data.month;all_data_temp.month(idx_all)];
            all_data.day = [all_data.day;all_data_temp.day(idx_all)];
            all_data.lat = [all_data.lat;all_data_temp.lat(idx_all)];
            all_data.lon = [all_data.lon;all_data_temp.lon(idx_all)];
            all_data.depth = [all_data.depth;all_data_temp.depth(idx_all)];
            all_data.Oxygen = [all_data.Oxygen;all_data_temp.DOprofile_stddepth(idx_all)];
            all_data.Temperature = [all_data.Temperature;all_data_temp.CODC_temperature_interp(idx_all)];
            all_data.Salinity = [all_data.Salinity;all_data_temp.CODC_salinity_interp(idx_all)];
            all_data.type = [all_data.type;all_data_temp.type(idx_all)];
            % all_data.mode = [all_data.mode;all_data_temp.(vars{v})(idx_all)];

        end
    end
end

% % test plot
% depths = unique(all_data.depth);
% for d = 1:length(depths)
%     idx = all_data.depth == depths(d);
%     o2(d) = mean(all_data.Oxygen(idx));
% end
% figure; plot(depths,o2);

% % remove real-time data
% idx_mode = logical(idx_mode);
% vars = fieldnames(all_data);
% for v = 1:length(vars)
%     all_data.(vars{v}) = all_data.(vars{v})(idx_mode);
% end

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

% % calculate float oxygen saturation
% Ts = log((298.15 - all_data.Temperature)./(273.15 + all_data.Temperature));
% % The coefficents below are from the second column of Table 1 of Garcia and
% % Gordon (1992)
% a0 =  5.80871; a1 =  3.20291; a2 =  4.17887; a3 =  5.10006;
% a4 = -9.86643e-2; a5 =  3.80369; b0 = -7.01577e-3; b1 = -7.70028e-3;
% b2 = -1.13864e-2; b3 = -9.51519e-3; c0 = -2.75915e-7;
% all_data.Oxygen_sat = exp(a0 + Ts.*(a1 + Ts.*(a2 + Ts.*(a3 + Ts.*(a4 + a5*Ts)))) + ...
%     all_data.Salinity.*(b0 + Ts.*(b1 + Ts.*(b2 + b3*Ts)) + c0*all_data.Salinity));
% all_data.Oxygen_sat_per = 100.*(all_data.Oxygen./all_data.Oxygen_sat);
% all_data.Oxygen_uncorrected = all_data.Oxygen;

% % index to float data
% idx_float = all_data.type == 3;
% idx_deep = all_data.pressure >= max_delta_pressure;
% idx_shal = all_data.pressure < max_delta_pressure;

% % adjust float data according to GOBAI-O2 deep and shallow float correction (June 2026)
% load(fullfile('/data2/GOBAI/o2/Data','float_corr_Jun-2026_D_A'));
% corr_fac = nan(size(all_data.Oxygen_sat));
% corr_fac(idx_float & idx_deep) = all_data.pressure(idx_float & idx_deep).*slp_deep + int_deep;
% corr_fac(idx_float & idx_shal) = all_data.pressure(idx_float & idx_shal).*slp_shal + int_shal;
% all_data.Oxygen(idx_float) = all_data.Oxygen(idx_float) - corr_fac;

% remove float data for test with ship data only
% vars = fieldnames(all_data);
% for v = 1:length(vars)
%     all_data.(vars{v})(idx_float) = [];
% end

% save data
save(godip_data_path('O2','Data',[vrs '_data_' num2str(y1) '_' num2str(y2)]),...
    'all_data','-v7.3');