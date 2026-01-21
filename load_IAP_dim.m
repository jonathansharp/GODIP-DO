% load_IAP_dim
%
% DESCRIPTION:
% This function is used to load the 
% spatiotemporal dimensions of the IAP
% temperature and salinity data product.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/20/2026

function TS = load_IAP_dim(fpath,y1,y2)

% load temperature and salinity
TS.Longitude = ncread([fpath 'IAPv4_Temp_monthly_1_6000m_year_2025_month_09.nc'],'lon');
TS.Latitude = ncread([fpath 'IAPv4_Temp_monthly_1_6000m_year_2025_month_09.nc'],'lat');
TS.Depth = ncread([fpath 'IAPv4_Temp_monthly_1_6000m_year_2025_month_09.nc'],'depth_std');
% compute dimensions
TS.xdim = length(TS.Longitude);
TS.ydim = length(TS.Latitude);
TS.zdim = length(TS.Depth);
% determine number of monthly timesteps
files = dir(fpath);
cnt = 1;
for n = 1:length(files)
    if contains(files(n).name,'.nc')
        date = cell2mat(extractBetween(files(n).name,'year_','.nc'));
        TS.years(cnt) = str2double(date(1:4));
        TS.months(cnt) = str2double(date(12:13));
        cnt = cnt+1;
    end
end
TS.years = TS.years';
TS.months = TS.months';
TS.Time = datenum(TS.years,TS.months,15);
TS.tdim = length(TS.Time);