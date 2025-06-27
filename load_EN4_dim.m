% load_EN4_dim
%
% DESCRIPTION:
% This function is used to load the 
% spatiotemporal dimensions of the EN4
% temperature and salinity data product.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 5/20/2025

function TS = load_EN4_dim(fpath,y1,y2)

% load temperature and salinity
TS.Longitude = ncread([fpath 'EN.4.2.2.f.analysis.c14.200001.nc'],'lon');
TS.Latitude = ncread([fpath 'EN.4.2.2.f.analysis.c14.200001.nc'],'lat');
TS.Depth = ncread([fpath 'EN.4.2.2.f.analysis.c14.200001.nc'],'depth');
% compute dimensions
TS.xdim = length(TS.Longitude);
TS.ydim = length(TS.Latitude);
TS.zdim = length(TS.Depth);
% determine number of monthly timesteps
files = dir(fpath);
cnt = 1;
for n = 1:length(files)
    if contains(files(n).name,'.nc')
        date = cell2mat(extractBetween(files(n).name,'c14.','.nc'));
        TS.years(cnt) = str2double(date(1:4));
        TS.months(cnt) = str2double(date(5:6));
        cnt = cnt+1;
    end
end
TS.years = TS.years';
TS.months = TS.months';
TS.Time = datenum(TS.years,TS.months,15);
TS.tdim = length(TS.Time);