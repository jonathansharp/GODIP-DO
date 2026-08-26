% load_model_dim
%
% DESCRIPTION:
% This function is used to load the 
% spatiotemporal dimensions of CMIP
% model output.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 9/16/2025

function [TS,months,weeks,timesteps] = load_model_dim(fpath)

% model info
inf = ncinfo(fpath);
% load dimensions
TS.Latitude = ncread(fpath,'y');
TS.Longitude = ncread(fpath,'x');
TS.Depth = ncread(fpath,'lev');
TS.Depth_Bounds = ncread(fpath,'lev_bnds');
% compute dimensions
TS.xdim = length(TS.Longitude);
TS.ydim = length(TS.Latitude);
TS.zdim = length(TS.Depth);
% determine number of monthly timesteps
TS.Time = ncread(fpath,'time');
timesteps = length(TS.Time);
months = 1:timesteps;
weeks = ones(timesteps,1);
% add years and months
time_att_idx = strcmp({inf.Variables.Name},'time');
cal_att_idx = strcmp({inf.Variables(time_att_idx).Attributes.Name},'calendar');
if strcmp(inf.Variables(time_att_idx).Attributes(cal_att_idx).Value,'360_day')
    % time on 360-day calendar
    dates = datevec(datenum(1850,1,1+(365.2425/360).*double(TS.Time)));
else
    % time on normal calendar
    dates = datevec(datenum(1850,1,1+double(TS.Time)));
end
TS.years = dates(:,1);
TS.months = dates(:,2);
clear dates
