function [model_path,param_path,temp_path,sal_path] = path_config(base_grid,param)

% model path
model_path = '/fast7/sharp/model/';

% temperature and salinity paths
temp_path = '/fast7/EN4.2.2/';
sal_path = '/fast7/EN4.2.2/';

% parameter path
if strcmp(param,'o2')
    param_path = '/fast4/o2/domip/';
elseif strcmp(param,'no3')
    param_path = '/fast5/no3/domip/';
elseif strcmp(param,'dic')
    param_path = '/fast6/dic/domip/';
end
