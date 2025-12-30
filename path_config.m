function fpaths = path_config(base_grid,param)

% model path
fpaths.model_path = '/fast7/model/';

% temperature and salinity paths
fpaths.temp_path = '/fast7/EN4.2.2/';
fpaths.sal_path = '/fast7/EN4.2.2/';

% parameter path
if strcmp(param,'o2')
    fpaths.param_path = '/fast4/o2/domip/';
elseif strcmp(param,'no3')
    fpaths.param_path = '/fast5/no3/domip/';
elseif strcmp(param,'dic')
    fpaths.param_path = '/fast6/dic/domip/';
end
