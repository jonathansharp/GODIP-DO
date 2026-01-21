function fpaths = path_config(base_grid,param)

% model path
fpaths.model_path = '/fast7/model/';

% temperature and salinity paths
if strcmp(base_grid,'EN4')
    fpaths.temp_path = '/fast7/EN4.2.2/';
    fpaths.sal_path = '/fast7/EN4.2.2/';
elseif strcmp(base_grid,'IAP')
    fpaths.temp_path = '/fast7/IAPv4_temp_monthly/';
    fpaths.sal_path = '/fast7/IAPv2_sal_monthly/';
end

% parameter path
if strcmp(param,'o2')
    fpaths.param_path = '/fast4/o2/GODIP-DO/';
elseif strcmp(param,'no3')
    fpaths.param_path = '/fast5/no3/GODIP-DO/';
elseif strcmp(param,'dic')
    fpaths.param_path = '/fast6/dic/GODIP-DO/';
end
