% script to download IAP data for GODIP-DO
url = 'http://www.ocean.iap.ac.cn/ftp/cheng/IAP_DO_data_19602025_July2026/IAP_DO_data_19602025_with_IAP_TS/';
for y = 1961:2025
    for m = 1:12
        try
            fname = ['CAS_DO_profiles_' num2str(y) '_' num2str(m) '.nc'];
            websave(['IAP_Profiles_Data/' fname],[url fname]);
        catch
            % do nothing
        end
    end
end
