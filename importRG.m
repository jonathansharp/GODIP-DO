%% Import Roemmich-Gilson Climatology
path = '/raid/sharp/matlab/O2_modelling/For_Hercules/RG_CLIM/RG_ArgoClim_Temperature_2019.nc';
RG.latitude = ncread(path,'LATITUDE');
RG.longitude = ncread(path,'LONGITUDE');
RG.pressure = ncread(path,'PRESSURE');
RG.time = ncread(path,'TIME');
RG.temp_anom = ncread(path,'ARGO_TEMPERATURE_ANOMALY');
RG.temp_mean = ncread(path,'ARGO_TEMPERATURE_MEAN');
path = '/raid/sharp/matlab/O2_modelling/For_Hercules/RG_CLIM/RG_ArgoClim_Salinity_2019.nc';
RG.sal_anom  = ncread(path,'ARGO_SALINITY_ANOMALY');
RG.sal_mean  = ncread(path,'ARGO_SALINITY_MEAN');

% Add temperature anomaly extensions
path = '/raid/sharp/matlab/O2_modelling/For_Hercules/RG_CLIM/';
RG.TEMP_ANOM_1901 = ncread([path 'RG_ArgoClim_201901_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1902 = ncread([path 'RG_ArgoClim_201902_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1903 = ncread([path 'RG_ArgoClim_201903_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1904 = ncread([path 'RG_ArgoClim_201904_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1905 = ncread([path 'RG_ArgoClim_201905_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1906 = ncread([path 'RG_ArgoClim_201906_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1907 = ncread([path 'RG_ArgoClim_201907_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1908 = ncread([path 'RG_ArgoClim_201908_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1909 = ncread([path 'RG_ArgoClim_201909_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1910 = ncread([path 'RG_ArgoClim_201910_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1911 = ncread([path 'RG_ArgoClim_201911_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_1912 = ncread([path 'RG_ArgoClim_201912_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2001 = ncread([path 'RG_ArgoClim_202001_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2002 = ncread([path 'RG_ArgoClim_202002_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2003 = ncread([path 'RG_ArgoClim_202003_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2004 = ncread([path 'RG_ArgoClim_202004_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2005 = ncread([path 'RG_ArgoClim_202005_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2006 = ncread([path 'RG_ArgoClim_202006_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2007 = ncread([path 'RG_ArgoClim_202007_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2008 = ncread([path 'RG_ArgoClim_202008_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2009 = ncread([path 'RG_ArgoClim_202009_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2010 = ncread([path 'RG_ArgoClim_202010_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2011 = ncread([path 'RG_ArgoClim_202011_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2012 = ncread([path 'RG_ArgoClim_202012_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2101 = ncread([path 'RG_ArgoClim_202101_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2102 = ncread([path 'RG_ArgoClim_202102_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2103 = ncread([path 'RG_ArgoClim_202103_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2104 = ncread([path 'RG_ArgoClim_202104_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2105 = ncread([path 'RG_ArgoClim_202105_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2106 = ncread([path 'RG_ArgoClim_202106_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2107 = ncread([path 'RG_ArgoClim_202107_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2108 = ncread([path 'RG_ArgoClim_202108_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2109 = ncread([path 'RG_ArgoClim_202109_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2110 = ncread([path 'RG_ArgoClim_202110_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2111 = ncread([path 'RG_ArgoClim_202111_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2112 = ncread([path 'RG_ArgoClim_202112_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2201 = ncread([path 'RG_ArgoClim_202201_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2202 = ncread([path 'RG_ArgoClim_202202_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2203 = ncread([path 'RG_ArgoClim_202203_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2204 = ncread([path 'RG_ArgoClim_202204_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2205 = ncread([path 'RG_ArgoClim_202205_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2206 = ncread([path 'RG_ArgoClim_202206_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2207 = ncread([path 'RG_ArgoClim_202207_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2208 = ncread([path 'RG_ArgoClim_202208_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2209 = ncread([path 'RG_ArgoClim_202209_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2210 = ncread([path 'RG_ArgoClim_202210_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2211 = ncread([path 'RG_ArgoClim_202211_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2212 = ncread([path 'RG_ArgoClim_202212_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2301 = ncread([path 'RG_ArgoClim_202301_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2302 = ncread([path 'RG_ArgoClim_202302_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2303 = ncread([path 'RG_ArgoClim_202303_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2304 = ncread([path 'RG_ArgoClim_202304_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2305 = ncread([path 'RG_ArgoClim_202305_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2306 = ncread([path 'RG_ArgoClim_202306_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2307 = ncread([path 'RG_ArgoClim_202307_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2308 = ncread([path 'RG_ArgoClim_202308_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2309 = ncread([path 'RG_ArgoClim_202309_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2310 = ncread([path 'RG_ArgoClim_202310_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2311 = ncread([path 'RG_ArgoClim_202311_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.TEMP_ANOM_2312 = ncread([path 'RG_ArgoClim_202312_2019.nc'],'ARGO_TEMPERATURE_ANOMALY');
RG.temp_anom = cat(4,RG.temp_anom,RG.TEMP_ANOM_1901,RG.TEMP_ANOM_1902,...
                     RG.TEMP_ANOM_1903,RG.TEMP_ANOM_1904,RG.TEMP_ANOM_1905,...
                     RG.TEMP_ANOM_1906,RG.TEMP_ANOM_1907,RG.TEMP_ANOM_1908,...
                     RG.TEMP_ANOM_1909,RG.TEMP_ANOM_1910,RG.TEMP_ANOM_1911,...
                     RG.TEMP_ANOM_1912,RG.TEMP_ANOM_2001,RG.TEMP_ANOM_2002,...
                     RG.TEMP_ANOM_2003,RG.TEMP_ANOM_2004,RG.TEMP_ANOM_2005,...
                     RG.TEMP_ANOM_2006,RG.TEMP_ANOM_2007,RG.TEMP_ANOM_2008,...
                     RG.TEMP_ANOM_2009,RG.TEMP_ANOM_2010,RG.TEMP_ANOM_2011,...
                     RG.TEMP_ANOM_2012,RG.TEMP_ANOM_2101,RG.TEMP_ANOM_2102,...
                     RG.TEMP_ANOM_2103,RG.TEMP_ANOM_2104,RG.TEMP_ANOM_2105,...
                     RG.TEMP_ANOM_2106,RG.TEMP_ANOM_2107,RG.TEMP_ANOM_2108,...
                     RG.TEMP_ANOM_2109,RG.TEMP_ANOM_2110,RG.TEMP_ANOM_2111,...
                     RG.TEMP_ANOM_2112,RG.TEMP_ANOM_2201,RG.TEMP_ANOM_2202,...
                     RG.TEMP_ANOM_2203,RG.TEMP_ANOM_2204,RG.TEMP_ANOM_2205,...
                     RG.TEMP_ANOM_2206,RG.TEMP_ANOM_2207,RG.TEMP_ANOM_2208,...
                     RG.TEMP_ANOM_2209,RG.TEMP_ANOM_2210,RG.TEMP_ANOM_2211,...
                     RG.TEMP_ANOM_2212,RG.TEMP_ANOM_2301,RG.TEMP_ANOM_2302,...
                     RG.TEMP_ANOM_2303,RG.TEMP_ANOM_2304,RG.TEMP_ANOM_2305,...
                     RG.TEMP_ANOM_2306,RG.TEMP_ANOM_2307,RG.TEMP_ANOM_2308,...
                     RG.TEMP_ANOM_2309,RG.TEMP_ANOM_2310,RG.TEMP_ANOM_2311,...
                     RG.TEMP_ANOM_2312);
% Add extended months
RG.time = cat(1,RG.time,(RG.time(end)+1:1:RG.time(end)+60)');
% Determine absolute temperatures from mean and anomalies
RG.temp = nan(size(RG.temp_anom));
for n = 1:size(RG.time,1)
RG.temp(:,:,:,n) = RG.temp_anom(:,:,:,n) + RG.temp_mean;
end
         
% Add salinity anomaly extensions
RG.SAL_ANOM_1901 = ncread([path 'RG_ArgoClim_201901_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1902 = ncread([path 'RG_ArgoClim_201902_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1903 = ncread([path 'RG_ArgoClim_201903_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1904 = ncread([path 'RG_ArgoClim_201904_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1905 = ncread([path 'RG_ArgoClim_201905_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1906 = ncread([path 'RG_ArgoClim_201906_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1907 = ncread([path 'RG_ArgoClim_201907_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1908 = ncread([path 'RG_ArgoClim_201908_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1909 = ncread([path 'RG_ArgoClim_201909_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1910 = ncread([path 'RG_ArgoClim_201910_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1911 = ncread([path 'RG_ArgoClim_201911_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_1912 = ncread([path 'RG_ArgoClim_201912_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2001 = ncread([path 'RG_ArgoClim_202001_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2002 = ncread([path 'RG_ArgoClim_202002_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2003 = ncread([path 'RG_ArgoClim_202003_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2004 = ncread([path 'RG_ArgoClim_202004_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2005 = ncread([path 'RG_ArgoClim_202005_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2006 = ncread([path 'RG_ArgoClim_202006_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2007 = ncread([path 'RG_ArgoClim_202007_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2008 = ncread([path 'RG_ArgoClim_202008_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2009 = ncread([path 'RG_ArgoClim_202009_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2010 = ncread([path 'RG_ArgoClim_202010_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2011 = ncread([path 'RG_ArgoClim_202011_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2012 = ncread([path 'RG_ArgoClim_202012_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2101 = ncread([path 'RG_ArgoClim_202101_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2102 = ncread([path 'RG_ArgoClim_202102_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2103 = ncread([path 'RG_ArgoClim_202103_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2104 = ncread([path 'RG_ArgoClim_202104_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2105 = ncread([path 'RG_ArgoClim_202105_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2106 = ncread([path 'RG_ArgoClim_202106_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2107 = ncread([path 'RG_ArgoClim_202107_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2108 = ncread([path 'RG_ArgoClim_202108_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2109 = ncread([path 'RG_ArgoClim_202109_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2110 = ncread([path 'RG_ArgoClim_202110_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2111 = ncread([path 'RG_ArgoClim_202111_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2112 = ncread([path 'RG_ArgoClim_202112_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2201 = ncread([path 'RG_ArgoClim_202201_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2202 = ncread([path 'RG_ArgoClim_202202_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2203 = ncread([path 'RG_ArgoClim_202203_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2204 = ncread([path 'RG_ArgoClim_202204_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2205 = ncread([path 'RG_ArgoClim_202205_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2206 = ncread([path 'RG_ArgoClim_202206_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2207 = ncread([path 'RG_ArgoClim_202207_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2208 = ncread([path 'RG_ArgoClim_202208_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2209 = ncread([path 'RG_ArgoClim_202209_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2210 = ncread([path 'RG_ArgoClim_202210_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2211 = ncread([path 'RG_ArgoClim_202211_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2212 = ncread([path 'RG_ArgoClim_202212_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2301 = ncread([path 'RG_ArgoClim_202301_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2302 = ncread([path 'RG_ArgoClim_202302_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2303 = ncread([path 'RG_ArgoClim_202303_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2304 = ncread([path 'RG_ArgoClim_202304_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2305 = ncread([path 'RG_ArgoClim_202305_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2306 = ncread([path 'RG_ArgoClim_202306_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2307 = ncread([path 'RG_ArgoClim_202307_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2308 = ncread([path 'RG_ArgoClim_202308_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2309 = ncread([path 'RG_ArgoClim_202309_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2310 = ncread([path 'RG_ArgoClim_202310_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2311 = ncread([path 'RG_ArgoClim_202311_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.SAL_ANOM_2312 = ncread([path 'RG_ArgoClim_202312_2019.nc'],'ARGO_SALINITY_ANOMALY');
RG.sal_anom = cat(4,RG.sal_anom,RG.SAL_ANOM_1901,RG.SAL_ANOM_1902,...
                     RG.SAL_ANOM_1903,RG.SAL_ANOM_1904,RG.SAL_ANOM_1905,...
                     RG.SAL_ANOM_1906,RG.SAL_ANOM_1907,RG.SAL_ANOM_1908,...
                     RG.SAL_ANOM_1909,RG.SAL_ANOM_1910,RG.SAL_ANOM_1911,...
                     RG.SAL_ANOM_1912,RG.SAL_ANOM_2001,RG.SAL_ANOM_2002,...
                     RG.SAL_ANOM_2003,RG.SAL_ANOM_2004,RG.SAL_ANOM_2005,...
                     RG.SAL_ANOM_2006,RG.SAL_ANOM_2007,RG.SAL_ANOM_2008,...
                     RG.SAL_ANOM_2009,RG.SAL_ANOM_2010,RG.SAL_ANOM_2011,...
                     RG.SAL_ANOM_2012,RG.SAL_ANOM_2101,RG.SAL_ANOM_2102,...
                     RG.SAL_ANOM_2103,RG.SAL_ANOM_2104,RG.SAL_ANOM_2105,...
                     RG.SAL_ANOM_2106,RG.SAL_ANOM_2107,RG.SAL_ANOM_2108,...
                     RG.SAL_ANOM_2109,RG.SAL_ANOM_2110,RG.SAL_ANOM_2111,...
                     RG.SAL_ANOM_2112,RG.SAL_ANOM_2201,RG.SAL_ANOM_2202,...
                     RG.SAL_ANOM_2203,RG.SAL_ANOM_2204,RG.SAL_ANOM_2205,...
                     RG.SAL_ANOM_2206,RG.SAL_ANOM_2207,RG.SAL_ANOM_2208,...
                     RG.SAL_ANOM_2209,RG.SAL_ANOM_2210,RG.SAL_ANOM_2211,...
                     RG.SAL_ANOM_2212,RG.SAL_ANOM_2301,RG.SAL_ANOM_2302,...
                     RG.SAL_ANOM_2303,RG.SAL_ANOM_2304,RG.SAL_ANOM_2305,...
                     RG.SAL_ANOM_2306,RG.SAL_ANOM_2307,RG.SAL_ANOM_2308,...
                     RG.SAL_ANOM_2309,RG.SAL_ANOM_2310,RG.SAL_ANOM_2311,...
                     RG.SAL_ANOM_2312);
% Determine absolute salinities from mean and anomalies
RG.sal = nan(size(RG.sal_anom));
for n = 1:size(RG.time,1)
RG.sal(:,:,:,n) = RG.sal_anom(:,:,:,n) + RG.sal_mean;
end

%% Remove concatenated fields from structure
fields = fieldnames(RG);
idx = contains(fields,'_ANOM_');
RG = rmfield(RG,fields(idx));
clear path n idx fields

% %% Plot T at 20 dbar
% figure; worldmap([-90 90],[20 380]);
% title('Annual mean at 20 dbar (RG09)','fontsize',16)
% set(gcf,'Position',[617, 599, 820, 420])
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',12);
% pcolorm(double(RG.latitude),double(RG.longitude),mean(RG.temp(:,:,3,:),4)');
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar; caxis([0 30]);
% colormap(cmocean('thermal',12));
% c.Label.String = ['Temperature ' char(176) 'C'];
% c.FontSize = 12;
% c.TickLength = 0;
% mlabel off; plabel off;
% if ~isfolder('Figures'); mkdir('Figures'); end
% if ~isfolder('Figures/Surface_Plots'); mkdir('Figures/Surface_Plots'); end
% exportgraphics(gcf,'Figures/Surface_Plots/temp_20_dbar_RG.png');
% close
% 
% %% Plot S at 20 dbar
% figure; worldmap([-90 90],[20 380]);
% title('Annual mean at 20 dbar (RG09)','fontsize',16)
% set(gcf,'Position',[617, 599, 820, 420])
% setm(gca,'ffacecolor','w');
% setm(gca,'fontsize',12);
% pcolorm(double(RG.latitude),double(RG.longitude),mean(RG.sal(:,:,3,:),4)');
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(land,'FaceColor',rgb('grey'));
% c=colorbar; caxis([32 38]);
% colormap(cmocean('haline',12));
% c.Label.String = 'Salinity';
% c.FontSize = 12;
% mlabel off; plabel off;
% if ~isfolder('Figures'); mkdir('Figures'); end
% if ~isfolder('Figures/Surface_Plots'); mkdir('Figures/Surface_Plots'); end
% exportgraphics(gcf,'Figures/Surface_Plots/sal_20_dbar_RG.png');
% close
