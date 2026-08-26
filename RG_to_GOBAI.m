% convert output from UW server to user-friendly GOBAI (.nc)

% Load RG file
load('Data/RG')
correct_RG;
% Load Uncertainty file
load('Data/Uncertainty')
% create filename
fname = ['Data/GOBAI-O2-' date '-v2.2'];

%% .nc (oxygen)

% create file and variables
if exist([fname '.nc'],'file'); delete([fname '.nc']); end
nccreate([fname '.nc'],'oxy','Format','64bit',...
    'Dimensions',{'lon',length(RG.longitude),'lat',...
    length(RG.latitude),'pres',length(RG.pressure),'time',...
    length(RG.time)},'Datatype','single');
nccreate([fname '.nc'],'lon','Dimensions',{'lon'},'Datatype','single');
nccreate([fname '.nc'],'lat','Dimensions',{'lat'},'Datatype','single');
nccreate([fname '.nc'],'pres','Dimensions',{'pres'},'Datatype','single');
nccreate([fname '.nc'],'time','Dimensions',{'time'},'Datatype','single');
nccreate([fname '.nc'],'temp','Dimensions',{'lon',...
    'lat','pres','time'},'Datatype','single');
nccreate([fname '.nc'],'sal','Dimensions',{'lon',...
    'lat','pres','time'},'Datatype','single');
nccreate([fname '.nc'],'uncer','Dimensions',{'lon',...
    'lat','pres','time'},'Datatype','single');

% write variables and attributes

ncwrite([fname '.nc'],'lon',RG.longitude);
ncwriteatt([fname '.nc'],'lon','long_name','Longitude');
ncwriteatt([fname '.nc'],'lon','units','degrees east');
ncwriteatt([fname '.nc'],'lon','source','grid');

ncwrite([fname '.nc'],'lat',RG.latitude);
ncwriteatt([fname '.nc'],'lat','long_name','Latitude');
ncwriteatt([fname '.nc'],'lat','units','degrees north');
ncwriteatt([fname '.nc'],'lat','source','grid');

ncwrite([fname '.nc'],'pres',RG.pressure);
ncwriteatt([fname '.nc'],'pres','long_name','Pressure');
ncwriteatt([fname '.nc'],'pres','units','decibars');
ncwriteatt([fname '.nc'],'pres','source','grid');

ncwrite([fname '.nc'],'time',RG.time-datenum(1950,1,1));
ncwriteatt([fname '.nc'],'time','long_name','Time');
ncwriteatt([fname '.nc'],'time','units','days since 1950-01-01');
ncwriteatt([fname '.nc'],'time','source','grid');

ncwrite([fname '.nc'],'oxy',RG.oxy_ENS);
ncwriteatt([fname '.nc'],'oxy','long_name','Dissolved Oxygen Concentration');
ncwriteatt([fname '.nc'],'oxy','units','micromoles per kilogram');
ncwriteatt([fname '.nc'],'oxy','source','machine learning scheme decribed in Sharp et al. (submitted)');

ncwrite([fname '.nc'],'temp',RG.temp);
ncwriteatt([fname '.nc'],'temp','long_name','Temperature');
ncwriteatt([fname '.nc'],'temp','units','degrees Celcius');
ncwriteatt([fname '.nc'],'temp','source','Roemmich and Gilson, 2009 (https://sio-argo.ucsd.edu/RG_Climatology.html)');

ncwrite([fname '.nc'],'sal',RG.sal);
ncwriteatt([fname '.nc'],'sal','long_name','Salinity');
ncwriteatt([fname '.nc'],'sal','units','n/a');
ncwriteatt([fname '.nc'],'sal','source','Roemmich and Gilson, 2009 (https://sio-argo.ucsd.edu/RG_Climatology.html)');

ncwrite([fname '.nc'],'uncer',Uncer.tot);
ncwriteatt([fname '.nc'],'uncer','long_name','Total Uncertainty');
ncwriteatt([fname '.nc'],'uncer','units','micromoles per kilogram');
ncwriteatt([fname '.nc'],'uncer','source','produced by combining three separate uncertainty sources in quadrature, as decribed in Sharp et al. (submitted)');

% clean up
clear all
