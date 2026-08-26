% this function defines the types of basins and data splits for each data type

function [basins,data_types] = basin_split(type)

if strcmp(type,'Data')
    basins = {'Atl' 'Pac' 'Ind' 'Arc' 'Med' 'NSou' 'SSou'};
    data_types = {'train' 'test' 'all'};
elseif strcmp(type,'All')
    basins = {'Atl' 'Pac' 'Ind' 'Arc' 'Med' 'NSou' 'SSou'};
elseif strcmp(type,'WOA')
    basins = {'Atl' 'Pac' 'Ind' 'Arc' 'Med' 'NSou' 'SSou'};
    data_types = {'train' 'test'};
elseif strcmp(type,'BSOSE')
    basins = {'Atl' 'Pac' 'Ind' 'NSou' 'SSou'};
    data_types = {'train' 'test'};
elseif strcmp(type,'ESM4')
    basins = {'Atl' 'Pac' 'Ind' 'Arc' 'Med' 'NSou' 'SSou'};
    data_types = {'train' 'test'};
end