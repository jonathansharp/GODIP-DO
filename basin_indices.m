function basin_index = basin_indices(type,OXY,Poly)

% define basins and data split types
[basins,data_types] = basin_split(type);

% split data into test and training
for b = 1:length(basins)
    for t = 3
        basin_index.(data_types{t}).(basins{b}) = zeros(size(OXY.(data_types{t}),1),1);
        if strcmp(basins{b},'Pac')
            basin_index.(data_types{t}).(basins{b}) = ...
                inpolygon(OXY.(data_types{t}).longitude,OXY.(data_types{t}).latitude,...
                Poly.([basins{b} '1'])(:,1),Poly.([basins{b} '1'])(:,2)) | ...
                inpolygon(OXY.(data_types{t}).longitude,OXY.(data_types{t}).latitude,...
                Poly.([basins{b} '2'])(:,1),Poly.([basins{b} '2'])(:,2));
        else
            basin_index.(data_types{t}).(basins{b}) = ...
                inpolygon(OXY.(data_types{t}).longitude,OXY.(data_types{t}).latitude,...
                Poly.(basins{b})(:,1),Poly.(basins{b})(:,2));
        end
    end
end