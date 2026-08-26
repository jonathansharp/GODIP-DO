% plot land areas
function plot_land(type,clr)

if strcmp(type,'xy')

    land = shaperead('landareas', 'UseGeoCoords',false);
    land(2).X(land(2).X>-150) = land(2).X(land(2).X>-150)-360;
    land(163).X(land(163).X>-150) = land(163).X(land(163).X>-150)-360;
    mapshow(land,'FaceColor',clr);

elseif strcmp(type,'map')

    land = shaperead('landareas', 'UseGeoCoords',true);
    geoshow(land,'FaceColor',clr);

end

clear land