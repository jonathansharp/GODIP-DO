function mean_map_plot(h,x,y,z,climits,ttl,cbr)
    h=worldmap('world');
    title(h,ttl)
    pcolorm(x,y,z);
    clim(climits);
    colormap(cbr);
    plot_land('map');
    plabel off;
    mlabel off;
end