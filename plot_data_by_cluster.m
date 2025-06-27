% plot_data_by_cluster
%
% DESCRIPTION:
% This function plots data points according to their most fitting cluster.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 6/10/2025

function plot_data_by_cluster(param_props,base_grid,file_date,...
    float_file_ext,num_clusters,numWorkers_predict,start_year,end_year)

%% plot data points by cluster
% load combined data
load(['O2/Data/wod_data_' num2str(start_year) '_' num2str(end_year) '.mat'],'all_data');
% load cluster data
load([param_props.dir_name '/Data/all_data_clusters_' base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');
% define pressure axis
pressures = sort(unique(all_data.depth));
% open parallel pool
tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;
% make plots
parfor p = 1:length(pressures)
    % plot data by cluster
    figure('visible','off','Position',[100 100 800 400]); hold on;
    set(gca,'fontsize',12);
    idx = all_data.depth == pressures(p);
    m_proj('robinson','lon',[20 380]);
    m_coast('patch',rgb('grey'));
    m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
    lon_temp = convert_lon(convert_lon(all_data.lon),'0-360');
    lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
    % use depth for CMIP models
    title(['Data by Cluster at ' num2str(pressures(p)) ' m (' base_grid ')'],'fontsize',16);
    m_scatter(lon_temp(idx),all_data.lat(idx),3,all_data_clusters.clusters(idx),'filled');
    mycolormap = [1,1,1;flipud(jet(num_clusters))]; % white then jet
    colormap(mycolormap); % white then jet
    clim([-0.5 num_clusters+0.5]);
    c=colorbar;
    c.Limits = [0.5 num_clusters+0.5];
    c.Label.String = 'Cluster';
    c.FontSize = 12;
    c.TickLength = 0;
    % save figure
    dname = [param_props.dir_name '/Figures/Clusters/' base_grid '_c' num2str(num_clusters)];
    if ~isfolder([pwd '/' dname]); mkdir(dname); end
    export_fig(gcf,[dname '/clustered_data_' num2str(pressures(p)) '.png'],...
        '-transparent','-silent');
    close
    % plot data by cluster probability
    for clst = 1:num_clusters
        figure('visible','off','Position',[100 100 800 400]); hold on;
        set(gca,'fontsize',12);
        idx = all_data.depth == pressures(p);
        m_proj('robinson','lon',[20 380]);
        m_coast('patch',rgb('grey'));
        m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
        lon_temp = convert_lon(convert_lon(all_data.lon),'0-360');
        lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
        % use depth for CMIP models
        title(['Data by Cluster at ' num2str(pressures(p)) ' m (' base_grid ')'],'fontsize',16);
        m_scatter(lon_temp(idx),all_data.latitude(idx),3,...
            all_data_clusters.(['c' num2str(clst)])(idx),'filled');
        colormap(customcolormap([0 1],[mycolormap(clst+1,:); 1 1 1]));
        clim([0 1]);
        c=colorbar;
        c.Limits = [0 1];
        c.Label.String = ['Cluster #' num2str(clst) ' Probability'];
        c.FontSize = 12;
        c.TickLength = 0;
        % save figure
        dname = [param_props.dir_name '/Figures/Clusters/' base_grid '_c' num2str(num_clusters)];
        if ~isfolder([pwd '/' dname]); mkdir(dname); end
        export_fig([dname '/clustered_data_probability_c' ...
            num2str(clst) '_' num2str(pressures(p)) '.png'],...
            '-transparent','-silent');
        close
    end
end

% end parallel session
delete(gcp('nocreate'));
