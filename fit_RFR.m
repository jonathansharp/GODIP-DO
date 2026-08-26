function train_sum = fit_RFR(type,OXY,basin_index,all_oxy_RF,importance)

% define basins
[basins] = basin_split(type);

% fit model
for b = 1:length(basins)
    
    % basin index for training
    if strcmp(type,'All')
        idx_train = basin_index.all.(basins{b});
    else
        idx_train = basin_index.train.(basins{b});
    end

    train_sum.(basins{b}) = sum(idx_train(:));
    
    % fit random forest
    numtrees = 600;
    minLeafSize = 10;
    NumPredictors = 6;
    if importance == 1
        oob_imp = 'on';
    else
        oob_imp = 'off';
    end
    % set up parallel pool
    % parpool;
    paroptions = statset('UseParallel',true);
    rng(2); % for reproducibility
    tic % start timing fit
    if strcmp(type,'All')
        model_oxy_RF = ...
            TreeBagger(numtrees,OXY.all(idx_train,all_oxy_RF),'oxygen',...
            'Method','regression',...
            'OOBPrediction','off',...
            'OOBPredictorImportance',oob_imp,...
            'MinLeafSize',minLeafSize,...
            'NumPredictorsToSample',NumPredictors);
    else
        model_oxy_RF = ...
            TreeBagger(numtrees,OXY.train(idx_train,all_oxy_RF),'oxygen',...
            'Method','regression',...
            'OOBPrediction','off',...
            'OOBPredictorImportance','off',...
            'MinLeafSize',minLeafSize,...
            'NumPredictorsToSample',NumPredictors);
    end
    toc % stop timing fit
    
    % end parallel session
    % delete(gcp('nocreate'))
    
    % Plot Out-of-Bag RMSE based on tree number
    % figure;
    % plot(sqrt(oobError(model_oxy_RF)),'k','LineWidth',2);
    % xlabel('Number of Grown Trees');
    % ylabel('Out-of-Bag Root Mean Squared Error');
    % close
    
    % Plot importance of each predictor
    if importance == 1
        figure;
        bar(model_oxy_RF.OOBPermutedPredictorDeltaError);
        xlabel('Predictor') ;
        ylabel('Out-of-Bag Feature Importance');
        xticklabels(OXY.all.Properties.VariableNames(all_oxy_RF(1:end-1)));
        exportgraphics(gcf,['Figures/RFR_importance_' type '_' basins{b} '.png']);
        close
    end
    
    % save and close model
    if ~exist('Models','dir'); mkdir('Models'); end
    save(['Models/model_oxy_RF_' type '_' basins{b}],'model_oxy_RF','-v7.3');
    clear model_oxy_RF
    
end