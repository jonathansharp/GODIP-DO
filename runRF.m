% run random forest model for specific basin
function oxy_RF = ...
    runRF(type,main_basin,sub_basin,oxy_data,idx_test,preds_oxy_RF,idx_sub)

if strcmp(main_basin,sub_basin)
    load(['Models/model_oxy_RF_' type '_' main_basin])
    % break into chunks if trying to pass a lot of data through the model
    if size(oxy_data,1) > 5e6
        chunk_size = round(size(oxy_data,1)./5e6,0);
        chunk = round(size(oxy_data,1)/chunk_size,0);
        oxy_RF = nan(sum(idx_test),1);
        for n = 1:chunk_size
            idx1 = chunk.*(n-1)+1;
            idx2 = chunk.*(n-1)+chunk;
            oxy_test_temp = oxy_data(idx1:idx2,:);
            idx_test_temp = idx_test(idx1:idx2);
            idx1 = find(isnan(oxy_RF),1,'first');
            idx2 = idx1+sum(idx_test_temp)-1;
            oxy_RF(idx1:idx2) = predict(model_oxy_RF,...
                oxy_test_temp(idx_test_temp,preds_oxy_RF));
        end
    else
        oxy_RF = predict(model_oxy_RF,oxy_data(idx_test,preds_oxy_RF));
    end
else
    if sum(idx_sub) > 0
        load(['Models/model_oxy_RF_' type '_' sub_basin])
        if size(oxy_data,1) > 5e6
            chunk_size = round(size(oxy_data,1)./5e6,0);
            chunk = round(size(oxy_data,1)/chunk_size,0);
            oxy_RF = zeros(sum(idx_test),1);
            oxy_RF_sub = nan(sum(idx_sub),1);
            for n = 1:chunk_size
                idx1 = chunk.*(n-1)+1;
                idx2 = chunk.*(n-1)+chunk;
                oxy_sub_temp = oxy_data(idx1:idx2,:);
                idx_sub_temp = idx_sub(idx1:idx2);
                idx1 = find(isnan(oxy_RF_sub),1,'first');
                idx2 = idx1+sum(idx_sub_temp)-1;
                oxy_RF_sub(idx1:idx2) = predict(model_oxy_RF,...
                    oxy_sub_temp(idx_sub_temp,preds_oxy_RF));
            end
            idx_sub_temp = idx_sub(idx_test);
            oxy_RF(idx_sub_temp) = oxy_RF_sub;
        else
            oxy_RF = zeros(size(idx_test));
            oxy_RF(idx_sub) = ...
                predict(model_oxy_RF,oxy_data(idx_sub,preds_oxy_RF));
            oxy_RF = oxy_RF(idx_test);
        end
    else
        oxy_RF = zeros(sum(idx_test),1);
    end
end

clear model_oxy_RF_val