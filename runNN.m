% run neural network for specific basin
function oxy_NN = ...
    runNN(type,main_basin,sub_basin,oxy_data,idx_test,preds_oxy_NN,idx_sub)

if strcmp(main_basin,sub_basin)
    load(['Models/model_oxy_NN_1_' type '_' main_basin])
    oxy_NN_1 = model_oxy_NN_1(table2array(oxy_data(idx_test,preds_oxy_NN))')';
    clear model_oxy_NN_1
    load(['Models/model_oxy_NN_2_' type '_' main_basin])
    oxy_NN_2 = model_oxy_NN_2(table2array(oxy_data(idx_test,preds_oxy_NN))')';
    clear model_oxy_NN_2
    load(['Models/model_oxy_NN_3_' type '_' main_basin])
    oxy_NN_3 = model_oxy_NN_3(table2array(oxy_data(idx_test,preds_oxy_NN))')';
    clear model_oxy_NN_3
    oxy_NN = mean([oxy_NN_1,oxy_NN_2,oxy_NN_3],2);
    clear oxy_NN_1 oxy_NN_2 oxy_NN_3
else
    if sum(idx_sub) > 0
        load(['Models/model_oxy_NN_1_' type '_' sub_basin])
        oxy_NN_1 = model_oxy_NN_1(table2array(oxy_data(idx_sub,preds_oxy_NN))')';
        clear model_oxy_NN_1
        load(['Models/model_oxy_NN_2_' type '_' sub_basin])
        oxy_NN_2 = model_oxy_NN_2(table2array(oxy_data(idx_sub,preds_oxy_NN))')';
        clear model_oxy_NN_2
        load(['Models/model_oxy_NN_3_' type '_' sub_basin])
        oxy_NN_3 = model_oxy_NN_3(table2array(oxy_data(idx_sub,preds_oxy_NN))')';
        clear model_oxy_NN_3
        oxy_NN = zeros(size(idx_test));
        oxy_NN(idx_sub) = mean([oxy_NN_1,oxy_NN_2,oxy_NN_3],2);
        clear oxy_NN_1 oxy_NN_2 oxy_NN_3
        oxy_NN = oxy_NN(idx_test);
    else
        oxy_NN = zeros(sum(idx_test),1);
    end
end
