function fit_NN(type,OXY,basin_index,preds_oxy_NN)

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
    
    % network sizes
    sizes1 = [20 15 10];
    sizes2 = [10 15 20];

    for n = 1:3

        rng(3); % for reproducibility
        % Neural network
        trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation
        hiddenLayer1Size = sizes1(n); % size of first hidden layer
        hiddenLayer2Size = sizes2(n); % size of second hidden layer
        % Create a Fitting Network
        if n < 4
            net = feedforwardnet([hiddenLayer1Size hiddenLayer2Size],trainFcn);
        else
            net = feedforwardnet(hiddenLayer1Size,trainFcn);
        end
        % Set training parameter criteria
        net.trainParam.max_fail = 6; % default: 6
        net.trainParam.mu_max = 1e+10; % default: 1e+10
        net.trainParam.min_grad = 1e-7; % default: 1e-7
        net.trainParam.epochs = 500; % default: 1000
        % Setup Division of Data for Training, Validation, Testing
        net.divideParam.trainRatio = 70/100;
        net.divideParam.valRatio = 15/100;
        net.divideParam.testRatio = 15/100;
        % Train the Network
        if strcmp(type,'All')
            model_oxy_NN = train(net,table2array(OXY.all(idx_train,preds_oxy_NN))',...
                table2array(OXY.all(idx_train,end))');
        else
            model_oxy_NN = train(net,table2array(OXY.train(idx_train,preds_oxy_NN))',...
                table2array(OXY.train(idx_train,end))');
        end
        % weird error was fixed by removing cpt_map toolbox from path 

        % save networks
        if n == 1
            % save and close model
            model_oxy_NN_1 = model_oxy_NN;
            if ~exist('Models','dir'); mkdir('Models'); end
            save(['Models/model_oxy_NN_1_' type '_' basins{b}],'model_oxy_NN_1','-v7.3');
            clear model_oxy_NN_1 model_oxy_NN
        elseif n == 2
            % save and close model
            model_oxy_NN_2 = model_oxy_NN;
            if ~exist('Models','dir'); mkdir('Models'); end
            save(['Models/model_oxy_NN_2_' type '_' basins{b}],'model_oxy_NN_2','-v7.3');
            clear model_oxy_NN_2 model_oxy_NN
        elseif n == 3
            % save and close model
            model_oxy_NN_3 = model_oxy_NN;
            if ~exist('Models','dir'); mkdir('Models'); end
            save(['Models/model_oxy_NN_3_' type '_' basins{b}],'model_oxy_NN_3','-v7.3');
            clear model_oxy_NN_3 model_oxy_NN
        end

    end

end