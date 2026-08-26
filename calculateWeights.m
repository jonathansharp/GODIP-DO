% calculate weights for an ancillary basin index from overall weights
function weights_basin = calculateWeights(weights,idx_basin)

weights_basin = zeros(size(weights));
weights_basin(idx_basin) = 1 - weights(idx_basin);