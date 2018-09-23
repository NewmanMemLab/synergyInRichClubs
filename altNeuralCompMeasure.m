% This function implements the method described by Chichilnisky, 2001.
%
% Inputs:
%   raster: #neurons X #time bins (in ms) matrix with 1's in bins with spikes, and 
%   0's elsewhere
%
%   duration: width of time window before spike in neuron that you want to
%   consider (in bins)
%
% Outputs:
%   ratio: #neurons X 1 vector of ratio of sum of squared errors of nonlinear
%   fit to sum of squared errors of linear fit. Smaller values correspond to
%   a better nonlinear fit. Larger values correspond to a better linear
%   fit. For node classification, one can perform a median split, for
%   example.
%
% Written by John Beggs, Sam Faber, and Ehren Newman, 2016, Indiana
% University

function ratio = altNeuralCompMeasure(raster,duration)

% pre-allocate
ratio = nan(size(raster,1),1);

raster = single(raster);
[nNeurons, ~] = size(raster);

for iNode = 1:nNeurons % for each neuron
    
    % Step 1: Compute the spike-triggered average (STA) of the input spikes.
    % The sum of the inputs occuring in the appropriate preceding time bins
    % divided by the total number of spikes.
    
    J = find(raster(iNode, :) ~= 0);
    [Inds] = find(2*duration < J);
    J = J(Inds);
    NumNeuron1Spikes = length(J);
    STA = zeros(nNeurons, duration+1);
    TIMERASTERn1 = raster(iNode, :);
    raster(iNode, :) = 0; % exclude spikes from the neuron itself
    for ii=1:(duration+1)
        STA(:, ii) = sum(raster(:, J-duration-1+ii)')'/NumNeuron1Spikes; %probability of firing after iNode fired 
    end
    
    % Step 2: Compute the generator signal -- the STA multiplied by the sum of
    % the input spikes in the appropriate preceding time bins
  
    generator = nan(size(raster));
    for ii=1:nNeurons
        generator(ii,:) = conv(raster(ii,:),STA(ii,:),'same'); 
    end
    generator = sum(generator,1);
    
    % Step 3: Estimate the average spike count in time bins (average number of spikes/bin) with nearly equal
    % generator signals. Repeat this for many/all values of g.
    
    nBins = 100;
    bins = linspace(min(generator),max(generator),nBins);
    [nSamp, gsBin] = histc(generator,bins);
    nSamp = nSamp ./ size(generator,2);
    
    % add neuron of interest back into timeraster to include its
    % own spiking in the prediction
    raster(iNode,:) = TIMERASTERn1;
    
    for bin = 1:nBins
        n_g(bin) = mean(raster(iNode,gsBin==bin));
    end
    
    % remove instances of very few samples
    bins(nSamp<0.005) = [];
    n_g(nSamp<0.005) = [];
    
    % Step 4: Find linear and sigmoid fits to the transfer function (n_g vs bins) 
    % for each neuron and calculate the sum of squared errors for each fit.
    
    lm = fitlm(bins,n_g);
    sumSqLError = sum((lm.predict-n_g').^2);
    
    [~,stat] = sigm_fit(bins,n_g,[],[],0);
    sumSqNlError  = sum((stat.ypred-n_g').^2);
    
    % ratio of the two sums of squared errors 
    ratio(iNode) = sumSqNlError/sumSqLError;
    
end

