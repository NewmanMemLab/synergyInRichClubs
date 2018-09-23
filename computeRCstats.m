% This function computes weighted, normalized rich club coefficients, as
% well as there p-values and z-scores (computed based on 500 randomly
% rewired networks).

% Input:
%   CIJ: A weighted, directed adjacency matrix
% Outputs: 
%   RCcoeffs: Vector of weighted, normalized rich club coefficients
%   RCpValues: Vector of p-values for each coefficient
%   RCzSCore: Vector of z-scores for each coefficient

% Written by Sam Faber
% 2016, Indiana University

function [RCcoeffs,RCpValues,RCzScore] = computeRCstats(CIJ)

% pre-allocate
randRCcoeffs = nan(500,size(CIJ,1));

% generate 500 null models 
for ii = 1:500
    % randomly rewire the network 
    CIJrand = dir_generate_srand_SamVer(CIJ); 
    % calculate the rich club coefficients for shuffled networks
    randRCcoeffs(ii,:) = rich_club_wd_SamVer(CIJrand); % each row is a random RCcoeff vector
end

% get RC coeffs for observed network
RCcoeffs = rich_club_wd_SamVer(CIJ);

% normalize RC coeffs
RCcoeffs = RCcoeffs./nanmean(randRCcoeffs);

% get p-values 
RCpValues = sum(randRCcoeffs-repmat(RCcoeffs,500,1) > eps)./500;

% get z-scores
RCzScore = (RCcoeffs - nanmean(randRCcoeffs))./nanstd(randRCcoeffs);
