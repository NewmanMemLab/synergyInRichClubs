% This function calculates the amount of synergy in the rich club for all possible rich clubs. 
% Synergy occurs in the rich club when the receiver is a member.
%
% Inputs:
%   network: weighted, directed network (adjacency matrix)
%
%   synVals: Nx4 matrix of synergy data where each row is a triad. Column 1 
%   is the receiver node, column 2 is the first transmitter, column 3 is the second
%   transmitter, and column 4 is the synergy value.
%
% Outputs:
%   meanSyn: Nx2 vector of mean synergy per triad values for each of N
%   possible rich club levels, sorted from strongest to weakest. 
%   Column 1 is mean synergy inside the rich club.
%   Column 2 is mean synergy outside the rich club. 
%
%   propSyn: Nx2 vector of proportion of network synergy values for each of
%   N possible rich club levels, sorted from strongest to weakest. 
%   Column 1 is proportion of synergy inside the rich club. 
%   Column 2 is proportion of synergy outside the rich club.
%   
% Written by Sam Faber, 2016, Indiana University

function [meanSyn,propSyn] = synInRCAllThresh(network, synVals)

% get node weights, and sort them
weightedDeg = sum(network,1)'+sum(network,2);
[sortedwDeg, sortedInds] = sort(weightedDeg,'descend'); 

% initialize
meanSyn = nan(size(find(sortedwDeg),1),2);
propSyn = nan(size(find(sortedwDeg),1),2);

% for each possible threshold of the network (node weight is threshold)
for iNode = 1:size(find(sortedwDeg),1) 
    nodes = sortedInds(1:iNode);
    
    % mean synergy per triad in the rich club
    meanSyn(iNode,1) = nanmean(synVals(ismember(synVals(:,1),nodes),4));
    
    % mean synergy per triad outside the rich club
    meanSyn(iNode,2) = nanmean(synVals(~ismember(synVals(:,1),nodes),4));
    
    % proportion of network synergy in the rich club
    propSyn(iNode,1) = (sum(synVals(ismember(synVals(:,1),nodes),4))/...
        sum(synVals(:,4)));
    
    % proportion of network synergy outside the rich club
    propSyn(iNode,2) = 1-propSyn(iNode,1); 
        
end

