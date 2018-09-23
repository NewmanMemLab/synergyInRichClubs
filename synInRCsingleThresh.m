% This function calculates the amount of synergy in the rich club, for each of 6 possible interactions
% between triads and the rich club, for a single (significant) rich club. 
%
% Inputs:
%   network: weighted, directed network (adjacency matrix)
%
%   synVals: Nx4 matrix of synergy data where each row is a triad. Column 1 
%   is the receiver node, column 2 is the first transmitter, column 3 is the second
%   transmitter, and column 4 is the synergy value.
%   
%   RCcoeffs: vector of weighted, normalized rich club coefficients for the
%   corresponding 'network'
%
%   RCpValues: vector of p-values corresponding to rich club coefficients
%   in 'RCcoeffs'
%
%   pValThresh: The p-value threshold of your choice. Default is p=0.01.
%
% Outputs:
%   meanSyn: 6x1 vector of mean synergy per triad values for each configuration of
%   triads interacting with the rich club. 
%
%   propSyn: 6x1 vector of proportion of network synergy values for
%   each configuration of triads interacting with the rich club.
%   
%   Triad interactions with the rich club include:
%   1: No triad members in the rich club.
%   2: A single transmitter in the rich club.
%   3: Both transmitters in the rich club.
%   4: Just the receiver in the rich club.
%   5: The receiver and a single transmitter in the rich club.
%   6: All triad members in the rich club.
%
%   Cases 1-3 correspond to synergy occurring outside the rich club.
%   Cases 4-6 correspond to synergy occurring inside the rich club.
%
% Written by Sam Faber, 2016, Indiana University

function [meanSyn,propSyn] = synInRCsingleThresh(network,synVals,RCcoeffs,RCpValues,pValThresh)

if nargin < 5
    pValThresh = 0.01;
end

% initialize outputs
meanSyn = nan(6,1);
propSyn = nan(6,1);

% calculate node weight and sort nodes
weight = sum(network,1)' + sum(network,2);
[~,  sortedNodes]  = sort(weight,'ascend');

% Randomly select a rich club with which to perform the analysis & identify nodes in the rich club
rcClassification = zeros(length(network),1);
if size(find(RCcoeffs>1 & ...
        RCpValues < pValThresh & RCpValues >= 0),2) ~= 0
    cutOff = randsample(find(RCcoeffs>1 & ...
        RCpValues < pValThresh & RCpValues >= 0),1);
    rcClassification(cutOff:end) = 1;
else
    error('Could not find significant rich club.')
end
       
% get receiver and transmitter node ID's
rec = synVals(:,1);
transA = synVals(:,2);
transB = synVals(:,3);

% get inds to rc nodes in the network
RCMatrix = zeros(length(network),1);
RCNodes = sortedNodes(logical(rcClassification));
RCMatrix(RCNodes) = 1;

% create configuration matrix to easily identify which of the 6 cases (of
% triad interaction with the rich club) we are in
RCSynNodeMatrix = [RCMatrix(rec), RCMatrix(transA), RCMatrix(transB)];

% get synergy values
syn = synVals(:,4);


% CASE 1 - triad completely outside RC
if  size(find(sum(RCSynNodeMatrix,2) == 0),1)~=0  
    meanSyn(1) = sum(syn(sum(RCSynNodeMatrix,2) == 0))/size(find(sum(RCSynNodeMatrix,2) == 0),1);
    propSyn(1) = sum(syn(sum(RCSynNodeMatrix,2) == 0)) / sum(syn);
end

% CASE 2 - Single Transmitter in RC
if size(find(ismember(RCSynNodeMatrix,[[0 1 0];[0 0 1]],'rows')),1)~=0
    meanSyn(2) = sum(syn(ismember(RCSynNodeMatrix,[[0 1 0];[0 0 1]],'rows') ))/...
        size(find(ismember(RCSynNodeMatrix,[[0 1 0];[0 0 1]],'rows')),1);
    propSyn(2) = sum(syn(ismember(RCSynNodeMatrix,[[0 1 0];[0 0 1]],'rows') )) / sum(syn);
end

% CASE 3 - Both Transmitters in RC
if size(find(ismember(RCSynNodeMatrix,[0 1 1],'rows')),1)~=0
    meanSyn(3) = sum(syn(ismember(RCSynNodeMatrix,[0 1 1],'rows') ))/...
        size(find(ismember(RCSynNodeMatrix,[0 1 1],'rows')),1);
    propSyn(3) = sum(syn(ismember(RCSynNodeMatrix,[0 1 1],'rows') )) / sum(syn);
end

% CASE 4 - Only Receiver in RC
if size(find(ismember(RCSynNodeMatrix,[1 0 0],'rows')),1)~=0
    meanSyn(4) = sum(syn(ismember(RCSynNodeMatrix,[1 0 0],'rows') ))/...
        size(find(ismember(RCSynNodeMatrix,[1 0 0],'rows')),1);
    propSyn(4) = sum(syn(ismember(RCSynNodeMatrix,[1 0 0],'rows') )) / sum(syn);
end

% CASE 5 - Receiver and Single Transmitter in RC
if size(find(ismember(RCSynNodeMatrix,[[1 1 0];[1 0 1]],'rows')),1)~=0
    meanSyn(5) = sum(syn(ismember(RCSynNodeMatrix,[[1 1 0];[1 0 1]],'rows') ))/...
        size(find(ismember(RCSynNodeMatrix,[[1 1 0];[1 0 1]],'rows')),1);
    propSyn(5) = sum(syn(ismember(RCSynNodeMatrix,[[1 1 0];[1 0 1]],'rows') )) / sum(syn);
end

% CASE 6 - Receiver and both Transmitters in RC
if  size(find(sum(RCSynNodeMatrix,2) == 3),1)~=0
    meanSyn(6) = sum(syn(sum(RCSynNodeMatrix,2) == 3 ))/...
        size(find(sum(RCSynNodeMatrix,2) == 3),1);
    propSyn(6) = sum(syn(sum(RCSynNodeMatrix,2) == 3)) / sum(syn);
end
