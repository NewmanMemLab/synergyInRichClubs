function   Rw = rich_club_wd_SamVer(CIJ)
%RICH_CLUB_WD	Rich club coefficients curve (weighted directed graph)
%
%   The weighted rich club coefficient, Rw, at level k is the fraction of
%   edge weights that connect nodes of weight k or higher out of the
%   maximum edge weights that such nodes might share.
%
%    Input:
%      CIJ:        weighted directed connection matrix
%
%   Output:
%       Rw:         rich-club curve
%
%
%   References:
%       T Opsahl et al. Phys Rev Lett, 2008, 101(16)
%       M van den Heuvel, O Sporns, J Neurosci 2011 31(44)
%
%   Martijn van den Heuvel, University Medical Center Utrecht, 2011

%   Modification History:
%   2011: Original
%   2015: Expanded documentation (Mika Rubinov)
%   2016: Adapted by Sam Faber, Indiana University, to use node weights
%   instead of degree.


% get node weights, rank them
weights = sum(CIJ,1)'+sum(CIJ,2);
weightsSorted = sort(weights);

% wrank contains the ranked weights of the network, with strongest connections on top
wrank = sort(CIJ(:), 'descend');


% loop over all possible rk-levels (richness parameter levels) 
for kk = 1:length(weightsSorted)

    SmallNodes=find(weights  < weightsSorted(kk)); 

    if isempty(SmallNodes);
        Rw(kk)=NaN;         %#ok<*AGROW>
        continue
    end
    
    % remove small nodes with weights < weightsSorted(kk)
    CutoutCIJ=CIJ;
    CutoutCIJ(SmallNodes,:)=[];
    CutoutCIJ(:,SmallNodes)=[];
    
    % total weight of connections in subset E>r
    Wr(kk) = sum(CutoutCIJ(:));

    % total number of connections in subset E>r
    Er = length(find(CutoutCIJ~=0));

    % E>r number of connections with max weight in network
    wrank_r = wrank(1:1:Er);
    
    % weighted rich-club coefficient
    Rw(kk)=Wr(kk) ./ sum(wrank_r);
    
end
 
    
    
    
    
        
