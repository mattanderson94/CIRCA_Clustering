function [gs, ngcount, counter, R, converged] = CIRCA_Clustering_gif(hedges,ng,gs,maxits) %#codegen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate Ascent clustering method:
% By Matt Anderson & James H Elder. For bug reports, please contact
% matt.anderson@soton.ac.uk

% INPUT ARGUMENTS:
% hedges = n x n symmetric similarity/adjacency matrix of pairwise similarity
% judgements. Missing values should be imputed as NaN. Values should be
% normalized to be within 0 and 1.
% ng = integer specifying number of clusters.
% maxits = number of coordinate ascent iterations. set to inf to ensure convergence.

% OUTPUT ARGUMENTS:
% gs = An nx1 vector of integers from 1-ng, which gives the cluster 
% assignments for every row/column (i.e., image/item) in hedges.
% ngcount = number of moves made via coordinate ascent.
% counter = number of proposals made. 
% R = Rand Index between similarity matrix and optimized clustering.
% converged = true if the algorithm reaches a stationary point (and no
% moves that improve the rand index are possible)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial checks:
if ~isequaln(hedges,hedges.')
    error('Similarity matrix must be symmetric')
end

if min(hedges,[],'all') < 0 || max(hedges,[],'all') > 1
    error('All similarity values must be in the 0->1 range')
end

if ~(rem(ng,1) == 0)
    error('Number of clusters must be a real, positive integer')
end

if ng >= size(hedges,1)
    error('Number of clusters must be smaller than number of images in similarity matrix')
end

nim = size(hedges,1); 
neg_hedges = 1-hedges; 
num_p = sum(~isnan(hedges))-1; % number of pairings per image

% randomly allocate images, ensuring each category appears at least once
bestgs = gs;

% initial similarity matrix for model
m0edges = bsxfun(@eq,gs,gs');

% initial match between human and model 
matches = m0edges.*hedges + (1-m0edges).*(neg_hedges);
matches = matches - diag(diag(matches));

% Starting Rand index for each image
Rims = nansum(matches)./(num_p);

% initialize loop vars
searching = true;
converged = false;
ngcount = 0; % number of deviations from coordinate ascent
numsearches = 0;
counter = 1; % proposal number

getRandIndex = @(x,y)nanmean((x == x').*y + (x ~= x').*(1-y),'all');

while searching

    numsearches = numsearches+1;
    
    % Enumerate all exhaustive combinations of images and categories & randomize order
    imi = 1:nim;
    imi = imi(randperm(nim));
    [~,gi] = sort(rand(nim,ng),2);
    
    for i = 1:nim
        
        % Propose to move image im to a different category. Try all
        % categories, and accept the first proposal that improves the rand index
        im = imi(i); 
        
        for j = 1:ng
            
            g = gi(i,j);
            
            % ignore non-moves
            if gs(im) == g
                continue
            end 
            
            % keep track of number of proposals
            counter = counter + 1;
            
            % compute rand index for image im, if moved to category g
            match = (gs==g).*hedges(im,:) + (gs~=g).*(neg_hedges(im,:));
            Rim = nansum(match)/(num_p(im));
            dR = (Rim - Rims(im));
                        
            % make sure we stick with the initialized number of clusters
            if dR > 0 && sum(gs == gs(im)) > 1
                Rims(im) = Rims(im) + dR;
                gs(im) = g;
                ngcount = ngcount + 1;
                break
            end
            
        end
        
    end
    
    % if there is no change in the rand index after attempting all
    % combinations, then we have reached a stationary point
    if isequal(bestgs, gs)
        R = getRandIndex(gs,hedges);
        converged = true;
        searching = false;
    end
    
    bestgs = gs;
    
    % exit early if max iteration is met
    if numsearches >= maxits
        R = getRandIndex(gs,hedges);
        searching = false;
    end
    
end
