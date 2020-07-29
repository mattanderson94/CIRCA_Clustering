function [R, ngcount, counter, gs, converged] = CoordinateAscent_RI2(hedges,ng,maxits) %#codegen

% Coordinate Ascent clustering method:

% hedges = nim x nim symmetric similarity matrix of probabilites that pairs
% of images belong to the same category/cluster. Missing values should be imputed as 'NaN's
% ng = number of clusters
% maxits = number of coordinate ascent iterations. set to inf to ensure convergence

nim = size(hedges,1); 
neg_hedges = 1-hedges; 
num_p = sum(~isnan(hedges))-1; % number of pairings per image

% randomly allocate images, ensuring each category appears at least once
x = 1:ng;
gs = unidrnd(ng,1,nim-ng); 
gs = [x,gs];
gs = gs(randperm(nim));
bestgs = gs;

% initial similarity matrix for model
m0edges = bsxfun(@eq,gs,gs');

% initial match between human and model 
matches = m0edges.*hedges + (1-m0edges).*(neg_hedges);
matches = matches - diag(diag(matches));

% Starting Rand index for each image
Rims = nansum(matches)./(num_p);

% Initial rand index for all images
R = nansum(matches(:))/nansum(num_p); 

% initialize loop vars
searching = true;
converged = false;
ngcount = 0; % number of deviations from coordinate ascent
numsearches = 0;
counter = 1; % proposal number

while searching

    numsearches = numsearches+1;
    % disp(['Starting Round ', num2str(numsearches), ' of coordinate ascent']);
    
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
            Rim = nansum(match)/num_p(im);
            dR = (Rim - Rims(im));
            Rprop = R + dR*(num_p(im)/sum(num_p));
                        
            % make sure we stick with the initial number of clusters
            if dR > 0 && sum(gs == gs(im)) > 1
                % update clustering
                ngcount = ngcount + 1; 
                R = Rprop;
                gs(im) = g;
                Rims(im) = Rim;
                break
            end
            
        end
        
    end
    
    % if there is no change in the rand index after attempting all
    % combinations, then we have reached a stationary point
    if isequal(bestgs, gs)
        converged = true;
        searching = false;
    end
    
    bestgs = gs;
    
    % exit early if max iteration is met
    if numsearches >= maxits
        break;
    end
    
end

% codegen CoordinateAscent_RI2 -args {coder.typeof(double(1), [1000 1000]), double(1), double(1)}