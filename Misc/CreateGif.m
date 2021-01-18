clear; clc; close all;

% simulate some random clustering data
ng = 15; % number of clusters/groups
nim = 250; % number of images

% ensure each group has at least one image in it
grps = 1:ng;
rgrps = unidrnd(ng,1,nim-ng); 
grps = [rgrps,grps];
grps = grps(randperm(nim));
hedges = double(grps == grps');

% we only care about the upper/lower triangle
hedges = triu(hedges,1);

% add a bit of noise
A = ones(nim);
idxsu = find(triu(A,1));
noise = rand(length(idxsu),1)-.5;
hedges(idxsu) = hedges(idxsu) + noise;

% sparsify the matrix a bit (to simulate missing data)
n_nans = 1000;
nanidxs = randperm(length(idxsu));
nanidxs = idxsu(nanidxs(1:n_nans));
hedges(nanidxs) = NaN;

% normalize between 0 and 1
minv = min(hedges(idxsu));
maxv = max(hedges(idxsu));
hedges(idxsu) = (hedges(idxsu) - minv)./(maxv-minv);

% make matrix symmetric
hedges = hedges+hedges';
hedges = hedges+diag(diag(A));

% now generate a clustering using the coordinate ascent clustering method

x = 1:ng;
gs = unidrnd(ng,1,nim-ng); 
gs = [x,gs];
gs = gs(randperm(nim));

h = figure;
axis tight manual;
axis square;
filename = 'CIRCAInAction.gif';

for i = 1:10
    [gs, ngcount, counter, R, converged] = CIRCA_Clustering_gif(hedges,ng,gs,1); 

    [~,sidx] = sort(gs);
    imagesc(hedges(sidx,sidx))
    % colorbar;
    axis off
    colormap copper
    set(gca,'LooseInset',get(gca,'TightInset'));
    
    drawnow
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    
    if i == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end
end

close;

