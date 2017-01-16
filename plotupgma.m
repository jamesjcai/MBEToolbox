function [topology] = plotupgma(dist,labels,fontsize)
%PLOTUPGMA - Performs UPGMA on distance matrix and produces plot of dendrogram
%
% Syntax: [topology] = plotupgma(dist,labels,fontsize)
%
% Inputs:
%    dist       - [n x n] symmetric distance matrix.
%    labels     - optional [n x q] matrix of group labels for dendrogram.
%    fontsize   - optional font size for labels [default = 10].
%
% Outputs:
%    topology   - [(n-1) x 4] matrix summarizing dendrogram topology:
%                 col 1 = 1st OTU/cluster being grouped at current step
%                 col 2 = 2nd OTU/cluster
%                 col 3 = ID of cluster being produced
%                 col 4 = distance at node
%
% See also: plotnjtree

% RE Strauss, 5/27/96
%   9/7/99 - miscellaneous changes for Matlab v5.
%   9/24/01 - check diagonal elements against eps rather than zero.

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



  if nargin<2, labels=[]; end; % #ok<SEPEX>
  if nargin<3, fontsize=[]; end;

  [n,p] = size(dist);
  if (n~=p || any(diag(dist)>eps))
    dist
    error('  UPGMA: input matrix is not a distance matrix');
  end;

  if (~isempty(labels))
    if (size(labels,1)~=n)
      error('  UPGMA: numbers of taxa and taxon labels do not match');
    end;
  end;


  suprt = 0;
  if (nargout > 1)
    suprt = 1;
  end;

  clstsize = ones(1,n);               % Number of elements in clusters/otus
  id = 1:n;                           % Cluster IDs
  topology = zeros(n-1,4);            % Output dendrogram-topology matrix

  plug = 10e6;
  dist = dist + eye(n)*plug;          % Replace diagonal with plugs

  for step = 1:(n-1)                  % Clustering steps
    min_dist = min(dist(:));            % Find minimum pairwise distance
    [ii,jj] = find(dist==min_dist);     % Find location of minimum
    k = 1;                              % Use first identified minimum
    while (ii(k)>jj(k))                 %   for which i<j
      k = k+1;
    end;
    i = ii(k);
    j = jj(k);
    if (id(i)<id(j))
      topology(step,:) = [id(i) id(j) n+step min_dist];
    else
      topology(step,:) = [id(j) id(i) n+step min_dist];
    end;
    id(i) = n+step;
    dist(i,j) = plug;
    dist(j,i) = plug;

    new_clstsize = clstsize(i) + clstsize(j);
    alpha_i = clstsize(i) / new_clstsize;
    alpha_j = clstsize(j) / new_clstsize;
    clstsize(i) = new_clstsize;

    for k = 1:n                         % For all other clusters/OTUs,
      if (k~=i && k~=j)                %   adjust distances to new cluster
        dist(k,i) = alpha_i * dist(k,i) + alpha_j * dist(k,j);
        dist(i,k) = alpha_i * dist(i,k) + alpha_j * dist(j,k);
        dist(k,j) = plug;
        dist(j,k) = plug;
      end;
    end;
  end; % for step


    plot_dendrogram(topology,labels,fontsize);

  return;