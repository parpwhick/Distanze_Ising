L = 150;
N = L^2;


[x,y]=ind2sub([L,L],1:N);
%draw AFM spin background
s = reshape(1/2 * (-1) .^ (x+y),L,L);
%draw a FM ellipse in the middle of the spin configuration
center = [L/2,L/2];
ellipse_axes = [1,1];
radius = L/3.5;
s(1/ellipse_axes(1)*(x-center(1)).^2+1/ellipse_axes(2)*(y-center(2)).^2 <= radius^2)=0.5;

%parameters of the Hamiltonian
t = 1;
lambda = 20;

%nearest neighbors to fill up the matrix
nnr = reshape(1:N,L,L);
nnd = circshift(nnr,-1);
nnu = circshift(nnr,1);
nnl = (circshift(nnr',1))';
nnr = (circshift(nnr',-1))';
colonne = repmat((1:N)',5,1);
righe = [nnl(:);nnr(:);nnd(:);nnu(:);(1:N)'];
%sum of the neighboring spins for every site
neigh_sum = sum(s([nnl(:),nnr(:),nnd(:),nnu(:)]),2);

% dopon spin up
values=([t*ones(4*N,1);lambda * (+1/2*s(:)+1/4) + 0.5*neigh_sum]);
H_up= sparse(righe,colonne,values);

% dopon spin down
values=([t*ones(4*N,1);lambda * (-1/2*s(:)+1/4) - 0.5*neigh_sum]);
H_down= sparse(righe,colonne,values);

opts.issym = true;