function [adiacenza,NN]=adiacenza_ising(lato)

N=lato*lato;

nnr = reshape(1:N,lato,lato);
nnd = circshift(nnr,-1);
nnu = circshift(nnr,1);
nnl = (circshift(nnr',1))';
nnr = (circshift(nnr',-1))';

NN=[nnu(:),nnd(:),nnl(:),nnr(:)];

[indici_riga,~,indici_colonna]=find(NN);
adiacenza=sparse(indici_riga,indici_colonna,1,N,N,N*4);

end