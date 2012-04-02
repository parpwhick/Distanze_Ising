function [catena,medie,E]=ising1d(N,beta,n,catena)
if (nargin <3)
    N=2000;
    beta=30;
    n=500;

end

exp_beta=exp(-beta*[-4 0 4]');

if (nargin < 4)
    catena=sign((rand(1,N)*2-1));
end

catena=reshape(catena,1,N);
medie=zeros(1,n);
E=zeros(1,n);

nnl=circshift(1:N,[1,1]);
nnr=circshift(1:N,[1,-1]);
steps=100*N;

for L=1:n   

    for i=1:steps
        j=ceil(rand()*N);
        dH = 2*catena(j).*(catena(nnl(j))+catena(nnr(j)));
        if(rand(1)<exp(-beta*dH))
             catena(j)=-catena(j);
        end
    end
    
    medie(L)=mean(catena);
    E(L)=sum(catena.*circshift(catena,[1,-1]));

end
