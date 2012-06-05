function [L] = MetroByNN(NN,beta,T,L0)
%function rho = MetroIsing2d(N,bt,T,L0)


% if nargin < 4
%     cicle=30;
% end

[N,zmax]=size(NN);

Z=sum(NN>0,2);

if (nargin == 4 && length(L0) == N)
    L=L0;
else
    N12=ceil(sqrt(N));
    N1=N12^2;
    L = 2*round(rand(N1,1))-1;
    L((N+1):N1) = 0;
    L=reshape(L,N12,N12);
end;

NN(NN(:)==0)=N+1;

% elenco=reshape(ceil((1:N^2)./N),N,N);
% elenco=elenco+elenco';
% dispari=find(mod(elenco,2));
% pari=find(~mod(elenco,2));

for t = 1:T,
    %    beta=bt+0.4*sin(2*pi/200*t);
    %exp_beta=exp(-beta*[-8 -4 0 4 8]');
    
    N2 = randperm(N);
    quanti=ceil(N/200);
    for k = 0:(floor(N/quanti)-1);
        
        i=N2((k*quanti+1):((k+1)*quanti));
        %j=NN(i,zetas(k));
        %nn=NN(i,NN(i,:)>0);
        %dH = 2*L(i)*sum(L(NN(i,1:Z(i))));
        dH = 2*L(i)'.*sum(L(NN(i,:)),2);
        chosen= (dH < 0) | (rand(quanti,1) < exp(-beta*dH));
        L(i(chosen))=-L(i(chosen));
    end;  
    
    %imagesc(L);
    %pause(0.01);
    
end;

end

