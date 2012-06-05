function [L] = MetroByAdj(adj,beta,T,L0)
%function rho = MetroIsing2d(N,bt,T,L0)


% if nargin < 4
%     cicle=30;
% end

[N,~]=size(adj);


if (nargin == 4 && length(L0) == N)
    L=L0;
else
%     N12=ceil(sqrt(N));
%     N1=N12^2;
%     L = 2*round(rand(N1,1))-1;
%     L((N+1):N1) = 0;
%     L=reshape(L,N12,N12);
     L = 2*round(rand(1,N))-1;
end;

% elenco=reshape(ceil((1:N^2)./N),N,N);
% elenco=elenco+elenco';
% dispari=find(mod(elenco,2));
% pari=find(~mod(elenco,2));

for t = 1:T    
    N2 = randperm(N);
    quanti=ceil(N/200);
    for k = 0:(floor(N/quanti)-1);
        
        i=N2((k*quanti+1):((k+1)*quanti));
        
        %versione scalare:
        %espressione ottimale per sfruttare la sparsita'
        %i=N2(k);
        %dH = 2*L(i)* (L * adj(:,i));
        %if (dH < 0) || (rand < exp(-beta*dH))
        %    L(i) = - L(i);
        %end;
        
        %vettoriale, 'quanti' elementi per volta
        dH = 2*L(i) .* ( L * adj(:,i));
        chosen= (dH < 0) | (rand(1,quanti) < exp(-beta*dH));
        L(i(chosen))=-L(i(chosen));
    end;  
    
    %imagesc(L);
    %pause(0.01);
    
end;

end

