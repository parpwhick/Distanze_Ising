lambda = 20;
t=1;
L= 10;
N = L^2;
H = zeros(N);

for i=1:(N)
H([nnu(i),nnl(i),nnr(i),nnd(i)],i)=t;
sum_nn = sum(s([nnu(i),nnl(i),nnr(i),nnd(i)]));
%sum_nn=0;
H(i,i) = lambda * (-1/2*s(i)+1/4) -1/2*sum_nn;
end
E_minus = eig(H);

H = zeros(N);
for i=1:(N)
H([nnu(i),nnl(i),nnr(i),nnd(i)],i)=t;
sum_nn = sum(s([nnu(i),nnl(i),nnr(i),nnd(i)]));
%sum_nn=0;
H(i,i) = lambda * (+1/2*s(i)+1/4) +1/2*sum_nn;
end
E_plus = eig(H);

plot(1:length(E_minus),E_minus,'.-b',1:length(E_plus),E_plus,'+-r');
legend('E down','E up')
title('Eigenvalues of the hopping term+constraint, over an AFM configuration')
ylabel('Eigenvalues');
xlabel('index')

