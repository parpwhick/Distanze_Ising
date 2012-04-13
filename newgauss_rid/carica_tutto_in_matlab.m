xy = @(beta) loglog(beta(:,1),beta(:,2),'-');
i=0;

beta100=load('beta_1.00');
i=i+1;
dati(i).beta=1.00;
dati(i).x=beta100(:,1);
dati(i).y=beta100(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1000=load('beta_10.00');
i=i+1;
dati(i).beta=10.00;
dati(i).x=beta1000(:,1);
dati(i).y=beta1000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta10000=load('beta_100.00');
i=i+1;
dati(i).beta=100.00;
dati(i).x=beta10000(:,1);
dati(i).y=beta10000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1500=load('beta_15.00');
i=i+1;
dati(i).beta=15.00;
dati(i).x=beta1500(:,1);
dati(i).y=beta1500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta2000=load('beta_20.00');
i=i+1;
dati(i).beta=20.00;
dati(i).x=beta2000(:,1);
dati(i).y=beta2000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta2500=load('beta_25.00');
i=i+1;
dati(i).beta=25.00;
dati(i).x=beta2500(:,1);
dati(i).y=beta2500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta3000=load('beta_30.00');
i=i+1;
dati(i).beta=30.00;
dati(i).x=beta3000(:,1);
dati(i).y=beta3000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta3500=load('beta_35.00');
i=i+1;
dati(i).beta=35.00;
dati(i).x=beta3500(:,1);
dati(i).y=beta3500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta4000=load('beta_40.00');
i=i+1;
dati(i).beta=40.00;
dati(i).x=beta4000(:,1);
dati(i).y=beta4000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta4500=load('beta_45.00');
i=i+1;
dati(i).beta=45.00;
dati(i).x=beta4500(:,1);
dati(i).y=beta4500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta500=load('beta_5.00');
i=i+1;
dati(i).beta=5.00;
dati(i).x=beta500(:,1);
dati(i).y=beta500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta5000=load('beta_50.00');
i=i+1;
dati(i).beta=50.00;
dati(i).x=beta5000(:,1);
dati(i).y=beta5000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta50000=load('beta_500.00');
i=i+1;
dati(i).beta=500.00;
dati(i).x=beta50000(:,1);
dati(i).y=beta50000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta5500=load('beta_55.00');
i=i+1;
dati(i).beta=55.00;
dati(i).x=beta5500(:,1);
dati(i).y=beta5500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta6000=load('beta_60.00');
i=i+1;
dati(i).beta=60.00;
dati(i).x=beta6000(:,1);
dati(i).y=beta6000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta7000=load('beta_70.00');
i=i+1;
dati(i).beta=70.00;
dati(i).x=beta7000(:,1);
dati(i).y=beta7000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta8000=load('beta_80.00');
i=i+1;
dati(i).beta=80.00;
dati(i).x=beta8000(:,1);
dati(i).y=beta8000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta9000=load('beta_90.00');
i=i+1;
dati(i).beta=90.00;
dati(i).x=beta9000(:,1);
dati(i).y=beta9000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

clear target;

figure(1)
semilogx(dati(1).x,dati(1).y);
hold all
for j=2:i
    semilogx(dati(j).x,dati(j).y);
end

legende=cell(1,i);
for j=1:i
    legende(j)={['\beta=',num2str(dati(j).beta)]};
end
%legend(legende);

figure(2);
%calcoliamo tutti i punti centrali!!!!
risultato=1:i;
beta=risultato;
for j=1:i
dati(j).xtarget=fitta_andamento(dati(j));
	disp(dati(j).xtarget);
	risultato(j)=dati(j).xtarget;
	beta(j)=dati(j).beta;
end

[beta,p]=sort(beta);
risultato=risultato(p);

figure(2);
indice = round(i/2);
risultato_riscalato = risultato/(risultato(indice)*-log(tanh(beta(indice))));
semilogy(beta,-1./log(tanh(beta)),beta,risultato_riscalato,'o');
hold on
plot(beta(indice),-1/log(tanh(beta(indice))),'s','MarkerSize',7,'MarkerFaceColor','red')
legend('log tanh \beta','lunghezza di riscalamento')
hold off
xlabel('\beta');

figure(3)
%semilogx(dati(1).x/dati(1).xtarget,dati(1).y);
%hold all
%for j=2:i
%    semilogx(dati(j).x/dati(j).xtarget,dati(j).y);
%end

loglog(beta,risultato,'-o');
legend('lunghezza di riscalamento')
xlabel('\beta');
ylabel('lunghezza di riscalamento');
