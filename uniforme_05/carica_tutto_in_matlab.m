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

beta100000=load('beta_1000.00');
i=i+1;
dati(i).beta=1000.00;
dati(i).x=beta100000(:,1);
dati(i).y=beta100000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1000000=load('beta_10000.00');
i=i+1;
dati(i).beta=10000.00;
dati(i).x=beta1000000(:,1);
dati(i).y=beta1000000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta10000000=load('beta_100000.00');
i=i+1;
dati(i).beta=100000.00;
dati(i).x=beta10000000(:,1);
dati(i).y=beta10000000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta12500=load('beta_125.00');
i=i+1;
dati(i).beta=125.00;
dati(i).x=beta12500(:,1);
dati(i).y=beta12500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1500=load('beta_15.00');
i=i+1;
dati(i).beta=15.00;
dati(i).x=beta1500(:,1);
dati(i).y=beta1500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta15000=load('beta_150.00');
i=i+1;
dati(i).beta=150.00;
dati(i).x=beta15000(:,1);
dati(i).y=beta15000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta17500=load('beta_175.00');
i=i+1;
dati(i).beta=175.00;
dati(i).x=beta17500(:,1);
dati(i).y=beta17500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta200=load('beta_2.00');
i=i+1;
dati(i).beta=2.00;
dati(i).x=beta200(:,1);
dati(i).y=beta200(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta2000=load('beta_20.00');
i=i+1;
dati(i).beta=20.00;
dati(i).x=beta2000(:,1);
dati(i).y=beta2000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta20000=load('beta_200.00');
i=i+1;
dati(i).beta=200.00;
dati(i).x=beta20000(:,1);
dati(i).y=beta20000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta200000=load('beta_2000.00');
i=i+1;
dati(i).beta=2000.00;
dati(i).x=beta200000(:,1);
dati(i).y=beta200000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta2000000=load('beta_20000.00');
i=i+1;
dati(i).beta=20000.00;
dati(i).x=beta2000000(:,1);
dati(i).y=beta2000000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta2500=load('beta_25.00');
i=i+1;
dati(i).beta=25.00;
dati(i).x=beta2500(:,1);
dati(i).y=beta2500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta300=load('beta_3.00');
i=i+1;
dati(i).beta=3.00;
dati(i).x=beta300(:,1);
dati(i).y=beta300(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta3000=load('beta_30.00');
i=i+1;
dati(i).beta=30.00;
dati(i).x=beta3000(:,1);
dati(i).y=beta3000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta30000=load('beta_300.00');
i=i+1;
dati(i).beta=300.00;
dati(i).x=beta30000(:,1);
dati(i).y=beta30000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta300000=load('beta_3000.00');
i=i+1;
dati(i).beta=3000.00;
dati(i).x=beta300000(:,1);
dati(i).y=beta300000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta3500=load('beta_35.00');
i=i+1;
dati(i).beta=35.00;
dati(i).x=beta3500(:,1);
dati(i).y=beta3500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta400=load('beta_4.00');
i=i+1;
dati(i).beta=4.00;
dati(i).x=beta400(:,1);
dati(i).y=beta400(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta4000=load('beta_40.00');
i=i+1;
dati(i).beta=40.00;
dati(i).x=beta4000(:,1);
dati(i).y=beta4000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta40000=load('beta_400.00');
i=i+1;
dati(i).beta=400.00;
dati(i).x=beta40000(:,1);
dati(i).y=beta40000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta400000=load('beta_4000.00');
i=i+1;
dati(i).beta=4000.00;
dati(i).x=beta400000(:,1);
dati(i).y=beta400000(:,2);
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

beta500000=load('beta_5000.00');
i=i+1;
dati(i).beta=5000.00;
dati(i).x=beta500000(:,1);
dati(i).y=beta500000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta600=load('beta_6.00');
i=i+1;
dati(i).beta=6.00;
dati(i).x=beta600(:,1);
dati(i).y=beta600(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta6000=load('beta_60.00');
i=i+1;
dati(i).beta=60.00;
dati(i).x=beta6000(:,1);
dati(i).y=beta6000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta60000=load('beta_600.00');
i=i+1;
dati(i).beta=600.00;
dati(i).x=beta60000(:,1);
dati(i).y=beta60000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta700=load('beta_7.00');
i=i+1;
dati(i).beta=7.00;
dati(i).x=beta700(:,1);
dati(i).y=beta700(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta7000=load('beta_70.00');
i=i+1;
dati(i).beta=70.00;
dati(i).x=beta7000(:,1);
dati(i).y=beta7000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta800=load('beta_8.00');
i=i+1;
dati(i).beta=8.00;
dati(i).x=beta800(:,1);
dati(i).y=beta800(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta8000=load('beta_80.00');
i=i+1;
dati(i).beta=80.00;
dati(i).x=beta8000(:,1);
dati(i).y=beta8000(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta900=load('beta_9.00');
i=i+1;
dati(i).beta=9.00;
dati(i).x=beta900(:,1);
dati(i).y=beta900(:,2);
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
