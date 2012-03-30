xy = @(beta) loglog(beta(:,1),beta(:,2),'-');
i=0;

beta100=load('beta_1.00');
i=i+1;
dati(i).beta=1.00;
dati(i).x=beta100(:,1);
dati(i).y=beta100(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta150=load('beta_1.50');
i=i+1;
dati(i).beta=1.50;
dati(i).x=beta150(:,1);
dati(i).y=beta150(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta200=load('beta_2.00');
i=i+1;
dati(i).beta=2.00;
dati(i).x=beta200(:,1);
dati(i).y=beta200(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta250=load('beta_2.50');
i=i+1;
dati(i).beta=2.50;
dati(i).x=beta250(:,1);
dati(i).y=beta250(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta300=load('beta_3.00');
i=i+1;
dati(i).beta=3.00;
dati(i).x=beta300(:,1);
dati(i).y=beta300(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta350=load('beta_3.50');
i=i+1;
dati(i).beta=3.50;
dati(i).x=beta350(:,1);
dati(i).y=beta350(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta400=load('beta_4.00');
i=i+1;
dati(i).beta=4.00;
dati(i).x=beta400(:,1);
dati(i).y=beta400(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta500=load('beta_5.00');
i=i+1;
dati(i).beta=5.00;
dati(i).x=beta500(:,1);
dati(i).y=beta500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta700=load('beta_7.00');
i=i+1;
dati(i).beta=7.00;
dati(i).x=beta700(:,1);
dati(i).y=beta700(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta800=load('beta_8.00');
i=i+1;
dati(i).beta=8.00;
dati(i).x=beta800(:,1);
dati(i).y=beta800(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta900=load('beta_9.00');
i=i+1;
dati(i).beta=9.00;
dati(i).x=beta900(:,1);
dati(i).y=beta900(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1200=load('beta_12.00');
i=i+1;
dati(i).beta=12.00;
dati(i).x=beta1200(:,1);
dati(i).y=beta1200(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1300=load('beta_13.00');
i=i+1;
dati(i).beta=13.00;
dati(i).x=beta1300(:,1);
dati(i).y=beta1300(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1900=load('beta_19.00');
i=i+1;
dati(i).beta=19.00;
dati(i).x=beta1900(:,1);
dati(i).y=beta1900(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta2800=load('beta_28.00');
i=i+1;
dati(i).beta=28.00;
dati(i).x=beta2800(:,1);
dati(i).y=beta2800(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta4200=load('beta_42.00');
i=i+1;
dati(i).beta=42.00;
dati(i).x=beta4200(:,1);
dati(i).y=beta4200(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta6300=load('beta_63.00');
i=i+1;
dati(i).beta=63.00;
dati(i).x=beta6300(:,1);
dati(i).y=beta6300(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta9400=load('beta_94.00');
i=i+1;
dati(i).beta=94.00;
dati(i).x=beta9400(:,1);
dati(i).y=beta9400(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta14100=load('beta_141.00');
i=i+1;
dati(i).beta=141.00;
dati(i).x=beta14100(:,1);
dati(i).y=beta14100(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta21100=load('beta_211.00');
i=i+1;
dati(i).beta=211.00;
dati(i).x=beta21100(:,1);
dati(i).y=beta21100(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta31600=load('beta_316.00');
i=i+1;
dati(i).beta=316.00;
dati(i).x=beta31600(:,1);
dati(i).y=beta31600(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta47400=load('beta_474.00');
i=i+1;
dati(i).beta=474.00;
dati(i).x=beta47400(:,1);
dati(i).y=beta47400(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta71100=load('beta_711.00');
i=i+1;
dati(i).beta=711.00;
dati(i).x=beta71100(:,1);
dati(i).y=beta71100(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta100000=load('beta_1000.00');
i=i+1;
dati(i).beta=1000.00;
dati(i).x=beta100000(:,1);
dati(i).y=beta100000(:,2);
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
