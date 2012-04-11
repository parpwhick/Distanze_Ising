xy = @(beta) loglog(beta(:,1),beta(:,2),'-');
i=0;

beta050=load('beta_0.50');
i=i+1;
dati(i).beta=0.50;
dati(i).x=beta050(:,1);
dati(i).y=beta050(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

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

beta1142=load('beta_11.42');
i=i+1;
dati(i).beta=11.42;
dati(i).x=beta1142(:,1);
dati(i).y=beta1142(:,2);
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

beta225=load('beta_2.25');
i=i+1;
dati(i).beta=2.25;
dati(i).x=beta225(:,1);
dati(i).y=beta225(:,2);
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

beta338=load('beta_3.38');
i=i+1;
dati(i).beta=3.38;
dati(i).x=beta338(:,1);
dati(i).y=beta338(:,2);
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

beta450=load('beta_4.50');
i=i+1;
dati(i).beta=4.50;
dati(i).x=beta450(:,1);
dati(i).y=beta450(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta500=load('beta_5.00');
i=i+1;
dati(i).beta=5.00;
dati(i).x=beta500(:,1);
dati(i).y=beta500(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta507=load('beta_5.07');
i=i+1;
dati(i).beta=5.07;
dati(i).x=beta507(:,1);
dati(i).y=beta507(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta600=load('beta_6.00');
i=i+1;
dati(i).beta=6.00;
dati(i).x=beta600(:,1);
dati(i).y=beta600(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta700=load('beta_7.00');
i=i+1;
dati(i).beta=7.00;
dati(i).x=beta700(:,1);
dati(i).y=beta700(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta761=load('beta_7.61');
i=i+1;
dati(i).beta=7.61;
dati(i).x=beta761(:,1);
dati(i).y=beta761(:,2);
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
