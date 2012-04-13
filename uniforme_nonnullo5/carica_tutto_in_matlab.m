xy = @(beta) loglog(beta(:,1),beta(:,2),'-');
i=0;

beta025=load('beta_0.25');
i=i+1;
dati(i).beta=0.25;
dati(i).x=beta025(:,1);
dati(i).y=beta025(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta050=load('beta_0.50');
i=i+1;
dati(i).beta=0.50;
dati(i).x=beta050(:,1);
dati(i).y=beta050(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta075=load('beta_0.75');
i=i+1;
dati(i).beta=0.75;
dati(i).x=beta075(:,1);
dati(i).y=beta075(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta100=load('beta_1.00');
i=i+1;
dati(i).beta=1.00;
dati(i).x=beta100(:,1);
dati(i).y=beta100(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta125=load('beta_1.25');
i=i+1;
dati(i).beta=1.25;
dati(i).x=beta125(:,1);
dati(i).y=beta125(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta150=load('beta_1.50');
i=i+1;
dati(i).beta=1.50;
dati(i).x=beta150(:,1);
dati(i).y=beta150(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta175=load('beta_1.75');
i=i+1;
dati(i).beta=1.75;
dati(i).x=beta175(:,1);
dati(i).y=beta175(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta200=load('beta_2.00');
i=i+1;
dati(i).beta=2.00;
dati(i).x=beta200(:,1);
dati(i).y=beta200(:,2);
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
