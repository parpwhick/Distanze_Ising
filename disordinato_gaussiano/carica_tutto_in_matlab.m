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

beta225=load('beta_2.25');
i=i+1;
dati(i).beta=2.25;
dati(i).x=beta225(:,1);
dati(i).y=beta225(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta338=load('beta_3.38');
i=i+1;
dati(i).beta=3.38;
dati(i).x=beta338(:,1);
dati(i).y=beta338(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta507=load('beta_5.07');
i=i+1;
dati(i).beta=5.07;
dati(i).x=beta507(:,1);
dati(i).y=beta507(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta761=load('beta_7.61');
i=i+1;
dati(i).beta=7.61;
dati(i).x=beta761(:,1);
dati(i).y=beta761(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1142=load('beta_11.42');
i=i+1;
dati(i).beta=11.42;
dati(i).x=beta1142(:,1);
dati(i).y=beta1142(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1713=load('beta_17.13');
i=i+1;
dati(i).beta=17.13;
dati(i).x=beta1713(:,1);
dati(i).y=beta1713(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta2570=load('beta_25.70');
i=i+1;
dati(i).beta=25.70;
dati(i).x=beta2570(:,1);
dati(i).y=beta2570(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta3855=load('beta_38.55');
i=i+1;
dati(i).beta=38.55;
dati(i).x=beta3855(:,1);
dati(i).y=beta3855(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta5782=load('beta_57.82');
i=i+1;
dati(i).beta=57.82;
dati(i).x=beta5782(:,1);
dati(i).y=beta5782(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta8673=load('beta_86.73');
i=i+1;
dati(i).beta=86.73;
dati(i).x=beta8673(:,1);
dati(i).y=beta8673(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta13009=load('beta_130.09');
i=i+1;
dati(i).beta=130.09;
dati(i).x=beta13009(:,1);
dati(i).y=beta13009(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta19513=load('beta_195.13');
i=i+1;
dati(i).beta=195.13;
dati(i).x=beta19513(:,1);
dati(i).y=beta19513(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta29269=load('beta_292.69');
i=i+1;
dati(i).beta=292.69;
dati(i).x=beta29269(:,1);
dati(i).y=beta29269(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta43903=load('beta_439.03');
i=i+1;
dati(i).beta=439.03;
dati(i).x=beta43903(:,1);
dati(i).y=beta43903(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta65854=load('beta_658.54');
i=i+1;
dati(i).beta=658.54;
dati(i).x=beta65854(:,1);
dati(i).y=beta65854(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta98781=load('beta_987.81');
i=i+1;
dati(i).beta=987.81;
dati(i).x=beta98781(:,1);
dati(i).y=beta98781(:,2);
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
