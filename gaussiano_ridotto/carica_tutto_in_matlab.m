xy = @(beta) loglog(beta(:,1),beta(:,2),'-');
i=0;

beta100=load('beta_1.00');
i=i+1;
dati(i).beta=1.00;
dati(i).x=beta100(:,1);
dati(i).y=beta100(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta130=load('beta_1.30');
i=i+1;
dati(i).beta=1.30;
dati(i).x=beta130(:,1);
dati(i).y=beta130(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta169=load('beta_1.69');
i=i+1;
dati(i).beta=1.69;
dati(i).x=beta169(:,1);
dati(i).y=beta169(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta220=load('beta_2.20');
i=i+1;
dati(i).beta=2.20;
dati(i).x=beta220(:,1);
dati(i).y=beta220(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta286=load('beta_2.86');
i=i+1;
dati(i).beta=2.86;
dati(i).x=beta286(:,1);
dati(i).y=beta286(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta372=load('beta_3.72');
i=i+1;
dati(i).beta=3.72;
dati(i).x=beta372(:,1);
dati(i).y=beta372(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta484=load('beta_4.84');
i=i+1;
dati(i).beta=4.84;
dati(i).x=beta484(:,1);
dati(i).y=beta484(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta629=load('beta_6.29');
i=i+1;
dati(i).beta=6.29;
dati(i).x=beta629(:,1);
dati(i).y=beta629(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta818=load('beta_8.18');
i=i+1;
dati(i).beta=8.18;
dati(i).x=beta818(:,1);
dati(i).y=beta818(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1063=load('beta_10.63');
i=i+1;
dati(i).beta=10.63;
dati(i).x=beta1063(:,1);
dati(i).y=beta1063(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1382=load('beta_13.82');
i=i+1;
dati(i).beta=13.82;
dati(i).x=beta1382(:,1);
dati(i).y=beta1382(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta1797=load('beta_17.97');
i=i+1;
dati(i).beta=17.97;
dati(i).x=beta1797(:,1);
dati(i).y=beta1797(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta2336=load('beta_23.36');
i=i+1;
dati(i).beta=23.36;
dati(i).x=beta2336(:,1);
dati(i).y=beta2336(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta3037=load('beta_30.37');
i=i+1;
dati(i).beta=30.37;
dati(i).x=beta3037(:,1);
dati(i).y=beta3037(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta3948=load('beta_39.48');
i=i+1;
dati(i).beta=39.48;
dati(i).x=beta3948(:,1);
dati(i).y=beta3948(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta5132=load('beta_51.32');
i=i+1;
dati(i).beta=51.32;
dati(i).x=beta5132(:,1);
dati(i).y=beta5132(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta6672=load('beta_66.72');
i=i+1;
dati(i).beta=66.72;
dati(i).x=beta6672(:,1);
dati(i).y=beta6672(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta8674=load('beta_86.74');
i=i+1;
dati(i).beta=86.74;
dati(i).x=beta8674(:,1);
dati(i).y=beta8674(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta11276=load('beta_112.76');
i=i+1;
dati(i).beta=112.76;
dati(i).x=beta11276(:,1);
dati(i).y=beta11276(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta14659=load('beta_146.59');
i=i+1;
dati(i).beta=146.59;
dati(i).x=beta14659(:,1);
dati(i).y=beta14659(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta19057=load('beta_190.57');
i=i+1;
dati(i).beta=190.57;
dati(i).x=beta19057(:,1);
dati(i).y=beta19057(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta24774=load('beta_247.74');
i=i+1;
dati(i).beta=247.74;
dati(i).x=beta24774(:,1);
dati(i).y=beta24774(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta32206=load('beta_322.06');
i=i+1;
dati(i).beta=322.06;
dati(i).x=beta32206(:,1);
dati(i).y=beta32206(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta41868=load('beta_418.68');
i=i+1;
dati(i).beta=418.68;
dati(i).x=beta41868(:,1);
dati(i).y=beta41868(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta54428=load('beta_544.28');
i=i+1;
dati(i).beta=544.28;
dati(i).x=beta54428(:,1);
dati(i).y=beta54428(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta70756=load('beta_707.56');
i=i+1;
dati(i).beta=707.56;
dati(i).x=beta70756(:,1);
dati(i).y=beta70756(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta91983=load('beta_919.83');
i=i+1;
dati(i).beta=919.83;
dati(i).x=beta91983(:,1);
dati(i).y=beta91983(:,2);
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
