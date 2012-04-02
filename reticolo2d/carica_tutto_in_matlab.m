xy = @(beta) loglog(beta(:,1),beta(:,2),'-');
i=0;

beta010=load('beta_0.10');
i=i+1;
dati(i).beta=0.10;
dati(i).x=beta010(:,1);
dati(i).y=beta010(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta015=load('beta_0.15');
i=i+1;
dati(i).beta=0.15;
dati(i).x=beta015(:,1);
dati(i).y=beta015(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta020=load('beta_0.20');
i=i+1;
dati(i).beta=0.20;
dati(i).x=beta020(:,1);
dati(i).y=beta020(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta025=load('beta_0.25');
i=i+1;
dati(i).beta=0.25;
dati(i).x=beta025(:,1);
dati(i).y=beta025(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta030=load('beta_0.30');
i=i+1;
dati(i).beta=0.30;
dati(i).x=beta030(:,1);
dati(i).y=beta030(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta035=load('beta_0.35');
i=i+1;
dati(i).beta=0.35;
dati(i).x=beta035(:,1);
dati(i).y=beta035(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta039=load('beta_0.39');
i=i+1;
dati(i).beta=0.39;
dati(i).x=beta039(:,1);
dati(i).y=beta039(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta040=load('beta_0.40');
i=i+1;
dati(i).beta=0.40;
dati(i).x=beta040(:,1);
dati(i).y=beta040(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta041=load('beta_0.41');
i=i+1;
dati(i).beta=0.41;
dati(i).x=beta041(:,1);
dati(i).y=beta041(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta042=load('beta_0.42');
i=i+1;
dati(i).beta=0.42;
dati(i).x=beta042(:,1);
dati(i).y=beta042(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta043=load('beta_0.43');
i=i+1;
dati(i).beta=0.43;
dati(i).x=beta043(:,1);
dati(i).y=beta043(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta044=load('beta_0.44');
i=i+1;
dati(i).beta=0.44;
dati(i).x=beta044(:,1);
dati(i).y=beta044(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta045=load('beta_0.45');
i=i+1;
dati(i).beta=0.45;
dati(i).x=beta045(:,1);
dati(i).y=beta045(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta046=load('beta_0.46');
i=i+1;
dati(i).beta=0.46;
dati(i).x=beta046(:,1);
dati(i).y=beta046(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta050=load('beta_0.50');
i=i+1;
dati(i).beta=0.50;
dati(i).x=beta050(:,1);
dati(i).y=beta050(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta055=load('beta_0.55');
i=i+1;
dati(i).beta=0.55;
dati(i).x=beta055(:,1);
dati(i).y=beta055(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta060=load('beta_0.60');
i=i+1;
dati(i).beta=0.60;
dati(i).x=beta060(:,1);
dati(i).y=beta060(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta065=load('beta_0.65');
i=i+1;
dati(i).beta=0.65;
dati(i).x=beta065(:,1);
dati(i).y=beta065(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta070=load('beta_0.70');
i=i+1;
dati(i).beta=0.70;
dati(i).x=beta070(:,1);
dati(i).y=beta070(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta075=load('beta_0.75');
i=i+1;
dati(i).beta=0.75;
dati(i).x=beta075(:,1);
dati(i).y=beta075(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta080=load('beta_0.80');
i=i+1;
dati(i).beta=0.80;
dati(i).x=beta080(:,1);
dati(i).y=beta080(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta085=load('beta_0.85');
i=i+1;
dati(i).beta=0.85;
dati(i).x=beta085(:,1);
dati(i).y=beta085(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta090=load('beta_0.90');
i=i+1;
dati(i).beta=0.90;
dati(i).x=beta090(:,1);
dati(i).y=beta090(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta095=load('beta_0.95');
i=i+1;
dati(i).beta=0.95;
dati(i).x=beta095(:,1);
dati(i).y=beta095(:,2);
dati(i).target=0.65;
dati(i).n=length(dati(i).x);

beta100=load('beta_1.00');
i=i+1;
dati(i).beta=1.00;
dati(i).x=beta100(:,1);
dati(i).y=beta100(:,2);
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
