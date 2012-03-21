xy = @(beta) loglog(beta(:,1),beta(:,2),'-');
i=0;

beta02=load('beta_0.2');
i=i+1;
dati(i).beta=0.2;
dati(i).x=beta02(:,1);
dati(i).y=beta02(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta06=load('beta_0.6');
i=i+1;
dati(i).beta=0.6;
dati(i).x=beta06(:,1);
dati(i).y=beta06(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta08=load('beta_0.8');
i=i+1;
dati(i).beta=0.8;
dati(i).x=beta08(:,1);
dati(i).y=beta08(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta10=load('beta_1.0');
i=i+1;
dati(i).beta=1.0;
dati(i).x=beta10(:,1);
dati(i).y=beta10(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta12=load('beta_1.2');
i=i+1;
dati(i).beta=1.2;
dati(i).x=beta12(:,1);
dati(i).y=beta12(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta13=load('beta_1.3');
i=i+1;
dati(i).beta=1.3;
dati(i).x=beta13(:,1);
dati(i).y=beta13(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta15=load('beta_1.5');
i=i+1;
dati(i).beta=1.5;
dati(i).x=beta15(:,1);
dati(i).y=beta15(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta16=load('beta_1.6');
i=i+1;
dati(i).beta=1.6;
dati(i).x=beta16(:,1);
dati(i).y=beta16(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta18=load('beta_1.8');
i=i+1;
dati(i).beta=1.8;
dati(i).x=beta18(:,1);
dati(i).y=beta18(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta19=load('beta_1.9');
i=i+1;
dati(i).beta=1.9;
dati(i).x=beta19(:,1);
dati(i).y=beta19(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta20=load('beta_2.0');
i=i+1;
dati(i).beta=2.0;
dati(i).x=beta20(:,1);
dati(i).y=beta20(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta21=load('beta_2.1');
i=i+1;
dati(i).beta=2.1;
dati(i).x=beta21(:,1);
dati(i).y=beta21(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta23=load('beta_2.3');
i=i+1;
dati(i).beta=2.3;
dati(i).x=beta23(:,1);
dati(i).y=beta23(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta24=load('beta_2.4');
i=i+1;
dati(i).beta=2.4;
dati(i).x=beta24(:,1);
dati(i).y=beta24(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta25=load('beta_2.5');
i=i+1;
dati(i).beta=2.5;
dati(i).x=beta25(:,1);
dati(i).y=beta25(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta27=load('beta_2.7');
i=i+1;
dati(i).beta=2.7;
dati(i).x=beta27(:,1);
dati(i).y=beta27(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta28=load('beta_2.8');
i=i+1;
dati(i).beta=2.8;
dati(i).x=beta28(:,1);
dati(i).y=beta28(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta30=load('beta_3.0');
i=i+1;
dati(i).beta=3.0;
dati(i).x=beta30(:,1);
dati(i).y=beta30(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta31=load('beta_3.1');
i=i+1;
dati(i).beta=3.1;
dati(i).x=beta31(:,1);
dati(i).y=beta31(:,2);
dati(i).target=dati(i).y(end)/2;
dati(i).n=length(dati(i).x);

beta32=load('beta_3.2');
i=i+1;
dati(i).beta=3.2;
dati(i).x=beta32(:,1);
dati(i).y=beta32(:,2);
dati(i).target=dati(i).y(end)/2;
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
legend('log tanh \beta','fattore di scaling sperimentale')
hold off
xlabel('\beta');

figure(3)
semilogx(dati(1).x/dati(1).xtarget,dati(1).y);
hold all
for j=2:i
    semilogx(dati(j).x/dati(j).xtarget,dati(j).y);
end
