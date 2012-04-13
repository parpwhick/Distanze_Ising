$header=<<EOF;
xy = @(beta) loglog(beta(:,1),beta(:,2),'-');
EOF

print $header;

$i=0;
print "i=0;\n\n";


while(<>){
	/beta_(\d+).(\d+)/; 
	print "beta$1$2=load('beta_$1.$2');\n";	
	print "i=i+1;\n";
	print "dati(i).beta=$1.$2;\n";
	print "dati(i).x=beta$1$2(:,1);\n";
	print "dati(i).y=beta$1$2(:,2);\n";
	print "dati(i).target=0.65;\n";
	#print "dati(i).target=dati(i).y(end)/2;\n";
	print "dati(i).n=length(dati(i).x);\n";
	print "\n";
}

$tail=<<EOF;
clear target;

figure(1)
semilogx(dati(1).x,dati(1).y);
hold all
for j=2:i
    semilogx(dati(j).x,dati(j).y);
end

legende=cell(1,i);
for j=1:i
    legende(j)={['\\beta=',num2str(dati(j).beta)]};
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
legend('log tanh \\beta','lunghezza di riscalamento')
hold off
xlabel('\\beta');

figure(3)
%semilogx(dati(1).x/dati(1).xtarget,dati(1).y);
%hold all
%for j=2:i
%    semilogx(dati(j).x/dati(j).xtarget,dati(j).y);
%end

loglog(beta,risultato,'-o');
legend('lunghezza di riscalamento')
xlabel('\\beta');
ylabel('lunghezza di riscalamento');
EOF

print $tail;
