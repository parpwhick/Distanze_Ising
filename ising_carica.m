clear dist*;
fclose('all');
clear fid*;

bins=50;

have_fuzzy = ~isempty(dir('output-fuzzy*'));

if(have_fuzzy)
    fid6=fopen('output-fuzzy.bin','r');
    fid7=fopen('output-fuzzyt.bin','r');
    
    %distanze fuzzy
    dist_f=fread(fid6,[n,n],'double');
    %distanze fuzzy topologiche
    dist_f_t=fread(fid7,[n,n],'double');
    
    [hist_f,labels6]=hist(dist_f(dist_f>0),100);
    [hist_f_t,labels7]=hist(dist_f_t(dist_f_t>0),100);
    
    dist_f=dist_f+dist_f';
    dist_f_t=dist_f_t+dist_f_t';
    
    plot(labels6,hist_f,'-*r',labels7,hist_f_t,'-*b');
    legend('Salto','Salto top.');
    xlabel('Distanza');
    ylabel('Frequenza');
    title('Partizioni con salto');
    
end