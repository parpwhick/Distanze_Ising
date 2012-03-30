imagesc(L');
colormap([1,1,1;0.7,0.7,0.7]);
hold on;

for i=1:20
    for j=1:20
        if(L(i,j)==1)
            plot(i,j,'^k');
        else
            plot(i,j,'+k');
        end
    end
end

for i=1:20
    for j=1:20
        if((i>1) && L(i,j)==L(i-1,j))
            line([i,i-1],[j,j],'LineWidth',3,'Color','k');
        end
        if((j>1) && L(i,j)==L(i,j-1))
            line([i,i],[j,j-1],'LineWidth',3,'Color','k');
        end
    end
end
