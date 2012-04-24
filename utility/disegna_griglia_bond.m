imagesc(L');
colormap([1,1,1;0.5,0.5,0.5]);
hold on;

set(gca,'YTick',[])
set(gca,'XTick',[])

for i=1:10
    for j=1:10
        if(L(i,j)==1)
            plot(i,j,'sk','MarkerSize',5,'MarkerFaceColor','b');
        else
            plot(i,j,'ok','MarkerSize',5,'MarkerFaceColor','r');
        end
    end
end

for i=1:10
    for j=1:10
        if((i>1) && L(i,j)==L(i-1,j))
            line([i,i-1],[j,j],'LineWidth',3,'Color','k');
        end
        if((j>1) && L(i,j)==L(i,j-1))
            line([i,i],[j,j-1],'LineWidth',3,'Color','k');
        end
    end
end
