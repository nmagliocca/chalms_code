%%%% relocation analyses

[starty,startx]=ind2sub([NLENGTH NWIDTH],reloc_stats{2,t});
[endy,endx]=ind2sub([NLENGTH NWIDTH],reloc_stats{3,t});

h1=figure;
set(h1,'Color','white')
plot(startx,starty,'.b','MarkerSize',15)
hold on
plot(endx,endy,'dr','MarkerSize',5)
warning('off')
legend('startpos','endpos','Location','NorthWest')
% reloc_xvec=sort([startx(1) endx(1)],'ascend');
% reloc_yvec=sort([starty(1) endy(1)],'ascend');
% plot(linspace(reloc_xvec(1),reloc_xvec(2),10),linspace(reloc_yvec(1),reloc_yvec(2),10),'-')
plot(linspace(startx(1),endx(1),3),linspace(starty(1),endy(1),3),'-')

% hold on
for i=2:length(reloc_stats{1,t})
%     reloc_xvec=sort([startx(i) endx(i)],'ascend');
%     reloc_yvec=sort([starty(i) endy(i)],'ascend');
%     plot(linspace(reloc_xvec(1),reloc_xvec(2),10),linspace(reloc_yvec(1),reloc_yvec(2),10),'-')
    plot(linspace(startx(i),endx(i),3),linspace(starty(i),endy(i),3),'-')
end
xlim([1 NWIDTH])
ylim([1 NLENGTH])
axis ij
title(sprintf('Relocations, t=%d',t))
set(gca,'xticklabel',[],'yticklabel',[])

%overlaid plots with LOTTYPE
h2=figure;
set(h2,'Color','white')
imagesc(reshape(LOTTYPE(:,t),NLENGTH,NWIDTH))
CMAP=colormap;
CMAP(1,:)=[1 1 1];
colormap(h2,CMAP);
colorbar
hold on
% plot(startx,starty,'.k','MarkerSize',15)
plot(startx,starty,'ok','MarkerSize',5)
% plot(endx,endy,'dr','MarkerSize',5,'MarkerFaceColor','red')
plot(endx,endy,'dr','MarkerSize',5)
warning('off')
legend('startpos','endpos','Location','NorthWest')
% reloc_xvec=sort([startx(1) endx(1)],'ascend');
% reloc_yvec=sort([starty(1) endy(1)],'ascend');
% plot(linspace(reloc_xvec(1),reloc_xvec(2),10),linspace(reloc_yvec(1),reloc_yvec(2),10),'-')
plot(linspace(startx(1),endx(1),3),linspace(starty(1),endy(1),3),'-k')

% hold on
for i=2:length(reloc_stats{1,t})
%     reloc_xvec=sort([startx(i) endx(i)],'ascend');
%     reloc_yvec=sort([starty(i) endy(i)],'ascend');
%     plot(linspace(reloc_xvec(1),reloc_xvec(2),10),linspace(reloc_yvec(1),reloc_yvec(2),10),'-')
    plot(linspace(startx(i),endx(i),3),linspace(starty(i),endy(i),3),'-k')
end
% xlim([1 NWIDTH])
% ylim([1 NLENGTH])
% axis ij
title(sprintf('Relocations, t=%d',t))
set(gca,'xticklabel',[],'yticklabel',[])
