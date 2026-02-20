colors = wes_palettes('FantasticFox2');

close all
figure
set(gcf,'color','white')
nbars = size(colors,1);
hold on
for i = 1:nbars
    bar(i,nbars+0.5,1,'FaceColor',colors(i,:),'EdgeColor',[1 1 1]*1);
end
set(gca,'xtick',[],'ytick',[])
xlim([0.5 nbars+0.5])
ylim([0.5 nbars+0.5])
axis square