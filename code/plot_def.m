function plot_def()

box off
set(gca,'TickDir','out')
set(gca,'TickLength',[.02;.02])
set(gca,'LabelFontSizeMultiplier',1.25)

t=get(gca,'Children');
for i=1:length(t)
    if strcmp(t(i).Type,'line')
        set(t(i),'LineWidth',2)
    end
    if strcmp(t(i).Type,'histogram')
        set(t(i),'EdgeColor','none')
    end
end