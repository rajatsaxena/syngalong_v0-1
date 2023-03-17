function plotSubnet(params,subnet,Tlist,type_est)

if nargin<4,
    type_est=subnet*0+1;
end
cmap = lines(max(type_est));

% [J,I]=meshgrid(1:size(result.crosscorrelogram,1),1:size(result.crosscorrelogram,1));
% ijlist = [I(sidx) J(sidx)];

maxc = zeros(length(subnet));

for i=1:length(subnet)
    for j=i:length(subnet)
                
        subplot(length(subnet),length(subnet),(i-1)*length(subnet)+j)
        
        scale = 1;
        corr=0;
%         scale = 1/length(Tlist{subnet(i)});
%         scale = 1/(((params.trange(2)-params.trange(1)))/(params.tn-2))/max(cellfun(@max,Tlist));
%         scale = 1/max(cellfun(@max,Tlist));
        dt = (((params.trange(2)-params.trange(1)))/(params.tn-2));
%         scale = 1/max(cellfun(@max,Tlist))/dt;
%         bins = max(cellfun(@max,Tlist))/dt;
        
        scale = 1/length(Tlist{subnet(i)})/dt;

%         scale = 1/bins/(length(Tlist{j})/bins*dt);
%         corr = length(Tlist{i})/bins/dt;
        
        
%         c = squeeze(result.crosscorrelogram(subnet(i),subnet(j),:));
        [c,~] = corr_fast_v3(Tlist{subnet(i)},Tlist{subnet(j)},params.trange(1),params.trange(2),params.tn);
        c=c(1:end-1);
%         if result.detected_cnx(subnet(i),subnet(j)) || result.detected_cnx(subnet(j),subnet(i))
%             bar(params.t,c*scale-corr,1,'EdgeColor','none')
%         elseif result.ssi(subnet(i),subnet(j))>params.ssi_max
%             bar(params.t,c*scale-corr,1,'EdgeColor','none','FaceAlpha',.5,'EdgeAlpha',.5,'FaceColor','r')
%         else
%             bar(params.t,c*scale-corr,1,'EdgeColor','none','FaceAlpha',.5,'EdgeAlpha',.5)
%         end
        bar(params.t,c*scale-corr,1,'EdgeColor','none','FaceColor',cmap(type_est(i),:))
        if i==j
            [tmp,~] = corr_fast_v3(Tlist{subnet(i)},Tlist{subnet(j)},params.trange(1),params.trange(2),params.tn);
            tmp(tmp==max(tmp))=NaN;
            c = tmp(1:end-1);
            bar(params.t,c*scale-corr,1,'EdgeColor','none','FaceColor',cmap(type_est(i),:))
        end
        maxc(i,j)=max(c*scale-corr);
%         hold on
% %         line([0 0],ylim(),'Color','k')
%         hold off
        
%         ylim([0 150])

        if j==i
%             ylabel(num2str(subnet(i)))
            xlabel(num2str(subnet(j)))
        end
        if i~=length(subnet)
            set(gca,'XTickLabel',[])
        end
        plot_def()
    end    
end

for i=1:length(subnet)
    for j=i:length(subnet)
        subplot(length(subnet),length(subnet),(i-1)*length(subnet)+j)
        ylim([0 max(maxc(:,j))])
        if i~=j
            set(gca,'YTickLabel',[])
        end
    end
end
