function plotPairs(result,params,sidx,offset)

if nargin<3, sidx=1:size(result.crosscorrelogram,1); end
if nargin<4, offset=1; end

figure(1)
clf
if size(sidx,2)>1
    ijlist = sidx;
    ijlist = ijlist(offset:end,:);
    for c=1:min(size(ijlist,1),20)
        subplot(5,4,c)
        i=ijlist(c,1);
        j=ijlist(c,2);
        
        bar(params.t,squeeze(result.crosscorrelogram(i,j,:)),1,'EdgeColor','none')
        hold on
        % plot(y(k,:))
        line([0 0],ylim(),'Color','k')
%         [m,s,p] = conv_jit(squeeze(result.crosscorrelogram(i,j,:)),params.jit_conv_window,true);
%         plot(params.t,m)
%         plot(params.t,m+s)
%         plot(params.t,m-s)
        title(sprintf('ss%.02f--zz%.02f',result.ssi(i,j),result.zmat(i,j)))
        ylabel([num2str(i) '->' num2str(j)])
        
        plot_def()
        hold off
    end
else
    
    
    [J,I]=meshgrid(1:size(result.crosscorrelogram,1),1:size(result.crosscorrelogram,1));
    ijlist = [I(sidx) J(sidx)];
    ijlist=ijlist(offset:end,:);
    
    for c=1:min(size(ijlist,1),20)
        subplot(5,4,c)
        
        i=ijlist(c,1);
        j=ijlist(c,2);
        
        bar(params.t,squeeze(result.crosscorrelogram(i,j,:)),1,'EdgeColor','none')
        hold on
        % plot(y(k,:))
        line([0 0],ylim(),'Color','k')
        [m,s,p] = conv_jit(squeeze(result.crosscorrelogram(i,j,:)),params.jit_conv_window,true);
        plot(params.t,m)
        plot(params.t,m+s)
        plot(params.t,m-s)
        title(sprintf('ss%.02f--zz%.02f',result.ssi(i,j),result.zmat(i,j)))
        ylabel([num2str(i) '->' num2str(j)])
        
        plot_def()
        hold off
    end
end