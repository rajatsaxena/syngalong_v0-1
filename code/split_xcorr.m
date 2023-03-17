
function [C,bins,n] = split_xcorr(ij,Tlist,params,splitby,splitN,doPlot)

n=[];
lower(splitby)
switch lower(splitby)
    case 'isi'
        isi = [NaN; diff(Tlist{ij(1)})];
        qisi = prctile(isi,linspace(0,100,splitN+1));
        [~,bins]=histc(isi,qisi);
        
        for i=1:splitN
            id = isi>qisi(i) & isi<=qisi(i+1);
            [c,deltaT] = spk_xcorr(Tlist{ij(1)}(id),Tlist{ij(2)},params.interval(1),params.interval(2),params.num_bins);
            C(i,:)=c;
            n(i)=sum(id);
            
            if doPlot
                subplot(2,splitN,i)
                bar(params.t,c,1,'EdgeColor','none')
                title(num2str(floor(qisi(i)*1000)))
                subplot(2,splitN,i+splitN)
                a = spk_xcorr(Tlist{ij(1)}(id),Tlist{ij(1)},params.interval(1),params.interval(2),params.num_bins);
                a(a==max(a))=NaN;
                A(i,:)=a;
                bar(params.t,a,1,'EdgeColor','none')
            end
        end
    case 'isi_log'
        isi = [NaN; diff(Tlist{ij(1)})];
        %         qisi = prctile(isi,linspace(0,100,splitN+1));
        qisi = logspace(log10(min(isi)),log10(max(isi)),splitN+1);
        [~,bins]=histc(isi,qisi);
        
        for i=1:splitN
            id = isi>qisi(i) & isi<=qisi(i+1);
            if sum(id)>0
                [c,deltaT] = spk_xcorr(Tlist{ij(1)}(id),Tlist{ij(2)},params.interval(1),params.interval(2),params.num_bins);
                C(i,:)=c;
            else
                C(i,:)=NaN;
                c=NaN;
            end
            n(i)=sum(id);
            
            if doPlot
                subplot(2,splitN,i)
                bar(params.t,c,1,'EdgeColor','none')
                title(num2str(floor(qisi(i)*1000)))
                subplot(2,splitN,i+splitN)
                a = spk_xcorr(Tlist{ij(1)}(id),Tlist{ij(1)},params.interval(1),params.interval(2),params.num_bins);
                a(a==max(a))=NaN;
                A(i,:)=a;
                bar(params.t,a,1,'EdgeColor','none')
            end
        end
    case 'time'
        tvec = linspace(0,max(cellfun(@max,Tlist)),splitN+1);
        
        bins=zeros(length(Tlist{ij(1)}),1);
        for i=1:splitN
            id = Tlist{ij(1)}>tvec(i) & Tlist{ij(1)}<=tvec(i+1);
            bins(id) = i;
            if sum(id)>0
                [c,deltaT] = spk_xcorr(Tlist{ij(1)}(id),Tlist{ij(2)},params.interval(1),params.interval(2),params.num_bins);
                c = c/sum(id);
                C(i,:)=c;
                n(i)=sum(id);
                if doPlot
                    subplot(2,splitN,i)
                    bar(params.t,c,1,'EdgeColor','none')
                    subplot(2,splitN,i+splitN)
                    a = spk_xcorr(Tlist{ij(1)}(id),Tlist{ij(1)},params.interval(1),params.interval(2),params.num_bins);
                    a(a==max(a))=NaN;
                    A(i,:)=a;
                    bar(params.t,a,1,'EdgeColor','none')
                end
            end
        end
    case 'epoch'
        tmat=splitN;
        splitN=size(tmat,1);
        bins=zeros(length(Tlist{ij(1)}),1);
        for i=1:splitN
            id = Tlist{ij(1)}>tmat(i,1) & Tlist{ij(1)}<=tmat(i,2);
            bins(id) = i;
            if sum(id)>0
                [c,deltaT] = spk_xcorr(Tlist{ij(1)}(id),Tlist{ij(2)},params.interval(1),params.interval(2),params.num_bins);
                c = c/sum(id);
                C(i,:)=c;
                n(i)=sum(id);
                if doPlot
                    subplot(2,splitN,i)
                    bar(params.t,c,1,'EdgeColor','none')
                    subplot(2,splitN,i+splitN)
                    a = spk_xcorr(Tlist{ij(1)}(id),Tlist{ij(1)},params.interval(1),params.interval(2),params.num_bins);
                    a(a==max(a))=NaN;
                    a = a/sum(id);
                    A(i,:)=a;
                    bar(params.t,a,1,'EdgeColor','none')
                end
            end
        end
    case 'lastspk'
        isi = getLastSpk(Tlist{ij(1)},Tlist{ij(2)});
        qisi = prctile(isi,linspace(0,100,splitN+1));
        [~,bins]=histc(isi,qisi);
        
        for i=1:splitN
            id = isi>qisi(i) & isi<=qisi(i+1);
            [c,deltaT] = spk_xcorr(Tlist{ij(1)}(id),Tlist{ij(2)},params.interval(1),params.interval(2),params.num_bins);
            C(i,:)=c;
            
            if doPlot
                subplot(2,splitN,i)
                bar(params.t,c,1,'EdgeColor','none')
                subplot(2,splitN,i+splitN)
                a = spk_xcorr(Tlist{ij(1)}(id),Tlist{ij(1)},params.interval(1),params.interval(2),params.num_bins);
                a(a==max(a))=NaN;
                A(i,:)=a;
                bar(params.t,a,1,'EdgeColor','none')
            end
        end
    case 'iso'
        isi1 = [NaN; diff(Tlist{ij(1)})];
        isi2 = [diff(Tlist{ij(1)}); NaN];
        isi = nanmin([isi1 isi2],[],2);
        qisi = prctile(isi,linspace(0,100,splitN+1));
        [~,bins]=histc(isi,qisi);
        
        for i=1:splitN
            id = isi>qisi(i) & isi<=qisi(i+1);
            [c,deltaT] = corr_fast_v3(Tlist{ij(1)}(id),Tlist{ij(2)},params.interval(1),params.interval(2),params.num_bins);
            C(i,:)=c;
            
            if doPlot
                %                 subplot(1,splitN,i)
                %                 bar(params.t,c,1,'EdgeColor','none')
                
                subplot(2,splitN,i)
                bar(params.t,c,1,'EdgeColor','none')
                title(num2str(floor(qisi(i)*1000)))
                
                
                subplot(2,splitN,i+splitN)
                a = corr_fast_v3(Tlist{ij(1)}(id),Tlist{ij(1)},params.interval(1),params.interval(2),params.num_bins);
                a(a==max(a))=NaN;
                A(i,:)=a;
                bar(params.t,a,1,'EdgeColor','none')
            end
        end
    case 'iso_sign'
        isi1 = [NaN; diff(Tlist{ij(1)})];
        isi2 = [diff(Tlist{ij(1)}); NaN];
        isi = nanmin([isi1 isi2],[],2);
        isi(isi==isi1)=-isi(isi==isi1);
        qisi = prctile(isi,linspace(0,100,splitN+1));
        [~,bins]=histc(isi,qisi);
        
        for i=1:splitN
            id = isi>qisi(i) & isi<=qisi(i+1);
            [c,deltaT] = corr_fast_v3(Tlist{ij(1)}(id),Tlist{ij(2)},params.interval(1),params.interval(2),params.num_bins);
            C(i,:)=c;
            
            if doPlot
                %                 subplot(1,splitN,i)
                %                 bar(params.t,c,1,'EdgeColor','none')
                
                subplot(2,splitN,i)
                bar(params.t,c,1,'EdgeColor','none')
                subplot(2,splitN,i+splitN)
                a = corr_fast_v3(Tlist{ij(1)}(id),Tlist{ij(1)},params.interval(1),params.interval(2),params.num_bins);
                a(a==max(a))=NaN;
                A(i,:)=a;
                bar(params.t,a,1,'EdgeColor','none')
            end
        end
    otherwise
        print('problem')
end

if doPlot
    for i=1:splitN
        %         subplot(1,splitN,i)
        subplot(2,splitN,i)
        ylim([0 max(C(:))*1.1])
        box off; set(gca,'TickDir','out');
        subplot(2,splitN,i+splitN)
        ylim([0 max(A(:))])
        box off; set(gca,'TickDir','out');
    end
end