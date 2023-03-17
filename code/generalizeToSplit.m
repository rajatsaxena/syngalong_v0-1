
function generalizeToSplit(ij,Tlist,params,fit,model_name,bins,refit_baseline)

if nargin<7,
    refit_baseline = true;
end

if params.nknots<2
    params.nknots=-1;
end

for i=1:max(bins)
    % Cross correlogram
    [c,~] = spk_xcorr(Tlist{ij(1)}(bins==i),Tlist{ij(2)},params.interval(1),params.interval(2),params.num_bins);
    N = sum(bins==i);

    % Auto correlograms (for convolutional models)
    a1 = spk_xcorr(Tlist{ij(1)}(bins==i),Tlist{ij(1)}(bins==i),2*params.interval(1)+params.dt/2,2*params.interval(2)-params.dt/2,2*(params.num_bins-1)+1);
    a2 = spk_xcorr(Tlist{ij(2)},Tlist{ij(2)},2*params.interval(1)+params.dt/2,2*params.interval(2)-params.dt/2,2*(params.num_bins-1)+1);

    [mdl,atoms] = composite_model_constructor(model_name,params,a1,a2);

    Xout = mdl.cov_construct(fit.feat_params);

    if refit_baseline
        geff = Xout*fit.glm_params;
        bb = glmfit(Xout(:,1:(params.nknots+2)),[c repmat(N,size(c))],'binomial','constant','off','offset',Xout(:,(params.nknots+3):end)*fit.glm_params((params.nknots+3):end));
        bout2 = fit.glm_params;
        bout2(1:(params.nknots+2))=bb;
    else
        bb=0;
        bout2=fit.glm_params;
    end

    yhat = 1./(1+exp(-Xout*bout2));
    lamsmoo = 1./(1+exp(-Xout(:,1:(params.nknots+2))*bout2(1:(params.nknots+2))));

    subplot(2,max(bins),2*i-1)
    [ap,~] = spk_xcorr(Tlist{ij(1)}(bins==i),Tlist{ij(1)},params.interval(1),params.interval(2),params.num_bins);
    ap(params.t==0)=NaN;
    bar(params.t,ap,1,'EdgeColor','flat')
    
    subplot(2,max(bins),2*i)
    hold on
    bar(params.t,c,1,'EdgeColor','flat')
    %         hold on
    plot(params.t,yhat*N,'LineWidth',2)
    plot(params.t,lamsmoo*N)
    hold off
    drawnow
end