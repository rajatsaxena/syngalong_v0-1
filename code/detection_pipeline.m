function result = detection_pipeline(Tlist,params,result)

result.zmat = zeros(length(Tlist),length(Tlist));
result.ssi = zeros(length(Tlist),length(Tlist));
result.pall = zeros(length(Tlist),length(Tlist),sum(params.t>params.jit_min_t))*NaN;
if nargin<3,
    result.crosscorrelogram = zeros(length(Tlist),length(Tlist),length(params.t));
    result.spikes_pre = zeros(length(Tlist),length(Tlist));
end

tic
for i=1:length(Tlist)
    fprintf('Neuron %03i...',i)
    for j=(i+1):length(Tlist)
        if nargin<3,
            [tmp,~] = corr_fast_v3(Tlist{i},Tlist{j},params.trange(1),params.trange(2),params.tn);

            result.crosscorrelogram(i,j,:) = tmp(1:(end-1)); % last bin of histc contains data == the last bin
            result.crosscorrelogram(j,i,:) = flipud(tmp(1:(end-1)));
            result.spikes_pre(i,j) = length(Tlist{i});
            result.spikes_pre(j,i) = length(Tlist{j});
        end
        
        % Approximate jitter method statistics
        if isfield(params,'mask')
            m=params.t'*0+NaN; s=params.t'*0+NaN; p=params.t'*0+NaN;
            [m(params.mask),s(params.mask),p(params.mask)] = conv_jit(squeeze(result.crosscorrelogram(i,j,params.mask)),params.jit_conv_window,true);
        else
            [m,s,p] = conv_jit(squeeze(result.crosscorrelogram(i,j,:)),params.jit_conv_window,true);
        end

        [p_ext(i,j),p1_loc] = min(min([p(params.t>params.jit_min_t) 1-p(params.t>params.jit_min_t)],[],2));
        p_ext_bin(i,j)=sum(params.t<=params.jit_min_t)+p1_loc;
        p_sign(i,j)=(p(p_ext_bin(i,j))>.5)*2-1;
%         
        result.zmat(i,j) = (result.crosscorrelogram(i,j,p_ext_bin(i,j))-m(p_ext_bin(i,j)))/s(p_ext_bin(i,j));
        [p_ext(j,i),p_ext_bin(j,i)] = min(min([p(params.t<-params.jit_min_t) 1-p(params.t<-params.jit_min_t)],[],2));
        p_sign(j,i)=(p(p_ext_bin(j,i))>.5)*2-1;
        result.zmat(j,i) = (result.crosscorrelogram(i,j,p_ext_bin(j,i))-m(p_ext_bin(j,i)))/s(p_ext_bin(j,i));
        
        result.pall(i,j,:) = p(params.t>params.jit_min_t);
        result.pall(j,i,:) = p(params.t<-params.jit_min_t);

        
        % Likelihood ratio relative to the saturated model...
        y=squeeze(result.crosscorrelogram(i,j,:));
        yij = y(params.t>params.jit_min_t); 
        mij = m(params.t>params.jit_min_t);
        N = result.spikes_pre(i,j);
        result.llr(i,j) = ((nansum(yij.*log(yij/N)+(N-yij).*log(1-yij/N))-nansum(yij.*log(mij/N)+(N-yij).*log(1-mij/N))))/length(yij)/log(2);
        llrt = (((yij.*log(yij/N+(yij==0))+(N-yij).*log(1-yij/N))-(yij.*log(mij/N+(mij==0))+(N-yij).*log(1-mij/N))))/log(2);
        result.llrmax(i,j) = max(llrt);
        result.llrmax2(i,j) = max(llrt(1:end-1)+llrt(2:end));
        result.llrmax3(i,j) = max(llrt(1:end-2)+llrt(2:end-1)+llrt(3:end));
%         result.pall(i,j,:) = chi2cdf(result.pall(i,j,:),1);
        yji = y(params.t<-params.jit_min_t); 
        mji = m(params.t<-params.jit_min_t);
        N = result.spikes_pre(j,i);
        result.llr(j,i) = ((nansum(yji.*log(yji/N)+(N-yji).*log(1-yji/N))-nansum(yji.*log(mji/N)+(N-yji).*log(1-mji/N))))/length(yji)/log(2);
        llrt = (((yji.*log(yji/N+(yji==0))+(N-yji).*log(1-yji/N))-(yji.*log(mji/N+(mji==0))+(N-yji).*log(1-mji/N))))/log(2);
        result.llrmax(j,i) = max(llrt);
%         result.pall(j,i,:) = chi2cdf(result.pall(j,i,:),1);
        result.llrmax2(j,i) = max(llrt(1:end-1)+llrt(2:end));
        result.llrmax3(j,i) = max(llrt(1:end-2)+llrt(2:end-1)+llrt(3:end));
        
        % Spike shadowing index
        mc = mean(result.crosscorrelogram(i,j,abs(params.t)<params.ssi_min_t));
        mleft = mean(result.crosscorrelogram(i,j,params.t<-params.ssi_min_t & params.t>-2*params.ssi_min_t));
        mrigt = mean(result.crosscorrelogram(i,j,params.t>params.ssi_min_t & params.t<2*params.ssi_min_t));
        result.ssi(i,j) = min(abs(mc-mleft)/mleft,abs(mc-mrigt)/mrigt);
        result.ssi(j,i) = result.ssi(i,j);
        
        % Exclude pairs with spike shadowing
        if result.ssi(i,j)>params.ssi_max || ~isfinite(result.ssi(i,j))
            result.pall(i,j,:) = NaN;
            result.pall(j,i,:) = NaN;
        end
    end
    toc
end


% fdr
pall_one_sided = result.pall;
pall_one_sided(pall_one_sided>0.5)=1-pall_one_sided(pall_one_sided>0.5);
pall_sorted = sort(pall_one_sided(:));
m = sum(isfinite(pall_sorted));
pall_sorted_k = find(pall_sorted>([1:length(pall_sorted)]'/m*params.alph),1);
result.detected_cnx = sum(pall_one_sided<pall_sorted(pall_sorted_k),3)>0;
% 
% % llr-based
% Z = result.llr(:);
% Z(result.ssi>params.ssi_max | ~isfinite(result.ssi))=0;
% [zs,sidx] = sort(abs(Z),'descend');
% sidx=sidx(zs>1 & result.llrmax3(sidx)>15);
% result.detected_cnx = result.detected_cnx*0;
% result.detected_cnx(sidx)=1;