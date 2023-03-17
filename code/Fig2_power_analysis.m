
%% cartoon - fig 2a
figure(1)
clf
N=100
e=0.15
p=0.2;
bcrit = binoinv([0.025 0.975],N,p);
barwid=0.8;
cmap=lines(1);

subplot(3,1,1)
bar(0:N,binopdf(0:N,N,p-e),barwid,'EdgeColor','none','FaceColor',cmap(1,:));
yl=ylim();
line([1 1]*bcrit(1),ylim())
plot_def()

subplot(3,1,2)
bar(0:N,binopdf(0:N,N,p),barwid,'EdgeColor','none','FaceColor',cmap(1,:));
ylim(yl)
line([1 1]*bcrit(1),ylim())
line([1 1]*bcrit(2),ylim())
plot_def()

subplot(3,1,3)
bar(0:N,binopdf(0:N,N,p+e),barwid,'EdgeColor','none','FaceColor',cmap(1,:));
ylim(yl)
line([1 1]*bcrit(2),ylim())
plot_def()

%% power calculation - fig 2b

nvec = [100 1000];
pvec = [0.2 0.1 0.05];
e80=[];
e80a=[];
zcrit = @(x) sqrt(2)*erfinv(2*x-1);
e = linspace(-1,1,2049);

power=[]; eq=[];

for j=1:length(pvec)
    for i=1:length(nvec)
    
        N=nvec(i);
        p = pvec(j);
        q = max(min(p+e,1),0);
        m = [0:N]';
        bcrit = p*N+[-1.96*sqrt(N*p*(1-p)) 1.96*sqrt(N*p*(1-p))];
        power(i,j,:) = normcdf(bcrit(1),N*q,sqrt(N.*q.*(1-q))) + normcdf(bcrit(2),N*q,sqrt(N.*q.*(1-q)),'upper');      
        
        eq(i,j,:) = 1+normcdf((sqrt(N*e.^2)+sqrt(p*(1-p))*-1.96)./sqrt(q.*(1-q)))-normcdf((sqrt(N*e.^2)+sqrt(p*(1-p))*1.96)./sqrt(q.*(1-q)));
        
    end
end

figure(2)
clf
subplot(1,2,1)
plot(e,squeeze(power(1,:,:)))

xl=[-1 1]*.2;
xlim(xl)
ylim([0 1])
line([0 0],ylim())
plot_def()
subplot(1,2,2)
plot(e,squeeze(power(2,:,:)))

xlim(xl)
ylim([0 1])
line([0 0],ylim())
plot_def()





%% approximation - fig 2c

e = linspace(-1,1,2049);
nvec = ceil(logspace(1,4,50));
pvec = [0.2 0.1 0.05];
e80=[];
e80a=[];
zcrit = @(x) sqrt(2)*erfinv(2*x-1);
est=[];

for j=1:length(pvec)
    for i=1:length(nvec)
    
        N=nvec(i);
        p = pvec(j);
        q = max(min(p+e,1),0);
        m = [0:N]';
        bcrit = p*N+[-1.96*sqrt(N*p*(1-p)) 1.96*sqrt(N*p*(1-p))];
        power = normcdf(bcrit(1),N*q,sqrt(N.*q.*(1-q))) + normcdf(bcrit(2),N*q,sqrt(N.*q.*(1-q)),'upper');
        e80a(i,j,1) = e(find(power<.8,1));
        if e80a(i,j,1)==-1, e80a(i,j,1)=NaN; end
        
        e80a(i,j,2) = e(find(power<.8 & q<1,1,'last'));
        if e80a(i,j,2)==1, e80a(i,j,2)=NaN; end
        
        
    end
    est(:,j,1) = p*(1-p)*(zcrit(0.05/2)-zcrit(1-0.2))./(sqrt(nvec*p*(1-p))-(p-1/2)*zcrit(1-0.2));
    est(:,j,2) = p*(1-p)*(zcrit(1-0.05/2)-zcrit(0.2))./(sqrt(nvec*p*(1-p))-(p-1/2)*zcrit(0.2));
end

p=0.5;
est0=p*(1-p)*(zcrit(0.8)-zcrit(0.05/2))./(sqrt(nvec*p*(1-p))+(-p+1/2)*zcrit(0.8));

figure(11)
clf
subplot(1,2,1)
semilogx(nvec,-squeeze(e80a(:,:,1)),'o')
hold on
plot(nvec,-squeeze(est(:,:,1)),'r','LineWidth',2)
plot(nvec,est0,'k','LineWidth',2)
hold off
plot_def()

subplot(1,2,2)
semilogx(nvec,squeeze(e80a(:,:,2)),'o')
hold on
plot(nvec,squeeze(est(:,:,2)),'r','LineWidth',2)
plot(nvec,est0,'k','LineWidth',2)
hold off
plot_def()

