clear;close all;clc;

mu  = 0.8;
D   = 0.3;
dx  = 0.1;
dt  = 0.01;
z   = -20:dx:20;
t   = 0:dt:10;

p   = exp(-z.^2/2) / sqrt(2*pi);
p(1) = 0; p(end) = 0;
p   = p / sum(p);
pent = p;
pmax = max(p);
zmax = max(z);

psave = p;

z = z/zmax;
p = p/pmax;

N = 7;
c2 = 2;
xmin = 0.2;
xmax = 0.4;

inds = find(z>=xmin & z<=xmax);
zred = z(inds);
pred = p(inds);

inds21 = find(z<=xmin);
zblack1 = z(inds21);
pblack1 = p(inds21);

inds22 = find(z>=xmax);
zblack2 = z(inds22);
pblack2 = p(inds22);

h = figure;
h.Position = [723 1 302 536];
subplot(2,1,1)
plot3(zeros(1,numel(zblack1)),zblack1,pblack1,'Color',[0 0 0 0.3],'LineWidth',2)
hold on
plot3(zeros(1,numel(zblack2)),zblack2,pblack2,'Color',[0 0 0 0.3],'LineWidth',2)
plot3(zeros(1,numel(zred)),zred,pred,'Color',[1 0 0 0.3],'LineWidth',3)
grid on
axis([0 142 -0.1 0.8 0 1])
view(119,28)

logp        = log2(pent);
logp(pent==0)  = 0;
sig          = -pent.*logp;
sigpred1     = sig;
sigpred2     = sig;
sigpred3     = sig;
Sroi        = sum(sig(find(z==xmin):find(z==xmax)));
Sroi_pred_full   = Sroi;
Sroi_pred_cons_only   = Sroi;
Sroi_pred3_noncons_only   = Sroi;

c = 0;
for k = 2:length(t)
    dpdx = [0, diff(p)/dx];
    d2pdx2 = [0, diff(p, 2)/dx^2, 0];
    p = p + dt * (-mu*dpdx + D*d2pdx2);
    p(1) = 0; p(end) = 0;
    p = p / sum(p);
    pent = p;
    p = p/pmax;

    psave(k,:) = p;

    pred = p(inds);
    pblack1 = p(inds21);
    pblack2 = p(inds22);

    if numel(regexp(num2str(k/N/c2),'\.','split')) == 1
        c = c + c2;
        plot3(c*ones(1,numel(zblack1)),zblack1,pblack1,'Color',[0 0 0 0.3],'LineWidth',2)
        plot3(c*ones(1,numel(zblack2)),zblack2,pblack2,'Color',[0 0 0 0.3],'LineWidth',2)
        plot3(c*ones(1,numel(zred)),zred,pred,'Color',[1 0 0 0.3],'LineWidth',3)
        axis([0 142 -0.1 0.8 0 1])
        drawnow
    end

    logp        = log2(pent);
    logp(pent==0)  = 0;
    sig         = -pent.*logp;
    Sroi(k)     = sum(sig(find(z==xmin):find(z==xmax)));
    dsigdz      = [0, diff(sig)/dx];
    d2sigdz2    = [0, diff(sig,2)/dx^2, 0];
    dlogpdz     = [0, diff(logp)/dx];
    sigpred1     = sigpred1 + dt*(D*d2sigdz2 - mu*dsigdz + D*pent.*(dlogpdz).^2);
    sigpred2     = sigpred2 + dt*(D*d2sigdz2 - mu*dsigdz);
    sigpred3     = sigpred3 + dt*(D*pent.*(dlogpdz).^2);
    Sroi_pred_full(k)   = sum(sigpred1(find(z==xmin):find(z==xmax)));
    Sroi_pred_cons_only(k)   = sum(sigpred2(find(z==xmin):find(z==xmax)));
    Sroi_pred3_noncons_only(k)   = sum(sigpred3(find(z==xmin):find(z==xmax)));

end
hold off

subplot(2,1,2)
plot(Sroi(2:end),'k','Linewidth',3)
hold on
plot(Sroi_pred_full(2:end),'g','Linewidth',3)
plot(Sroi_pred_cons_only(2:end),'c','Linewidth',3)
plot(Sroi_pred3_noncons_only(2:end),'r','Linewidth',3)
legend('ground truth','full model','cons. only','non cons. only','Position',[0.65 0.2 0.05 0.1])
hold off