function [out] = plot_flxdn_magnetizingFluxLinkage(pathname,filename)

if nargin()<2
    [filename, pathname, ~] = uigetfile(checkPathSyntax([cd '\.mat']), 'LOAD DATA');
end

load([pathname filename])

out.SOL.fMagFund  = nan(size(out.SOL.th));
out.SOL.fMagHarm  = nan(size(out.SOL.th));
out.SOL.fMagStat  = nan(size(out.SOL.th));

[kw,~] = calcKwTh0(geo);

if geo.ps==1
    x = out.SOL.Bg(:,1);
    for ii=1:length(out.SOL.th)
        % sin and cos reference for the alignment during rotation
        fundRef.sin = sin((x-out.SOL.th(ii)/geo.p+out.SOL.th(1)/geo.p)*pi/180*geo.p);
        fundRef.cos = cos((x-out.SOL.th(ii)/geo.p+out.SOL.th(1)/geo.p)*pi/180*geo.p);
        y = out.SOL.Bg(:,ii+1);
        % total with harmonics
        tmpSin = mean(y.*sign(fundRef.sin))*pi*(geo.r+geo.g/2)/1000/geo.p*geo.l/1000*geo.win.Ns*kw(1);
        tmpCos = mean(y.*sign(fundRef.cos))*pi*(geo.r+geo.g/2)/1000/geo.p*geo.l/1000*geo.win.Ns*kw(1);
        out.SOL.fMagHarm(ii) = (tmpSin.^2+tmpCos.^2).^0.5;
        % fundamental
        tmpSin = trapz(x/(x(end)),y.*fundRef.sin)*pi*(geo.r+geo.g/2)/1000/geo.p*geo.l/1000*geo.win.Ns*kw(1);
        tmpCos = trapz(x/(x(end)),y.*fundRef.cos)*pi*(geo.r+geo.g/2)/1000/geo.p*geo.l/1000*geo.win.Ns*kw(1);
        out.SOL.fMagFund(ii) = (tmpSin.^2+tmpCos.^2).^0.5;
        % total with stator harmonics
        ns = 6*geo.q;
        yTot = [y;-y];
        a = fft(yTot);
        harmStat = [1 (1:1:1000)*ns+1 (1:1:1000)*ns-1];
        harmStat = sort(harmStat);
        harmStatOK = harmStat(harmStat<numel(a)/2-1);
        aHarm = a(2:end-1);
        aHarmNew = zeros(size(aHarm));
        for jj=1:length(harmStatOK)
            aHarmNew(harmStatOK(jj)) = aHarm(harmStatOK(jj));
            aHarmNew(end-harmStatOK(jj)+1) = aHarm(end-harmStatOK(jj)+1);
        end
        aNew = [0; aHarmNew; 0];
        yNew = ifft(aNew,'symmetric');
        yNewHalf = yNew(1:end/2);
        tmpSin = mean(yNewHalf.*sign(fundRef.sin))*pi*(geo.r+geo.g/2)/1000/geo.p*geo.l/1000*geo.win.Ns;
        tmpCos = mean(yNewHalf.*sign(fundRef.cos))*pi*(geo.r+geo.g/2)/1000/geo.p*geo.l/1000*geo.win.Ns;
        out.SOL.fMagStat(ii) = (tmpSin.^2+tmpCos.^2).^0.5;

        % old
        %out.SOL.fMagHarm(ii) = mean(y.*sign(fundRef))*pi*geo.r/1000/geo.p*geo.l/1000*geo.win.Ns;
        %out.SOL.fMagFund(ii) = trapz(x/(x(end)),y.*fundRef)*pi*geo.r/1000/geo.p*geo.l/1000*geo.win.Ns;
    end
else
    error('Work in progress...')
end

% multi-3-phase computation
% if geo.win.n3phase>1
for ii=1:geo.win.n3phase
    th = out.SOL.th;

    fa = out.SOL.fa(ii,:);
    fb = out.SOL.fb(ii,:);
    fc = out.SOL.fc(ii,:);
    fdq = abc2dq(fa,fb,fc,(th+geo.th0(ii)-geo.th0(1))*pi/180);
    out.SOL.sets(ii).fa = fa;
    out.SOL.sets(ii).fb = fb;
    out.SOL.sets(ii).fc = fc;
    out.SOL.sets(ii).fd = fdq(1,:);
    out.SOL.sets(ii).fq = fdq(2,:);
    out.SOL.sets(ii).f0 = (fa+fb+fc)/3;

    ia = out.SOL.ia(ii,:);
    ib = out.SOL.ib(ii,:);
    ic = out.SOL.ic(ii,:);
    idq = abc2dq(ia,ib,ic,(th+geo.th0(ii)-geo.th0(1))*pi/180);
    out.SOL.sets(ii).ia = ia;
    out.SOL.sets(ii).ib = ib;
    out.SOL.sets(ii).ic = ic;
    out.SOL.sets(ii).id = idq(1,:);
    out.SOL.sets(ii).iq = idq(2,:);
    out.SOL.sets(ii).i0 = (ia+ib+ic)/3;
end
% end



% average computations
out.fMagFund = mean(out.SOL.fMagFund);
out.fMagHarm = mean(out.SOL.fMagHarm);
out.fMagStat = mean(out.SOL.fMagStat);
if isfield(out.SOL,'sets')
    for ii=1:length(out.SOL.sets)
        out.(['fd' int2str(ii)]) = mean(out.SOL.sets(ii).fd);
        out.(['fq' int2str(ii)]) = mean(out.SOL.sets(ii).fq);
    end
end

% wf - refers to complete waveform
nRep = 360/per.delta_sim_singt; % number of repetition needed

wf.fd       = repmat(out.SOL.fd,1,nRep);            
wf.fq       = repmat(out.SOL.fq,1,nRep);            
wf.fMagFund = repmat(out.SOL.fMagFund,1,nRep);
wf.fMagHarm = repmat(out.SOL.fMagHarm,1,nRep);
wf.fMagStat = repmat(out.SOL.fMagStat,1,nRep);
wf.T        = repmat(out.SOL.T,1,nRep);              
wf.th       = linspace(0,360,length(wf.fd)+1);
wf.th       = wf.th(1:end-1);

fph = phaseQuantityDecoding(out.SOL.fa,out.SOL.fb,out.SOL.fc,per.delta_sim_singt);
wf.fa = [fph.a fph.a(:,1)];
wf.fb = [fph.b fph.b(:,1)];
wf.fc = [fph.c fph.c(:,1)];

if isfield(out.SOL,'sets')
    for ii=1:length(out.SOL.sets)
        wf.sets(ii).fd = repmat(out.SOL.sets(ii).fd,1,nRep);
        wf.sets(ii).fq = repmat(out.SOL.sets(ii).fq,1,nRep);
        wf.sets(ii).id = repmat(out.SOL.sets(ii).id,1,nRep);
        wf.sets(ii).iq = repmat(out.SOL.sets(ii).iq,1,nRep);
        
        tmp = phaseQuantityDecoding(out.SOL.sets(ii).fa,out.SOL.sets(ii).fb,out.SOL.sets(ii).fc,per.delta_sim_singt);
        wf.sets(ii).fa = tmp.a;
        wf.sets(ii).fb = tmp.b;
        wf.sets(ii).fc = tmp.c;
        wf.sets(ii).f0 = (tmp.a+tmp.b+tmp.c)/3;

        tmp = phaseQuantityDecoding(out.SOL.sets(ii).ia,out.SOL.sets(ii).ib,out.SOL.sets(ii).ic,per.delta_sim_singt);
        wf.sets(ii).ia = tmp.a;
        wf.sets(ii).ib = tmp.b;
        wf.sets(ii).ic = tmp.c;
        wf.sets(ii).i0 = (tmp.a+tmp.b+tmp.c)/3;
    end
end

% plots

nRep = 360/per.delta_sim_singt; % number of repetition needed

fd       = [repmat(out.SOL.fd,1,nRep) out.SOL.fd(1)];  % last point added for plot
fq       = [repmat(out.SOL.fq,1,nRep) out.SOL.fq(1)];  % last point added for plot
fMagFund = [repmat(out.SOL.fMagFund,1,nRep) out.SOL.fMagFund(1)];  % last point added for plot
fMagHarm = [repmat(out.SOL.fMagHarm,1,nRep) out.SOL.fMagHarm(1)];  % last point added for plot
fMagStat = [repmat(out.SOL.fMagStat,1,nRep) out.SOL.fMagStat(1)];  % last point added for plot

fph = phaseQuantityDecoding(out.SOL.fa,out.SOL.fb,out.SOL.fc,per.delta_sim_singt);
fa = [fph.a fph.a(:,1)];
fb = [fph.b fph.b(:,1)];
fc = [fph.c fph.c(:,1)];

% if isfield(out.SOL,'sets')
%     for ii=1:length(out.SOL.sets)
%         sets(ii).fd = [repmat(out.SOL.sets(ii).fd,1,nRep) out.SOL.sets(ii).fd(1)];  % last point added for plot
%         sets(ii).fq = [repmat(out.SOL.sets(ii).fq,1,nRep) out.SOL.sets(ii).fq(1)];  % last point added for plot
%     end
% end

th = linspace(0,360,length(fd));

fMax = max([fd(:);fq(:);fMagHarm(:);fMagFund(:);fMagStat(:);fa(:);fb(:);fc(:)]);

for ii=1:2
    hfig(ii) = figure();
    figSetting()
    hax(ii) = axes('OuterPosition',[0 0 1 1]);
    xlabel('$\theta$ [elt deg]')
    ylabel('[Vs]')
    set(gca,'XLim',[0 360],'XTick',0:60:360);
    set(gca,'YLim',fMax*1.1*[-1 1]);
    plot([0 360],[0 0],'-k','LineWidth',1,'HandleVisibility','off')
    set(hax(ii),'ColorOrderIndex',1)
    switch ii
        case 1
            set(hfig(ii),'FileName',[pathname 'MagnetizingVSdq.fig'])
        case 2
            set(hfig(ii),'FileName',[pathname 'MagnetizingVSabc.fig'])
    end
    hleg(ii) = legend(hax(ii),'show','Location','northeast');
end

plot(hax(1),th,fd,'-','DisplayName','$\lambda_d$')
plot(hax(1),th,fq,'-','DisplayName','$\lambda_q$')
plot(hax(1),th,abs(fd+j*fq),'--','DisplayName','$|\lambda_{dq}|$')
plot(hax(1),th,fMagFund,'-k','DisplayName','$\lambda_{magnetizing}^{fundamental}$')
plot(hax(1),th,fMagHarm,':k','DisplayName','$\lambda_{magnetizing}^{total}$')
plot(hax(1),th,fMagStat,'--k','DisplayName','$\lambda_{magnetizing}^{stator}$')
if exist('sets','var')
    for ii=1:length(sets)
        set(hax(1),'ColorOrderIndex',ii+3)
        plot(hax(1),th,sets(ii).fd,':','DisplayName',['$\lambda_{d' int2str(ii) '}$'])
        set(hax(1),'ColorOrderIndex',ii+3)
        plot(hax(1),th,sets(ii).fq,'--','DisplayName',['$\lambda_{q' int2str(ii) '}$'])
    end
end

plot(hax(2),th,fa,'-','DisplayName','$\lambda_a$')
plot(hax(2),th,fb,'-','DisplayName','$\lambda_b$')
plot(hax(2),th,fc,'-','DisplayName','$\lambda_c$')
plot(hax(2),th,fMagFund,'-k','DisplayName','$\lambda_{magnetizing}^{fundamental}$')
plot(hax(2),th,fMagHarm,':k','DisplayName','$\lambda_{magnetizing}^{total}$')
plot(hax(2),th,fMagStat,'--k','DisplayName','$\lambda_{magnetizing}^{stator}$')


save([pathname 'magnetizingFluxLinkage.mat'],'out','geo','per','mat','wf')
for ii=1:length(hfig)
    savePrintFigure(hfig(ii));
end

if nargout()==0
    clear out
end


