function [normOut] = calc_normalizedPMSMparametrs(motorModel,figFlag)

%% calcola i parametri del modello di macchina semplificato (Ld,Lq,Fm) in punti di lavoro definiti: velocità base, velocità di massima potenza, limite UGO

if nargin()==0
    [filename,pathname,~] = uigetfile([cs '\*,mat'],'Select a motor model');
    load([pathname filename],'motorModel')
    figFlag = 1;
elseif nargin()==1
    figFlag = 1;
end


% translate to PM axis
if strcmp(motorModel.data.axisType,'SR')
    motorModel = MMM_sr2pm(motorModel);
end

% compute the inductance maps

motorModel.AppInductanceMap_dq = MMM_eval_appInductanceMap(motorModel);
motorModel.IncInductanceMap_dq = MMM_eval_inductanceMap(motorModel);

% compute the operating limit
Plim = OpLimEval(motorModel,motorModel.data.Imax,motorModel.data.Vdc);

% AOC and ASC points

fdfq = motorModel.FluxMap_dq;
AppInd = motorModel.AppInductanceMap_dq;
IncInd = motorModel.IncInductanceMap_dq;


fm = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,0,0);

% Tmax @ base speed
baseSpeed.n   = Plim.n_A;
baseSpeed.idq = Plim.id_A+j*Plim .iq_A;
baseSpeed.fdq = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,real(baseSpeed.idq),imag(baseSpeed.idq))+j*interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,real(baseSpeed.idq),imag(baseSpeed.idq));
baseSpeed.T   = Plim.T_A;

% UGO limit
ugoSpeed.n   = baseSpeed.n*abs(baseSpeed.fdq)/fm;
ugoSpeed.n   = motorModel.data.Vdc/sqrt(3)/fm/motorModel.data.p*30/pi;
ugoSpeed.idq = interp1(Plim.n,Plim.id,ugoSpeed.n)+j*interp1(Plim.n,Plim.iq,ugoSpeed.n);
ugoSpeed.fdq = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,real(ugoSpeed.idq),imag(ugoSpeed.idq))+j*interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,real(ugoSpeed.idq),imag(ugoSpeed.idq));
ugoSpeed.T   = interp1(Plim.n,Plim.T,ugoSpeed.n);

% Pmax
[~,index] = max(Plim.T.*Plim.n*pi/30);
pmaxSpeed.n   = Plim.n(index);
pmaxSpeed.idq = Plim.id(index)+j*Plim.iq(index);
pmaxSpeed.fdq = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,real(pmaxSpeed.idq),imag(pmaxSpeed.idq))+j*interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,real(pmaxSpeed.idq),imag(pmaxSpeed.idq));
pmaxSpeed.T   = Plim.T(index);

% CPSR
Ptmp = Plim.P-Plim.T_A*Plim.n_A*pi/30;
ntmp = Plim.n;
if Ptmp(2)<Ptmp(end)
    cpsrSpeed.n = Inf;
else
    [~,index] = max(Ptmp);
    ntmp2 = ntmp(index:end);
    Ptmp2 = Ptmp(index:end);
    cpsrSpeed.n = interp1(Ptmp2,ntmp2,0);
end

if isfinite(cpsrSpeed.n)
    cpsrSpeed.idq = interp1(Plim.n,Plim.id,cpsrSpeed.n)+j*interp1(Plim.n,Plim.iq,cpsrSpeed.n);
    cpsrSpeed.fdq = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,real(cpsrSpeed.idq),imag(cpsrSpeed.idq))+j*interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,real(cpsrSpeed.idq),imag(cpsrSpeed.idq));
    cpsrSpeed.T   = interp1(Plim.n,Plim.T,cpsrSpeed.n);
else
    cpsrSpeed.idq = NaN;
    cpsrSpeed.fdq = NaN;
    cpsrSpeed.T   = NaN;
end

% profiles
limit.n   = Plim.n;
limit.idq = Plim.id+j*Plim.iq;
limit.fdq = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,real(limit.idq),imag(limit.idq))+j*interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,real(limit.idq),imag(limit.idq));
limit.T   = Plim.T;


%% calcolo SSSC e TSC currents (ch e HWC)

iPlot = unique(fdfq.Id);
fPlot = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,iPlot,zeros(size(iPlot)));

ich = -interp1(fPlot,iPlot,0);
baseSpeed.iHWC = -interp1(fPlot,iPlot,-abs(baseSpeed.fdq));
ugoSpeed.iHWC  = -interp1(fPlot,iPlot,-abs(ugoSpeed.fdq));
pmaxSpeed.iHWC = -interp1(fPlot,iPlot,-abs(pmaxSpeed.fdq));
cpsrSpeed.iHWC = -interp1(fPlot,iPlot,-abs(cpsrSpeed.fdq));
limit.iHWC     = -interp1(fPlot,iPlot,-abs(limit.fdq));

%% calcolo induttanze apparenti ed incrementali

baseSpeed.Lapp.Ld = interp2(AppInd.Id,AppInd.Iq,AppInd.Ld,real(baseSpeed.idq),imag(baseSpeed.idq));
baseSpeed.Lapp.Lq = interp2(AppInd.Id,AppInd.Iq,AppInd.Lq,real(baseSpeed.idq),imag(baseSpeed.idq));
baseSpeed.Lapp.Fm = interp2(AppInd.Id,AppInd.Iq,AppInd.Fm,real(baseSpeed.idq),imag(baseSpeed.idq));
baseSpeed.Linc.Ld = interp2(IncInd.Id,IncInd.Iq,IncInd.Ldd,real(baseSpeed.idq),imag(baseSpeed.idq));
baseSpeed.Linc.Lq = interp2(IncInd.Id,IncInd.Iq,IncInd.Lqq,real(baseSpeed.idq),imag(baseSpeed.idq));

ugoSpeed.Lapp.Ld = interp2(AppInd.Id,AppInd.Iq,AppInd.Ld,real(ugoSpeed.idq),imag(ugoSpeed.idq));
ugoSpeed.Lapp.Lq = interp2(AppInd.Id,AppInd.Iq,AppInd.Lq,real(ugoSpeed.idq),imag(ugoSpeed.idq));
ugoSpeed.Lapp.Fm = interp2(AppInd.Id,AppInd.Iq,AppInd.Fm,real(ugoSpeed.idq),imag(ugoSpeed.idq));
ugoSpeed.Linc.Ld = interp2(IncInd.Id,IncInd.Iq,IncInd.Ldd,real(ugoSpeed.idq),imag(ugoSpeed.idq));
ugoSpeed.Linc.Lq = interp2(IncInd.Id,IncInd.Iq,IncInd.Lqq,real(ugoSpeed.idq),imag(ugoSpeed.idq));

pmaxSpeed.Lapp.Ld = interp2(AppInd.Id,AppInd.Iq,AppInd.Ld,real(pmaxSpeed.idq),imag(pmaxSpeed.idq));
pmaxSpeed.Lapp.Lq = interp2(AppInd.Id,AppInd.Iq,AppInd.Lq,real(pmaxSpeed.idq),imag(pmaxSpeed.idq));
pmaxSpeed.Lapp.Fm = interp2(AppInd.Id,AppInd.Iq,AppInd.Fm,real(pmaxSpeed.idq),imag(pmaxSpeed.idq));
pmaxSpeed.Linc.Ld = interp2(IncInd.Id,IncInd.Iq,IncInd.Ldd,real(pmaxSpeed.idq),imag(pmaxSpeed.idq));
pmaxSpeed.Linc.Lq = interp2(IncInd.Id,IncInd.Iq,IncInd.Lqq,real(pmaxSpeed.idq),imag(pmaxSpeed.idq));

cpsrSpeed.Lapp.Ld = interp2(AppInd.Id,AppInd.Iq,AppInd.Ld,real(cpsrSpeed.idq),imag(cpsrSpeed.idq));
cpsrSpeed.Lapp.Lq = interp2(AppInd.Id,AppInd.Iq,AppInd.Lq,real(cpsrSpeed.idq),imag(cpsrSpeed.idq));
cpsrSpeed.Lapp.Fm = interp2(AppInd.Id,AppInd.Iq,AppInd.Fm,real(cpsrSpeed.idq),imag(cpsrSpeed.idq));
cpsrSpeed.Linc.Ld = interp2(IncInd.Id,IncInd.Iq,IncInd.Ldd,real(cpsrSpeed.idq),imag(cpsrSpeed.idq));
cpsrSpeed.Linc.Lq = interp2(IncInd.Id,IncInd.Iq,IncInd.Lqq,real(cpsrSpeed.idq),imag(cpsrSpeed.idq));

limit.Lapp.Ld = interp2(AppInd.Id,AppInd.Iq,AppInd.Ld,real(limit.idq),imag(limit.idq));
limit.Lapp.Lq = interp2(AppInd.Id,AppInd.Iq,AppInd.Lq,real(limit.idq),imag(limit.idq));
limit.Lapp.Fm = interp2(AppInd.Id,AppInd.Iq,AppInd.Fm,real(limit.idq),imag(limit.idq));
limit.Linc.Ld = interp2(IncInd.Id,IncInd.Iq,IncInd.Ldd,real(limit.idq),imag(limit.idq));
limit.Linc.Lq = interp2(IncInd.Id,IncInd.Iq,IncInd.Lqq,real(limit.idq),imag(limit.idq));


%% output data

normOut.motorModel = motorModel;
normOut.Plim       = Plim;
normOut.PM.fm      = fm;
normOut.PM.ich     = ich;
normOut.PM.Lapp.Ld = interp2(AppInd.Id,AppInd.Iq,AppInd.Ld,-ich,0);
normOut.PM.Lapp.Lq = interp2(AppInd.Id,AppInd.Iq,AppInd.Lq,-ich,0);
normOut.PM.Lapp.Fm = interp2(AppInd.Id,AppInd.Iq,AppInd.Fm,-ich,0);
normOut.PM.Linc.Ld = interp2(IncInd.Id,IncInd.Iq,IncInd.Ldd,-ich,0);
normOut.PM.Linc.Lq = interp2(IncInd.Id,IncInd.Iq,IncInd.Lqq,-ich,0);
normOut.baseSpeed  = baseSpeed;
normOut.ugoSpeed   = ugoSpeed;
normOut.pmaxSpeed  = pmaxSpeed;
normOut.cpsrSpeed  = cpsrSpeed;
normOut.limit      = limit;

%% figures

if figFlag

    hfig(1) = figure();
    figSetting(16,10)
    hax(1) = axes('OuterPosition',[0 0 1 1]);
    xlabel('$n$ (rpm)')
    ylabel('$T$ (Nm)')
    hleg(1) = legend('show','Location','northeast');
    plot(limit.n,limit.T,'-b','DisplayName','curve')
    plot(baseSpeed.n,baseSpeed.T,'bo','DisplayName','base point')
    plot(ugoSpeed.n,ugoSpeed.T,'bd','DisplayName','UGO speed')
    plot(cpsrSpeed.n,cpsrSpeed.T,'bs','DisplayName','CPSR speed')

    colors = get(gca,'ColorOrder');

    hfig(2) = figure();
    figSetting(16,10)
    hax(2) = axes('OuterPosition',[0 0 1 1]);
    xlabel('$n$ (rpm)')
    ylabel('$P$ (W)')
    hleg(2) = legend('show','Location','northeast');
    plot(limit.n,limit.T.*limit.n*pi/30,'-b','DisplayName','curve')
    plot(baseSpeed.n,baseSpeed.T.*baseSpeed.n*pi/30,'bo','DisplayName','base point')
    plot(ugoSpeed.n,ugoSpeed.T.*ugoSpeed.n*pi/30,'bd','DisplayName','UGO speed')
    plot(cpsrSpeed.n,cpsrSpeed.T.*cpsrSpeed.n*pi/30,'bs','DisplayName','CPSR speed')

    hfig(3) = figure();
    figSetting(16,10)
    hax(3) = axes('OuterPosition',[0 0 1 1]);
    xlabel('$n$ (rpm)')
    hleg(3) = legend('show','Location','northeast');
    set(hleg(3),'NumColumns',2);
    yyaxis('left')
    set(gca,'YColor','k')
    set(gca,'YLim',[0 1],'YTick',0:0.1:1)
    ylabel('$\lambda_m/\lambda_{max}$')
    plot(baseSpeed.n*[1,1],[0 10],'-k','DisplayName','$n_b$')
    plot(ugoSpeed.n*[1,1],[0 10],'--k','DisplayName','$n_{UGO}$')
    plot(cpsrSpeed.n*[1,1],[0 10],':k','DisplayName','$n_{CPSR}$')

    plot(limit.n,limit.Lapp.Fm./abs(baseSpeed.fdq),'-','Color',colors(1,:),'DisplayName','$\frac{\lambda_m}{\lambda_{max}}$')
    yyaxis('right')
    set(gca,'YColor','k')
    set(gca,'YLim',[0 10],'YTick',0:1:10)
    ylabel('$\xi$')
    plot(limit.n,limit.Lapp.Lq./limit.Lapp.Ld,'-','Color',colors(2,:),'DisplayName','$\xi_{apparent}$')
    plot(limit.n,limit.Linc.Lq./limit.Linc.Ld,'-','Color',colors(3,:),'DisplayName','$\xi_{incremental}$')

end









