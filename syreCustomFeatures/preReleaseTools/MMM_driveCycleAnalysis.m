% Copyright 2020
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [driveCycleResults] = MMM_driveCycleAnalysis(motorModel,setup)

if nargin()==1
    [filename,pathname] = uigetfile(checkPathSyntax([cd '\']),'Select a drive cycle file');
    tmp = load([pathname filename]);
    setup.n = tmp.n;
    setup.T = tmp.T;
    setup.t = tmp.t;
    setup.figFlag = 0;
end

% if setup.figFlag
resFolder = checkPathSyntax([motorModel.data.pathname motorModel.data.motorName '_results\MMM results\driveCycle_' datestr(now,30) '\']);

for ii=1:6
    hfig(ii) = figure();
    figSetting();
    hax(ii) = axes(...
        'OuterPosition',[0 0 1 1],...
        'XLim',[min(setup.t) max(setup.t)]);
    xlabel('$t$ [s]')
    switch ii
        case 1
            ylabel('$T$ [Nm]')
            set(hfig(ii),'FileName',[resFolder 'torqueVStime.fig'])
        case 2
            ylabel('$n$ [rpm]')
            set(hfig(ii),'FileName',[resFolder 'speedVStime.fig'])
        case 3
            ylabel('$P$ [W]')
            set(hfig(ii),'FileName',[resFolder 'powerVStime.fig'])
        case 4
            ylabel('$I_{ph}$ [Apk]')
            set(hfig(ii),'FileName',[resFolder 'currentVStime.fig'])
        case 5
            ylabel('$V_{ph}$ [Vpk]')
            set(hfig(ii),'FileName',[resFolder 'voltageVStime.fig'])
        case 6
            ylabel('$\lambda$ [Vs]')
            set(hfig(ii),'FileName',[resFolder 'fluxVStime.fig'])
    end
    hleg(ii) = legend(hax(ii),'show','Location','northeast');
end

plot(hax(1),setup.t,setup.T,'-r','DisplayName','reference')
hT = plot(hax(1),0,0,'-g','DisplayName','output');

plot(hax(2),setup.t,setup.n,'-r','DisplayName','reference')
hn = plot(hax(2),0,0,'-g','DisplayName','output');

hP = plot(hax(3),0,0,'DisplayName','output power');
hL = plot(hax(3),0,0,'-r','DisplayName','loss');

plot(hax(4),[min(setup.t) max(setup.t)],motorModel.data.Imax*[1,1],'-r','DisplayName','limit');
hI = plot(hax(4),0,0,'-g','DisplayName','actual');

plot(hax(5),[min(setup.t) max(setup.t)],motorModel.data.Vdc/sqrt(3)*[1,1],'-r','DisplayName','limit');
hV = plot(hax(5),0,0,'-g','DisplayName','actual');

hF = plot(hax(6),0,0,'-g','DisplayName','actual');


% end


driveCycleResults.n     = nan(size(setup.t));
driveCycleResults.t     = nan(size(setup.t));
driveCycleResults.T     = nan(size(setup.t));
driveCycleResults.Id    = nan(size(setup.t));
driveCycleResults.Iq    = nan(size(setup.t));
driveCycleResults.Fd    = nan(size(setup.t));
driveCycleResults.Fq    = nan(size(setup.t));
driveCycleResults.Tem   = nan(size(setup.t));
driveCycleResults.Vo    = nan(size(setup.t));
driveCycleResults.Io    = nan(size(setup.t));
driveCycleResults.Im    = nan(size(setup.t));
driveCycleResults.Iph   = nan(size(setup.t));
driveCycleResults.Vph   = nan(size(setup.t));
driveCycleResults.PF    = nan(size(setup.t));
driveCycleResults.P     = nan(size(setup.t));
driveCycleResults.Ploss = nan(size(setup.t));
driveCycleResults.Pjs   = nan(size(setup.t));
driveCycleResults.PjDC  = nan(size(setup.t));
driveCycleResults.PjAC  = nan(size(setup.t));
driveCycleResults.Pfe   = nan(size(setup.t));
driveCycleResults.Pfes  = nan(size(setup.t));
driveCycleResults.Pfer  = nan(size(setup.t));
driveCycleResults.Ppm   = nan(size(setup.t));
driveCycleResults.Pjr   = nan(size(setup.t));
driveCycleResults.Pmech = nan(size(setup.t));
driveCycleResults.Eo    = nan(size(setup.t));
driveCycleResults.Ife   = nan(size(setup.t));
driveCycleResults.slip  = nan(size(setup.t));
driveCycleResults.Ir    = nan(size(setup.t));
driveCycleResults.Rs    = nan(size(setup.t));
driveCycleResults.eff   = nan(size(setup.t));

disp('Drive cycle computation...')
fprintf(' %06.2f%%',0)

for ii=1:length(setup.t)
    tmp = calcTnPoint(motorModel,setup.T(ii),setup.n(ii));
    driveCycleResults.t(ii)     = setup.t(ii);
    driveCycleResults.n(ii)     = setup.n(ii);
    driveCycleResults.T(ii)     = tmp.T;
    driveCycleResults.Id(ii)    = tmp.Id;
    driveCycleResults.Iq(ii)    = tmp.Iq;
    driveCycleResults.Fd(ii)    = tmp.Fd;
    driveCycleResults.Fq(ii)    = tmp.Fq;
    driveCycleResults.Tem(ii)   = tmp.Tem;
    driveCycleResults.Vo(ii)    = tmp.Vo;
    driveCycleResults.Io(ii)    = tmp.Io;
    driveCycleResults.Im(ii)    = tmp.Im;
    driveCycleResults.Iph(ii)   = tmp.Iph;
    driveCycleResults.Vph(ii)   = tmp.Vph;
    driveCycleResults.PF(ii)    = tmp.PF;
    driveCycleResults.P(ii)     = tmp.P;
    driveCycleResults.Ploss(ii) = tmp.Ploss;
    driveCycleResults.Pjs(ii)   = tmp.Pjs;
    driveCycleResults.PjDC(ii)  = tmp.PjDC;
    driveCycleResults.PjAC(ii)  = tmp.PjAC;
    driveCycleResults.Pfe(ii)   = tmp.Pfe;
    driveCycleResults.Pfes(ii)  = tmp.Pfes;
    driveCycleResults.Pfer(ii)  = tmp.Pfer;
    driveCycleResults.Ppm(ii)   = tmp.Ppm;
    driveCycleResults.Pjr(ii)   = tmp.Pjr;
    driveCycleResults.Pmech(ii) = tmp.Pmech;
    driveCycleResults.Eo(ii)    = tmp.Eo;
    driveCycleResults.Ife(ii)   = tmp.Ife;
    driveCycleResults.slip(ii)  = tmp.slip;
    driveCycleResults.Ir(ii)    = tmp.Ir;
    driveCycleResults.Rs(ii)    = tmp.Rs;
    driveCycleResults.eff(ii)   = tmp.eff;
    
    if (setup.figFlag || ii==length(setup.t))
        set(hT,'XData',driveCycleResults.t,'YData',driveCycleResults.T);
        set(hn,'xData',driveCycleResults.t,'YData',driveCycleResults.n);
        set(hP,'xData',driveCycleResults.t,'YData',driveCycleResults.P);
        set(hL,'xData',driveCycleResults.t,'YData',driveCycleResults.Ploss);
        set(hI,'xData',driveCycleResults.t,'YData',abs(driveCycleResults.Iph));
        set(hV,'xData',driveCycleResults.t,'YData',abs(driveCycleResults.Vph));
        set(hF,'xData',driveCycleResults.t,'YData',abs(driveCycleResults.Fd+j*driveCycleResults.Fq));
        drawnow()
    end

%     disp(['Point ' int2str(ii)])
    fprintf('\b\b\b\b\b\b\b\b')
    fprintf(' %06.2f%%',ii/length(setup.t)*100)
end

disp(' ')
disp('Drive cycle computed!')

mkdir(resFolder)
save([resFolder 'driveCycleResults.mat'],'setup','motorModel','driveCycleResults');
for ii=1:length(hfig)
    savePrintFigure(hfig(ii));
end

if nargout()==0
    clear driveCycleResults;
end

