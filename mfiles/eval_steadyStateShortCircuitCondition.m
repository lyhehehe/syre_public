% Copyright 2025
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

function [outSSSC,resFolder] = eval_steadyStateShortCircuitCondition(dataIn)

pathname=dataIn.currentpathname;
filemot = strrep(dataIn.currentfilename,'.mat','.fem');

load([dataIn.currentpathname dataIn.currentfilename]);

RatedCurrent = dataIn.RatedCurrent;
CurrLoPP = dataIn.CurrLoPP;
%SimulatedCurrent = dataIn.SimulatedCurrent;
SimulatedCurrent = RatedCurrent*CurrLoPP;
GammaPP  = dataIn.GammaPP;
BrPP = dataIn.BrPP;
NumOfRotPosPP = dataIn.NumOfRotPosPP;
AngularSpanPP = dataIn.AngularSpanPP;
per.flag3phaseSet = dataIn.Active3PhaseSets;
n3ph = dataIn.Num3PhaseCircuit;
flag3phaseSet = per.flag3phaseSet;


per.EvalSpeed = dataIn.EvalSpeed;

clc;

if ~isfield(geo,'axisType')
    if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype')||strcmp(geo.RotType,'Spoke-type')
        geo.axisType = 'PM';
    else
        geo.axisType = 'SR';
    end
end

if ~strcmp(geo.axisType,dataIn.axisType)
    %geo.axisType = dataIn.axisType;
    if strcmp(dataIn.axisType,'PM')
        geo.th0 = geo.th0 - 90;
    else
        geo.th0 = geo.th0 + 90;
    end
end


eval_type = dataIn.EvalType;

per.overload=CurrLoPP;
per.i0 = RatedCurrent;
per.BrPP=BrPP;
per.gamma = GammaPP;


per.nsim_singt = NumOfRotPosPP;       % # simulated positions
per.delta_sim_singt = AngularSpanPP;  % angular span of simulation
per.if = dataIn.FieldCurrent;         % Field Current EESM

motname=[pathname filemot(1:end-4) '.fem'];

per.custom_act = 0;

if strcmp(dataSet.axisType,'PM')
    gammaFault=-180-5;
else
    gammaFault=90+5;
end



MaxIter = 20;
i0 = per.i0;
fTolPU = 0.1;

per.flag3phaseSet = ones(1,n3ph);

disp('Starting FEMM simulations...')

if sum(per.flag3phaseSet)==0
    per.overload = 0;
    [~,~,~,out,~] = FEMMfitness([],geo,per,mat,'singt',motname);
    fTol = fTolPU*abs(out.fd+j*out.fq);
else
    per.overload = CurrLoPP.*flag3phaseSet;
    [~,~,~,out,~] = FEMMfitness([],geo,per,mat,'singt',motname);
    fTol = nan(1,geo.win.n3phase);
    for ss=1:geo.win.n3phase
        if flag3phaseSet(ss)==0
            fTol(ss) = abs(out.sets(ss).fd+j*out.sets(ss).fq)*fTolPU;
        end
    end
end

disp(' Reference simulation done')

done=0;
ii=1;

idqIter = nan(MaxIter,geo.win.n3phase);
fdqIter = nan(MaxIter,geo.win.n3phase);
TIter   = nan(MaxIter,1);

while ~done
    if ii==1
        per.overload = CurrLoPP;
        idqIter(ii,:) = i0.*per.overload.*cosd(GammaPP)+j*i0.*per.overload.*sind(GammaPP);
    elseif ii==2
        idqIter(ii,:) = i0.*per.overload.*flag3phaseSet.*cosd(GammaPP)+j*i0.*per.overload.*flag3phaseSet.*sind(GammaPP);
    elseif ii==3
        for ff=1:n3ph
            if flag3phaseSet(ff)==0
                idqIter(ii,ff) = i0*cosd(gammaFault)+j*i0*sind(gammaFault);
            else
                idqIter(ii,ff) = idqIter(ii-1,ff);
            end
        end
    else
        for ff=1:n3ph
            if todo(ff)==1
                fd = real(fdqIter(2:ii-1,ff));
                fq = imag(fdqIter(2:ii-1,ff));
                id = real(idqIter(2:ii-1,ff));
                iq = imag(idqIter(2:ii-1,ff));
                [idOK,index] = unique(id);
                fdOK = fd(index);
                [iqOK,index] = unique(iq);
                fqOK = fq(index);
                idNew = interp1(fdOK,idOK,0,'linear','extrap');
                iqNew = interp1(fqOK,iqOK,0,'linear','extrap');
                idqIter(ii,ff) = idNew+j*iqNew;
            else
                idqIter(ii,ff) = idqIter(ii-1,ff);
            end
        end
    end

    per.overload = abs(idqIter(ii,:))/i0;
    per.gamma    = angle(idqIter(ii,:))*180/pi;

    [~,~,~,out,~] = FEMMfitness([],geo,per,mat,'singt',motname);

    todo = zeros(1,n3ph);
    for ff=1:n3ph
        fdqIter(ii,ff) = out.sets(ff).fd+j*out.sets(ff).fq;
        if abs(fdqIter(ii,ff))>fTol(ff)
            todo(ff)=1;
        end
    end
    TIter(ii) = out.T;

    disp([' Iteration ' int2str(ii) ': done'])

    if sum(todo)==0
        done = 1;
    else
        done = 0;
        ii = ii+1;
    end

    if ii>MaxIter
        done=1;
        ii=ii-1;
    end
end


disp('FEMM simulations done!!!')

flagSC = zeros(size(flag3phaseSet));
flagSC(flag3phaseSet==0) = 1;

outSSSC.idq     = idqIter(ii,:);
outSSSC.fdq     = fdqIter(ii,:);
outSSSC.idqIter = idqIter(1:ii,:);
outSSSC.fdqIter = fdqIter(1:ii,:);
outSSSC.TIter   = TIter(1:ii);
outSSSC.flagSC  = flagSC;
outSSSC.out     = out;


resFolder = checkPathSyntax([pathname filemot(1:end-4) '_results\FEA results\steadyStateShortCircuit_' datestr(now,30) '\']);
mkdir(resFolder)

save([resFolder 'outSSSC.mat'],'outSSSC','dataSet')

iMax = max(abs(outSSSC.idqIter(:)));
fMax = max(abs(outSSSC.fdqIter(:)));


figNames{1} = 'CurrentsDQ';
figNames{2} = 'FluxLinkagesDQ';

for ff=1:length(figNames)
    hfig(ff) = figure();
    figSetting(12,12)
    hax(ff) = axes('OuterPosition',[0 0 1 1]);
    set(hfig(ff),'Name',figNames{ff},'FileName',[resFolder figNames{ff} '.fig']);
    hleg(ff) = legend(hax(ff),'show','Location','southeast');
    switch ff
        case 1
            xlabel('$i_d$ (A)')
            ylabel('$i_q$ (A)')
            set(gca,'XLim',iMax*[-1 1],'YLim',iMax*[-1 1])
            plot(iMax*[-1 1],[0 0],'-k','LineWidth',1,'HandleVisibility','off');
            plot([0 0],iMax*[-1 1],'-k','LineWidth',1,'HandleVisibility','off');
        case 2
            xlabel('$\lambda_d$ (Vs)')
            ylabel('$\lambda_q$ (Vs)')
            set(gca,'XLim',fMax*[-1 1],'YLim',fMax*[-1 1])
            plot(fMax*[-1 1],[0 0],'-k','LineWidth',1,'HandleVisibility','off');
            plot([0 0],fMax*[-1 1],'-k','LineWidth',1,'HandleVisibility','off');
    end
    set(hax(ff),'ColorOrderIndex',1)
end

colors = get(hax(1),'ColorOrder');


for ii=1:n3ph
    plot(hax(1),real(outSSSC.idq(ii)),imag(outSSSC.idq(ii)),'x','Color',colors(ii,:),'DisplayName',['set ' int2str(ii) ' - final'])
    plot(hax(1),real(outSSSC.idqIter(1,ii)),imag(outSSSC.idqIter(1,ii)),'o','Color',colors(ii,:),'DisplayName',['set ' int2str(ii) ' - initial'])
    plot(hax(2),real(outSSSC.fdq(ii)),imag(outSSSC.fdq(ii)),'x','Color',colors(ii,:),'DisplayName',['set ' int2str(ii) ' - final'])
    plot(hax(2),real(outSSSC.fdqIter(1,ii)),imag(outSSSC.fdqIter(1,ii)),'o','Color',colors(ii,:),'DisplayName',['set ' int2str(ii) ' - initial'])
end

for ff=1:length(hfig)
    savePrintFigure(hfig(ff));
end

if nargout()==0
    clear outSSSC resFolder
end