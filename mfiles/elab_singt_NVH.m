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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NVHsources] = elab_singt_NVH(geo,per,out,forceOut,NewDir)

figNames{1}  = 'NVH_WF_torque';
figNames{2}  = 'NVH_FFT_torque';
figNames{3}  = 'NVH_WF_forceR_tooth';
figNames{4}  = 'NVH_WF_forceT_tooth';
figNames{5}  = 'NVH_FFT_forceR_tooth';
figNames{6}  = 'NVH_FFT_forceT_tooth';

for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting();
    hax(ii) = axes('OuterPosition',[0 0 1 1]);
    set(gcf,'Name',figNames{ii},'FileName',[NewDir figNames{ii} '.fig'])
    switch ii
        case 1
            xlabel('$\theta_{rot}$ (elt deg)')
            ylabel('$T$ (Nm)')
            set(gca,'XLim',[0 360],'XTick',0:60:360)
            plot([0 360],[0 0],'-k','LineWidth',1,'HandleVisibility','off')
        case 2
            xlabel('elt harmonic order')
            ylabel('$T$ (Nm)')
            set(gca,'XLim',[-1 49],'XTick',[0 6:6:96])
        case 3
            xlabel('$\theta_{rot}$ (elt deg)')
            ylabel('$F_r$ (N)')
            set(gca,'XLim',[0 360],'XTick',0:360/(2*geo.p):360)
            plot([0 360],[0 0],'-k','LineWidth',1,'HandleVisibility','off')
        case 4
            xlabel('$\theta_{rot}$ (elt deg)')
            ylabel('$F_t$ (N)')
            set(gca,'XLim',[0 360],'XTick',0:360/(2*geo.p):360)
            plot([0 360],[0 0],'-k','LineWidth',1,'HandleVisibility','off')
        case 5
            xlabel('elt harmonic order')
            ylabel('$F_r$ (N)')
            set(gca,'XLim',[-1 49],'XTick',[0 6:6:96])
        case 6
            xlabel('elt harmonic order')
            ylabel('$F_t$ (N)')
            set(gca,'XLim',[-1 49],'XTick',[0 6:6:96])
    end
    set(hax(ii),'ColorOrderIndex',1)
end




% torque
thRot = out.SOL.th-out.SOL.th(1);
xdeg  = per.delta_sim_singt;
T     = out.SOL.T;
% time symmetry (half period --> full period)
if xdeg==360
    disp('Simulation done on 360 elt deg rotation. OK')
elseif xdeg==180
    disp('Simulation done on 180 elt deg rotation. Symmetry applied')
    thRot = [thRot thRot+180];
    T     = [T T];
else
    disp(['Simulation done on ' int2str(xdeg) ' elt deg rotation. Symmetry not possible']);
end

plot(hax(1),thRot,T,'-b');
a = fft(T);
Tcont = abs(a(1))/length(a);
Tharm = 2*abs(a(2:end))/length(a);
bar(hax(2),0,Tcont,'FaceColor','b','EdgeColor',[0 0 0.5])
bar(hax(2),1:1:length(Tharm),Tharm,'FaceColor','b','EdgeColor',[0 0 0.5])

% force
th = forceOut.thRot;
Fr = forceOut.FrTooth(1,:);
Ft = forceOut.FtTooth(1,:);
plot(hax(3),th,Fr,'-b')
plot(hax(4),th,Ft,'-b')
a = fft(Fr);
FRcont = abs(a(1))/length(a);
FRharm = 2*abs(a(2:end))/length(a);
bar(hax(5),0,FRcont,'FaceColor','b','EdgeColor',[0 0 0.5])
bar(hax(5),1:1:length(FRharm),FRharm,'FaceColor','b','EdgeColor',[0 0 0.5])
a = fft(Ft);
FTcont = abs(a(1))/length(a);
FTharm = 2*abs(a(2:end))/length(a);
bar(hax(6),0,FTcont,'FaceColor','b','EdgeColor',[0 0 0.5])
bar(hax(6),1:1:length(FTharm),FTharm,'FaceColor','b','EdgeColor',[0 0 0.5])







% output
NVHsources   = forceOut;
NVHsources.T = T;

save([NewDir 'NVHsources.mat'],'NVHsources')

for ii=1:length(hfig)
    savePrintFigure(hfig(ii));
end


if nargout()==0
    clear NVHsources
end

