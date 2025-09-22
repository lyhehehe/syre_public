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

function eval_operatingPointCOMSOL(dataIn)

pathname = dataIn.currentpathname;
filename = dataIn.currentfilename;
% filemot = strrep(dataSet.currentfilename, '.mat');              %,'.mph');
load([dataIn.currentpathname dataIn.currentfilename]);

RatedCurrent = dataIn.RatedCurrent;
CurrLoPP = dataIn.CurrLoPP;
SimulatedCurrent = dataIn.SimulatedCurrent;
GammaPP  = dataIn.GammaPP;
BrPP = dataIn.BrPP;
NumOfRotPosPP = dataIn.NumOfRotPosPP;
AngularSpanPP = dataIn.AngularSpanPP;
NumGrid = dataIn.NumGrid;

per.EvalSpeed = dataIn.EvalSpeed;

clc;

overload_temp =  CurrLoPP;   % current to be simulated
gamma_temp = GammaPP;        % current phase angle
Br = BrPP;                   % remanence of all barriers magnets

eval_type = dataIn.EvalType;

per.overload=CurrLoPP;
per.i0 = RatedCurrent;
per.BrPP=BrPP;

per.nsim_singt = NumOfRotPosPP;       % # simulated positions
per.delta_sim_singt = AngularSpanPP;  % angular span of simulation

offset = zeros(1,length(CurrLoPP));

geometry = [];
output = [];
tempDirName = [];
% fileMotWithPath=[pathname filemot];

geo0 = geo;
mat0 = mat;
performance = per;
geoTmp = geo0;
perTmp = performance;
matTmp = mat0;

%[~,geometry,~,output,tempDirName] = 
 
[geo,mat,out] = COMSOLfitness([],geoTmp,perTmp,matTmp,eval_type,pathname,filename);  %,fileMotWithPath);


%%---------Post-processing (Creazione Grafici)----------%%


% Format current
iStr = num2str(RatedCurrent*CurrLoPP, 3);  % 3 significant digits
iStr = strrep(iStr, '.', 'A');        % Replace dot with 'A'

% Format gamma
gammaStr = num2str(GammaPP, 4);       % 4 significant digits
gammaStr = strrep(gammaStr, '.', 'd'); % Replace dot with 'd'

if ~contains(gammaStr, 'd')
    gammaStr = [gammaStr 'd'];
end

% Final folder name
folder_name = ['T_eval_' iStr '_' gammaStr '_' int2str(dataIn.tempPP) 'deg'];
pathname_solved = fullfile(pathname, [filename(1:end-4) '_results'], 'COMSOL', folder_name);

% Create COMSOL folder if it doesn’t exist
if ~exist(pathname_solved, 'dir')
    mkdir(pathname_solved);
end

% Plot Torque - IPF
figure();
set(gcf,'color','w');

subplot(2,1,1);
hold on;
plot((out.th - out.th(1)), out.T, 'b', 'LineWidth', 1.5);
mean_Tor = (mean(out.T));
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('[Nm]', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
title(['Mean Torque = ', num2str(mean_Tor), ' Nm'], 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
box on;

subplot(2,1,2);
hold on;
plot((out.th - out.th(1)), out.IPF, 'b', 'LineWidth', 1.5);
mean_IPF = mean(out.IPF);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('IPF', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
title(['Mean IPF = ', num2str(mean_IPF)], 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
box on;

% Save plot
%saveas(gcf, fullfile(pathname_solved, [filename(1:end-4) '_T_gamma.png']));
%exportgraphics(gcf, fullfile(pathname_solved, [filename(1:end-4) '_T_gamma.png']), 'Resolution', 300);
savefig(gcf, fullfile(pathname_solved, [filename(1:end-4) '_T_gamma.fig']));

% Plot concatenated flux dq
figure();
set(gcf,'color','w');

subplot(2,1,1);                    % Plot lambda_d
hold on;
plot((out.th - out.th(1)), out.fd, 'b', 'LineWidth', 1.5);
mean_lambda_d = mean(out.fd);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('$\lambda_d$ [Vs]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);      
title(['Mean $\lambda_d$ = ', num2str(mean_lambda_d), ' Vs'], 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
box on;

subplot(2,1,2);                    % Plot lambda_q
hold on;
plot((out.th - out.th(1)), out.fq, 'b', 'LineWidth', 1.5);
mean_lambda_q = mean(out.fq);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('$\lambda_q$ [Vs]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
title(['Mean $\lambda_q$ = ', num2str(mean_lambda_q), 'Vs'], 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
box on;

% Save plot
%saveas(gcf, fullfilepathname_solved, [filename(1:end-4) '_plot_flux.png']));
%exportgraphics(gcf, fullfile(pathname_solved, [filename(1:end-4) '_plot_flux.png']), 'Resolution', 300);
savefig(gcf, fullfile(pathname_solved, [filename(1:end-4) '_plot_flux.fig']));

% Plot concatenated flux abc - currents abc
figure();
set(gcf,'color','w');

subplot(2,1,1);                                    % Plot lambda abc
plot((out.th - out.th(1)), out.fa, 'b', 'LineWidth', 1.5);
hold on;
plot((out.th - out.th(1)), out.fb, 'r', 'LineWidth', 1.5);
hold on;
plot((out.th - out.th(1)), out.fc, 'g', 'LineWidth', 1.5);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('$\lambda_{abc}$ [Vs]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
title('Phase flux linkages', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12);
box on;

subplot(2,1,2);                                    % Plot i abc
plot((out.th - out.th(1)), out.ia, 'b', 'LineWidth', 1.5);
hold on;
plot((out.th - out.th(1)), out.ib, 'r', 'LineWidth', 1.5);
hold on;
plot((out.th - out.th(1)), out.ic, 'g', 'LineWidth', 1.5);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('$i_{abc}$ [A]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
title('Phase currents', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12);
box on;
legend({'$i_a$', '$i_b$', '$i_c$'}, 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);

% Save plot
%saveas(gcf, fullfile(pathname_solved, [filename(1:end-4) '_plot_phase.png']));
%exportgraphics(gcf, fullfile(pathname_solved, [filename(1:end-4) '_plot_phase.png']), 'Resolution', 300);
savefig(gcf, fullfile(pathname_solved, [filename(1:end-4) '_plot_phase.fig']));