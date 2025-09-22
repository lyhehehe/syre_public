% Copyright 2024
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [SOL] = simulate_xdeg_COMSOL(geo,per,eval_type,pathname,filename,pathnameIn)

% Connessione a COMSOL e Apertura modello impostato
import com.comsol.model.*
import com.comsol.model.util.*
model = mphopen([pathname filename(1:end-4) '.mph']);

load(fullfile(pathnameIn, [filename(1:end-4), '.mat']));

% Output directory definition for CSV files
csv_dir = fullfile (pathnameIn, [filename(1:end-4) '_Comsol\']);
if ~exist(csv_dir, 'dir')
    mkdir(csv_dir);
end

% Creazione dei percorsi completi per i file CSV
torque_csv = fullfile(csv_dir, [filename(1:end-4), '_Torque.csv']);
flux_csv = fullfile(csv_dir, [filename(1:end-4), '_Flux.csv']);
current_csv = fullfile(csv_dir, [filename(1:end-4), '_Current.csv']);
Pfe_r_csv = fullfile(csv_dir, [filename(1:end-4), '_Pfe_r.csv']);
Pfe_s_csv = fullfile(csv_dir, [filename(1:end-4), '_Pfe_s.csv']);
PM_losses_csv = fullfile(csv_dir, [filename(1:end-4), '_PM_losses.csv']);
PjrBar_csv = fullfile(csv_dir, [filename(1:end-4), '_PjrBar.csv']);

eval_type = 'singt';

% ============== Verifica geometria del rotore ============== %

    % if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype')
    %     geo.axisType = 'PM';
    %     geo.th0 = geo.th0 - 90;
    % else
    %     geo.axisType = 'SR';
    %     geo.th0 = geo.th0 + 90;
    % end

% Impostazione: discretizzation potential magnetic vector  (linear = 1, quadratic = 2, cubic = 3)
model.component('comp1').physics('rmm').prop('ShapeProperty').set('order_magneticvectorpotential', 1);


%theta_e = per.delta_sim_singt*(per.nsim_singt-1)/per.nsim_singt;    % intervallo angolare elettrico [deg]
%theta = theta_e*pi/180/p;                                           % angolo meccanico [rad]
% ts = theta / w;                                                    % delta time steps [s]
% t = ts*(per.nsim_singt-1);                                         % time step [s]
p = geo.p;                                                           % paia poli macchina
w = per.EvalSpeed*pi/30;                                             % velocità di rotazione [rad/s]
t = (per.delta_sim_singt/(p*w))*pi/180;                                  % time x simulation                 *(per.nsim_singt-1)/per.nsim_singt; 
%ts = t/(per.nsim_singt-1);                                               % number x function RANGE
range = linspace(0,t,per.nsim_singt);
freq = w*p/2/pi;                                                     % frequenza di alimentazione [Hz]
theta_m_tot = w*t*180/pi;                                            % angolo meccanico di simulazione [deg]
theta_e_tot = 2*pi*freq*t*180/pi;                                    % angolo elettrico di simulazione [deg]

% nsim = per.nsim_singt;                                             % numero di posizioni di rotore simulate 
% xdeg = per.delta_sim_singt;                                        % angolo elettrico di simulazione [deg]
% gamma = dataSet.GammaPP;                                           % angolo vettore corrente Idq dataIn [deg]
% th0 = geo.th0;                                                     % offset [deg]

if mod(geo.ps, 2)==0
    N_sim_sectors = geo.p;                %number of mech. sectors simulated
else
    N_sim_sectors = geo.p*2;
end 

% ============== Definizione studio ============== %


model.param().set('t', '0 [s]');                                          % impostazione tempo 0 s per la simulazione 
model.param().set('l', [num2str(geo.l) ' [mm]']);                         % impostazione lunghezza assiale l 
model.study('std1').setGenIntermediatePlots(true);
% Time dependent
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('tlist', num2str(range));                                 %['range(0,' num2str(ts) ',' num2str(t) ')']);
model.study('std1').feature('time').set('plot', true);
% TimetoFrequency Losses
model.study('std1').create('emloss', 'TimeToFrequencyLosses');
model.study('std1').feature('emloss').set('endtmethod', 'userdef');
% w_loss = 3*(w/2/pi);
model.study('std1').feature('emloss').set('lossstarttime', '0');
model.study('std1').feature('emloss').set('lossendtime', num2str(t));
% % Studio Circuito                                     %if "CIR" is activated again, the following row have to be restored  
% model.study('std1').feature('stat').setSolveFor('/physics/cir', true);
% model.study('std1').feature('time').setSolveFor('/physics/cir', true);
% model.study('std1').feature('emloss').setSolveFor('/physics/cir', true);
model.study('std1').feature('emloss').setSolveFor('/physics/rmm', true);
model.study('std1').feature('emloss').set('fftendtime', num2str(t));      %è ridondante, non so dove va (per ora)

% ============== Solver Settings ============== %

model.sol().create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 100);
model.sol('sol1').feature('s1').feature('fc1').set('ntolfact', 1);
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 100);
model.sol('sol1').feature('s1').feature('fc1').set('ntolfact', 1);
model.sol('sol1').feature('s1').feature().remove('fcDef');
model.sol('sol1').create('su1', 'StoreSolution');
model.sol('sol1').create('st2', 'StudyStep');
model.sol('sol1').feature('st2').set('study', 'std1');
model.sol('sol1').feature('st2').set('studystep', 'time');
model.sol('sol1').create('v2', 'Variables');
model.sol('sol1').feature('v2').set('initmethod', 'sol');
model.sol('sol1').feature('v2').set('initsol', 'sol1');
model.sol('sol1').feature('v2').set('initsoluse', 'sol2');
model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
model.sol('sol1').feature('v2').set('notsol', 'sol1');
model.sol('sol1').feature('v2').set('notsoluse', 'sol2');
model.sol('sol1').feature('v2').set('control', 'time');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', num2str(range));    %['range(0,' num2str(ts) ',' num2str(t) ')']
model.sol('sol1').feature('t1').set('plot', 'on');
model.sol('sol1').feature('t1').set('plotgroup', 'Default');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('tout', 'tstepsclosest');
model.sol('sol1').feature('t1').set('atolglobalvaluemethod', 'factor');
model.sol('sol1').feature('t1').set('reacf', true);
model.sol('sol1').feature('t1').set('storeudot', true);
model.sol('sol1').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol1').feature('t1').set('endtimeinterpolation', true);
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 25);
model.sol('sol1').feature('t1').feature('fc1').set('ntolfact', 0.2);
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 25);
model.sol('sol1').feature('t1').feature('fc1').set('ntolfact', 0.2);
model.sol('sol1').feature('t1').feature().remove('fcDef');
model.sol('sol1').create('su2', 'StoreSolution');
model.sol('sol1').create('st3', 'StudyStep');
model.sol('sol1').feature('st3').set('study', 'std1');
model.sol('sol1').feature('st3').set('studystep', 'emloss');
model.sol('sol1').create('v3', 'Variables');
model.sol('sol1').feature('v3').set('initmethod', 'sol');
model.sol('sol1').feature('v3').set('initsol', 'sol1');
model.sol('sol1').feature('v3').set('initsoluse', 'sol3');
model.sol('sol1').feature('v3').set('notsolmethod', 'sol');
model.sol('sol1').feature('v3').set('notsol', 'sol1');
model.sol('sol1').feature('v3').set('control', 'emloss');
model.sol('sol1').create('fft1', 'FFT');
model.sol('sol1').feature('fft1').set('control', 'emloss');
model.study('std1').feature('emloss').set('fftstarttime', '0');
model.study('std1').feature('emloss').set('fftendtime', num2str(t));           %%%%%%%%%%%%%%%%%
model.study('std1').feature('emloss').set('fftmaxfreq', ['6/(' num2str(t) '-0)']);
model.sol('sol1').create('su3', 'StoreSolution');
model.sol('sol1').create('cms1', 'CombineSolution');
model.sol('sol1').feature('cms1').set('soloper', 'gensum');
model.sol('sol1').feature('cms1').setIndex('gensumexpressionactive', 'on', 7);
model.sol('sol1').feature('cms1').setIndex('gensumexpression', 'comp1.rmm.Qloss', 7);
model.sol('sol1').feature('cms1').set('control', 'emloss');
model.sol('sol1').feature('v3').set('notsolnum', 'auto');
model.sol('sol1').feature('v3').set('notsolvertype', 'solnum');
model.sol('sol1').attach('std1');
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', false);
model.sol('sol1').feature('t1').feature('fc1').set('dtech', 'auto');
model.sol('sol1').feature('t1').feature('fc1').set('rstep', 4);
                                                                                      %if "CIR" is activated again, the following row have to be restored 
% % Impostazione bobine nel circuit interface (NB: non spostare i comandi da qui!!)  
% model.component('comp1').physics('cir').feature('termI1').set('V_src', 'root.comp1.rmm.VCoil_1');
% model.component('comp1').physics('cir').feature('termI2').set('V_src', 'root.comp1.rmm.VCoil_2'); 
% model.component('comp1').physics('cir').feature('termI3').set('V_src', 'root.comp1.rmm.VCoil_3');

mphrun(model);                  % run intero modello (con barra di progresso)

% model.study('std1').run();      % run singolo studio (senza barra di progresso)

pause(0.1);

%% Post-processing (Creazione Dataset)

% Coppia Motrice (Arkkio)
model.result().numerical().create('gev1', 'EvalGlobal');
model.result().table().create('tbl1', 'Table');
model.result().export().create('tbl1', 'Table');

% Flussi concatenati dq
model.result().numerical().create('gev2', 'EvalGlobal');
model.result().table().create('tbl2', 'Table');
model.result().export().create('tbl2', 'Table');

% Correnti abc
model.result().numerical().create('gev3', 'EvalGlobal');
model.result().table().create('tbl3', 'Table');
model.result().export().create('tbl3', 'Table');

% Pfe rotore
model.result().numerical().create('int1', 'IntSurface');
model.result().table().create('tbl4', 'Table');
model.result().export().create('tbl4', 'Table');

% Pfe statore
model.result().numerical().create('int2', 'IntSurface');
model.result().table().create('tbl5', 'Table');
model.result().export().create('tbl5', 'Table');

% PM losses
model.result().numerical().create('int3', 'IntSurface');
model.result().table().create('tbl6', 'Table');
model.result().export().create('tbl6', 'Table');

% PM losses
model.result().numerical().create('int4', 'IntSurface');
model.result().table().create('tbl7', 'Table');
model.result().export().create('tbl7', 'Table');

%% Post-processing (Salvataggio Risultati)

% Coppia Motrice (Arkkio)
model.result().numerical('gev1').set('data', 'dset3');
model.result().numerical('gev1').set('expr', {'rmm.Tark_1'});
model.result().numerical('gev1').set('descr', {'Axial Torque'});
model.result().table('tbl1').comments('Global Evaluation 1');
model.result().numerical('gev1').set('table', 'tbl1');
model.result().numerical('gev1').setResult();
model.result().export('tbl1').set('filename', torque_csv);
model.result().export('tbl1').set('table', 'tbl1');
model.result().export('tbl1').run();

% Flussi concatenati dq
model.result().numerical('gev2').set('data', 'dset3');
% model.result().numerical('gev2').set('expr', {'rmm.PhiCoil_1', 'rmm.PhiCoil_2', 'rmm.PhiCoil_3'});  %%%% here- questa vale per il circuito
model.result().numerical('gev2').set('expr', {'rmm.PhiWinding_1_aPh1', 'rmm.PhiWinding_1_aPh2', 'rmm.PhiWinding_1_aPh3'});
model.result().numerical('gev2').set('descr', {'Lambda 1', 'Lambda 2', 'Lambda 3'});
model.result().table('tbl2').comments('Global Evaluation 2');
model.result().numerical('gev2').set('table', 'tbl2');
model.result().numerical('gev2').setResult();
model.result().export('tbl2').set('filename', flux_csv);
model.result().export('tbl2').set('table', 'tbl2');
model.result().export('tbl2').run();

% Correnti abc
model.result().numerical('gev3').set('data', 'dset3');
%model.result().numerical('gev3').set('expr', {'rmm.ICoil_1', 'rmm.ICoil_2', 'rmm.ICoil_3'});  
model.result().numerical('gev3').set('expr', {'rmm.IWinding_1_aPh1', 'rmm.IWinding_1_aPh2', 'rmm.IWinding_1_aPh3'}); 
model.result().numerical('gev3').set('descr', {'Coil current', 'Coil current', 'Coil current'});
model.result().table('tbl3').comments('Global Evaluation 3');
model.result().numerical('gev3').set('table', 'tbl3');
model.result().numerical('gev3').setResult();
model.result().export('tbl3').set('filename', current_csv);
model.result().export('tbl3').set('table', 'tbl3');
model.result().export('tbl3').run();

% Pfe rotore
model.result().numerical('int1').set('intvolume', true);
model.result().numerical('int1').set('expr', 'rmm.Qh');
model.result().numerical('int1').set('descr', 'Volumetric loss density, electromagnetic');
model.result().numerical('int1').setIndex('expr', 'rmm.Qh*l', 0);
model.result().numerical('int1').selection().set(geo.BC_fe_r);
model.result().table('tbl4').comments('Surface Integration 1');
model.result().numerical('int1').set('table', 'tbl4');
model.result().numerical('int1').setResult();
model.result().export('tbl4').set('filename', Pfe_r_csv);
model.result().export('tbl4').set('table', 'tbl4');
model.result().export('tbl4').run();

% Pfe statore
model.result().numerical('int2').set('intvolume', true);
model.result().numerical('int2').set('expr', 'rmm.Qh');
model.result().numerical('int2').set('descr', 'Volumetric loss density, electromagnetic');
model.result().numerical('int2').setIndex('expr', 'rmm.Qh*l', 0);
model.result().numerical('int2').selection().set(geo.BC_fe_s);
model.result().table('tbl5').comments('Surface Integration 1');
model.result().numerical('int2').set('table', 'tbl5');
model.result().numerical('int2').setResult();
model.result().export('tbl5').set('filename', Pfe_s_csv);
model.result().export('tbl5').set('table', 'tbl5');
model.result().export('tbl5').run();

% PM losses
model.result().numerical('int3').set('intvolume', true);
model.result().numerical('int3').set('expr', 'rmm.Qh');
model.result().numerical('int3').set('descr', 'Volumetric loss density, electromagnetic');
model.result().numerical('int3').setIndex('expr', 'rmm.Qh*l', 0);
model.result().numerical('int3').selection().set(geo.BC_pm);
model.result().table('tbl6').comments('Surface Integration 1');
model.result().numerical('int3').set('table', 'tbl6');
model.result().numerical('int3').setResult();
model.result().export('tbl6').set('filename', PM_losses_csv);
model.result().export('tbl6').set('table', 'tbl6');
model.result().export('tbl6').run();

% Pj barriere di flusso rotore
model.result().numerical('int4').set('intvolume', true);
model.result().numerical('int4').set('expr', 'rmm.Qh');
model.result().numerical('int4').set('descr', 'Volumetric loss density, electromagnetic');
model.result().numerical('int4').setIndex('expr', 'rmm.Qh*l', 0);
model.result().numerical('int4').selection().set(geo.Bar);
model.result().table('tbl7').comments('Surface Integration 1');
model.result().numerical('int4').set('table', 'tbl7');
model.result().numerical('int4').setResult();
model.result().export('tbl7').set('filename', PjrBar_csv);
model.result().export('tbl7').set('table', 'tbl7');
model.result().export('tbl7').run();

pause(0.1);

% Saving Solved Model
RatedCurrent = dataSet.RatedCurrent;
CurrLoPP = dataSet.CurrLoPP;
SimulatedCurrent = RatedCurrent*CurrLoPP;
GammaPP  = dataSet.GammaPP;
tempPP = dataSet.tempPP;

% Format current
iStr = num2str(SimulatedCurrent, 3);  % 3 significant digits
iStr = strrep(iStr, '.', 'A');        % Replace dot with 'A'

% Format gamma
gammaStr = num2str(GammaPP, 4);       % 4 significant digits
gammaStr = strrep(gammaStr, '.', 'd'); % Replace dot with 'd'
if ~contains(gammaStr, 'd')
    gammaStr = [gammaStr 'd'];
end

% Final folder name
folder_name = ['T_eval_' iStr '_' gammaStr '_' int2str(tempPP) 'deg'];

pathname_solved = fullfile(pathnameIn, [filename(1:end-4) '_results'], 'COMSOL', folder_name);
% Create COMSOL folder if it doesn’t exist
if ~exist(pathname_solved, 'dir')
    mkdir(pathname_solved);
end

mphsave(model, fullfile(pathname_solved, [filename(1:end-4) '_solved.mph']));

%% Post-processing (Impostazione struttura SOL)

% Coppia Motrice (Arkkio)
axial_torque = readmatrix(torque_csv, 'NumHeaderLines', 5);
tempo = [axial_torque(:, 1)]';
Tor = [axial_torque(:, 2)]';
Tor = repmat(Tor,1,360/per.delta_sim_singt);               %note: as electric simulation angle shpuld be used a real multiple of 360
% theta_elt = tempo*2*pi*freq+geo.th0*pi/180;%-pi/2;
% theta_elt_deg = theta_elt*180/pi;
theta_elt_deg = linspace(0,360,(per.nsim_singt)*360/per.delta_sim_singt);
theta_0 = geo.th0;

if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype') || strcmp(geo.RotType,'SPM-Halbach') || strcmp(geo.RotType, 'Seg')
    if geo.axisType == 'PM'
        theta_0 = theta_0-90;
    else
        theta_0 = (theta_0);
    end
elseif strcmp(geo.RotType,'Spoke-type')
    geo.axisType = 'PM';
    theta_0 = theta_0;
else
    geo.axisType = 'SR';
    theta_0 = theta_0;
end

% theta_abc2dq = (linspace(0,360-per.delta_sim_singt/per.nsim_singt,(per.nsim_singt)*360/per.delta_sim_singt) + (geo.th0-90) )*pi/180;
theta_abc2dq = (linspace(0,360-per.delta_sim_singt/per.nsim_singt,(per.nsim_singt)*360/per.delta_sim_singt) + (theta_0-90) )*pi/180;

% Flussi concatenati 
concatenated_flux = readmatrix(flux_csv, 'NumHeaderLines', 5);                                          
lambda_1 = [concatenated_flux(:, 2)]';                                                                  
lambda_2 = [concatenated_flux(:, 3)]';                                                                  
lambda_3 = [concatenated_flux(:, 4)]';                                                                  
lambda_360 = phaseQuantityDecoding(lambda_1,lambda_2,lambda_3,per.delta_sim_singt);                     
lambda_1 = lambda_360.a;
lambda_2 = lambda_360.b;
lambda_3 = lambda_360.c;

lambda_dq = abc2dq(lambda_1, lambda_2, lambda_3, theta_abc2dq);    %theta_elt);
lambda_d = lambda_dq(1, :);
lambda_q = lambda_dq(2, :);

lambda_mod = abs(lambda_d + 1i*lambda_q);         % modulo flusso [Vs]
lambda_arg = angle(lambda_mod(end));              % angolo flusso [rad]

% Correnti
currents = readmatrix(current_csv, 'NumHeaderLines', 5);
%tempo_i = [currents(:, 1)]';
i_1 = [currents(:, 2)]'; 
i_2 = [currents(:, 3)]'; 
i_3 = [currents(:, 4)]';
i_360 = phaseQuantityDecoding(i_1,i_2,i_3,per.delta_sim_singt);
i_1 = i_360.a;
i_2 = i_360.b;
i_3 = i_360.c;
%theta_i_deg = tempo_i*2*pi*freq/pi*180+geo.th0-90;

i_dq = abc2dq(i_1, i_2, i_3, theta_abc2dq);%theta_i_deg*pi/180);
id_pp = i_dq(1, :);
iq_pp = i_dq(2, :);

% Pfe rotore
Pfe_rotore = readmatrix(Pfe_r_csv, 'NumHeaderLines', 5);
Pfer = Pfe_rotore(1: 1);

% Pfe statore
Pfe_statore = readmatrix(Pfe_s_csv, 'NumHeaderLines', 5);
Pfes = Pfe_statore(1: 1);

% PM losses
Ppm_mat = readmatrix(PM_losses_csv, 'NumHeaderLines', 5);
Ppm = Ppm_mat(1: 1);

% Pj Barriere di flusso rotore
PjrBar_mat = readmatrix(PjrBar_csv, 'NumHeaderLines', 5);
PjrBar = PjrBar_mat(1: 1);

SOL.Pfes   = Pfes;
SOL.Pfer   = Pfer;
SOL.PjrBar = PjrBar;
SOL.Ppm    = Ppm;

SOL.th = theta_elt_deg;
SOL.id = id_pp;
SOL.iq = iq_pp;
SOL.fd = lambda_d*N_sim_sectors;
SOL.fq = lambda_q*N_sim_sectors;

SOL.T = Tor;

SOL.ia = i_1;
SOL.ib = i_2;
SOL.ic = i_3;

SOL.fa = lambda_1*N_sim_sectors;
SOL.fb = lambda_2*N_sim_sectors;
SOL.fc = lambda_3*N_sim_sectors;
