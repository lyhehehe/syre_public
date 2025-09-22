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

function eval_vonMisesStress_COMSOL(dataIn)

pathname = dataIn.currentpathname;
filename = dataIn.currentfilename;
load([pathname filename]);

 % Storing the speed
evalSpeed = dataIn.EvalSpeed; 
if evalSpeed==0
    error('Please set a speed higher than zero!!!')
end

% Creating Comsol directory
[~, fileNameNoExt, ~] = fileparts(filename);
Comsol_dir = fullfile(dataIn.currentpathname, [fileNameNoExt '_Comsol']);
if ~isfolder(Comsol_dir)
    mkdir(Comsol_dir);
end

% Creating DXF for structural
geo.rot_mec = [];
geo.flux_barr = [];
geo.mag = [];

% Rotor on geometry (without magnets)
%geo.rot_mec = geo.rotor(geo.rotor(:, end-1) ~= 6, :);
geo.rot_mec = geo.rotor(ismember(geo.rotor(:, end-1), [1, 5]), :);
syreToDxf(NaN, geo.rot_mec, Comsol_dir, [fileNameNoExt '_rot_mec.dxf']);

% Barriers on geometry
geo.flux_barr = geo.rotor(geo.rotor(:, end-1) == 1, :);
syreToDxf(NaN, geo.flux_barr, Comsol_dir, [fileNameNoExt '_flux_barr.dxf']);

% Magnets on geometry
geo.mag = geo.rotor(geo.rotor(:, end-1) == 6, :);
syreToDxf(NaN, geo.mag, Comsol_dir, [fileNameNoExt '_magn.dxf']);

%Importing DXF files
rot_mec_dxf  = fullfile(Comsol_dir, [filename(1:end-4), '_rot_mec.dxf']);
magn_dxf = fullfile(Comsol_dir, [filename(1:end-4), '_magn.dxf']);
flux_barr_dxf = fullfile(Comsol_dir, [filename(1:end-4), '_flux_barr.dxf']);

% Creating csv directory
csv_dir = fullfile(pathname, [filename(1:end-4) '_Comsol']);
if ~exist(csv_dir, 'dir')
    mkdir(csv_dir);
end

% Creating csv path
Stress_csv = fullfile(csv_dir, [filename(1:end-4), '_Stress.csv']);

%COMSOL connection
import com.comsol.model.*
import com.comsol.model.util.*

% Opening of the model saved on the temporary folder
%model = mphopen([pathname filename(1:end-4), '.mph']);
model = mphopen([pathname filename(1:end-4) '_Comsol\' filename(1:end-4), '.mph']);

% Setting component2
model.component().create('comp2', true);
model.component('comp2').geom().create('geom2', 2);
model.component('comp2').mesh().create('mesh2');
model.component('comp2').geom('geom2').lengthUnit('mm');

% Setting geometry
model.component('comp2').geom('geom2').create('imp1', 'Import');
model.component('comp2').geom('geom2').feature('imp1').set('filename', rot_mec_dxf);
model.component('comp2').geom('geom2').run('imp1');
model.component('comp2').geom('geom2').create('imp2', 'Import');
model.component('comp2').geom('geom2').feature('imp2').set('filename', flux_barr_dxf);
model.component('comp2').geom('geom2').run('imp2');
model.component('comp2').geom('geom2').create('dif1', 'Difference');
model.component('comp2').geom('geom2').feature('dif1').selection('input').set('imp1');
model.component('comp2').geom('geom2').feature('dif1').selection('input2').set('imp2');
model.component('comp2').geom('geom2').run('dif1');
model.component('comp2').geom('geom2').create('imp3', 'Import');
model.component('comp2').geom('geom2').feature('imp3').set('filename', magn_dxf);
model.component('comp2').geom('geom2').run('imp3');
model.component('comp2').geom('geom2').feature('fin').set('action', 'assembly');
model.component('comp2').geom('geom2').run('fin');

% Setting selections
model.component('comp2').selection().create('disk3', 'Disk');
disk3 = model.component('comp2').selection('disk3');
disk3 = model.component('comp2').selection('disk3').set('entitydim', 2);

model.component('comp2').selection().create('disk4', 'Disk');
disk4 = model.component('comp2').selection('disk4');
disk4 = model.component('comp2').selection('disk4').set('entitydim', 1);
disk4 = model.component('comp2').selection('disk4').set('r', 0.1);

% Defining study
model.study().create('std2');
model.study('std2').create('stat', 'Stationary');
model.study('std2').feature('stat').setSolveFor('/physics/rmm', true);
model.study('std2').feature('stat').setEntry('activate', 'rmm', false);
model.component('comp2').physics().create('solid', 'SolidMechanics', 'geom2');
model.component('comp2').physics('solid').prop('d').set('d', 0.001);
model.study('std2').feature('stat').setSolveFor('/physics/solid', true);


% definizione materiali

%BH Curve Ferro
BH_curve = [mat.Stator.BH(:,2),mat.Stator.BH(:,1)];
%Bf_max = max(BH_curve(:,2));
BH_curve_cell = arrayfun(@(x, y) {num2str(x), num2str(y)}, BH_curve(:, 1), BH_curve(:, 2), 'UniformOutput', false);    % Converte la matrice BH_curve in una cella di stringhe con due colonne
BH_curve_table = vertcat(BH_curve_cell{:});                                                                            % Converte la cella di stringhe in una matrice di stringhe

% PM
PMName = mat.LayerMag.MatName;
Br = per.BrPP;
mu_PM = mat.LayerMag.mu;
rho_PM = mat.LayerMag.kgm3;
model.component('comp2').material().create('mat5', 'Common');
model.component('comp2').material('mat5').propertyGroup().create('RemanentFluxDensity', 'Remanent flux density');
model.component('comp2').material('mat5').label(PMName);
model.component('comp2').material('mat5').set('family', 'chrome');
model.component('comp2').material('mat5').propertyGroup('def').set('electricconductivity', {'1/1.4[uohm*m]', '0', '0', '0', '1/1.4[uohm*m]', '0', '0', '0', '1/1.4[uohm*m]'});
model.component('comp2').material('mat5').propertyGroup('def').set('relpermittivity', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp2').material('mat5').propertyGroup('RemanentFluxDensity').set('murec', {num2str(mu_PM), '0', '0', '0', num2str(mu_PM), '0', '0', '0', num2str(mu_PM)});
model.component('comp2').material('mat5').propertyGroup('RemanentFluxDensity').set('normBr', {num2str(Br)});
model.component('comp2').material('mat5').propertyGroup().create('Enu', 'Youngs_modulus_and_Poissons_ratio');
model.component('comp2').material('mat5').propertyGroup('Enu').set('E', '175e9'); % manual input
model.component('comp2').material('mat5').propertyGroup('Enu').set('nu', '0.24'); % manual input
model.component('comp2').material('mat5').propertyGroup('def').set('density', num2str(rho_PM));

% Silicon Iron
IronName = mat.Rotor.MatName;
rho_Iron = mat.Rotor.kgm3;
model.component('comp2').material().create('mat6', 'Common');
model.component('comp2').material('mat6').propertyGroup().create('BHCurve', 'B-H Curve');
model.component('comp2').material('mat6').propertyGroup('BHCurve').func().create('BH', 'Interpolation');
model.component('comp2').material('mat6').label(IronName);
model.component('comp2').material('mat6').propertyGroup('def').set('electricconductivity', {'1.851852[MS/m]', '0', '0', '0', '1.851852[MS/m]', '0', '0', '0', '1.851852[MS/m]'});
model.component('comp2').material('mat6').propertyGroup('def').set('relpermittivity', {'1[1]', '0', '0', '0', '1[1]', '0', '0', '0', '1[1]'});
model.component('comp2').material('mat6').propertyGroup('BHCurve').label('B-H Curve');
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').label('Interpolation 1');
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('table', BH_curve_table);
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('extrap', 'linear');
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('fununit', {'T'});
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('argunit', {'A/m'});
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('defineinv', true);
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('defineprimfun', true);
model.component('comp2').material('mat6').propertyGroup('BHCurve').set('normB', 'BH(normHin)');
model.component('comp2').material('mat6').propertyGroup('BHCurve').set('normH', 'BH_inv(normBin)');
model.component('comp2').material('mat6').propertyGroup('BHCurve').set('Wpm', 'BH_prim(normHin)');
model.component('comp2').material('mat6').propertyGroup('BHCurve').descr('normHin', 'Magnetic field norm');
model.component('comp2').material('mat6').propertyGroup('BHCurve').descr('normBin', 'Magnetic flux density norm');
model.component('comp2').material('mat6').propertyGroup('BHCurve').addInput('magneticfield');
model.component('comp2').material('mat6').propertyGroup('BHCurve').addInput('magneticfluxdensity');
model.component('comp2').material('mat6').propertyGroup().create('Enu', 'Youngs_modulus_and_Poissons_ratio');
model.component('comp2').material('mat6').propertyGroup('Enu').set('E', '200e9'); % manual input
model.component('comp2').material('mat6').propertyGroup('Enu').set('nu', '0.3'); % manual input
model.component('comp2').material('mat6').propertyGroup('def').set('density', rho_Iron);

% Material assignment
selNumber = [];
fer = [];
pm = [];

for kk = 1:size(geo.BLKLABELS.rotore.xy, 1)
    x = geo.BLKLABELS.rotore.xy(kk, 1);
    y = geo.BLKLABELS.rotore.xy(kk, 2);
    disk3.set('posx', x);
    disk3.set('posy', y);
    selNumber = disk3.entities();
    if geo.BLKLABELS.rotore.xy(kk, 3) == 5
       fer = horzcat(fer, selNumber);
    elseif geo.BLKLABELS.rotore.xy(kk, 3) == 6
       pm = horzcat(pm, selNumber);
    end
end

%model.component('comp2').material('mat6').selection().set(fer);
% Assign to COMSOL materials
model.component('comp2').material('mat6').selection().set(unique(fer));  % iron
model.component('comp2').material('mat5').selection().set(unique(pm));   % PM

% Assegnazione boundary conditions
% setting rotating frame
model.component('comp2').physics('solid').create('rotf1', 'RotatingFrame', 2);
model.component('comp2').physics('solid').feature('rotf1').set('SpinSoftening', false);
model.component('comp2').physics('solid').feature('rotf1').set('Ovm', evalSpeed/30*pi);


% setting sector symmetry
model.component('comp2').physics('solid').create('sym1', 'SymmetrySolid', 1);

simm_bound = [];

for kk = 1:size(geo.BLKLABELS.rotore.boundary, 1)
    x = geo.BLKLABELS.rotore.boundary(kk, 1);
    y = geo.BLKLABELS.rotore.boundary(kk, 2);
    disk4.set('posx', x);
    disk4.set('posy', y);
    selNumber = disk4.entities();
    simm_bound = horzcat(simm_bound, selNumber);
end

model.component('comp2').physics('solid').feature('sym1').selection().set(simm_bound);

selNumber_shaft = [];

raggio = geo.Ar;
alpha_mecc = pi/geo.p/2;
x_mecc_m = raggio*cos(alpha_mecc);
y_mecc_m = raggio*sin(alpha_mecc);
disk4.set('posx', x_mecc_m);
disk4.set('posy', y_mecc_m);
selNumber_shaft = disk4.entities();

model.component('comp2').physics('solid').create('fix1', 'Fixed', 1);
model.component('comp2').physics('solid').feature('fix1').selection().set(selNumber_shaft);

% mphrun(model)

model.study('std2').run();      % run singolo studio (senza barra di progresso)

%% ====== Post-processing ====== %%

% nel caso in cui si apra un file già risolto, bisogna cambiare i numeri
% dei plotgroup (pg), pg1 diventa pg5 e così via e verificare che i dataset
% siano quelli corretti

% Final folder name
folder_name = ['structural_' int2str(evalSpeed) 'rpm'];

pathname_solved = fullfile(pathname, [filename(1:end-4) '_results'], 'COMSOL', folder_name);

% Create COMSOL folder if it doesn’t exist
if ~exist(pathname_solved, 'dir')
    mkdir(pathname_solved);
end

% Width in cm
w    = 12;
% Height in cm
h   = 10;

% von Mises deformation plot
model.result().create('pg1', 'PlotGroup2D');
model.result('pg1').set('data', 'dset2');
model.result('pg1').set('defaultPlotID', 'stress');
model.result('pg1').label('Stress (solid)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').set('edges', false);
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('expr', {'solid.misesGp'});
model.result('pg1').feature('surf1').set('threshold', 'manual');
model.result('pg1').feature('surf1').set('thresholdvalue', 0.2);
model.result('pg1').feature('surf1').set('colortable', 'Rainbow');
model.result('pg1').feature('surf1').set('colortabletrans', 'none');
model.result('pg1').feature('surf1').set('colorscalemode', 'linear');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg1').feature('surf1').set('colortable', 'Prism');
model.result('pg1').feature('surf1').create('def1', 'Deform');
model.result('pg1').feature('surf1').feature('def1').set('scaleactive', true);
model.result('pg1').feature('surf1').feature('def1').set('scale', 100);
model.result('pg1').feature('surf1').feature('def1').set('expr', {'u', 'v'});
model.result('pg1').feature('surf1').feature('def1').set('descr', 'Displacement field');
model.result('pg1').feature('surf1').set('unit', 'MPa');
model.result('pg1').feature('surf1').set('colortable', 'RainbowLightClassic');

% Points exceeding the yield limit
model.result().export().create('plot1', 'pg1', 'surf1', 'Plot');
model.result().export('plot1').set('filename', Stress_csv);
model.result().export('plot1').run();

Stress = readmatrix(fullfile(Stress_csv));

[out] = eval_maxStress_COMSOL(Stress,mat);
% Count how many points exceeded the material yield limit
numOver = size(out.stressrot, 1);

if numOver > 0
    % Creating figure
    figure();
    figSetting();  % custom formatting
    axis equal;
    axis tight;
    set(gca, 'FontName', 'Times', 'FontSize', 12, ...
         'TickLabelInterpreter', 'latex', ...
         'GridLineStyle', ':');

    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [5 5 w h]);
    set(gcf, 'Color', [1 1 1]);
    
    xlabel('[mm]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
    ylabel('[mm]', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
    title('Von Mises Stress Limit Exceeded', ...
        'FontWeight', 'normal', 'FontName', 'Times', ...
        'FontSize', 12, 'Interpreter', 'latex');

    % Plot rotor geometry and overstressed points
    GUI_Plot_Machine(gca, geo.rotor(1:end-3, :));
    plot(out.stressrot(:,1), out.stressrot(:,2), 'ro', ...
         'MarkerSize', 4, 'LineWidth', 1.5);

    % Save figure
    savefig(gcf, fullfile(pathname_solved, ...
        [filename(1:end-4) '_stress_limit_exceeded.fig']));
    
    % Results message
    disp([num2str(numOver), ' nodes exceed the maximum material stress.']);
else
    disp('No nodes exceed the maximum material stress.');
end


% Creating figure
figDeform = figure();
figSetting()
mphplot(model, 'pg1', 'rangenum', 1);
title('von Mises Stress [MPa] - deformation scale = 100', 'Interpreter', 'latex', 'FontWeight', 'normal');

% Saving figure
savefig(figDeform, fullfile(pathname_solved, [filename(1:end-4) 'Deformation.fig']));
%-------------%

% Displacement plot
model.result().create('pg2', 'PlotGroup2D');
model.result('pg2').set('data', 'dset2');
model.result('pg2').label('Displacement');
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').set('edges', false);
model.result('pg2').create('surf1', 'Surface');
model.result('pg2').feature('surf1').create('def1', 'Deform');
model.result('pg2').feature('surf1').feature('def1').set('scaleactive', true);
model.result('pg2').feature('surf1').feature('def1').set('scale', 1);
model.result('pg2').feature('surf1').set('colortable', 'RainbowLightClassic');

% Creating figure
figDips = figure();
figSetting()
mphplot(model, 'pg2', 'rangenum', 1);
title('Displacement [mm]', 'Interpreter', 'latex', 'FontWeight', 'normal');


% Saving figure
savefig(figDips, fullfile(pathname_solved, [filename(1:end-4) 'Displacement.fig']));
%-------------%

% von Mises Stress plot
model.result().create("pg3", "PlotGroup2D");
model.result("pg3").set("data", "dset2");
model.result("pg3").create("surf1", "Surface");
model.result('pg3').set('edges', false);
model.result("pg3").feature("surf1").set("expr", "solid.misesGp");
model.result("pg3").feature("surf1").set("descr", "von Mises stress");
model.result("pg3").feature("surf1").set("unit", "MPa");
model.result("pg3").feature("surf1").set("colortable", "RainbowLightClassic");

% Creating figure
figStress = figure();
figSetting()
mphplot(model, 'pg3', 'rangenum', 1);
title('von Mises Stress [MPa]', 'Interpreter', 'latex', 'FontWeight', 'normal');

% Saving figure
savefig(figStress, fullfile(pathname_solved, [filename(1:end-4) 'von_Mises_stress.fig']));
%-------------%

%Saving solved model
mphsave(model, fullfile(pathname_solved, [filename(1:end-4) '_solved.mph']));
