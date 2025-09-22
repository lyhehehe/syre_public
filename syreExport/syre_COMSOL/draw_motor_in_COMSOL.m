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

function [geo,mat,dataSet] = draw_motor_in_COMSOL(geo,mat,pathIn,nameIn,dataSet)

% Connessione a COMSOL e Apertura modello default
import com.comsol.model.*
import com.comsol.model.util.*
model = mphopen('Test_auto.mph');

if ~isfile([pathIn nameIn(1:end-4),'_dxf'])
    button='Yes';
else
    button = questdlg('DXF existing. Replace it?','SELECT','Yes','No','Yes');
end

% if strcmp(button,'Yes')
%     syreToDxf(geo.stator, geo.rotor,pathname, filename);
% end

% Estrazione del nome del file senza estensione
[~, file_name, ~] = fileparts(nameIn);

%Caricamento dei file
load([pathIn nameIn(1:end-4),'.mat']);

comsol_stator = [geo.stator
                geo.r+geo.g/2 0 geo.r+geo.g 0 NaN NaN 0 2 max(geo.stator(:,end)+1);
                0 0 geo.r+geo.g/2 0 (geo.r+geo.g/2)*cos(pi/geo.p*geo.ps) (geo.r+geo.g/2)*sin(pi/geo.p*geo.ps) 1 2 max(geo.stator(:,end)+1);
                (geo.r+geo.g/2)*cos(pi/geo.p*geo.ps) (geo.r+geo.g/2)*sin(pi/geo.p*geo.ps) (geo.r+geo.g)*cos(pi/geo.p*geo.ps) (geo.r+geo.g)*sin(pi/geo.p*geo.ps) NaN NaN 0 2 max(geo.stator(:,end)+1)];
comsol_rotor = [geo.rotor
                geo.r 0 geo.r+geo.g/2 0 NaN NaN 0 1 max(geo.rotor(:,end)+1);
                0 0 geo.r+geo.g/2 0 (geo.r+geo.g/2)*cos(pi/geo.p*geo.ps) (geo.r+geo.g/2)*sin(pi/geo.p*geo.ps) 1 1 max(geo.rotor(:,end)+1);
                geo.r*cos(pi/geo.p*geo.ps) geo.r*sin(pi/geo.p*geo.ps) (geo.r+geo.g/2)*cos(pi/geo.p*geo.ps) (geo.r+geo.g/2)*sin(pi/geo.p*geo.ps) NaN NaN 0 1 max(geo.rotor(:,end)+1)];
Comsol_dir = strcat(pathIn,strrep(nameIn,'.mph','_Comsol\'));
syreToDxf(comsol_stator,NaN,Comsol_dir,strcat(file_name,'_stat.dxf'));
syreToDxf(NaN,comsol_rotor,Comsol_dir,strcat(file_name,'_rot.dxf'));

% Caricamento dei file DXF
rot_dxf = fullfile(Comsol_dir, [file_name, '_rot.dxf']);
stat_dxf = fullfile(Comsol_dir, [file_name, '_stat.dxf']);

% ============== Import e Costruzione della Geometria ============== %

geom = model.component('comp1').geom('geom1');
model.component('comp1').geom('geom1').create('imp1', 'Import');
model.component('comp1').geom('geom1').feature('imp1').set('filename', rot_dxf); 
model.component('comp1').geom('geom1').feature('imp1').importData();
model.component('comp1').geom('geom1').create('imp2', 'Import');
model.component('comp1').geom('geom1').feature('imp2').set('filename', stat_dxf);
model.component('comp1').geom('geom1').feature('imp2').importData();
model.component('comp1').geom('geom1').feature('fin').set('action', 'assembly');
model.component('comp1').geom('geom1').run('fin');
model.component('comp1').geom('geom1').lengthUnit('mm');

% Implementazione lunghezza assiale
model.component('comp1').physics('rmm').prop('d').set('d', 'l');

% Implementazione Air Gap
%rotore
xre2_n = geo.r + geo.g/2;
yre2_n = 0;
xre3_n = (geo.r + geo.g/2)*cos(pi/geo.p*geo.ps);
yre3_n = (geo.r + geo.g/2)*sin(pi/geo.p*geo.ps);

xra2_n = geo.r;
yra2_n = 0;
xra3_n = geo.r*cos(pi/geo.p*geo.ps);
yra3_n = geo.r*sin(pi/geo.p*geo.ps);

index_el_r = max(geo.rotor(:,9)) + 1;
air_r = 1;

%statore
xse1_n = geo.r + geo.g;
yse1_n = 0;
xse2_n = (geo.r + geo.g)*cos(pi/geo.p*geo.ps);
yse2_n = (geo.r + geo.g)*sin(pi/geo.p*geo.ps);

xsi1_n = geo.r + geo.g/2;
ysi1_n = 0;
xsi2_n = (geo.r + geo.g/2)*cos(pi/geo.p*geo.ps);
ysi2_n = (geo.r + geo.g/2)*sin(pi/geo.p*geo.ps);

index_el_s = max(geo.stator(:,9)) + 1;
air_s = 2;

%assegnazione coordinate
geo_stator = [];
geo_rotor = [];
geo.stator_n = [xsi1_n, ysi1_n, xse1_n, yse1_n, NaN, NaN, 0, air_s, index_el_s; 
                0, 0, xse1_n, yse1_n, xse2_n, yse2_n, 1, air_s, index_el_s;
                xse2_n, yse2_n, xsi2_n, ysi2_n, NaN, NaN, 0, air_s, index_el_s
                ];
geo.rotor_n = [xra2_n, yra2_n, xre2_n, yre2_n, NaN, NaN, 0, air_r, index_el_r;
               0, 0, xre2_n, yre2_n, xre3_n, yre3_n, 1, air_r, index_el_r;
               xre3_n, yre3_n, xra3_n, yra3_n, NaN, NaN, 0, air_r, index_el_r;];

geo_stator = [geo.stator; geo.stator_n];
geo_rotor = [geo.rotor; geo.rotor_n];

% Ndomains = max(geo.stator(:,end))+(geo.Qs)+(geo.Qs-1)+max(geo.rotor(:,end));
% Ndomains_n = max(geo_stator(:,end))+(geo.Qs)+(geo.Qs-1)+max(geo_rotor(:,end));
Ndomains_n = geom.getNDomains();  

if strcmp(dataSet.TypeOfRotor,'EESM')
    Ndomains_n = Ndomains_n - 1;
end

%Ndomains_n = 55   x Thor but why?
% Correzione matrici boundaries
stat_boundary = [];
rot_boundary = [];

%rotore
x_bnd_r_1 = (xra2_n + xre2_n)/2;
y_bnd_r_1 = (yra2_n + yre2_n)/2;
x_bnd_r_2 = (xre3_n + xra3_n)/2;
y_bnd_r_2 = (yre3_n + yra3_n)/2;

%statore
x_bnd_s_1 = (xsi1_n + xse1_n)/2;
y_bnd_s_1 = (ysi1_n + yse1_n)/2;
x_bnd_s_2 = (xse2_n + xsi2_n)/2;
y_bnd_s_2 = (yse2_n + ysi2_n)/2;

%assegnazione coordinate
gm.s.boundary = [x_bnd_s_1, y_bnd_s_1, 10;
                 x_bnd_s_2, y_bnd_s_2, 10
                 ];
gm.r.boundary = [x_bnd_r_1, y_bnd_r_1, 10;
                 x_bnd_r_2, y_bnd_r_2, 10
                 ];

stat_boundary = [geo.BLKLABELS.statore.boundary; gm.s.boundary];
rot_boundary = [geo.BLKLABELS.rotore.boundary; gm.r.boundary];

% Correzione matrici materiali (BLKLABELS)
stat_mat = [];
rot_mat = [];

%definizione punto medio
xm_arc = geo.r*cos(pi/geo.p); 
ym_arc = geo.r*sin(pi/geo.p);

%rotore
xm_rot = xm_arc + geo.g/4;
ym_rot = ym_arc + geo.g/4;

%statore
xm_stat = xm_arc + (3*geo.g/4);
ym_stat = ym_arc + (3*geo.g/4);

stat_mat = [geo.BLKLABELS.statore.xy; xm_stat, ym_stat, 2, 1.6667, 1];
rot_mat = [geo.BLKLABELS.rotore.xy; xm_rot, ym_rot, 1, 1.6667, 1, 0, 0, 0];

% ============== Assegnazione dei Materiali ============== %

%BH Curve M270-35A
% BH_curve_invertita = mat.Stator.BH(:,:);    
% BH_curve = zeros(200,2);                    
% BH_curve(:,1) = BH_curve_invertita(:,2);
% BH_curve(:,2) = BH_curve_invertita(:,1);             %note: this code requires the BH_tab with inverted columns
BH_curve = [mat.Stator.BH(:,2),mat.Stator.BH(:,1)];
Bf_max = max(BH_curve(:,2));
BH_curve_cell = arrayfun(@(x, y) {num2str(x), num2str(y)}, BH_curve(:, 1), BH_curve(:, 2), 'UniformOutput', false);    % Converte la matrice BH_curve in una cella di stringhe con due colonne
BH_curve_table = vertcat(BH_curve_cell{:});                                                                            % Converte la cella di stringhe in una matrice di stringhe

% Air
model.component('comp1').material().create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup('def').func().create('eta', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func().create('Cp', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func().create('rho', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func().create('k', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func().create('cs', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func().create('an1', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func().create('an2', 'Analytic');
model.component('comp1').material('mat1').propertyGroup().create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat1').propertyGroup().create('NonlinearModel', 'Nonlinear model');
model.component('comp1').material('mat1').propertyGroup().create('idealGas', 'Ideal gas');
model.component('comp1').material('mat1').propertyGroup('idealGas').func().create('Cp_IG', 'Piecewise');
model.component('comp1').material('mat1').label('Air');
model.component('comp1').material('mat1').set('family', 'air');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0', '1600.0', '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0', '1600.0', '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/R_const[K*mol/J]/T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('args', {'pA', 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argunit', {'Pa', 'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('plotargs', {'pA', '101325', '101325'; 'T', '273.15', '293.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0', '1600.0', '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*R_const[K*mol/J]/0.02897*T)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('fununit', 'm/s');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('plotargs', {'T', '273.15', '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('expr', '-1/rho(pA,T)*d(rho(pA,T),T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('args', {'pA', 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('fununit', '1/K');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('argunit', {'Pa', 'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotargs', {'pA', '101325', '101325'; 'T', '273.15', '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('funcname', 'muB');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('expr', '0.6*eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotargs', {'T', '200', '1600'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', '');
model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', '');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(pA,T)', '0', '0', '0', 'alpha_p(pA,T)', '0', '0', '0', 'alpha_p(pA,T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '0.02897[kg/mol]');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', 'muB(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0[S/m]', '0', '0', '0', '0[S/m]', '0', '0', '0', '0[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(pA,T)');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T)', '0', '0', '0', 'k(T)', '0', '0', '0', 'k(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T)');
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('def').addInput('pressure');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat1').propertyGroup('NonlinearModel').set('BA', '(def.gamma+1)/2');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp_IG').label('Piecewise 2');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp_IG').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp_IG').set('pieces', {'200.0', '1600.0', '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp_IG').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp_IG').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('Rs', 'R_const/Mn');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('heatcapacity', 'Cp_IG(T)');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('molarmass', '0.02897');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('pressure');
model.component('comp1').material('mat1').materialType('nonSolid');
model.component('comp1').material('mat1').set('family', 'air');

% Silicon Iron (SPLIT STATOR AND ROTOR?)
IronName = mat.Stator.MatName;
model.component('comp1').material().create('mat2', 'Common');
model.component('comp1').material('mat2').propertyGroup().create('BHCurve', 'B-H Curve');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func().create('BH', 'Interpolation');
model.component('comp1').material('mat2').label(IronName);  
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'0'}); 
model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', {'1[1]', '0', '0', '0', '1[1]', '0', '0', '0', '1[1]'});
model.component('comp1').material('mat2').propertyGroup('BHCurve').label('B-H Curve');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').label('Interpolation 1');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('table', BH_curve_table);
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('extrap', 'linear');  %'const' instead of 'linear' - eh no, suggestions of chatgtp not always are good 
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('fununit', 'T');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('argunit', 'A/m');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('defineinv', true);
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('defineprimfun', true);
model.component('comp1').material('mat2').propertyGroup('BHCurve').set('normB', 'BH(normHin)');
model.component('comp1').material('mat2').propertyGroup('BHCurve').set('normH', 'BH_inv(normBin)');
model.component('comp1').material('mat2').propertyGroup('BHCurve').set('Wpm', 'BH_prim(normHin)');
model.component('comp1').material('mat2').propertyGroup('BHCurve').descr('normHin', 'Magnetic field norm');
model.component('comp1').material('mat2').propertyGroup('BHCurve').descr('normBin', 'Magnetic flux density norm');
model.component('comp1').material('mat2').propertyGroup('BHCurve').addInput('magneticfield');
model.component('comp1').material('mat2').propertyGroup('BHCurve').addInput('magneticfluxdensity');
model.component('comp1').material('mat2').set('family', 'plastic');

% Copper
density_copper = mat.SlotCond.kgm3;               %[kg/m^3]
model.component('comp1').material().create('mat3', 'Common');
model.component('comp1').material('mat3').propertyGroup().create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat3').propertyGroup().create('linzRes', 'Linearized resistivity');
model.component('comp1').material('mat3').label('Copper');
model.component('comp1').material('mat3').set('family', 'copper');
model.component('comp1').material('mat3').propertyGroup('def').set('relpermeability', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat3').propertyGroup('def').set('electricconductivity', {'5.998e7[S/m]', '0', '0', '0', '5.998e7[S/m]', '0', '0', '0', '5.998e7[S/m]'});
model.component('comp1').material('mat3').propertyGroup('def').set('heatcapacity', '385[J/(kg*K)]');
model.component('comp1').material('mat3').propertyGroup('def').set('relpermittivity', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat3').propertyGroup('def').set('emissivity', '0.5');
model.component('comp1').material('mat3').propertyGroup('def').set('density', {num2str(density_copper)} );
model.component('comp1').material('mat3').propertyGroup('def').set('thermalconductivity', {'400[W/(m*K)]', '0', '0', '0', '400[W/(m*K)]', '0', '0', '0', '400[W/(m*K)]'});
model.component('comp1').material('mat3').propertyGroup('Enu').set('E', '126e9[Pa]');
model.component('comp1').material('mat3').propertyGroup('Enu').set('nu', '0.34');
model.component('comp1').material('mat3').propertyGroup('linzRes').set('rho0', '1.667e-8[ohm*m]');
model.component('comp1').material('mat3').propertyGroup('linzRes').set('alpha', '3.862e-3[1/K]');
model.component('comp1').material('mat3').propertyGroup('linzRes').set('Tref', '293.15[K]');
model.component('comp1').material('mat3').propertyGroup('linzRes').addInput('temperature');
model.component('comp1').material('mat3').set('family', 'copper');

% N52 PM
Br = per.BrPP;
mu = mat.LayerMag.mu; 
model.component('comp1').material().create('mat4', 'Common');
model.component('comp1').material('mat4').propertyGroup().create('RemanentFluxDensity', 'Remanent flux density');
model.component('comp1').material('mat4').label('N52 (Sintered NdFeB)');
model.component('comp1').material('mat4').set('family', 'chrome');
model.component('comp1').material('mat4').propertyGroup('def').set('electricconductivity', {'1/1.4[uohm*m]', '0', '0', '0', '1/1.4[uohm*m]', '0', '0', '0', '1/1.4[uohm*m]'});
model.component('comp1').material('mat4').propertyGroup('def').set('relpermittivity', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat4').propertyGroup('RemanentFluxDensity').set('murec', {num2str(mu), '0', '0', '0', num2str(mu), '0', '0', '0', num2str(mu)});
model.component('comp1').material('mat4').propertyGroup('RemanentFluxDensity').set('normBr', {num2str(Br)});
g = model.component('comp1').material('mat4').propertyGroup('def');   
g.set('density','7500[kg/m^3]');
g.set('heatcapacity','400[J/(kg*K)]');
g.set('thermalconductivity',{'7[W/(m*K)]','0','0','0','7[W/(m*K)]','0','0','0','7[W/(m*K)]'});
model.component('comp1').material('mat4').set('family', 'chrome');

% Definizione selection inspector
model.component('comp1').selection().create('disk1', 'Disk');
model.component('comp1').selection().create('disk2', 'Disk');

% Assegnazione default Aria 
sel = 1:Ndomains_n;
model.component('comp1').material('mat1').selection().set(sel);
model.component('comp1').physics('rmm').feature('al1').label("Ampere's Law - Air");
model.component('comp1').physics('rmm').feature('al1').create('loss1', 'LossCalculation', 2);

% Vettori per assegnazione materiali
tmp = [];
tmp_rot = rot_mat(:, 1:3);
tmp_stat = stat_mat(:, 1:3);
tmp = [tmp_rot; tmp_stat]; 
disk1 = model.component('comp1').selection('disk1');

% Vettori per settare i domini
A = []; 
B = [];
C = [];

% Assegnazione materiali
for kk = 1:Ndomains_n
    x = tmp(kk,1);
    y = tmp(kk,2);
    disk1.set('posx', x);
    disk1.set('posy', y);
    selNumber = disk1.entities();
    switch tmp(kk, 3)
        case {3, 8}
            A = horzcat(A, selNumber);
        case {4, 5}
            B = horzcat(B, selNumber);
        case 6
            C = horzcat(C, selNumber);
    end
end

model.component('comp1').material('mat3').selection().set(A);
model.component('comp1').material('mat2').selection().set(B);
model.component('comp1').material('mat4').selection().set(C);

% ============== Fitting Steinmetz ============== %

% Calcolo tempo di simulazione

p = geo.p;                                  % paia poli macchina
w = per.EvalSpeed*pi/30;                    % velocità di rotazione [rad/s]
freq = w*p/2/pi;                            % frequenza di alimentazione [Hz]

ff = linspace(50, freq, 51);
Bf = linspace(0, Bf_max, 51);
kh = mat.Rotor.kh;
ke = mat.Rotor.ke;
alpha = mat.Rotor.alpha;
beta = mat.Rotor.beta;
pfe_f = kh .*(ff .^alpha) .*(Bf .^beta) +  ke .*(ff .^2) .*(Bf .^2);

[fitresult, gof] = createFit(ff, Bf, pfe_f);
KH = fitresult.KH;
ALPHA = fitresult.ALPHA;
BETA = fitresult.BETA;

% ============== Assegnazione Condizioni a Contorno ============== %

% Definizione boundary condition PM 
tmp_rot_righe = size(rot_mat, 1);
AM = [];
BC_pm = [];

for kk = 1:tmp_rot_righe            
    x = tmp(kk,1);
    y = tmp(kk,2);
    disk1.set('posx', x);
    disk1.set('posy', y);
    selNumber = disk1.entities();
    if rot_mat(kk, 3)==6
        BC_pm = horzcat(BC_pm, selNumber);
        xm = rot_mat(kk, 6);
        ym = rot_mat(kk, 7);
        zm = rot_mat(kk, 8);
        AM = [xm, ym, zm];
        alnumber_m = ['al_m' num2str(kk+1)];
        model.component('comp1').physics('rmm').create(alnumber_m, 'AmperesLaw', 2);
        model.component('comp1').physics('rmm').feature(alnumber_m).selection().set(selNumber);
        model.component('comp1').physics('rmm').feature(alnumber_m).set('ConstitutiveRelationBH', 'RemanentFluxDensity');
        model.component('comp1').physics('rmm').feature(alnumber_m).set('e_crel_BH_RemanentFluxDensity', AM);
        model.component('comp1').physics('rmm').feature(alnumber_m).create('loss1', 'LossCalculation', 2);
        % model.component("comp1").physics("rmm").feature("cmag2").selection().set(selNumber);
        % model.component("comp1").physics("rmm").feature("cmag2").feature("north1").selection().set(11, 19);
        % model.component("comp1").physics("rmm").feature("cmag2").feature("south1").selection().set(7, 15);
    end
end

% ============== Definizione conducting magnets ============== %  %not
%needed but useful for others pourposes (for example when the magnetization
%vector is unknow)
% model.component("comp1").physics("rmm").create("cmag2", "ConductingMagnet", 2);
% model.component("comp1").physics("rmm").feature().move("cmag2", 6);
% model.component("comp1").physics("rmm").feature("cmag2").label("Conducting Magnet");
% model.component("comp1").physics("rmm").feature("cmag2").set("sigma_mat", "userdef");
% model.component("comp1").physics("rmm").feature("cmag2").set("sigma_mat", "from_mat");
% model.component("comp1").physics("rmm").feature("cmag2").selection().set(5, 7);
% model.component("comp1").physics("rmm").feature("cmag2").feature("north1").selection().set(11, 19);
% model.component("comp1").physics("rmm").feature("cmag2").feature("south1").selection().set(7, 15);


% Memorizzazione domini barriere di flusso rotore             
Bar = [];

for kk = 1:tmp_rot_righe
    x = tmp(kk,1);
    y = tmp(kk,2);                           
    disk1.set('posx', x);
    disk1.set('posy', y);
    selNumber = disk1.entities();
    if rot_mat(kk, 3)==1
        Bar = horzcat(Bar, selNumber);
    end
end

% Definizione boundary condition su ferro di statore
BC_fe_s = [];

for kk = 1:size(tmp_stat, 1)
    x = tmp_stat(kk,1);
    y = tmp_stat(kk,2);
    disk1.set('posx', x);                         
    disk1.set('posy', y);                                
    selNumber = disk1.entities();                 
    if tmp_stat(kk, 3)==4                             
       BC_fe_s = horzcat(BC_fe_s, selNumber);          
    end                                                  
end                                                     

model.component('comp1').physics('rmm').create('al_fs', 'AmperesLaw', 2);
model.component('comp1').physics('rmm').feature('al_fs').selection().set(BC_fe_s);
model.component('comp1').physics('rmm').feature('al_fs').set('ConstitutiveRelationBH', 'BHCurve');
model.component('comp1').physics('rmm').feature('al_fs').create('loss_f', 'LossCalculation', 2);
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('LossModel', 'Steinmetz');
model.component('comp1').physics('rmm').feature('al_fs').label("Ampere's Law - Stator");
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('kh_steinmetz', KH);
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('alpha', ALPHA);
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('beta_steinmetz', BETA);

% Definizione boundary condition su ferro di rotore
BC_fe_r = [];

for kk = 1:size(tmp_rot, 1)
    x = tmp_rot(kk,1);
    y = tmp_rot(kk,2);
    disk1.set('posx', x);
    disk1.set('posy', y);
    selNumber = disk1.entities();
    if tmp_rot(kk,3)==5
       BC_fe_r = horzcat(BC_fe_r, selNumber);          %VA A DESTINAZIONE  |
    end                                                %                   |                    
end                                                    %                   |
                                                       %                   |
model.component('comp1').physics('rmm').create('al_fr', 'AmperesLaw', 2);% |
model.component('comp1').physics('rmm').feature('al_fr').selection().set(BC_fe_r);
model.component('comp1').physics('rmm').feature('al_fr').set('ConstitutiveRelationBH', 'BHCurve');
model.component('comp1').physics('rmm').feature('al_fr').create('loss_f', 'LossCalculation', 2);
model.component('comp1').physics('rmm').feature('al_fr').feature('loss_f').set('LossModel', 'Steinmetz');
model.component('comp1').physics('rmm').feature('al_fr').label("Ampere's Law - Rotor");
model.component('comp1').physics('rmm').feature('al_fr').feature('loss_f').set('kh_steinmetz', KH);
model.component('comp1').physics('rmm').feature('al_fr').feature('loss_f').set('alpha', ALPHA);
model.component('comp1').physics('rmm').feature('al_fr').feature('loss_f').set('beta_steinmetz', BETA);

% ============== Definizione multiphase_winding ============== % 
avv = geo.win.avv;
[num_righe_avv, num_colonne_avv] = size(avv);
coil1 = [];
coil2 = [];
coil3 = [];
coil_1 = [];
coil_2 = [];
coil_3 = [];
Ac = reshape(A,2,[]);            %matrice ordinata dei domini dei coils
                                 %note!: the domains are ordinated 
                                 %counterclockwise starting from x_axis
                                 %(from left to right)


% model.component('comp1').physics('rmm').create('coil1', 'Coil', 2);
% model.component('comp1').physics('rmm').feature('coil1').label('Phase 1');
% model.component('comp1').physics('rmm').feature('coil1').set('ConductorModel', 'Multi');
% model.component('comp1').physics('rmm').feature('coil1').set('coilGroup', true);
% model.component('comp1').physics('rmm').feature('coil1').set('CoilExcitation', 'CircuitCurrent');
% model.component('comp1').physics('rmm').feature('coil1').set('N', {num2str(geo.win.Nbob*2*geo.p)});
% model.component('comp1').physics('rmm').feature('coil1').set('AreaFrom', 'FillingFactor');
% model.component('comp1').physics('rmm').feature('coil1').set('FillingFactor', {num2str(geo.win.kcu)});
% model.component('comp1').physics('rmm').create('coil2', 'Coil', 2);
% model.component('comp1').physics('rmm').feature('coil2').label('Phase 2');
% model.component('comp1').physics('rmm').feature('coil2').set('ConductorModel', 'Multi');
% model.component('comp1').physics('rmm').feature('coil2').set('coilGroup', true);
% model.component('comp1').physics('rmm').feature('coil2').set('CoilExcitation', 'CircuitCurrent');
% model.component('comp1').physics('rmm').feature('coil2').set('N', {num2str(geo.win.Nbob*2*geo.p)});
% model.component('comp1').physics('rmm').feature('coil2').set('AreaFrom', 'FillingFactor');
% model.component('comp1').physics('rmm').feature('coil2').set('FillingFactor', {num2str(geo.win.kcu)});
% model.component('comp1').physics('rmm').create('coil3', 'Coil', 2);
% model.component('comp1').physics('rmm').feature('coil3').label('Phase 3');
% model.component('comp1').physics('rmm').feature('coil3').set('ConductorModel', 'Multi');
% model.component('comp1').physics('rmm').feature('coil3').set('coilGroup', true);
% model.component('comp1').physics('rmm').feature('coil3').set('CoilExcitation', 'CircuitCurrent');
% model.component('comp1').physics('rmm').feature('coil3').set('N', {num2str(geo.win.Nbob*2*geo.p)});
% model.component('comp1').physics('rmm').feature('coil3').set('AreaFrom', 'FillingFactor');
% model.component('comp1').physics('rmm').feature('coil3').set('FillingFactor', {num2str(geo.win.kcu)});
% 
% loop implemented to assign to the model the right configuration (designed in geo.win.avv)
for i = 1:num_righe_avv
    for c = 1:num_colonne_avv                                 
        switch avv(i, c)
            case 1
                coil1 = horzcat(coil1, Ac(i, c));
            case 3
                coil3 = horzcat(coil3, Ac(i, c));
            case 2
                coil2 = horzcat(coil2, Ac(i, c));
            case -1
                coil_1 = horzcat(coil_1, Ac(i, c));   %the "underscore" stays for "minus" 
            case -3
                coil_3 = horzcat(coil_3, Ac(i, c));
            case -2
                coil_2 = horzcat(coil_2, Ac(i, c));
            

        end
    end
end

% % Assegna i domini alle bobine     - FEATURE: COIL - NOT USED ANYMORE 
% model.component('comp1').physics('rmm').feature('coil1').selection().set(coil1);
% model.component('comp1').physics('rmm').feature('coil2').selection().set(coil2);
% model.component('comp1').physics('rmm').feature('coil3').selection().set(coil3);
% 
% % Assegna i domini alla bobina con corrente inversa 
% model.component('comp1').physics('rmm').feature('coil3').create('rcd1', 'ReverseCoilGroupDomain', 2);
% model.component('comp1').physics('rmm').feature('coil3').feature('rcd1').selection().set(coil3);
% model.component('comp1').physics('rmm').feature('coil3').feature('rcd1').label('Reverse Current Phase 3');
% 
% %Assegna calcolo perdite alle bobine
% model.component('comp1').physics('rmm').feature('coil1').create('loss1', 'LossCalculation', 2);
% model.component('comp1').physics('rmm').feature('coil2').create('loss1', 'LossCalculation', 2);
% model.component('comp1').physics('rmm').feature('coil3').create('loss1', 'LossCalculation', 2); 

gamma = dataSet.GammaPP;           % angle vector I wrt to d_axis [deg]
teta_0 = geo.th0;                  % angle between dq_plane and alphabeta_plane [deg]

% [~,phase1_offset] = calcKwTh0(geo);
% phase1_offset = phase1_offset+360/(6*geo.p*geo.q*geo.win.n3phase)/2*geo.p;    %first slot in 360/(6pq)/2 position

if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype') || strcmp(geo.RotType,'SPM-Halbach') || strcmp(geo.RotType, 'Seg')
    if geo.axisType == 'PM'
        teta_0 = teta_0-90;
    else
        teta_0 = teta_0;
    end
elseif strcmp(geo.RotType,'Spoke-type')
    geo.axisType = 'PM';
    teta_0 = teta_0;
else
    geo.axisType = 'SR';
    teta_0 = teta_0;
end

%p = geo.p;                                  % pole pairs                                                  -already defined                 
%w = per.EvalSpeed*pi/30;                    % rotation speed - mechanic[rad/s]                            -already defined 
%freq_t = w*p/2/pi;                          % power frequency (frequenzadi alimentazione) [Hz]            -already defined
Ipk = per.i0*per.overload;
q = geo.q;                                    % number of slots/pole/phase
N_poles = p*2;
N_slots_simulated = geo.Qs;
N_turns_per_slot = dataSet.TurnsInSeries/p/q/2;  %NOTE: "/2" is necessary to consider the physical separation of a single slot in two parts
N_simulated_sectors = 360/per.delta_sim_singt;
N_conductors = dataSet.SlotConductorNumber;
Slot_filling_factor = dataSet.SlotFillFactor;


model.component("comp1").physics("rmm").create("wnd1", "MultiphaseWinding", 2);
model.component("comp1").physics("rmm").feature("wnd1").label("Multiphase Winding - prova");
model.component("comp1").physics("rmm").feature("wnd1").set("Ipk", {num2str(Ipk)});                            
model.component("comp1").physics("rmm").feature("wnd1").set("alpha_i", {num2str((teta_0 + gamma)*pi/180)});  %rad!  
model.component("comp1").physics("rmm").feature("wnd1").set("freq_t", {num2str(freq)});                     
model.component("comp1").physics("rmm").feature("wnd1").set("WindingLayout", "automatic");
model.component("comp1").physics("rmm").feature("wnd1").set("NoPoles", {num2str(N_poles)});                   
model.component("comp1").physics("rmm").feature("wnd1").set("NoSlots", {num2str(N_slots_simulated)});         
model.component("comp1").physics("rmm").feature("wnd1").set("NoCoilsPerSlot", {num2str(N_conductors)});  
model.component("comp1").physics("rmm").feature("wnd1").set("N", {num2str(N_turns_per_slot)});                          
model.component("comp1").physics("rmm").feature("wnd1").set("sigmaCoil", {num2str(mat.SlotCond.sigma)});            
model.component("comp1").physics("rmm").feature("wnd1").set("AreaFrom", "FillingFactor");
model.component("comp1").physics("rmm").feature("wnd1").set("FillingFactor", {num2str(Slot_filling_factor)});        
model.component("comp1").physics("rmm").feature("wnd1").set("SectorSettingsType", "UserDefined");
model.component("comp1").physics("rmm").feature("wnd1").set("nsectors", {num2str(N_simulated_sectors)});                
model.component("comp1").physics("rmm").feature("wnd1").selection().set([coil1, coil2, coil3, coil_1, coil_2, coil_3]);
model.component("comp1").physics("rmm").feature("wnd1").create("aPh1", "Phase");
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh1").label("Automatic Phase 1");
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh1").create("rcd1", "ReversedCurrentDirection", 2);
model.component("comp1").physics("rmm").feature("wnd1").create("aPh2", "Phase");
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh2").label("Automatic Phase 2");
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh2").create("rcd1", "ReversedCurrentDirection", 2);
model.component("comp1").physics("rmm").feature("wnd1").create("aPh3", "Phase");
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh3").label("Automatic Phase 3");
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh3").create("rcd1", "ReversedCurrentDirection", 2);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh1").selection().set([coil1, coil_1]);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh1").feature("rcd1").selection().set(coil_1);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh1").active(true);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh1").feature("rcd1").active(true);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh2").selection().set([coil2, coil_2]);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh2").feature("rcd1").selection().set(coil_2);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh2").active(true);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh2").feature("rcd1").active(true);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh3").selection().set([coil3, coil_3]);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh3").feature("rcd1").selection().set(coil_3);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh3").active(true);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh3").feature("rcd1").active(true);
model.component("comp1").physics("rmm").feature("wnd1").feature().move("aPh1", 3);
model.component("comp1").physics("rmm").feature("wnd1").feature().move("aPh2", 3);
model.component("comp1").physics("rmm").feature("wnd1").feature().move("aPh3", 3);
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh1").set("alpha_o", "0[deg] ");            %check the order
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh2").set("alpha_o", "-120[deg]");
model.component("comp1").physics("rmm").feature("wnd1").feature("aPh3").set("alpha_o", "-240[deg]");
model.component("comp1").physics("rmm").feature("wnd1").create("loss1", "LossCalculation", 2);


% Vettori assegnazione periodic condition sui bordi
tmp = [];
tmp_rot = rot_boundary(:, 1:3);
tmp_stat = stat_boundary(:, 1:3);
tmp_righe_s = size(tmp_stat, 1);
tmp_righe_r = size(tmp_rot, 1);
AP_s = [];
AP_r = [];  

disk2 = model.component('comp1').selection('disk2').set('entitydim', 1);
disk2 = model.component('comp1').selection('disk2').set('r', 0.01);


% Definizione periodic condition statore (continuità/antiperiodicità)
% note: ps = number of simmetric sectors electrically equivalent
for kk = 1:tmp_righe_s
    x = tmp_stat(kk,1);
    y = tmp_stat(kk,2);
    disk2.set('posx', x);
    disk2.set('posy', y);
    selNumber = disk2.entities();       
    if tmp_stat(kk,3)==10
        AP_s = horzcat(AP_s, selNumber);
    end
end

model.component('comp1').physics('rmm').create('pc1', 'PeriodicCondition', 1);
model.component('comp1').physics('rmm').feature('pc1').selection().set(AP_s);
if mod(geo.ps, 2)==0
model.component('comp1').physics('rmm').feature('pc1').set('PeriodicType', 'Continuity');
else
model.component('comp1').physics('rmm').feature('pc1').set('PeriodicType', 'AntiPeriodicity');
end
model.component('comp1').physics('rmm').feature('pc1').label("Periodic Condition - Stator");

% Definizione periodic condition rotor (continuità/antiperiodicità)
for kk = 1:tmp_righe_r
    x = tmp_rot(kk,1);
    y = tmp_rot(kk,2);
    %disp([x, y])
    if isfinite(x) && isfinite(y)            %this loop prevents the crash caused by unreadble variable like "nah"
        disk2.set('posx', x);                %i ignore the reason that leads to the creation of such variable in the mat file (.geo)
        disk2.set('posy', y);                % with this loop there will be not problems anymore
    else
        %warning('Coordinata non valida: x=%g, y=%g. Salto questa selezione.', x, y)
        continue;
    end
    % disk2.set('posx', abs(x));
    % disk2.set('posy', abs(y));
    selNumber = disk2.entities();       
    if tmp_rot(kk,3)==10
        AP_r = horzcat(AP_r, selNumber);
    end
end

model.component('comp1').physics('rmm').create('pc2', 'PeriodicCondition', 1);
model.component('comp1').physics('rmm').feature('pc2').selection().set(AP_r);
if mod(geo.ps, 2)==0               
model.component('comp1').physics('rmm').feature('pc2').set('PeriodicType', 'Continuity');
else
model.component('comp1').physics('rmm').feature('pc2').set('PeriodicType', 'AntiPeriodicity');
end
model.component('comp1').physics('rmm').feature('pc2').label("Periodic Condition - Rotor");

% Definizione di Sector Symmetry su Air Gap                   

model.component('comp1').physics('rmm').create('ssc1', 'SectorSymmetry', 1);
model.component('comp1').physics('rmm').feature('ssc1').set('pairs', 'ap1');
if mod(geo.ps, 2)==0
model.component('comp1').physics('rmm').feature('ssc1').set('nsector', geo.p);
model.component('comp1').physics('rmm').feature('ssc1').set('PeriodicType', 'Continuity');
else
model.component('comp1').physics('rmm').feature('ssc1').set('nsector', geo.p*2);
model.component('comp1').physics('rmm').feature('ssc1').set('PeriodicType', 'AntiPeriodicity');
end
model.component('comp1').physics('rmm').feature('ssc1').set('constraintOptions', 'weakConstraints');

% Definizione Arkkio Torque Calculation
model.component('comp1').physics('rmm').create('ark1', 'ArkkioTorqueCalculation', 2);

% ============== Definizione Moving Mesh ============== %           
T = [];
tmp = rot_mat;

for kk = 1:tmp_rot_righe        
    x = tmp(kk,1);
    y = tmp(kk,2);
    disk1.set('posx', x);
    disk1.set('posy', y);   
    selNumber = disk1.entities();
    T = horzcat(T, selNumber);
end

model.component('comp1').common().create('rot1', 'RotatingDomain');
model.component('comp1').common('rot1').selection().set(T);
model.component('comp1').common('rot1').set('rotationType', 'rotationalVelocity');
model.component('comp1').common('rot1').set('rotationalVelocityExpression', 'constantAngularVelocity');
model.component('comp1').common('rot1').set('angularVelocity', w);



% ============== Definizione Circuito ============== % 
% NOTE: feature CIRCUIT - not used anymore - x others pourposes can be implemented 
%substituted by MULTIPHASE COILS

% % import com.comsol.model.*
% % import com.comsol.model.util.*
% % model = mphload('syreDefaultMotor_solved.mph');
% % model.physics('rmm').feature

% model.component('comp1').physics().create('cir', 'Circuit', 'geom1');
% model.component('comp1').physics('cir').create('I1', 'CurrentSourceCircuit', -1);
% model.component('comp1').physics('cir').create('I2', 'CurrentSourceCircuit', -1);
% model.component('comp1').physics('cir').create('I3', 'CurrentSourceCircuit', -1);
% model.component('comp1').physics('cir').create('termI1', 'ModelTerminalIV', -1);
% model.component('comp1').physics('cir').create('termI2', 'ModelTerminalIV', -1);
% model.component('comp1').physics('cir').create('termI3', 'ModelTerminalIV', -1);
% model.component('comp1').physics('cir').feature('I2').setIndex('Connections', 1, 0, 0);
% model.component('comp1').physics('cir').feature('I2').setIndex('Connections', 3, 1, 0);
% model.component('comp1').physics('cir').feature('I3').setIndex('Connections', 1, 0, 0);
% model.component('comp1').physics('cir').feature('I3').setIndex('Connections', 4, 1, 0);
% model.component('comp1').physics('cir').feature('termI1').set('Connections', 2);
% model.component('comp1').physics('cir').feature('termI2').set('Connections', 3);
% model.component('comp1').physics('cir').feature('termI3').set('Connections', 4);
% model.component('comp1').physics('cir').create('R1', 'Resistor', -1);
% model.component('comp1').physics('cir').feature('R1').set('R', '100000 [Ω]');
% model.component('comp1').physics('cir').feature('R1').setIndex('Connections', 1, 0, 0);
% model.component('comp1').physics('cir').feature('R1').setIndex('Connections', 0, 1, 0);
% model.component('comp1').physics('cir').feature('I1').set('sourceType', 'SineSource');
% model.component('comp1').physics('cir').feature('I2').set('sourceType', 'SineSource');
% model.component('comp1').physics('cir').feature('I3').set('sourceType', 'SineSource');
% model.component('comp1').physics('cir').feature('I1').set('value', [num2str(Imod) ' [A]']);
% model.component('comp1').physics('cir').feature('I1').set('freq', [num2str(freq) ' [Hz]']);
% model.component('comp1').physics('cir').feature('I1').set('phase', [num2str(theta_i)]);
% model.component('comp1').physics('cir').feature('I2').set('value', [num2str(Imod) ' [A]']);
% model.component('comp1').physics('cir').feature('I2').set('freq', [num2str(freq) ' [Hz]']);
% model.component('comp1').physics('cir').feature('I2').set('phase', [num2str(theta_i) ' - 120*pi/180']);
% model.component('comp1').physics('cir').feature('I3').set('value', [num2str(Imod) ' [A]']);
% model.component('comp1').physics('cir').feature('I3').set('freq', [num2str(freq) ' [Hz]']);
% model.component('comp1').physics('cir').feature('I3').set('phase', [num2str(theta_i) ' + 120*pi/180']);

% ============== Costruzione Mesh ============== %
model.component("comp1").mesh("mesh1").automatic(true);
model.component('comp1').mesh('mesh1').autoMeshSize(3);    % mesh size (1-10) 1 = extremely fine; 10 = extremely coarse 
model.component('comp1').mesh('mesh1').run();

% ============== Salvataggio Modello Inizzializzato ============== %

geo.BC_fe_s = BC_fe_s;
geo.BC_fe_r = BC_fe_r;
geo.BC_pm = BC_pm;
geo.Bar = Bar;

mphsave(model, [Comsol_dir nameIn(1:end-4), '.mph']);
