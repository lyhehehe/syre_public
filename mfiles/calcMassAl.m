% Copyright 2016
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

function Mass = calcMassAl(geo,mat)
% 
% Evaluate the rotor conductor mass (estimation)


% Data
rho = mat.BarCond.kgm3;  % [kg/m^3]

if strcmp(geo.RotType,'IM')
    lt    = geo.IM.lt;
    wt    = geo.IM.wt;
    ttd   = geo.IM.ttd;
    r     = geo.r;
    l     = geo.l;
    Nbars = geo.IM.Nbars;
    ws    = 2*pi*(r-lt/2)/Nbars-wt;
    Mass = rho*(lt-ttd)/1000*ws/1000*l/1000*Nbars;
elseif strcmp(geo.RotType, 'EESM')
    % rotor = geo.rotor;
    p = geo.p;
    l = geo.l;
    % pos = 1;
    % for ii = 1:length(rotor)
    %     if rotor(ii,9) == 2
    %         xR(pos) = rotor(ii,1);
    %         yR(pos) = rotor(ii,2);
    %         pos = pos + 1;
    %     end
    % end
    % xR(pos) = xR(1);
    % yR(pos) = yR(1);
    
    Acu  = geo.Acoilf;
    kcuf = geo.win.kcuf;
    Mass = round(rho*Acu/1e6*2*2*p*l/1e3*kcuf,2);
else
    Mass = 0;
end
