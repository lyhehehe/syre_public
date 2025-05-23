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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MMM_save_dqtMap(motorModel)

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = checkPathSyntax([motName '_results\MMM results\' 'dqtMap Model - ' int2str(motorModel.data.tempPM) 'deg\']);

if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder])
end

Id     = motorModel.FluxMap_dq.Id;
Iq     = motorModel.FluxMap_dq.Iq;
Fd     = motorModel.FluxMap_dq.Fd;
Fq     = motorModel.FluxMap_dq.Fq;
T      = motorModel.FluxMap_dq.T;
dT     = motorModel.FluxMap_dq.dT;
dTpp   = motorModel.FluxMap_dq.dTpp;
dqtMap = motorModel.FluxMap_dqt;

per.tempPP = motorModel.data.tempPM;
dataSet.axisType = motorModel.data.axisType;


%save([pathname resFolder 'dqtMap.mat'],'dqtMap')
save([pathname resFolder 'fdfq_idiq_n256_dqt.mat'],'Id','Iq','Fd','Fq','T','dT','dTpp','dqtMap','per','dataSet')