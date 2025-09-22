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

function [fdfq,tempPM] = MMM_load_fdfq(data,p)

if isfield(data,'Id')
    fdfq.Id = data.Id;
    fdfq.Iq = data.Iq;
    fdfq.Fd = data.Fd;
    fdfq.Fq = data.Fq;
    if isfield(data,'T')
        fdfq.T = data.T;
        % check torque sign
        Tflux = fdfq.Fd.*fdfq.Iq-fdfq.Fq.*fdfq.Id;
        flagSign = sign(fdfq.T)~=sign(Tflux);
        lim = 0.4;
        if (sum(flagSign(:))>(lim*numel(flagSign)))
            disp('Torque sign corrected')
            fdfq.T = -fdfq.T;
        end
    else
        fdfq.T = 3/2*p*(fdfq.Fd.*fdfq.Iq-fdfq.Fq.*fdfq.Id);
    end
    if isfield(data,'dT')
        fdfq.dT = data.dT;
    else
        fdfq.dT = nan(size(fdfq.Id));
    end
    if isfield(data,'dTpp')
        fdfq.dTpp = data.dTpp;
    else
        fdfq.dTpp = nan(size(fdfq.Id));
    end
    if isfield(data,'per')
        tempPM = data.per.tempPP;
    else
        tempPM = NaN;
    end
    if isfield(data,'We')
        fdfq.We = data.We;
        fdfq.Wc = data.Wc;
    end
    
    if isfield(data,'IM')
        fdfq.IM = data.IM;
    end

    if length(size(data.Id)) == 3 % 3D Matrix in case of EESM topology
        if isfield(data,'Ir')
            fdfq.Ir = data.Ir;
        else
            fdfq.Ir = nan(size(fdfq.Id));
        end

        if isfield(data,'Fr')
            fdfq.Fr = data.Fr;
        else
            fdfq.Fr = nan(size(fdfq.Id));
        end
    end
%     % check axis style
%     if ~isfield(data.dataSet,'axisType')
%         if strcmp(data.dataSet.TypeOfRotor,'SPM')||strcmp(data.dataSet.TypeOfRotor,'Vtype')
%             data.dataSet.axisType = 'PM';
%         else
%             data.dataSet.axisType = 'SR';
%         end
%     end
%     if ~strcmp(data.dataSet.axisType,motorModel.data.axisType)
%         tmp.FluxMap_dq = fdfq;
%         if strcmp(dataSet.axisType,'SR')
%             tmp = MMM_sr2pm(tmp);
%         else
%             tmp = MMM_pm2sr(tmp);
%         end
%         fdfq = tmp.FluxMap_dq;
%     end
else
    fdfq = [];
    tempPM = [];
end






