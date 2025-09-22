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

function [AOA] = MMM_eval_AOA(motorModel,method)

%% Load data
fdfq = motorModel.FluxMap_dq;
axisType = motorModel.data.axisType;

kRefine = 5;
kRefine3D = 3;
if ~strcmp(motorModel.data.motorType,'EE')
    fdfq = mapsReInterpolation(fdfq,'Id','Iq',kRefine*size(fdfq.Id,1),'linear');
else
    fdfq = mapsReInterpolation3D(fdfq,'Id','Iq','Ir',kRefine3D*size(fdfq.Id,1),'linear');
end

Id   = fdfq.Id;
Iq   = fdfq.Iq;
Fd   = fdfq.Fd;
Fq   = fdfq.Fq;
T    = fdfq.T;
dTpp = fdfq.dTpp;
PF   = sin(atan2(Iq,Id)-atan2(Fq,Fd));


%% Extract curves

%MTPA
if ~strcmp(motorModel.data.motorType,'EE')
    [id,iq] = calcOptCtrl(Id,Iq,T,abs(Id+j*Iq),axisType);
else
    dim_3D = size(fdfq.Ir,3);
    max_opt = 501;
    id = NaN(dim_3D,max_opt);
    iq = NaN(dim_3D,max_opt);
    Ivect = linspace(0,max(abs(Id(:)+j*Iq(:))),max_opt);
    for ii = 1:dim_3D
        [id(ii,:),iq(ii,:)] = calcOptCtrl(Id(:,:,ii)',Iq(:,:,ii)',T(:,:,ii)',abs(Id(:,:,ii)'+j*Iq(:,:,ii)'),axisType,Ivect,1);
    end
end

% Check MTPA
if strcmp(motorModel.data.axisType,'SR')
    index = find(id==max(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    % index = find(id==min(Id,[],'all'));
    % id(index) = NaN;
    % iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq<0);
    id(index) = NaN;
    iq(index) = NaN;

else
    % index = find(id==max(Id,[],'all'));
    % id(index) = NaN;
    % iq(index) = NaN;
    index = find(id==min(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id>0);
    id(index) = NaN;
    iq(index) = NaN;
end

if ~strcmp(motorModel.data.motorType,'EE')
    index = ~isnan(id);
    id = id(index);
    iq = iq(index);
    
    index = ~isnan(iq);
    id = id(index);
    iq = iq(index);
    
    if ((id(1)~=0)&&(iq(1)~=0))
        id = [0 id];
        iq = [0 iq];
    end
end

MTPA.id   = id;
MTPA.iq   = iq;
if strcmp(motorModel.data.motorType,'EE')
    temp_ir   = fdfq.Ir(1,1,:);
    MTPA.ir_layers = temp_ir(:);
    MTPV.ir_layers = temp_ir(:);
end

% MTPA.fd   = interp2(Id,Iq,Fd,id,iq);
% MTPA.fq   = interp2(Id,Iq,Fq,id,iq);
% MTPA.T    = interp2(Id,Iq,T,id,iq);
% MTPA.dTpp = interp2(Id,Iq,dTpp,id,iq);

% MTPV
if ~strcmp(motorModel.data.motorType,'EE')
    [id,iq] = calcOptCtrl(Id,Iq,T,abs(Fd+j*Fq),axisType);
else
    dim_3D = size(fdfq.Ir,3);
    max_opt = 501;
    id = NaN(dim_3D,max_opt);
    iq = NaN(dim_3D,max_opt);
    Ivect = linspace(0,max(abs(Id(:)+j*Iq(:))),max_opt);
    for ii = 1:dim_3D
        [id(ii,:),iq(ii,:)] = calcOptCtrl(Id(:,:,ii)',Iq(:,:,ii)',T(:,:,ii)',abs(Fd(:,:,ii)'+j*Fq(:,:,ii)'),axisType,Ivect,1);
    end
end


% Check MTPV
if strcmp(motorModel.data.axisType,'SR')
    index = find(id==max(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id==min(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq<0);
    id(index) = NaN;
    iq(index) = NaN;
else
    index = find(id==max(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id==min(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id>0);
    id(index) = NaN;
    iq(index) = NaN;
end

if ~strcmp(motorModel.data.motorType,'EE')
    index = ~isnan(id);
    id = id(index);
    iq = iq(index);
    
    index = ~isnan(iq);
    id = id(index);
    iq = iq(index);
end

MTPV.id = id;
MTPV.iq = iq;

% MPFPA
if ~strcmp(motorModel.data.motorType,'EE')
    [id,iq] = calcOptCtrl(Id,Iq,PF,abs(Id+j*Iq),axisType);
else
    dim_3D = size(fdfq.Ir,3);
    max_opt = 501;
    id = NaN(dim_3D,max_opt);
    iq = NaN(dim_3D,max_opt);
    Ivect = linspace(0,max(abs(Id(:)+j*Iq(:))),max_opt);
    for ii = 1:dim_3D
        [id(ii,:),iq(ii,:)] = calcOptCtrl(Id(:,:,ii)',Iq(:,:,ii)',PF(:,:,ii)',abs(Id(:,:,ii)'+j*Iq(:,:,ii)'),axisType,Ivect,1);
    end
end


if strcmp(motorModel.data.axisType,'SR')
    index = find(id==max(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id==min(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
else
    index = find(id==max(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id==min(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
end

MPFPA.id = id;
MPFPA.iq = iq;

%% interpolate (if needed)

nPoints = 101;


if strcmp(method,'Fit')
    if strcmp(motorModel.data.motorType,'EE')
        warning('Fit option not implemented yet for EESM, proceeding with LUT instead!!!')
    else
        ft = fittype('poly8');
        if strcmp(axisType,'SR')
            opts = fitoptions(...
                'Method','LinearLeastSquares',...
                'Lower',[zeros(1,8) 0],...
                'Upper',[inf*ones(1,8),0]);
            [xData,yData] = prepareCurveData(MTPA.id,MTPA.iq);
            fitFun = fit(xData,yData,ft,opts);
            MTPA.id = linspace(0,max(MTPA.id),nPoints)';
            MTPA.iq = fitFun(MTPA.id);
            if ~isempty(MTPV.id)
                opts = fitoptions(...
                    'Method','LinearLeastSquares',...
                    'Lower',[zeros(1,8)],...
                    'Upper',[inf*ones(1,8)]);
                [xData,yData] = prepareCurveData(MTPV.id,MTPV.iq);
                fitFun = fit(xData,yData,ft,opts);
                MTPV.id = linspace(0,max(MTPV.id),nPoints)';
                MTPV.iq = fitFun(MTPV.id);
            end
        else
            opts = fitoptions(...
                'Method','LinearLeastSquares',...
                'Upper',[zeros(1,8) 0],...
                'Lower',[-inf*ones(1,8),0]);
            [xData,yData] = prepareCurveData(MTPA.iq,MTPA.id);
            fitFun = fit(xData,yData,ft,opts);
            MTPA.iq = linspace(0,max(MTPA.iq),nPoints)';
            MTPA.id = fitFun(MTPA.iq);
            if ~isempty(MTPV.iq)
                opts = fitoptions(...
                    'Method','LinearLeastSquares',...
                    'Upper',[zeros(1,8)],...
                    'Lower',[-inf*ones(1,8)]);
                [xData,yData] = prepareCurveData(MTPV.iq,MTPV.id);
                fitFun = fit(xData,yData,ft,opts);
                MTPV.iq = linspace(0,max(MTPV.iq),nPoints)';
                MTPV.id = fitFun(MTPV.iq);
            end
        end
        AOA.method = 'Fit';
        
        % filt for trajectories outside the domain
        index = find(MTPA.id>max(Id,[],'all'));
        MTPA.id(index) = NaN;
        MTPA.iq(index) = NaN;
        index = find(MTPA.id<min(Id,[],'all'));
        MTPA.id(index) = NaN;
        MTPA.iq(index) = NaN;
        index = find(MTPA.iq>max(Iq,[],'all'));
        MTPA.id(index) = NaN;
        MTPA.iq(index) = NaN;
        index = find(MTPA.iq<min(Iq,[],'all'));
        MTPA.id(index) = NaN;
        MTPA.iq(index) = NaN;
    
        index = ~isnan(MTPA.id);
        MTPA.id = MTPA.id(index);
        MTPA.iq = MTPA.iq(index);
        
        index = find(MTPV.id>max(Id,[],'all'));
        MTPV.id(index) = NaN;
        MTPV.iq(index) = NaN;
        index = find(MTPV.id<min(Id,[],'all'));
        MTPV.id(index) = NaN;
        MTPV.iq(index) = NaN;
        index = find(MTPV.iq>max(Iq,[],'all'));
        MTPV.id(index) = NaN;
        MTPV.iq(index) = NaN;
        index = find(MTPV.iq<min(Iq,[],'all'));
        MTPV.id(index) = NaN;
        MTPV.iq(index) = NaN;
    
        index = ~isnan(MTPV.id);
        MTPV.id = MTPV.id(index);
        MTPV.iq = MTPV.iq(index);
    end
else
    AOA.method = 'LUT';
end


%% Interp flux linkages and torque
if ~strcmp(motorModel.data.motorType,'EE')
    MTPA.fd    = interp2(Id,Iq,Fd,MTPA.id,MTPA.iq);
    MTPA.fq    = interp2(Id,Iq,Fq,MTPA.id,MTPA.iq);
    MTPA.T     = interp2(Id,Iq,T,MTPA.id,MTPA.iq);
    MTPA.dTpp  = interp2(Id,Iq,dTpp,MTPA.id,MTPA.iq);
    
    MTPV.fd    = interp2(Id,Iq,Fd,MTPV.id,MTPV.iq);
    MTPV.fq    = interp2(Id,Iq,Fq,MTPV.id,MTPV.iq);
    MTPV.T     = interp2(Id,Iq,T,MTPV.id,MTPV.iq);
    MTPV.dTpp  = interp2(Id,Iq,dTpp,MTPV.id,MTPV.iq);
    
    MPFPA.fd   = interp2(Id,Iq,Fd,MPFPA.id,MPFPA.iq);
    MPFPA.fq   = interp2(Id,Iq,Fq,MPFPA.id,MPFPA.iq);
    MPFPA.T    = interp2(Id,Iq,T,MPFPA.id,MPFPA.iq);
    MPFPA.dTpp = interp2(Id,Iq,dTpp,MPFPA.id,MPFPA.iq);
else
    dim_3D = size(fdfq.Ir,3);
    for ii=1:dim_3D
        Id_vect = Id(:,1,1);
        Iq_vect = Iq(1,:,1);
        [Id_mg, Iq_mg] = meshgrid(Id_vect(:), Iq_vect(:));

        MTPA.fd(ii,:)    = interp2(Id_mg,Iq_mg,Fd(:,:,ii),MTPA.id(1,:),MTPA.iq(1,:));
        MTPA.fq(ii,:)    = interp2(Id_mg,Iq_mg,Fq(:,:,ii),MTPA.id(1,:),MTPA.iq(1,:));
        MTPA.T(ii,:)     = interp2(Id_mg,Iq_mg,T(:,:,ii),MTPA.id(1,:),MTPA.iq(1,:));
        MTPA.dTpp(ii,:)  = interp2(Id_mg,Iq_mg,dTpp(:,:,ii),MTPA.id(1,:),MTPA.iq(1,:));
        
        MTPV.fd(ii,:)    = interp2(Id_mg,Iq_mg,Fd(:,:,ii),MTPV.id(1,:),MTPV.iq(1,:));
        MTPV.fq(ii,:)    = interp2(Id_mg,Iq_mg,Fq(:,:,ii),MTPV.id(1,:),MTPV.iq(1,:));
        MTPV.T(ii,:)     = interp2(Id_mg,Iq_mg,T(:,:,ii),MTPV.id(1,:),MTPV.iq(1,:));
        MTPV.dTpp(ii,:)  = interp2(Id_mg,Iq_mg,dTpp(:,:,ii),MTPV.id(1,:),MTPV.iq(1,:));
    
        MPFPA.fd(ii,:)   = interp2(Id_mg,Iq_mg,Fd(:,:,ii),MPFPA.id(1,:),MPFPA.iq(1,:));
        MPFPA.fq(ii,:)   = interp2(Id_mg,Iq_mg,Fq(:,:,ii),MPFPA.id(1,:),MPFPA.iq(1,:));
        MPFPA.T(ii,:)    = interp2(Id_mg,Iq_mg,T(:,:,ii),MPFPA.id(1,:),MPFPA.iq(1,:));
        MPFPA.dTpp(ii,:) = interp2(Id_mg,Iq_mg,dTpp(:,:,ii),MPFPA.id(1,:),MPFPA.iq(1,:));
    end
end

AOA.MTPA  = MTPA;
AOA.MTPV  = MTPV;
AOA.MPFPA = MPFPA;




















