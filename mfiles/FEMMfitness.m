% Copyright 2014
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

function [cost,geo,mat,out,pathname] = FEMMfitness(RQ,geo,per,mat,eval_type,filenameIn)

% FEMMfitness
% runs FEMM simulation from file (existing machine) or file+RQ (MODE, RQ is the set of inputs)
% - creates a temp dir
% - moves or draw the machine into the temp dir
% - simulate_xdeg(eval_type)
% - calc out from SOL
% - saves the .mat file into the temp dir

[~,filename,ext] = fileparts(filenameIn);
filename = [filename ext]; % fem file name

[~,pathname]=createTempDir();


if ~isempty(RQ)

    % MODE optimization (RQ geometry)
    RQ=RQ';
    geo.pathname=pwd();

    %     options.iteration=0;
    %     options.currentgen=1;
    %     options.PopulationSize=1;

    if strcmp(eval_type,'MO_OA')
%         RQ % debug .. when syre crashes it is useful to have visibility of last RQ
    end

    [geo,gamma,mat] = interpretRQ(RQ,geo,mat);
    per.gamma=gamma;

    [geo,mat] = draw_motor_in_FEMM(geo,mat,pathname,filename);

    [~,geo] = calc_endTurnLength(geo);
    [~,geo] = calc_endTurnFieldLength(geo);
 
    %     flag_OptCurrConst = 1;
    switch per.flag_OptCurrConst
        case 0 % constant thermal loading
            per.loss = NaN;
            per.J    = NaN;
        case 1 % constant current density
            per.kj   = NaN;
            per.Loss = NaN;
        case 2 % constant current
            per.kj   = NaN;
            per.Loss = NaN;
            per.J    = NaN;
    end
    per = calc_i0(geo,per,mat);
    % warning('Define the ratio between stator and rotor current density')
    per.Jf = per.J*per.JfPU;
    per = calc_if(geo,per,mat);
    per.if = per.if0;

    %     if any(strcmp(geo.OBJnames,'Fdq0'))
    %         per0 = per;
    %         per0.overload = 0;
    %         per0.gamma = 0;
    %         per0.nsim_singt = 1;
    %     end

else
    % post proc or FEMM simulation (existing geometry)
    copyfile(filenameIn,[pathname filename]); % copy .fem in the temporary folder
end

mat.LayerMag.Br = per.BrPP;
mat.LayerMag.Hc = per.BrPP/(4e-7*pi*mat.LayerMag.mu);

flagSim = 1;
if ~isempty(RQ)
    if per.MechStressOptCheck
        simSetup.evalSpeed = geo.nmax;
        simSetup.meshSize  = 'coarse';
        simSetup.flagFull  = 0;
        simSetup.shaftBC   = 1;
        warning('off')
        %     [structModel] = syre2pde(geo,mat,simSetup);
        simSetup.filename = filename;
        simSetup.pathname = pathname;
        [out.structModel] = femm2pde(geo,mat,simSetup);
        [out.sVonMises,~,out.structModel] = calcVonMisesStress(out.structModel);
        [outMech] = eval_maxStress(out.structModel,out.sVonMises,geo,mat);
        figure
        figSetting
        pdeplot(out.structModel)
        saveas(gcf,[pathname 'mechMesh.fig']);
        close
        warning('on')
        
        if sum([outMech.nPointOverRad outMech.nPointOverTan])>100
            flagSim=0;
        end
        %     if any(outMech.stress_T>mat.Rotor.sigma_max*10^6)  || any(outMech.stress_R>mat.Rotor.sigma_max*10^6)
        %         flagSim = 0;
        %     end
        %         if any(outMech.stress_R>mat.Rotor.sigma_max*10^6)
        %             flagSim = 0;
        %         end
%         if outMech.MaxStress>mat.Rotor.sigma_max*10^6
%             flagSim=0;
%         end
    end

    if geo.statorYokeDivision
        % GalFerContest - structural simulation
        simSetup.evalSpeed = geo.nmax;
        simSetup.meshSize  = 'FEMM original';
        simSetup.flagFull  = 0;
        simSetup.shaftBC   = 1;
        simSetup.meshShaft = 0;
        warning('off')
        simSetup.filename = filename;
        simSetup.pathname = pathname;
        warning('on')
        [structModel] = femm2pde(geo,mat,simSetup);
        [sVonMises,~,structModel] = calcVonMisesStress(structModel,0);
        [outStruct] = eval_maxStress(structModel,sVonMises,geo,mat);
        % GalFerContest - 3D thermal simulation
        if exist([pathname filename(1:end-4) '.fem'],'file')&&~exist([pathname filename(1:end-4) '.ans'],'file')
            openfemm(1);
            opendocument([pathname filename(1:end-4) '.fem']);
            mi_createmesh;
            mi_analyze(1);
            mi_loadsolution;
            closefemm();
        end
        [tempCuAvg,~,tempCuMax] = solve_therm(geo,per,mat,[pathname filename(1:end-4)]);
    end

end

if (strcmp(eval_type,'idemag')||strcmp(eval_type,'idemagmap')||strcmp(eval_type,'demagArea'))
    flagSim = 0;
    flagDemag = 1;
else
    flagDemag = 0;
end
% tic
if flagSim
    if strcmp(geo.RotType,'IM')
        [SOL] = simulate_FOC_IM(geo,per,mat,eval_type,pathname,filename);
    else
        [SOL] = simulate_xdeg(geo,per,mat,eval_type,pathname,filename);
    end
    % standard results
    out.id     = mean(SOL.id(:));                                   % [A]
    out.iq     = mean(SOL.iq(:));                                   % [A]
    out.fd     = mean(SOL.fd);                                      % [Vs]
    out.fq     = mean(SOL.fq);                                      % [Vs]
    out.T      = mean(SOL.T);                                       % [Nm]
    out.dT     = std(SOL.T);                                        % [Nm]
    out.dTpu   = std(SOL.T)/out.T;                                  % [pu]
    out.dTpp   = max(SOL.T)-min(SOL.T);                             % [Nm]
    out.IPF    = sin(atan2(out.iq,out.id)-atan2(out.fq,out.fd));
    out.We     = mean(SOL.we);                                      % [J]
    out.Wc     = mean(SOL.wc);                                      % [J]
    out.SOL    = SOL;

    if isfield(out.SOL,'if')
       out.if = mean(SOL.if(:));
       out.ff = mean(SOL.ff(:));
    end

    if geo.win.n3phase>1
        for ff=1:geo.win.n3phase
            th = out.SOL.th;

            fa = out.SOL.fa(ff,:);
            fb = out.SOL.fb(ff,:);
            fc = out.SOL.fc(ff,:);
            fdq = abc2dq(fa,fb,fc,(th+geo.th0(ff)-geo.th0(1))*pi/180);
            out.SOL.sets(ff).fa = fa;
            out.SOL.sets(ff).fb = fb;
            out.SOL.sets(ff).fc = fc;
            out.SOL.sets(ff).fd = fdq(1,:);
            out.SOL.sets(ff).fq = fdq(2,:);
            out.SOL.sets(ff).f0 = (fa+fb+fc)/3;

            ia = out.SOL.ia(ff,:);
            ib = out.SOL.ib(ff,:);
            ic = out.SOL.ic(ff,:);
            idq = abc2dq(ia,ib,ic,(th+geo.th0(ff)-geo.th0(1))*pi/180);
            out.SOL.sets(ff).ia = ia;
            out.SOL.sets(ff).ib = ib;
            out.SOL.sets(ff).ic = ic;
            out.SOL.sets(ff).id = idq(1,:);
            out.SOL.sets(ff).iq = idq(2,:);
            out.SOL.sets(ff).i0 = (ia+ib+ic)/3;

            out.sets(ff).id = mean(out.SOL.sets(ff).id);
            out.sets(ff).iq = mean(out.SOL.sets(ff).iq);
            out.sets(ff).fd = mean(out.SOL.sets(ff).fd);
            out.sets(ff).fq = mean(out.SOL.sets(ff).fq);
        end
    end
    

    % check Torque sign
    if sign(out.T)~=sign(out.fd*out.iq-out.fq*out.id)
        out.T = -out.T;
        out.SOL.T = -out.SOL.T;
        % warning('Torque sign correction')
    end

    if isfield(SOL,'F')
        out.F=mean(SOL.F);
    end

    if isfield(SOL,'psh')
        out.Pfes_h        = sum(sum(SOL.psh))*(2*geo.p/geo.ps);
        out.Pfes_c        = sum(sum(SOL.psc))*(2*geo.p/geo.ps);
        out.Pfer_h        = sum(sum(SOL.prh))*(2*geo.p/geo.ps);
        out.Pfer_c        = sum(sum(SOL.prc))*(2*geo.p/geo.ps);
        out.Ppm           = sum(sum(SOL.ppm))*(2*geo.p/geo.ps);
        out.ppm_no3D      = sum(sum(SOL.ppm_no3D))*(2*geo.p/geo.ps);
        out.ppm_noRFno3D  = sum(sum(SOL.ppm_noRFno3D))*(2*geo.p/geo.ps);
        out.Ppm_breakdown = SOL.ppm_PM*(2*geo.p/geo.ps);
        out.Pfe           = out.Pfes_h + out.Pfes_c + out.Pfer_h + out.Pfer_c;
        out.velDim        = per.EvalSpeed;

        if strcmp(eval_type,'singmIron')
            % remove all the debug data from SOL, to avoid excessive data size
            SOL = rmfield(SOL,'psh');
            SOL = rmfield(SOL,'psc');
            SOL = rmfield(SOL,'prh');
            SOL = rmfield(SOL,'prc');
            SOL = rmfield(SOL,'ppm');
            %SOL = rmfield(SOL,'ppm_RF');
            %SOL = rmfield(SOL,'ppm_noRF');
            SOL = rmfield(SOL,'ppm_PM');
            SOL = rmfield(SOL,'freq');
            SOL = rmfield(SOL,'bs');
            SOL = rmfield(SOL,'br');
            SOL = rmfield(SOL,'am');
            SOL = rmfield(SOL,'Jm');
            SOL = rmfield(SOL,'pos');
            SOL = rmfield(SOL,'vol');
            SOL = rmfield(SOL,'groNo');
            out.SOL = SOL;
        end
    end

    if strcmp(geo.RotType,'IM')
        out.IM = SOL.IM;
    end

else
    out.T    = 0;
    out.dTpp = 10^50;
    out.IPF  = -10^50;
end
% toc

if flagDemag
    SOL = simulate_xdeg(geo,per,mat,eval_type,pathname,filename);
    
    out.id   = mean(SOL.id);
    out.iq   = mean(SOL.iq);
    out.SOL  = SOL;
    out.Bmin = min(SOL.Bmin);
    out.dPM  = max(SOL.dPM);
end

% Variables necessary for MTPA calculation
flagMTPA = 0;

maxIter = 4;
gammaStep = 2;
direction = 0;

if isfield(per,'if0')
    per.if = per.if0;
else
    per.if = 0;
end
if ~isempty(RQ)
    RQ(end) = 90;
end

if any(strcmp ('gamma', geo.RQnames))
    flagMTPA = 0;
else
    flagMTPA = 1;
end


if ~isempty(RQ)     % MODE optimization (RQ geometry)

    % Cost functions
    cost = zeros(1,length(geo.OBJnames));
    temp1 = 1;
    % Torque
    if strcmp(geo.OBJnames{temp1},'Torque')
        if ~flagMTPA
            cost(temp1) = -out.T;
        else
            % aggiungere ricerca MTPA
            gamma0 = RQ(end);
            jj = 1;
            done = 0;
            TVect    = nan(1,maxIter);
            gVect    = nan(1,maxIter);
            dTppVect = nan(1,maxIter);
            idqVect  = nan(1,maxIter);
            fdqVect  = nan(1,maxIter);

            perTmp = per;
            TmpSOL_old = [];
            TmpSOL_new = [];


             while ~done
                    if jj==1
                        gammaSim = gamma0;
                    elseif jj==2
                        gammaSim = gamma0+gammaStep;
                    elseif jj==3
                        gammaSim = gamma0-gammaStep;
                    else
                        gammaSim = gammaSim+direction*gammaStep;
                    end

                     RQ(end) = gammaSim;
                     TmpSOL = simulate_xdeg(geo,perTmp,mat,eval_type,pathname,filename);
                     TVect(jj)    = mean(TmpSOL.T);
                     gVect(jj)    = gammaSim;
                     dTppVect(jj) = max(TmpSOL.T) - min(TmpSOL.T);
                     idqVect(jj)  = mean(TmpSOL.id)+j*mean(TmpSOL.iq);
                     fdqVect(jj)  = mean(TmpSOL.fd)+j*mean(TmpSOL.fq);
                     
                    % To save the last 2 results of the iteractive process
                     TmpSOL_old = TmpSOL_new;
                     TmpSOL_new = TmpSOL;


                     if jj==3
                            [~,index] = max(TVect,[],'omitnan');
                            if index==1
                               done=1;
                            elseif index==2
                                   direction=+1;
                            else
                                   direction=-1;
                            end
                      elseif jj>3
                          if TVect(jj)<TVect(jj-1)
                             done=1;
                          end
                     end

                     if jj==maxIter
                          done=1;
                     end

                     jj = jj+1;
                     disp(['Simulation ' int2str(jj-1) ' done'])
             end

             [~,index] = max(TVect,[],'omitnan');

            

             % Output data
              OUT.geo   = geo;
              OUT.per   = perTmp;
              OUT.mat   = mat;
              OUT.T     = TVect(index);
              OUT.dTpp  = dTppVect(index);
              OUT.gamma = gVect(index);
              OUT.idq   = idqVect(index);
              OUT.fdq   = fdqVect(index);
              OUT.RQ    = RQ;
              %OUT.nFEMM = nFEMM;
              OUT.Pjs   = 3/2*per.Rs*abs(OUT.idq)^2;

              % Identify Rf value depending on rotor geometry
              if strcmp(geo.RotType, 'EESM')
                 Rf=per.Rf;
                 OUT.Pjf   = Rf*out.if;
              else
                  Rf=0;
                  OUT.Pjf = nan;
              end

              if mean(TmpSOL_new.T) > mean(TmpSOL_old.T)
                  OUT.SOL = TmpSOL_new;
              else
                  OUT.SOL = TmpSOL_old;
              end

             out.SOL = OUT.SOL;
             out.id     = real(OUT.idq);                                   % [A]
             out.iq     = imag(OUT.idq);                                   % [A]
             out.fd     = real(OUT.fdq);                                   % [Vs]
             out.fq     = imag(OUT.fdq);                                   % [Vs]
             out.T      = OUT.T;                                           % [Nm]
             out.dTpp   = OUT.dTpp;                                        % [Nm]
             out.gamma = OUT.gamma;                                        % [degrees] MTPA angle

             cost(temp1) = -mean(OUT.T);
             temp1 = temp1+1;
   
        end
        
    end

    % Torque Ripple
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'TorRip')
        %         cost(temp1) = out.dTpu*100;
        cost(temp1) = out.dTpp;
        temp1 = temp1+1;
    end
    % Copper Mass
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassCu')
        if flagSim
            cost(temp1) = calcMassCu(geo,mat);
            temp1=temp1+1;
        else
            cost(temp1) = 10^50;
        end
    end
    % PM Mass
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassPM')
        if flagSim
            cost(temp1) = calcMassPM(geo,mat);
            temp1=temp1+1;
        else
            cost(temp1) = 10^50;
        end
    end

    % Power Factor
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'PF')
        cost(temp1) = -out.IPF;
        temp1=temp1+1;
    end

    % No Load flux
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'Fdq0')
        if flagSim
            per0 = per;
            per0.overload = 0;
            per0.gamma = 0;
            per0.nsim_singt = 1;
            per0.nsim_MOOA = 1;
            [SOL0] = simulate_xdeg(geo,per0,mat,eval_type,pathname,filename);
            cost(temp1) = abs(SOL0.fd+j*SOL0.fq);
        else
            cost(temp1) = 10^50;
        end
        temp1 = temp1+1;
    end

     % Structural properties
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MechStress')
        cost(temp1) = max([outMech.sigmaRadMax outMech.sigmaTanMax])/1e6;
        temp1=temp1+1;
    end
    
    % Stator Joule loss
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'Pjs')
        if flagSim
            cost(temp1) = 3/2*per.Rs*abs(out.id+j*out.iq)^2;
        else
            cost(temp1) = 10^50;
        end
        temp1 = temp1+1;
    end

    % Field Joule loss
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'Pjf')
        if flagSim
            cost(temp1) = per.Rf*out.if^2;
        else
            cost(temp1) = 10^50;
        end
        temp1 = temp1+1;
    end


    % Structural GalFerContest
    if geo.statorYokeDivision
        cost = [cost outStruct.sigmaTotPrc/1e6 tempCuMax tempCuAvg];
    else
        % penalize weak solutions
        for ii = 1:length(cost)
            if cost(ii)>per.objs(ii,1) && per.objs(ii,3)==0
                if per.objs(ii,1)>0
                    cost(ii)=cost(ii)*10;  % minimization problem
                else
                    cost(ii)=cost(ii)*0.1; % maximization problem
                end
            end
    
            if ((cost(ii)<per.objs(ii,1)-per.objs(ii,3)) || (cost(ii)>per.objs(ii,1)+per.objs(ii,3))) && per.objs(ii,3)
                if per.objs(ii,1)>0
                    cost(ii)=cost(ii)*10;  % minimization problem
                else
                    cost(ii)=cost(ii)*0.1; % maximization problem
                end
            end
        end
    end

    %     dataSetPath = strcat(thisfilepath,'\dataSet.mat');    %OCT
    load('dataSet.mat');
    geo.RQ = RQ;

    [dataSet] = SaveInformation(geo,mat,dataSet);
    if isoctave()            %OCT
        save('-mat7-binary', strrep(filename,'.fem','.mat'),'geo','cost','per','dataSet','mat');
    else
        save([pathname strrep(filename,'.fem','.mat')],'geo','cost','per','dataSet','mat','out');
    end
else
    cost = [];
    save([pathname strrep(filename,'.fem','.mat')],'geo','out','mat','per');   % save geo and out
end

