function [motorModel] = MMM_build_constantParameterPMSM(setup)


if nargin()==0
    % setup.Ld   = 2.6450e-4;
    % setup.Lq   = 5.7130e-4;
    % setup.Fm   = 0.1134;
    % setup.Vdc  = 550;
    % setup.Imax = 700;
    % setup.p    = 3;
    % setup.nmax = 18000;
    prompt = {'Ld (H)','Lq (H)','Fm (Vs)','Vdc (V)','Imax (Apk)','p','nmax (rpm)'};
    name   = 'Linear model parameters';
    answer = {
        '2.6450e-4';
        '5.7130e-4';
        '0.1134';
        '550';
        '700';
        '3';
        '18000';
        };
    answer = inputdlg(prompt,name,1,answer);
    setup.Ld   = eval(answer{1});
    setup.Lq   = eval(answer{2});
    setup.Fm   = eval(answer{3});
    setup.Vdc  = eval(answer{4});
    setup.Imax = eval(answer{5});
    setup.p    = eval(answer{6});
    setup.nmax = eval(answer{7});
end


Ld   = setup.Ld;
Lq   = setup.Lq;
Fm   = setup.Fm;
Vdc  = setup.Vdc;
Imax = setup.Imax;
p    = setup.p;
nmax = setup.nmax;

Ich = Fm/Ld;

Itmp = max([Imax Ich]*1.1);

Id = linspace(-Itmp,Itmp,257);
Iq = linspace(0,Itmp,257);

[Id,Iq] = meshgrid(Id,Iq);
Fd = Ld*Id+Fm;
Fq = Lq*Iq;
T = 3/2*p*(Fd.*Iq-Fq.*Id);
dTpp = zeros(size(Id));
dT   = zeros(size(Id));


fdfq.Id   = Id;
fdfq.Iq   = Iq;
fdfq.Fd   = Fd;
fdfq.Fq   = Fq;
fdfq.T    = T;
fdfq.dTpp = dTpp;
fdfq.dT   = dT;


motorModel = MMM_defaultMotorModel();
motorModel = MMM_back_compatibility(motorModel,0);

motorModel.data.i0         = Imax;
motorModel.data.Imax       = Imax;
motorModel.data.Vdc        = Vdc;
motorModel.data.p          = p;
motorModel.data.nmax       = nmax;
motorModel.data.Rs         = 0;
motorModel.data.l          = 1;
motorModel.data.R          = 1;
motorModel.data.n3phase    = 1;
motorModel.data.axisType   = 'PM';
motorModel.data.Ns         = 1;
motorModel.data.nCurr      = 1;
motorModel.data.motorName  = 'ConstantParameterMotor';
motorModel.data.tempPM     = 20;
motorModel.data.tempVectPM = 20;

motorModel.TnSetup.nCurrent = 1;

motorModel.FluxMap_dq = fdfq;

motorModel = rmfield(motorModel,'dataSet');


if Fm==0
    motorModel.data.motorType = 'SR';
else
    motorModel.data.motorType = 'PM';
end