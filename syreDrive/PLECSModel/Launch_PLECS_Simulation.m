

% Add server for launching PLECS simulation 

addpath('matlab-jsonrpc-main\');
URL = 'http://localhost:1080';
proxy = jsonrpc(URL,'Timeout', 100);



load motorModel

MTPA    = motorModel.controlTrajectories.MTPA;
i0      = motorModel.data.i0;
idRef = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.id,i0);
iqRef = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.iq,i0);
nRef = motorModel.data.n0;


fileName = [motorModel.data.motorName '_Motor_ctrl'];

% Solver option
optStruct.SolverOpts.StartTime = 0.000;
optStruct.SolverOpts.TimeSpan  = 0.5;

optStruct.ModelVars.idRef = idRef;
optStruct.ModelVars.iqRef = iqRef;
optStruct.ModelVars.n_load  = n0;

% Launch simulation 
out = proxy.plecs.simulate(fileName,optStruct);



%% PLot Simulation Results

time = out.Time;
ia = out.Values(1,:);
ib = out.Values(2,:);
ic = out.Values(3,:);
fd = out.Values(4,:);
fq = out.Values(5,:);
id = out.Values(6,:);
iq = out.Values(7,:);
id_ref = out.Values(8,:);
iq_ref = out.Values(9,:);
Tm = out.Values(10,:);
theta_e = out.Values(11,:);



figure()
figSetting()
plot(time,ia,'DisplayName','$i_a$');
plot(time,ib,'DisplayName','$i_b$');
plot(time,ic,'DisplayName','$i_c$');
legend('Location','northeast');
xlabel('$t$ [s]');
ylabel('$i_{abc}$ [A]');


figure()
figSetting(12,12)
subplot(2,1,1)
plot(time,id,'DisplayName','$i_d$');
plot(time,id_ref,'DisplayName','$i_d^*$');
legend('Location','northeast');
xlabel('$t$ [s]');
ylabel('$i_{d}$ [A]');

subplot(2,1,2)
plot(time,iq,'DisplayName','$i_q$');
plot(time,iq_ref,'DisplayName','$i_q^*$');
legend('Location','northeast');
xlabel('$t$ [s]');
ylabel('$i_{q}$ [A]');


figure()
figSetting(12,15)
subplot(3,1,1)
plot(time,fd,'DisplayName','$\lambda_d$');
legend('Location','northeast');
xlabel('$t$ [s]');
ylabel('$\lambda_d$ [Vs]');

subplot(3,1,2)
plot(time,fq,'DisplayName','$\lambda_q$');
legend('Location','northeast');
xlabel('$t$ [s]');
ylabel('$\lambda_q$ [Vs]');

subplot(3,1,3)
plot(time,Tm,'DisplayName','$T_m$');
legend('Location','northeast');
xlabel('$t$ [s]');
ylabel('$T_m$ [Nm]');
















