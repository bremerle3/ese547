save_plots = 0; % Flag to bypass saving
rtd  = 180/pi;
dt = 0.002;
t = 0:dt:2;
w = logspace(-1,3,500);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);

q_idx = 25;

%**************************************************************************
% Aircraft Model Data
%**************************************************************************
% Model (taken from example 5.2, the aircraft pitch axis plant data)
% Aero constants used in the aircraft model:
Za_V = -1.3046;
Ma = 47.711;
Zd_V = -.2142;
Md = -104.83;
V_fps = 886.78; % (fps)
Za = V_fps*Za_V;
Zd = V_fps*Zd_V;
w_act = 2.*pi*11; % actuator natural frequency rps
z_act = 0.707; % actuator damping
grav = 32.174; % (fps2)
Mq = 0;

%**************************************************************************
% Plant model for analysis with actuator
%**************************************************************************
% States are AOA (rad) pitch rate (rps)
Ap = [Za_V 1. Zd_V 0.;
    Ma 0. Md 0.;
    0. 0. 0. 1.;
    0. 0. -w_act*w_act -2*z_act*w_act];
% Input is dele (rad)
Bp = [0.; 0.; 0.; w_act*w_act ];
% Outputs are Az (fps2) AOA (rad) pitch rate (rps) dele (rad) deledot (rps)
Cp = [ Za 0. Zd 0.;
    eye(4)];
Dp = 0.*Cp*Bp;

%*******************************************
% Design the RSLQR using 3rd order Aw and Bw
%*******************************************
% State Feedback Design
% states are [ int-err AOA q ]
% Aw = [ 0. V_fps*Za_V  0. ;
%     0.   Za_V  1. ;
%     0.   Ma    Mq ];
% % Input is dele-cmd
% Bw = [ V_fps*Zd_V ;
%     Zd_V;
%     Md];

Aw = [[0 V_fps*Za_V 0 0 0]; [zeros(4,1), Ap]];

Bw = [V_fps*Zd_V; Bp];


    
       

% Setup range of penalties for the LQR
Q=0.*Aw; %sets up a matrix Q that is the size of Aw and is all 0s
R=1;
xeig=[];
%qq is a vector which each element scales the LQR Q matrix
qq=logspace(-6,-1,50); %these are the varying values of q11
xopenloop=eig(Aw);
% Preallocate matrices for speed
az_st = 0.*ones(numel(qq),numel(t));
q_st = 0.*ones(numel(qq),numel(t));
del_st = 0.*ones(numel(qq),numel(t));
deldot_st = 0.*ones(numel(qq),numel(t));
Kx_lqr_st = (ones(numel(qq),5));
Kx_lqr_st(:,:) = complex(1,1);

% Loop LQR design increasing the LQR penealty (qq)
% and compute the time domain and frequency domain metrics as the penalty
% is varied
npts = numel(qq);
ii = 1;
while ii < npts,
    % for ii = 1:numel(qq),
    % The only penalty is the 1,1 element on the error state
    % Desing the gains
    Q(1,1)=qq(ii);
    [Kx_lqr,~,~]=lqr(Aw,Bw,Q,R);
    % populate the controller matrices
    Ac = 0.;
    Bc1 = [1. 0. 0. 0. 0.];
    Bc2 = -1;
    Cc = -Kx_lqr(1);
    Dc1 = [0. -Kx_lqr(2:3) 0. 0.];
    Dc2 = 0.;
    
    % Form the closed loop system
    Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
    Acl = [ (Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);
        (Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
    Bcl = [ Bp*Z*Dc2;
        (Bc2+Bc1*Dp*Z*Dc2)];
    Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
    Dcl =(Dp*Z*Dc2);
    sys_cl = ss(Acl,Bcl,Ccl,Dcl);
    % Compute closed loop eigenvalues for root locus plot
    xx=eig(Acl);
    xeig=[xeig;xx];
    
    % SS model of loop gain at the plant input Lu
    A_Lu = [ Ap 0.*Bp*Cc; Bc1*Cp Ac];
    B_Lu = [ Bp; Bc1*Dp];
    C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
    D_Lu = -[ Dc1*Dp];
    sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
    magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w))));
    wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
    sr = sigma(sys_Lu,w,3);
    srmin = min(abs(sr));
    rd = sigma(sys_Lu,w,2);
    rdmin = min(abs(rd));
    if rdmin <= 0.25,
        disp('Exit loop on rdmin')
        ii = npts;
    end
    
    y = step(sys_cl,t);
    az = y(:,1); % acceleration (fps2)
    q = y(:,3).*rtd; % pitch rate (dps)
    aze = abs(ones(size(az))-az); % error for az
    taur = 0.; taus= 0.; % rise time and settling time
    fv = aze(numel(aze)); % final value of the error
    e_n = aze - fv*ones(size(aze)) - 0.36*ones(size(aze));
    e_n1 = abs(e_n) + e_n;
    taur = crosst(e_n1,t); % rise time
    e_n = aze - fv*ones(size(aze)) - 0.05*ones(size(aze));
    e_n1 = abs(e_n) + e_n;
    taus = crosst(e_n1,t); % settling time
    azmin = abs(min(az))*100; % undershoot
    azmax = (abs(max(az))-1)*100; % overshoot
    dmax = max(abs(y(:,4)))*rtd*grav; % compute in per g commanded
    ddmax = max(abs(y(:,5)))*rtd*grav;
    metric=[qq(ii) rdmin srmin wc taur taus azmin azmax dmax ddmax];
    data(ii,:) = metric;
    az_st(ii,:) = az';
    q_st(ii,:) = q';
    del_st(ii,:) = rtd*y(:,4);
    deldot_st(ii,:) = rtd*y(:,5);
    Kx_lqr_st(ii,:) = Kx_lqr;
    ii = ii+1;
end

% figure; plot(data(:,4),data(:,2),'x'); title('rdmin');
% figure; plot(data(:,4),data(:,3),'x'); title('srmin');
% figure; plot(data(:,4),data(:,5),'x'); title('Rise Time');
% figure; plot(data(:,4),data(:,6),'x'); title('Settling Time');
% figure; plot(data(:,4),data(:,7),'x'); title('undershoot');
% figure; plot(data(:,4),data(:,8),'x'); title('overshoot');
% figure; plot(data(:,4),data(:,9),'x'); title('dmax');
% figure; plot(data(:,4),data(:,10),'x'); title('ddmax');


% figure; plot(t,az_st(q_idx,:),'x'); title('az (acceleration)');grid on;
% figure; plot(t,q_st(q_idx,:),'x'); title('Pitch Rate (deg/sec)');grid on;
% figure; plot(t,del_st(q_idx,:),'x'); title('Elevon (deg)');grid on;
% figure; plot(t,deldot_st(q_idx,:),'x'); title('Elevon rate (deg/sec)');grid on;

Kx_lqr = Kx_lqr_st(q_idx,:);
F = Aw - Bw*Kx_lqr;
[evecs, evals] = eig(F);
evals = sum(evals,1);
evals_sorted = sort(evals);
evec_idx = [];
for ii = 1:3
    evec_idx = [evec_idx, find(evals == evals_sorted(ii))];
end
dominant_vecs = evecs(:,evec_idx);
C_proj = [eye(3), repmat([0 0], 3,1)];
K_proj = Kx_lqr*dominant_vecs*(C_proj*dominant_vecs)^-1;

    Kx_lqr2 = K_proj;
    % populate the controller matrices
    Ac = 0.;
    Bc1 = [1. 0. 0. 0. 0.];
    Bc2 = -1;
    Cc = -Kx_lqr2(1);
    Dc1 = [0. -Kx_lqr2(2:3) 0. 0.];
    Dc2 = 0.;
    
    % Form the closed loop system
    Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
    Acl = [ (Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);
        (Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
    Bcl = [ Bp*Z*Dc2;
        (Bc2+Bc1*Dp*Z*Dc2)];
    Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
    Dcl =(Dp*Z*Dc2);
    sys_cl = ss(Acl,Bcl,Ccl,Dcl);
    % Compute closed loop eigenvalues for root locus plot
    xx=eig(Acl);
    xeig=[xeig;xx];
    
    % SS model of loop gain at the plant input Lu
    A_Lu = [ Ap 0.*Bp*Cc; Bc1*Cp Ac];
    B_Lu = [ Bp; Bc1*Dp];
    C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
    D_Lu = -[ Dc1*Dp];
    sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
    magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w))));
    wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
    sr = sigma(sys_Lu,w,3);
    srmin = min(abs(sr));
    rd = sigma(sys_Lu,w,2);
    rdmin = min(abs(rd));
    if rdmin <= 0.25,
        disp('Exit loop on rdmin')
        ii = npts;
    end
    
    y = step(sys_cl,t);
    az = y(:,1); % acceleration (fps2)
    q = y(:,3).*rtd; % pitch rate (dps)
    del =rtd*y(:,4)
    deldot= rtd*y(:,5);
    aze = abs(ones(size(az))-az); % error for az
    taur = 0.; taus= 0.; % rise time and settling time
    fv = aze(numel(aze)); % final value of the error
    e_n = aze - fv*ones(size(aze)) - 0.36*ones(size(aze));
    e_n1 = abs(e_n) + e_n;
    taur = crosst(e_n1,t); % rise time
    e_n = aze - fv*ones(size(aze)) - 0.05*ones(size(aze));
    e_n1 = abs(e_n) + e_n;
    taus = crosst(e_n1,t); % settling time
    azmin = abs(min(az))*100; % undershoot
    azmax = (abs(max(az))-1)*100; % overshoot
    dmax = max(abs(y(:,4)))*rtd*grav; % compute in per g commanded
    ddmax = max(abs(y(:,5)))*rtd*grav;

rt_5_2 = taur;save('../HW5/rt_5_2.mat', 'rt_5_2');
st_5_2 = taus;save('../HW5/st_5_2.mat', 'st_5_2');
    
figure; plot(t,az); title('az (acceleration)');grid on;xlabel('time (sec)');
az_5_2 = az; save('5_2_az.mat', 'az_5_2');
figure; plot(t,q); title('Pitch Rate (deg/sec)');grid on;xlabel('time (sec)');
q_5_2 = q; save('5_2_q.mat', 'q_5_2');
figure; plot(t,del); title('Elevon (deg)');grid on;xlabel('time (sec)');
del_5_2 = del; save('5_2_del.mat', 'del_5_2');
figure; plot(t,deldot); title('Elevon rate (deg/sec)');grid on;xlabel('time (sec)');
deldot_5_2 = deldot; save('5_2_deldot.mat', 'deldot_5_2');

pause;
close all;


Lu = freqresp(sys_Lu,w);
SRu_lqr = sigma(sys_Lu,w,3);
SRu_min = min(abs(SRu_lqr));
RDu_lqr = sigma(sys_Lu,w,2);
RDu_min = min(abs(RDu_lqr));


figure('Name','Nyquist Plot at Plant Input'),
plot(xx1,yy1,'r',real(squeeze(Lu)),imag(squeeze(Lu)),'k','LineWidth',2);grid
axis([-5 5 -5 5]);
legend('Unit Circle','RSLQR','Location','Best');
xlabel('Re(Lu)')
ylabel('Im(Lu)')
title('Nyquist Plot at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Bode Magnitude at Plant Input'),
semilogx(w,20*log10(abs(squeeze(Lu))),'k',w,20*log10(abs(squeeze(Lu))),'b','LineWidth',2);grid
xlabel('Frequency (rps)')
ylabel('Mag')
title('Bode at Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Return Difference at Plant Input'),
semilogx(w,20*log10(abs(RDu_lqr)),'k','LineWidth',2);grid
obltr_min = ['RD min = ' num2str(RDu_min)];
legend(obltr_min,'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Stability Robustness at Plant Input'),
semilogx(w,20*log10(abs(SRu_lqr)),'r--','LineWidth',2);grid
legend(['min(I+invLu) = ' num2str(SRu_min)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Stability Robustness at Plant Input2 with bode'),
semilogx(w,20*log10(abs(SRu_lqr)),'k','LineWidth',2);grid
obltr_min = ['SR min = ' num2str(SRu_min)];
legend(obltr_min,'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

disp('Classical Margins')
%allmargin(sys_Lu)

disp('  ')
disp('SV Margins')
RDu_nGM = 1/(1+RDu_min);
RDu_pGM = 1/(1-RDu_min);
RDu_Pha = 2*asin(RDu_min/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi ;
disp('RDu_nGM RDu_pGM RDu_Pha')
disp([num2str(RDu_nGM) ' ' num2str(RDu_pGM) ' ' num2str(RDu_Pha)])
disp([num2str(RDu_nGM_dB) ' ' num2str(RDu_pGM_dB) ' ' num2str(RDu_Pha_deg)])
SRu_nGM = 1-SRu_min;
SRu_pGM = 1+SRu_min;
SRu_Pha = 2*asin(SRu_min/2);
SRu_nGM_dB = 20*log10(SRu_nGM);
SRu_pGM_dB = 20*log10(SRu_pGM);
SRu_Pha_deg = 180*SRu_Pha/pi ;
disp('SRu_nGM SRu_pGM SRu_Pha')
disp([num2str(SRu_nGM) ' ' num2str(SRu_pGM) ' ' num2str(SRu_Pha)])
disp([num2str(SRu_nGM_dB) ' ' num2str(SRu_pGM_dB) ' ' num2str(SRu_Pha_deg)])
disp('  ')

margin(sys_Lu)