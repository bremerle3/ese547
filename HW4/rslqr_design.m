rtd  = 180/pi;
dt = 0.002;
t = 0:dt:2;
w = logspace(1,2,500);

%**************************************************************************
% Aircraft Model Data
%**************************************************************************
% Model (taken from example 5.2, the aircraft pitch axis plant data)
% Aero constants used in the aircraft model:
Za_V = -1.3046;
Ma = 47.711;
Zd_V = -.2142;
Md = -104.83;
V = 886.78; % (fps)
Za = V*Za_V;
Zd = V*Zd_V;
w_act = 11*rtd; % actuator natural frequency rps
z_act = 0.707; % actuator damping
grav = 32.174; % (fps2)
V_fps     =   329.127;
ZapV  =   -1.05273;
ZdpV  =   -0.0343;

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
Aw = [ 0. V_fps*ZapV  0. ;
    0.   ZapV  1. ;
    0.   Ma    0 ];
% Input is dele-cmd
Bw = [ V_fps*ZdpV ;
    ZdpV;
    Md];


% Setup range of penalties for the LQR
Q=0.*Aw; %sets up a matrix Q that is the size of Aw and is all 0s
R=1;
xeig=[];
%qq is a vector which each element scales the LQR Q matrix
qq=logspace(-6,-1,50); %these are the varying values of q11
xopenloop=eig(Aw);
% Preallocate matrices for speed
az_st = 0.*ones(numel(qq),numel(t));
del_st = 0.*ones(numel(qq),numel(t));
deldot_st = 0.*ones(numel(qq),numel(t));

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
    del_st(ii,:) = rtd*y(:,4);
    deldot_st(ii,:) = rtd*y(:,5);
    ii = ii+1;
    
end



% figure; plot(data(:,4),data(:,5),'x'); title('Rise Time');
% figure; plot(data(:,4),data(:,6),'x'); title('Settling Time');
% figure; plot(data(:,4),data(:,2),'x'); title('rdmin');
% figure; plot(data(:,4),data(:,3),'x'); title('srmin');
% figure; plot(data(:,4),data(:,7),'x'); title('Overshoot');
% figure; plot(data(:,4),data(:,6),'x'); title('Undershoot');

figure; plot(data(1:24,4),data(1:24,5),'x'); title('Rise Time');
figure; plot(data(1:24,4),data(1:24,6),'x'); title('Settling Time');
figure; plot(data(1:24,4),data(1:24,2),'x'); title('rdmin');
figure; plot(data(1:24,4),data(1:24,3),'x'); title('srmin');
figure; plot(data(1:24,4),data(1:24,7),'x'); title('Overshoot');
figure; plot(data(1:24,4),data(1:24,6),'x'); title('Undershoot');