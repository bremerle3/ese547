rtd  = 180/pi;
dt = 0.002;
t = 0:dt:2;
w = logspace(-1,3,500);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);



% Public release airframe parameters
Za = -1.05273;
Zd = -0.0343;
Ma = -2.3294;
Mq = -1.03341;
Md = -1.1684;
% Vehicle velocity (fps):
V = 329.127;
wa = 2*pi*13.; za = 0.6;
grav = 32.174; % (fps2)
%% Plant Model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plant model with Az state
% states = [Az, q, deltae, deltae_dot]
% controls = [deltaec]
% outputs = [Az, q, deltae, deltae_dot]
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ap = [ Za V*Za 0. V*Zd
Ma/(V*Za) Mq (Md-Ma*Zd/Za) 0.;
0. 0. 0. 1.;
0. 0. -wa^2 -2*za*wa];
Bp = [ 0 ; 0 ; 0 ; wa^2];
Cp = [1 0. 0. 0.;
     eye(4)
     ];
Dp = 0.*Cp*Bp;



Aw = [[0 1 0 0 0]; [zeros(4,1), Ap]];

Bw = [0; Bp];

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
    %Q(1,1)=qq(ii);
    Q(1,1)=0.2448;
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

Kx_lqr = Kx_lqr_st(20,:);
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

figure; plot(t,az); title('az (acceleration)');grid on;
figure; plot(t,q); title('Pitch Rate (deg/sec)');grid on;
figure; plot(t,del); title('Elevon (deg)');grid on;
figure; plot(t,deldot); title('Elevon rate (deg/sec)');grid on;
