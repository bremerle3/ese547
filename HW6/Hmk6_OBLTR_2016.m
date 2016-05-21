%
% Homework 6 OBLTR Design
% Kevin A. Wise
% ESE 547 Washington University
% Spring 2016
%
clear all
close all
clc
format short e

disp('****************************** Program Start ****************');
plot_file_name = 'Homework_6.ppt';
save_plots = 0; % Flag to bypass saving
w = logspace(-3,4,1000);
t = linspace(0,2.5,1000);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);
rtd = 180/pi;

%**************************************************************************
% Aircraft Model Data
%**************************************************************************
% Model (taken from example 5.2, the aircraft pitch axis plant data)
% Aero constants used in the aircraft model:
Za_V = -1.3046;
Ma   = 47.711;
Zd_V = -.2142;
Md   = -104.83;
V    = 886.78; % (fps)
Za   = V*Za_V;
Zd   = V*Zd_V;
w_act = 2.*pi*11; % actuator natural frequency rps
z_act = 0.707;    % actuator damping

grav = 32.174; % (fps2)

%**************************************************************************
% Plant model  
%**************************************************************************
% States are AOA (rad)  pitch rate (rps)
    Ap = [Za_V  1.  ; 
          Ma    0.  ];
% Input is dele (rad)
    Bp = [ Zd_V;
           Md ];
% Outputs are Az (fps2) AOA (rad) pitch rate (rps)
    Cp = [  Za   0. ; 
            eye(2)];
    Dp =  [ Zd ;
            0.;
            0.]; 
    
    disp('Open loop Plant Model')
    disp('Plant states = [ AOA q ]')
    disp('Plant inputs u = dele_cmd')
    disp('Plant Output y = [ Az AOA q]')
    Ap
    Bp
    Cp
    Dp

%**************************************************************************
% RSLQR Design Model Without The Actuator
%**************************************************************************
% design a RSLQR command to track Az using state feeback controller 
% form the wiggle system:
% state vector: [int(e) (fps), AOA (rad), pitch rate q (rps)]' 
Aw = [0.  Za     0.; 
      0.  Za_V   1.; 
      0.   Ma    0.];

Bw = [   Zd;
        Zd_V;
         Md]; 

% Setup range of penalties for the LQR
Q=0.*Aw; %sets up a matrix Q that is the size of Aw and is all 0s
R=1; 
xeig=[];
%qq is a vector which each element scales the LQR Q matrix
qq=logspace(-6,-1,50); %these are the varying values of q11 

xopenloop=eig(Aw);

% Use number 20 from Hmk 4
ip = 20;
Q(1,1)=qq(ip);

[Klqr,~,~] = lqr(Aw,Bw,Q,R);
   
% populate the controller matrices  
    Ac =  0.;
    Bc1 = [1. 0. 0.];
    Bc2 =  -1;
    Cc = -Klqr(1);
    Dc1 = [0. -Klqr(2:3)];
    Dc2 = 0.;

    disp('RSLQR State feedback Controller')
    disp('Controller states = int-err')
    disp('Controller inputs y = [ Az AOA q ]  r = [ Azcmd ]')
    disp('Controller Output u = dele_cmd')
    Ac
    Bc1
    Bc2
    Cc 
    Dc1 
    Dc2

% Form the closed loop system
     Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
     Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc);
         (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
     Bcl = [       Bp*Z*Dc2;
         (Bc2+Bc1*Dp*Z*Dc2)];
     Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
     Dcl =(Dp*Z*Dc2);
     sys_cl = ss(Acl,Bcl,Ccl,Dcl);

     % Time Domain Analysis
     y = step(sys_cl,t);
     az = y(:,1); %  acceleration (fps2)
     aze = abs(ones(size(az))-az);  % error for az
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
%      dmax = max(abs(y(:,4)))*rtd*grav; % compute in per g commanded
%      ddmax = max(abs(y(:,5)))*rtd*grav;
     
% SS model of loop gain at the plant input Lu
    Ain = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac];
    Bin = [ Bp; Bc1*Dp];
    Cin = -[ Dc1*Cp Cc];%change sign for loop gain
    Din = -[ Dc1*Dp];
    sys_Lu = ss(Ain,Bin,Cin,Din);
%SS model of loop gain L at the plant output
    Aout = [ Ap Bp*Cc;  0.*Bc1*Cp Ac];
    Bout = [ Bp*Dc1; Bc1];
    Cout = -[ Cp Dp*Cc];%change sign for loop gain
    Dout = -[ Dp*Dc1];
    sys_Ly = ss(Aout,Bout,Cout,Dout);

    Lu = freqresp(sys_Lu,w);
    Ly = freqresp(sys_Ly,w);
    T  = freqresp(sys_cl,w);
    S = 1 - T;
     
    magdb = 20*log10(abs(squeeze(Lu)));
    wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
    sr = sigma(sys_Lu,w,3);
    SRu_min = min(abs(sr));
    rd = sigma(sys_Lu,w,2);
    RDu_min = min(abs(rd));

% Save for comparison with OBLTR
    dmax  = 0.;
    ddmax = 0.;
    lqr_metric = [qq(ip) RDu_min SRu_min wc taur taus azmin azmax dmax ddmax];
    az_lqr  = az;
    Lu_lqr  = Lu;    
    SRu_lqr = sr;
    RDu_lqr = rd;
    T_lqr   = T;
    S_lqr   = S;
    
%*************************************************************************
%   OBLTR
%*************************************************************************
%  pstatnam = {'int-err' 'aoa' 'q'};
%  pinnam = {'dele' };
%  poutnam = {'int-err' 'q' };
    Cw = [ 1 0 0;
           0 0 1];
    Dw = [ 0;
            0];
% Plant disturbance and measurement noise covariance amtrices
    Q0 =  max(svd(Bw))*eye(3);
    R0 = eye(2);

% [L,P,E] = LQE(A,G,C,Q,R,N) 
% % Solve for the steady state Kf gain matrix
%disp('LQG Design');

    [K_lqg, P, E] = lqr(Aw', Cw', Q0, R0);

    disp('FARE P Matrix')
    P_lqg    = P
    disp('Observer Gain Matrix')
    L_lqg = K_lqg'
% FARE solution check%
% important to check this. This should be zero
    kf_eq = Aw*P + P*Aw' + Q0 - P*Cw'*inv(R0)*Cw*P;
    norm_kf_eq = norm(kf_eq);
    disp(['FARE Solution Tolerance is: ', num2str(norm_kf_eq)]);


% Must have a zero D- matrix to square up the system with this algorithm
    sys_1 = ss(Aw, Bw, Cw, Dw)
    sys_2 = square_system_inputs_pp(sys_1,-200)

% Compute the Zero Dynamics of the Squared-Up Extended Observer
% Writes them as transfer functions to the workspace
[ny,nu] = size(sys_2.d);
for ii = 1:ny,
    for jj = 1:nu,
        disp(['i_ny = ' num2str(ii) ' j_nu = ' num2str(jj) ])
        [n1,d1] = ss2tf(sys_2.a,sys_2.b(:,jj),sys_2.c(ii,:),sys_2.d(ii,jj));
        zpk(tf(n1,d1))
    end
end

nn = tzero(sys_2);
if isempty(nn)
    disp('No Finite Transmission Zeros in Squared-Up Extended Observer');
else
    disp('Extended Plant Transmission Zeros Bbar:');
    disp(['  Min Real Part = ', num2str(min(real(nn)))]);
    disp(['  Max Real Part = ', num2str(max(real(nn)))]);
end;

Bw_mod = sys_2.b;
Dw_mod = sys_2.d;

% Normalize squared-up columns to min SV of Bw
min_svd_Bw = min(svd(Bw));
mp = 1;
nCw_meas = 2;
for ii = mp+1:nCw_meas
    Bw_mod(:,ii) = Bw_mod(:,ii)./norm(Bw_mod(:,ii))*min_svd_Bw;
end;

Bw_mod 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Check squared-up C*B rank
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
min_svd_CB = min(svd(Cw*Bw_mod));

if min_svd_CB < 1e-05
  disp('Extended C*B is rank defficient');
end;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Check extended system transmission zeros
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sysw_zeros = tzero(Aw, Bw_mod, Cw, Dw_mod);

if isempty(sysw_zeros)
  disp('No Transmission Zeros in Extended System');
else
  disp('Extended System Transmission Zeros:');
  damp(sysw_zeros);
end;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Observer weights
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%rho = 1.
rho = 0.8

Qwe = Q0 +  (1/rho)*(Bw_mod)*(Bw_mod)';

Rwe = rho*R0;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Observer Gain (ARE olution)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[Kwf, Pwf, Ewf] = lqr(Aw', Cw', Qwe, Rwe);
Lv = Kwf'
% ARE filter solution check
lqrf_eq = Pwf*Aw' + Aw*Pwf + Qwe - Pwf*Cw'*inv(Rwe)*Cw*Pwf;
norm_lqrf_eq = norm(lqrf_eq);
disp(['Filter ARE Solution Tolerance is: ', num2str(norm_lqrf_eq)]);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% P*B directions for adaptive design
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SVD decompposition
[U, Lambda, V] = svd(Bw_mod'*Cw'*inv(sqrt(R0)));
W = (U*V')';

% PB = C testing
Pw = inv(Pwf);

% Coefficients on e_x' = (x_hat - x)' in state feedback adaptive laws
PwBw = Pw*Bw;

% Coefficients on e_y' = (y_hat - y)' in output feedback adaptive laws
Cw_computed = Cw'*inv(sqrt(R0))*W;
Cw_computed = Cw_computed(:,1:mp)';

% Adaptation directions: angle between Pw*Bw and Cw_computed'
angle_PBC_deg = acosd(Cw_computed*PwBw/norm(Cw_computed)/norm(PwBw));
disp([  'Adaptation Directions Angle: ', num2str(angle_PBC_deg), ' deg']);

F_r = [ -1; 0.; 0.];
Ac_lqg = Aw-Bw*Klqr-Lv*Cw;
Bc_lqg = [ Lv F_r];
Cc_lqg =  -Klqr;
Dc_lqg = [ 0. 0. 0.];


% Plant Model
disp('Modify the plant to only feedback int-err and q')
Cp_save = Cp;
Dp_save = Dp;
Ap
Bp
Cp = [Cp(1,:);Cp(3,:)]
Dp = [Dp(1,:);Dp(3,:)]

% % Model controller in state space form
Ac = [       0.*ones(1,4)      ;
    Lv(:,1)  (Aw-Bw*Klqr-Lv*Cw)];
Bc1 = [   1          0.;
    0.*Lv(:,1)  Lv(:,2)];
Bc2 = [ -1;
    -1;
    0.*ones(2,1)];
Cc = [0. -Klqr];
Dc1 = [ 0. 0.];
Dc2 = 0.;

disp('   ')
disp('LQG Compensator')
sys_c_zpk = ss(Ac,Bc1,Cc,Dc1);
zpk(sys_c_zpk)
disp('   ')


disp('OBLTR Controller');
Ac
Bc1
Bc2
Cc
Dc1
Dc2
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc);
    (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
Bcl = [       Bp*Z*Dc2;
    (Bc2+Bc1*Dp*Z*Dc2)];
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl = Dp*Z*Dc2;
Ccl_Az = Ccl(1,:);
Dcl_Az =Dcl(1,:);


disp('Closed Loop System Eigenvalues');
damp(Acl)
disp('RSLQR A-BK Eigenvalues');
damp(Aw-Bw*Klqr)
disp('Observer A-LC');
damp(Aw-Lv*Cw)
disp(' Compensator A-LC-BK Eigenvalues');
damp(Aw-Bw*Klqr-Lv*Cw)

% Step Response Time Domain Analysis - this confirms that the system is
% connected properly and is stable
[y_obltr,~] = step(Acl,Bcl,Ccl_Az,Dcl_Az,1,t);
sys_obltr   = ss(Acl,Bcl,Ccl_Az,Dcl_Az);

%plot step response
figure('Name','OBLTR Step Sim'),
plot(t,az_lqr,'ok',t,y_obltr,'xb','LineWidth',2);
lqr_min = ['RSLQR Tr Ts = ' num2str(lqr_metric(5)) ' ' num2str(lqr_metric(6))];
legend('RSLQR',lqr_min,'Location','best')
title('Acceleration Step Response');
xlabel('Time (sec)');ylabel('Az (fps2)');
title('Accel Time histories');
grid;
if(save_plots == 1) saveppt2(plot_file_name); end


%SS model of loop gain Lu at the plant input
Ain = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac];
Bin = [ Bp; Bc1*Dp];
Cin = -[ Dc1*Cp Cc];%change sign for loop gain
Din = -[ Dc1*Dp];
sys_Lu = ss(Ain,Bin,Cin,Din);

%SS model of loop gain L at the plant output
Aout = [ Ap Bp*Cc;  0.*Bc1*Cp Ac];
Bout = [ Bp*Dc1; Bc1];
Cout = -[ Cp Dp*Cc];%change sign for loop gain
Dout = -[ Dp*Dc1];
sys_Ly = ss(Aout,Bout,Cout,Dout);

Lu_obltr = freqresp(sys_Lu,w);
Ly_obltr = freqresp(sys_Ly,w);
T_obltr  = freqresp(sys_cl,w);
S_obltr = 1 - T_obltr;

    magdb = 20*log10(abs(squeeze(Lu_obltr)));
    wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
    SRu_obltr = sigma(sys_Lu,w,3);
    SRu_min = min(abs(SRu_obltr));
    RDu_obltr = sigma(sys_Lu,w,2);
    RDu_min = min(abs(RDu_obltr));
    
%
figure('Name','Nyquist Plot at Plant Input'),
plot(xx1,yy1,'r',real(squeeze(Lu_lqr)),imag(squeeze(Lu_lqr)),'k',...
    real(squeeze(Lu_obltr)),imag(squeeze(Lu_obltr)),'b','LineWidth',2);grid
axis([-5 5 -5 5]);
legend('Unit Circle','RSLQR','OBLTR','Location','Best');
xlabel('Re(Lu)')
ylabel('Im(Lu)')
title('Nyquist Plot at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Bode Magnitude at Plant Input'),
semilogx(w,20*log10(abs(squeeze(Lu_lqr))),'k',w,20*log10(abs(squeeze(Lu_obltr))),'b','LineWidth',2);grid
legend('RSLQR','OBLTR','Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag')
title('Bode at Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Return Difference at Plant Input'),
semilogx(w,20*log10(abs(RDu_lqr)),'k',w,20*log10(abs(RDu_obltr)),'b','LineWidth',2);grid
lqr_min = ['RSLQR RD min = ' num2str(lqr_metric(2))];
obltr_min = ['OBLTR RD min = ' num2str(RDu_min)];
legend(lqr_min,obltr_min,'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Stability Robustness at Plant Input'),
semilogx(w,20*log10(abs(SRu_lqr)),'k',w,20*log10(abs(SRu_obltr)),'b','LineWidth',2);grid
lqr_min = ['RSLQR SR min = ' num2str(lqr_metric(3))];
obltr_min = ['OBLTR SR min = ' num2str(SRu_min)];
legend(lqr_min,obltr_min,'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Comp Sensitivity at Output'),
Mag1 = 20*log10(abs(squeeze(T_lqr(1,1,:))));
Mag2 = 20*log10(abs(squeeze(T_obltr(1,1,:))));
semilogx(w,Mag1,'k',w,Mag2,'b','LineWidth',2);grid
v_min_lqr = max(max(abs(squeeze(T_lqr(1,1,:)))));
v_min_obltr = max(max(abs(squeeze(T_obltr(1,1,:)))));
lqr_min = ['RSLQR T max = ' num2str(v_min_lqr)];
obltr_min = ['OBLTR T max = ' num2str(v_min_obltr)];
legend(lqr_min,obltr_min,'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag (dB)')
title('Comp Sensitivity')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Sensitivity at Output'),
Mag1 = 20*log10(abs(squeeze(S_lqr(1,1,:))));
Mag2 = 20*log10(abs(squeeze(S_obltr(1,1,:))));
semilogx(w,Mag1,'k',w,Mag2,'b','LineWidth',2);grid
v_min_lqr = max(max(abs(squeeze(S_lqr(1,1,:)))));
v_min_obltr = max(max(abs(squeeze(S_obltr(1,1,:)))));
lqr_min = ['RSLQR S max = ' num2str(v_min_lqr)];
obltr_min = ['OBLTR S max = ' num2str(v_min_obltr)];
legend(lqr_min,obltr_min,'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag (dB)')
title('Sensitivity')
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
allmargin(sys_Lu)

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

   
   return

   
   