%
% Stability Analysis at plant Input and Output
% K.A. Wise  29 September 2015
%
% Modified:
%


clear all
close all

format short e

disp('****************************** Program Start ****************');
plot_file_name = 'stab_analysis1.ppt';
save_plots = 0; % Flag to bypass saving

rtd  = 180/pi;
dt = 0.002;
t = 0:dt:2;
w = logspace(-1,3,500);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);


%*******************************************
% Define aircraft Plant Model (Ap,Bp,Cp,Dp)
% xp = [ AOA q ]^T as states, u = dele the control
% Outputs are y = [ Az  q]^T.
%*******************************************
% Public release airframe parameters
ZapV  =   -1.05273;
ZdpV  =   -0.0343;
Ma    =   -2.3294;
Mq    =   -1.03341;
Md    =   -1.1684;
V_fps     =   329.127;

% Define the nominal plant matrices
Ap = [ZapV    1;
       Ma    Mq];
Bp = [ ZdpV
        Md];
Cp = [V_fps*ZapV 0. % Az
          1.   0.   % AOA
          0.   1.]; % q
Dp = [V_fps*ZdpV
          0.
          0.];
pstatnam = {'AOA-P rad' 'q-P rps'};
poutnam = {'Az_P ftps2' 'AOA-P rad' 'q-P rps' };
pinnam = {'dele-cmd-P rad' };

%*******************************************
% Design the RSLQR using 3rd order Aw and Bw
%*******************************************
% State Feedback Design
% states are [ int-err AOA q ]
Aw = [ 0. V_fps*ZapV  0. ;
    0.   ZapV  1. ;
    0.   Ma    Mq ];
% Input is dele-cmd
Bw = [ V_fps*ZdpV ;
    ZdpV;
    Md];

%LQR penalty matrix Q
Q_LQR = 0.*eye(3);
Q_LQR(1,1) = 0.2448;
%Q_LQR(1,1) = 0.1;
R_LQR = 1.;
% Solve for RSLQR state feedback gains
[Kx_lqr,Pw]=lqr(Aw,Bw,Q_LQR,R_LQR);

% LQR solution check
lqr_eq = Pw*Aw + Aw'*Pw + Q_LQR - Pw*Bw*inv(R_LQR)*Bw'*Pw;
norm_lqr_eq = norm(lqr_eq);
disp(['LQR ARE Solution Tolerance is: ', num2str(norm_lqr_eq)]);

% LQR State feedback Controller
% states are int-err 
% inputs are y = [ Az AOA q]
%            r = [ Azcmd ]
% Output is u = dele_cmd
Ac =  0.;
Bc1 = [1. 0. 0.];
Bc2 =  -1;
Cc = -Kx_lqr(1);
Dc1 = [0. -Kx_lqr(2:3)];
Dc2 = 0.;

%****************************************************
% Close the loop
% Plant form  xdot = Apx + Bpu;
%                      y = Cpx +Dpu
% Controller xcdot = Acxc + Bc1y + Bc2r
%                      u = Ccxc + Dc1y + Dc2r
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc);
    (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
Bcl = [       Bp*Z*Dc2;
    (Bc2+Bc1*Dp*Z*Dc2)];
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
Ccl_Az = Ccl(1,:);
Dcl_Az =Dcl(1,:);
sys_cl = ss(Acl,Bcl,Ccl,Dcl);

    
% Step Response Time Domain Analysis - this confirms that the system is
% connected properly and is stable
[y_lqr,x_lqr] = step(Acl,Bcl,Ccl_Az,Dcl_Az,1,t);
sys_lqr = ss(Acl,Bcl,Ccl_Az,Dcl_Az);

figure('Name','LQR Step Sim'),
plot(t,y_lqr,'LineWidth',2);grid
legend('Az','Location','Best');
xlabel('Time (sec)')
ylabel('Az fps2')
if(save_plots == 1) saveppt2(plot_file_name); end

%SS model of loop gain Lu at the plant input
Ain = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac];
Bin = [ Bp; Bc1*Dp];
Cin = -[ Dc1*Cp Cc];%change sign for loop gain
Din = -[ Dc1*Dp];
sys_u = ss(Ain,Bin,Cin,Din);

%SS model of loop gain L at the plant output
Aout = [ Ap Bp*Cc;  0.*Bc1*Cp Ac];
Bout = [ Bp*Dc1; Bc1];
Cout = -[ Cp Dp*Cc];%change sign for loop gain
Dout = -[ Dp*Dc1];
sys_y = ss(Aout,Bout,Cout,Dout);

Lu = freqresp(sys_u,w);
Ly = freqresp(sys_y,w);
T  = freqresp(sys_cl,w);
S = 1 - T;

[nCp,nAp] = size(Cp);
[~,nBp]   = size(Bp);

% Plant Input Freqeuncy Domain Analysis
for i=1:numel(w),
    s = sqrt(-1)*w(i);
    GG = Cp*inv(s*eye(size(Ap))-Ap)*Bp+Dp;
    KK = Cc*inv(s*eye(size(Ac))-Ac)*Bc1+Dc1;
    Lu_lqr(i)  = -KK*GG;
    RDu_lqr(i)  = 1. + Lu_lqr(i);
    SRu_lqr(i) = 1. + 1./Lu_lqr(i);
    Lin = Cin*inv(s*eye(size(Ain))-Ain)*Bin+Din;
    Lout = Cout*inv(s*eye(size(Aout))-Aout)*Bout+Dout;
    % This loops computes SISO freq response at plant input one loop open,
    % rest closed
    for jj = 1:nBp,
        Fyjj = eye(size(Lin));
        Fyjj(jj,jj) = 0.;
        Tujj(jj,i) = inv(eye(size(Lin)) + Lin*Fyjj)*Lin;
    end
    % This loops computes SISO freq response at plant output one loop open,
    % rest closed
    for jj = 1:nCp,
        Fyjj = eye(size(Lout));
        Fyjj(jj,jj) = 0.;
        Tyjj = inv(eye(size(Lout)) + Lout*Fyjj)*Lout;
        Tyj(jj,i) = Tyjj(jj,jj);
    end
    Sout = inv(eye(size(Lout))+Lout);
    Tout = Sout*Lout;
    Tout2 = Lout*Sout;
    sens_lqr(i) = max(Sout(1,1));
    compsens_lqr11(i) = max(Tout(1,1));
    compsens_lqr12(i) = max(Tout2(1,1));
end

% SISO loop gains at output
for k = 1:nCp
% Close other loops for SISO margins
idx_closed = setxor(k,1:nCp);
Ly_k = feedback(sys_y,eye(size(sys_y,1)-1),idx_closed,idx_closed);
[output_GM_design(k), output_PM_design(k), ~ , output_wcp_design(k)] = margin(Ly_k(k,k));
end
output_GM     = output_GM_design;
output_GM_design     = 20*log10(output_GM_design);
     
%
figure('Name','Nyquist Plot at Plant Input'),
plot(xx1,yy1,'r',real(squeeze(Lu)),imag(squeeze(Lu)),'k',real(Lu_lqr),imag(Lu_lqr),'r--',...
real(Tujj),imag(Tujj),'g--','LineWidth',2);grid
axis([-2 2 -2 2]);
legend('Unit Circle','Loop','Short Form','Location','Best');
xlabel('Re(L)')
ylabel('Im(L)')
title('Nyquist Plot at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Nyquist Plot at Plant Output'),
hold on
for jj = 1:nCp,
plot(xx1,yy1,'r',real(Tyj(jj,:)),imag(Tyj(jj,:)),'b--','LineWidth',2);grid
axis([-2 2 -2 2]);
xlabel('Re(L)')
ylabel('Im(L)')
title('Nyquist Tyj at Plant Output')
end
hold off
if(save_plots == 1) saveppt2(plot_file_name); end

for jj =1:nCp,
    figure('Name','Nyquist Plot at Plant Output'),
    plot(xx1,yy1,'r',real(Tyj(jj,:)),imag(Tyj(jj,:)),'r--','LineWidth',2);grid
    legend([num2str(jj) ' Nyquist of Tyj'],'Location','Best');
    axis([-2 2 -2 2]);
    xlabel('Re(L)')
    ylabel('Im(L)')
    title('Nyquist Tyj at Plant Output')
if(save_plots == 1) saveppt2(plot_file_name); end
end

for jj=1:nCp,
    figure('Name','RD at Plant Output'),
    semilogx(w,20*log10(abs(1.+Tyj(jj,:))),'b--','LineWidth',2);grid
    v_min = min(min(abs(1.+Tyj(jj,:))));
    legend([num2str(jj) ' min(1.+Tyj) = ' num2str(v_min)],'Location','Best');
    xlabel('Frequency (rps)')
    ylabel('Mag')
        title('RD 1.+Tyj at Plant Output')
if(save_plots == 1) saveppt2(plot_file_name); end
end



figure('Name','Bode Magnitude at Plant Input'),
semilogx(w,20*log10(abs(Lu_lqr)),'k',w,20*log10(abs(squeeze(Lu))),'r--','LineWidth',2);grid
legend('Loop','Short Form','Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag')
title('Bode at Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Return Difference at Plant Input'),
semilogx(w,20*log10(abs(RDu_lqr)),'r--','LineWidth',2);grid
RDu_min = min(min(abs(RDu_lqr)));
legend([num2str(jj) ' min(I+Lu) = ' num2str(RDu_min)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Stability Robustness at Plant Input'),
semilogx(w,20*log10(abs(SRu_lqr)),'r--','LineWidth',2);grid
SRu_min = min(min(abs(SRu_lqr)));
legend([num2str(jj) ' min(I+invLu) = ' num2str(SRu_min)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

for jj = 1:nCp,
    figure('Name','Comp Sensitivity at Output'),
    semilogx(w,20*log10(abs(squeeze(T(jj,1,:)))),'b','LineWidth',2);grid
    v_min = max(max(abs(squeeze(T(jj,1,:)))));
    legend([num2str(jj) ' T max = ' num2str(v_min)],'Location','Best');
    xlabel('Frequency (rps)')
    ylabel('Mag (dB)')
    title('Comp Sensitivity')
if(save_plots == 1) saveppt2(plot_file_name); end
end

for jj = 1:nCp,
    figure('Name','Sensitivity at Output'),
    semilogx(w,20*log10(abs(squeeze(S(jj,1,:)))),'b','LineWidth',2);grid
    v_min = max(max(abs(squeeze(S(jj,1,:)))));
    legend([num2str(jj) ' S max = ' num2str(v_min)],'Location','Best');
    xlabel('Frequency (rps)')
    ylabel('Mag (dB)')
    title('Sensitivity')
if(save_plots == 1) saveppt2(plot_file_name); end
end

% Return Difference at plant output
sv_rd_y = sigma(sys_y,w,2);
min_sv_rd_y = sv_rd_y(nCp,:);
min_rd_y = min(abs(min_sv_rd_y));

figure('Name','min SV of RDM at Plant Output'),
semilogx(w,20*log10(abs(min_sv_rd_y)),'k','LineWidth',2);grid
legend(['min(sv(I+Ly)) = ' num2str(min_rd_y)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

disp('Classical Margins')
allmargin(sys_u)

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

margin(sys_u)



