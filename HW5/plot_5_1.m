dt = 0.002;
t = 0:dt:2;

load('5_1_az.mat');
load('5_1_q.mat');
load('5_1_del.mat');
load('5_1_deldot.mat');
load('rt_5_1.mat');
load('st_5_1.mat');

load('hw4_az.mat');
load('hw4_q.mat');
load('hw4_del.mat');
load('hw4_deldot.mat');
load('../HW5/hw4_rt.mat');
load('../HW5/hw4_st.mat');

figure; plot(t,az_hw4,'ok', t,az_5_1,'xb'); 
title('az (acceleration)');grid on;xlabel('time (sec)');
str1 = ['No Actuator, Tr = ' num2str(hw4_rt) ' Ts = ' num2str(hw4_st)];
str2 = ['With Actuator, Tr = ' num2str(rt_5_1) ' Ts = ' num2str(st_5_1)];
legend({str1,str2}, 'Position',[0.39,0.45,0.25,0.1]);

figure; plot(t,q_hw4,'ok', t,q_5_1,'xb'); 
title('Pitch Rate (deg/sec)');grid on;xlabel('time (sec)');
legend('No Actuator','With Actuator');

figure; plot(t,del_hw4,'ok', t,del_5_1,'xb'); 
title('Elevon (deg)');grid on;xlabel('time (sec)'); 
legend('No Actuator','With Actuator');

figure; plot(t,deldot_hw4,'ok', t,deldot_5_1,'xb'); 
title('Elevon rate (deg/sec)');grid on;xlabel('time (sec)');
legend('No Actuator','With Actuator');

