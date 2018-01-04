clear all
close all 
clear global
clc


%% Define some constants
d=1;
%options for bodeplots (in Hz)
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';

%step_length = 30; 
Ts = 0.001;
k = [0:Ts:49];
N = length(k); %number of timesteps

zero = [0];
poles = [0.9 0.9];
K = 1;
sys_d = zpk(zero, poles, K, Ts);
figure;impulse(sys_d);
%p1 is inderdaad gelijk aan 1, zoals gezegd in Bristow
% ==> Asymptotic stability


%%%% Learning algorithm
Q = 1;
L = 0.5;

%% Task 
%% Define reference traject
nr_front = 5;
nr = 30;
nr_tot = 45;
ref = polytraj(1, Ts, nr, nr_tot);
ref = [zeros(nr_front,1); ref];
t = (0:1:(nr_front + nr_tot)-1)'*Ts;

figure(4);plot(t, ref);

%Define input
u0 = zeros(length(ref),1);

%Initial state of system
y0 = lsim(sys_d, u0, t);

%shift needed for the for loop which contains the ILC
y_new = y0;
up = u0;

%first error 
clear e;
e(:,1) = ref - y_new;

i=1;
while (norm(e(:,i))>1e-11 && i<500)
    %shift error with d (due to plant delay)
    e_shift = e(:,i);
    e_shift = [e_shift((1+d):end);e_shift(end)]; %put e_shift(end) as last element (better than just 0)
    
    %u_new = up + kp*e_shift;
    u_new = up + L*e_shift; %
   
    y_new = lsim(sys_d, u_new, t); %Same model as plant ==> perfect learning filter
%   y_new = lsim(Pp_d, u_new, t); %Without Q-filter and with plant mismodelling

%     figure(5);
%     subplot(3,1,1); plot(t, ref, t, y_new);    legend('ref', 'y_new');
%     subplot(3,1,2); plot(t, ref-y_new); title('plot of the error')
%     subplot(3,1,3); plot(t, up, t, u_new); title('plot of new and previous error'); legend('up', 'u_new');

    
    e(:,i+1) = ref - y_new;
    up = u_new;

    i = i+1;
end 


%Plot 2-norm of error signals
[~, c] = size(e);
for i=1:c
    e_norm(i) = norm(e(:,i));
end
figure(5);semilogy(1:1:i, e_norm, 'b+');title('logarithmic plot of 2-norm of error per trial');

%figure of error
figure(6); plot(t, ref-y_new); title('Tracking error of first iteration'); xlabel('time([s])');

%plot last output
figure(7); plot(t, ref, t, y_new); title('plot of the output of last trial'); legend('ref', 'y_new')

%% Check for AS

%% check for asymptotic stability

%Asymptotic stability condition: rho(Q*(I-L*P))<1 (with rho = spectral radius)
%Determine impulseresponses and take q0, l0 and p1
impP = impulse(sys_d);
p1 = impP(2)

Q = tf(1,1,Ts);
impQ = impulse(Q);
q0 = impQ(1)

L = tf(L, 1, Ts);
impL = impulse(L);
l0 = impL(1)

AS_crit = abs(q0*(1-l0*p1))


%% Check for transient behaviours

%No Q-filter

figure(8);bode((1-series(series(tf([1 0], [0 1],Ts), L), sys_d)), tf(1,1, Ts), opts); legend('1-zLP', '1/Q');
[mag, phase, wout] = bode(series(tf(1,1,Ts), (1-series(series(tf([1 0], [0 1],Ts), L), sys_d))));
MaxMonConv = max(mag)
if (max(mag)<1)
    disp('ILC should converge monotonically')
else 
    disp('Monotonic convergence criterium ILC not fulfilled')
end


return;