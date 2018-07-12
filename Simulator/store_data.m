function store_data (step, P_HMI_noMA_UA, xtrue, q_D, T_D)

global DATA XX PX

alpha= [-sin(XX(3)); cos(XX(3)); 0];
sig_hat= sqrt(alpha'*PX*alpha);

% Detector
DATA.q_D(step)= q_D;
DATA.T_D(step)= T_D;

% HMI -- with UA
DATA.P_HMI(step)= P_HMI_noMA_UA;

% Error
DATA.epsXX(step,:)= abs(xtrue - XX)';
DATA.stdXX(step,:)= 3*sqrt(diag(PX)');

% Error state of interest
DATA.eps(step)= (alpha'* (xtrue - XX) );
DATA.stdEps(step)= 3*sig_hat;

% Path
DATA.path(step,:)= XX(1:2)';

