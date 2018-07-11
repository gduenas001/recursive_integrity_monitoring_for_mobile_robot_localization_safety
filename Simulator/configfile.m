%% This file configures parameters and load/initialize variables

addpath('./utilities');
set(0,'DefaultFigureWindowStyle','normal')
format compact

% Global variables
global XX PX BIAS LM POSES  DATA PARAMS SWITCH IA

% creates POSES & LM
create_path_and_map();

% Number of time epochs to run the simulation
PARAMS.numEpochs= 400;

%% Controls
PARAMS.dt= 0.1;
PARAMS.m= 3;
PARAMS.at_waypoint= 5;
PARAMS.V= 10;
PARAMS.maxRateG= deg2rad(50);
PARAMS.maxG= deg2rad(30);
PARAMS.wheelbase= 1; % metres, vehicle wheel-base
PARAMS.veh= [0 -4    -4   0; 
              0 0.5  -0.5 0]; % vehicle animation

% Initial control inputs
iwp= 1; 
G= deg2rad(0);

% control noises
PARAMS.sigmaV= 0.1; % m/s ( default= 0.3 )
PARAMS.sigmaG= deg2rad(1); % radians ( default= (3.0*pi/180) )
PARAMS.Q= [PARAMS.sigmaV^2, 0; 0, PARAMS.sigmaG^2];


%% Measurements
PARAMS.maxRange= inf;%25
PARAMS.m_F= 2; % d.o.f. of one measurement

% observation noises
PARAMS.sigmaR= 0.3; % metres ( default 0.1 )
PARAMS.sigmaB= deg2rad(2); % radians ( default (1.0*pi/180) )
PARAMS.R= [PARAMS.sigmaR^2 0; 0 PARAMS.sigmaB^2];


%% Integrity 
PARAMS.I_REQ= 1e-5; % Integrity risk requirement
PARAMS.I_T= 0.001; % set threshold for the local NN
PARAMS.I_FOV= 1e-9;
PARAMS.alert_limit= 1; 
PARAMS.T2= chi2inv(1-PARAMS.I_T,PARAMS.m_F); % threshold for the local NN
PARAMS.P_IA_max= PARAMS.I_FOV*4; % for the landmark selection
PARAMS.C_REQ= 1e-5; % continuity risk allocation
PARAMS.Preceding_Horizon= 1; % epochs
PARAMS.P_UA= 10^-3; % assuming that it is constant for the whole landmarks.
PARAMS.P_ME= 0; % Misdetection probability
PARAMS.P_CA= 1; % Prior correct association probability

% switches
SWITCH.control_noise= 1; % if 0, velocity and gamma are perfect
SWITCH.sensor_noise= 1; % if 0, measurements are perfect
SWITCH.inflate_noise= 0; % if 1, the estimated Q and R are inflated (ie, add stabilising noise)
SWITCH.heading_known= 0; % if 1, the vehicle heading is observed directly at each iteration
SWITCH.seed_random= 1; % if not 0, seed the randn() with its value at beginning of simulation (for repeatability)
SWITCH.graphics= 1; % if 0, avoids plotting most animation data to maximise simulation speed
SWITCH.association= 1; % if 0, associations are given; if 1, they are estimated using the LNN
SWITCH.ME= 0; % If 0, no mis-extractions; if 1, there are mis-extractions.
if SWITCH.seed_random, rand('state',SWITCH.seed_random), randn('state',SWITCH.seed_random), end


%% Initializations 
 
% True & estimated state
xtrue= POSES(1,:)';
XX= xtrue; 
BIAS= zeros(PARAMS.m,1);

% Initial pose covariance to very small value
PX= diag( eps* ones(1,PARAMS.m) );

PARAMS.ftag= 1:size(LM,2);     % identifier for each landmark
da_table= zeros(1,size(LM,2)); % data association table 

% more initializations
epoch= 0; IA= 0;
DATA.epsXX= zeros(5000,3);
DATA.stdXX= zeros(5000,3);
DATA.eps= zeros(5000,1);
DATA.stdEps= zeros(5000,1);
DATA.P_HMI= zeros(5000,1);
DATA.path= zeros(5000,2);
DATA.P_MA_k= ones(5000,1);
DATA.P_MA_K= ones(5000,1);
DATA.IA= zeros(5000,1);
DATA.numAssoc= zeros(5000,1);

DATA.gamma= cell(5000,1);
DATA.stngamma= cell(5000,1);
DATA.PCAt= ones(5000,1);
DATA.realPCA= zeros(5000,1);
DATA.calcPCA= zeros(5000,1);
DATA.bias= zeros(5000,3);
DATA.bias_interest= zeros(5000,1);
DATA.T_D= zeros(5000,1);
DATA.q_D= zeros(5000,1);
DATA.lambda2= zeros(5000,1);
DATA.lambda2_current= zeros(5000,1);



