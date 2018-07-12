function [P_HMI_H,P_HMI,Fault_mag,Fault_dir,H_M,L_M,L_pp_M,Y_M]=...
    integrity_monitoring_fault (Phi_M ,H_M ,L_M ,L_pp_M ,A_k ,Y_M , alpha)

% PX : states prediction covarience matrix
% Phi_M : state tranision matrix over the horizon, including the current (concatenated)
% H_M : observation matrix over the horizon, excluding the current (concatenated)
% L_M : Kalman gain over the horizon, excluding the current (concatenated)
% L_p_M : L_prime over the horizon, excluding the current (concatenated)
% L_pp_M : L_prime_prime over the horizon, excluding the current (concatenated)
% A_k : () current
% Y_M : Innovation vector covarience matrix during the horizon

global PX PARAMS
syms Fault_magnitude

% For a cleaner notation
m= PARAMS.m;
M= PARAMS.M;
m_F= PARAMS.m_F;
[lm,n_L]= field_of_view_landmarks();
n= n_L*m_F;
n_M= n*(M+1); % measurments in the PH
nL_M= n_L*(M+1);

% TODO this depends on the miss-association probability
P_H= 1e-3;



% Initializa variables for current time
H_k= zeros(n, m); % Observation matrix at the current time step
h_k= zeros(n, 1); % Expected measurement

% models for the current time
for t= 1:n_L
    idx= ((t-1)*m_F)+1:t*m_F;
    [h_t,H_t]= compute_lm_model(lm(:,t));
    H_k( idx,:)= H_t;
    h_k( idx,:)= h_t;
end

V= kron(eye(n_L),PARAMS.R); % current measurements covarience matrix
Y_k= H_k*PX*H_k' + V; % Current innovations covarience matrix

% Update the innovation vector covarience matrix for the new PH
Y_M(n+1:end,n+1:end)= Y_M(1:end-n,1:end-n);
Y_M(1:n,1:n)= Y_k;


L_k= PX * H_k' / Y_k;
P_Hat= PX - L_k*H_k*PX;
Lk_pp= Phi_M(1:m,1:m) - L_k*H_k* Phi_M(1:m,1:m); % Kalman Gain prime-prime 

% Update the horizon increasing the size (TODO: not change sizes)
L_M= [L_k;L_M];
L_pp_M= [Lk_pp; L_pp_M];
H_M= [H_k; H_M];

if isempty(A_k)
    
    % Create A_k^M for the first time (later it will be done recursively)
    A_k= inf* ones( m, n_M + m );
    A_k(: , 1:n)= L_k;
    for i= 1:M
        if i == 1
            Dummy_Variable= L_pp_M(1:i*m,:);
        else
            Dummy_Variable= Dummy_Variable * L_pp_M( (i-1)*m+1:i*m, : );
        end
        A_k(:, n*i + 1 : n*(i+1) )= Dummy_Variable * L_M( i*m+1 : (i+1)*m, : );
    end
    A_k( :,n_M+1 : end )= Dummy_Variable * L_pp_M( M*m + 1 : (M+1)*m, : );    
    
else % The previous A_(k-1) is an input to the function -> A_k is computed recursively
    
    % Calculating A_k recursively
    A_k= [L_k, Lk_pp * A_k(:,1:n*M), Lk_pp * A_k(:,end-m+1:end) / L_pp_M(:,end-m+1:end)];
end

% Augmented B
B_bar= inf* ones( n_M , n_M+m );
A_prev= Lk_pp \ A_k( : , n + 1:end );
B_bar(1:n,:)= [eye(n), -H_k*Phi_M(1:m,:)*A_prev];

% Recursive computation of B
for i= 1:M
    A_prev= L_pp_M( i*m+1:i*m + m , :) \ A_prev(:,(n)+1:end);
    B= [eye(n), -H_M(n*i+1:n*(i+1),:)*Phi_M(m*i+1:m*(i+1),:)*A_prev];
    B_bar(i*n+1:(i+1)*n,1:n*i)= 0;
    B_bar(i*n+1:(i+1)*n, n*i+1:end)= B;
end
M_k= B_bar' / Y_M * B_bar;

% Removing the last element (TODO: do not change size)
H_M= H_M(1:end-(n),:);
L_M= L_M(1:end-m,:);
L_pp_M= L_pp_M(1:end-m,:);


% Initialization of fault slope for every hypothesis (concatenated)
% Organization:
% (only previous state faults ;;
% both previous state faults, and the whole combination of single measurement faults

n_H=  nL_M + 1;
Fault_dir= inf* ones( (n_M+m) * n_H, 1 );
Fault_mag= inf* ones( n_H, 1 );

% Detector threshold including the PH
T_D= chi2inv(1 - PARAMS.C_REQ, n_M);

% Initialization of P(HMI | H)
P_HMI_H= inf* ones( n_H, 1 );

%% Fault-free case done separately

% E matrix for only previous state faults
E= zeros( m, n_M+m );
E(:, end-m+1:end)= eye(m);

Fault_dir( 1 : n_M+m )= E' /(E*M_k*E') * E * A_k' * alpha;
Fault_dir( 1 : n_M+m )= Fault_dir( 1 : n_M+m ) / norm( Fault_dir( 1 : n_M+m ) );


sigma2_hat= alpha'*P_Hat*alpha;
fx_hat_dir= alpha'*A_k*Fault_dir( 1 : n_M+m );
M_dir= Fault_dir( 1 : n_M+m )' * M_k * Fault_dir( 1 : n_M+m );
[Fault_mag(1), P_HMI_H(1)]= fminbnd( @(Fault_magnitude)...
    -( (1-cdf('Normal',PARAMS.alert_limit, Fault_magnitude*fx_hat_dir, sigma2_hat) +...
    cdf('Normal',-PARAMS.alert_limit,Fault_magnitude* fx_hat_dir, sigma2_hat)) *...
    cdf('Noncentral Chi-square', T_D, n_M, Fault_magnitude^2* M_dir) ), -10, 10);
P_HMI_H(1)= -P_HMI_H(1);

% Initialize integrity
P_HMI= P_HMI_H(1) * P_H;

%% Loop over hypotheses with faulted measurements in the PH (only 1 fault)
for i= 1:nL_M 
    
    
    % E matrix
    E= zeros( m + m_F , n_M + m );
    E( end-m+1 : end , end-m+1:end )= eye(m); % previous bias
    E( 1:m_F , (i-1)*m_F+1 : i*m_F )= eye(m_F); % landmark i faulted
    
    % Worst-case fault direction
    Fault_dir( (n_M+m)*i+1 : (i+1)*(n_M+m) )= E' / (E*M_k*E') * E * A_k' * alpha;
    Fault_dir( (n_M+m)*i+1 : (i+1)*(n_M+m) )= ...
        Fault_dir( (n_M+m)*i+1 : (i+1)*(n_M+m) ) /norm( (n_M+m)*i+1 : (i+1)*(n_M+m) );
    
    % worst-case fault magnitude
    sigma2_hat= alpha'*P_Hat*alpha;
    fx_hat_dir= alpha'*A_k*Fault_dir( (n_M+m)*i+1 : (i+1)*(n_M+m) );
    M_dir= Fault_dir( (n_M+m)*i+1 : (i+1)*(n_M+m) )'*M_k*Fault_dir((n_M+m)*i+1 : (i+1)*(n_M+m));

    [Fault_mag(i+1), P_HMI_H(i+1)]= fminbnd( @(Fault_magnitude)...
        -((1-cdf('Normal',PARAMS.alert_limit, Fault_magnitude*...
        fx_hat_dir, sigma2_hat)+...
        cdf('Normal',-PARAMS.alert_limit,Fault_magnitude * fx_hat_dir, sigma2_hat))...
        * cdf('Noncentral Chi-square',T_D,m_F*nL_M,...
        Fault_magnitude^2 * M_dir )), -10, 10);
    P_HMI_H(i+1)= -P_HMI_H(i+1);
    
    % Add P(HMI | H) to the integrity risk
    P_HMI= P_HMI + P_HMI_H(i+1) * P_H;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,H]= compute_lm_model(lm)

global XX

dx= lm(1) - XX(1);
dy= lm(2) - XX(2);
d2= dx^2 + dy^2;
d= sqrt(d2);

% calculate h
h= [d;
    pi_to_pi( atan2(dy,dx) - XX(3) )];

% calculate H
H = [-dx/d, -dy/d,  0;
    dy/d2, -dx/d2, -1];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lm,nL]= field_of_view_landmarks ()

global XX PX LM PARAMS

% % calculate how much we need to include in the EFOV
%lambda_FV= sqrt( max(eig(Px(1:2,1:2))) );
%EFV= -sqrt(2) * lambda_FV * norminv(PARAMS.I_FOV/2,0,1);
% EFOV= sqrt(2) * lambda_FOV * sqrt( chi2inv(1 - PARAMS.I_FOV,1) ); % same as previous

% Get all visible landmarks, assuming no mis-extractions here
%idf= get_visible_landmarks(xx,PARAMS.maxRange+EFV, 0);

%lm= LM(:,idf);
%nL= length(idf);

% Assuming that LIDAR range is infinte, and there are no misextractions
lm= LM;
nL= size(LM,2);
