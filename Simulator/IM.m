function [P_HMI, H_M_cell, Y_M, A_k, L_M_cell, Lpp_M_cell]=...
    IM (Phi_M , H_M_cell, A_k ,Y_M , alpha, L_M_cell, Lpp_M_cell)

% PX : states prediction covarience matrix
% Phi_M : state tranision matrix over the horizon, including the current (concatenated)
% H_M : observation matrix over the horizon, excluding the current (concatenated)
% L_M : Kalman gain over the horizon, excluding the current (concatenated)
% L_pp_M : L_prime_prime over the horizon, excluding the current (concatenated)
% A_k : previous A_k for recursive computations
% Y_M : Innovation vector covarience matrix during the horizon

global PX PARAMS


% For a cleaner notation
m= PARAMS.m;
M= PARAMS.M;
m_F= PARAMS.m_F;
[lm,n_L]= field_of_view_landmarks();
n= n_L*m_F;
n_M= n*(M+1); % measurments in the PH
nL_M= n_L*(M+1);
n_H=  nL_M + 1; % number of hypotheses

% TODO this depends on the miss-association probability
P_H= 1e-2;

% Initializa variables for current time
Hk= zeros(n, m); % Observation matrix at the current time step
h_k= zeros(n, 1); % Expected measurement

% models for the current time
for t= 1:n_L
    idx= ((t-1)*m_F)+1:t*m_F;
    [h_t,H_t]= compute_lm_model(lm(:,t));
    Hk( idx,:)= H_t;
    h_k( idx,:)= h_t;
end

V= kron(eye(n_L),PARAMS.R); % current measurements covarience matrix
Y_k= Hk*PX*Hk' + V; % Current innovations covarience matrix

% Update the innovation vector covarience matrix for the new PH
Y_M(n+1:end,n+1:end)= Y_M(1:end-n,1:end-n);
Y_M(1:n,1:n)= Y_k;


Lk= PX * Hk' / Y_k;
P_Hat= PX - Lk*Hk*PX;
Lk_pp= Phi_M(1:m,1:m) - Lk*Hk* Phi_M(1:m,1:m); % Kalman Gain prime-prime 

% Add cells
H_M_cell= [ {Hk}, H_M_cell];
L_M_cell= [ {Lk} ,L_M_cell];
Lpp_M_cell= [ {Lk_pp} ,Lpp_M_cell];

if 0%~isempty(A_k)% The previous A_(k-1) is an input to the function -> A_k is computed recursively
    A_k= [Lk, Lk_pp*A_k];
    A_k(:,end-m-n+1:end-m)= [];
    A_k(:,end-m+1:end)= A_k(:,end-m+1:end) / Lpp_M_cell{end};
    
else % Create A_k^M for the first time (later it will be done recursively)
    A_k= inf* ones( m, n_M + m );
    A_k(: , 1:n)= Lk;
    for i= 1:M
        if i == 1
            Dummy_Variable= Lk_pp;
        else
            Dummy_Variable= Dummy_Variable * Lpp_M_cell{i};
        end
        A_k(:, n*i + 1 : n*(i+1) )= Dummy_Variable * L_M_cell{i+1};
    end
    A_k( :,n_M+1 : end )= Dummy_Variable * Lpp_M_cell{ M + 1 };
end

% Augmented B
B_bar= inf* ones( n_M , n_M+m );
A_prev= Lk_pp \ A_k( : , n + 1:end );
B_bar(1:n,:)= [eye(n), -Hk*Phi_M(1:m,:)*A_prev];

% Recursive computation of B
for i= 1:M
    A_prev= Lpp_M_cell{i+1} \ A_prev(:,(n)+1:end);
    B= [eye(n), -H_M_cell{i+1} * Phi_M(m*i+1:m*(i+1),:) * A_prev];
    B_bar(i*n+1:(i+1)*n,1:n*i)= 0;
    B_bar(i*n+1:(i+1)*n, n*i+1:end)= B;
end
M_k= B_bar' / Y_M * B_bar;

% Using Cells
H_M_cell(end)= [];
L_M_cell(end)= [];
Lpp_M_cell(end)= [];


% Detector threshold including the PH
T_D= chi2inv(1 - PARAMS.C_REQ, n_M);


%% Loop over hypotheses in the PH (only 1 fault)
P_HMI= 0;
for i= 1:n_H
    
    if i == 1 % E matrix for only previous state faults
        E= zeros( m, n_M+m );
        E(:, end-m+1:end)= eye(m);
    else % E matrix with faults in the PH
        E= zeros( m + m_F , n_M + m );
        E( end-m+1 : end , end-m+1:end )= eye(m); % previous bias
        E( 1:m_F , (i-2)*m_F+1 : (i-1)*m_F )= eye(m_F); % landmark i faulted 
    end
    
    % Worst-case fault direction
    f_M_dir= E' / (E*M_k*E') * E * A_k' * alpha;
    f_M_dir= f_M_dir / norm(f_M_dir);

    % worst-case fault magnitude
    sigma2_hat= alpha'*P_Hat*alpha;
    fx_hat_dir= alpha' * A_k * f_M_dir;
    M_dir= f_M_dir' * M_k * f_M_dir;

    [~, P_HMI_H]= fminbnd( @(f_M_mag)...
        -((1-cdf('Normal',PARAMS.alert_limit, f_M_mag * fx_hat_dir, sigma2_hat)+...
        cdf('Normal',-PARAMS.alert_limit,f_M_mag * fx_hat_dir, sigma2_hat))...
        * cdf('Noncentral Chi-square',T_D,m_F*nL_M, f_M_mag^2 * M_dir )), -10, 10);
    P_HMI_H= -P_HMI_H;
    
    % Add P(HMI | H) to the integrity risk
    P_HMI= P_HMI + P_HMI_H * P_H;
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
