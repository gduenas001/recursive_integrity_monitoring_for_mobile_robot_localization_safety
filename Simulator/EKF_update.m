function [q_D, T_D, gamma,gamma_M]= EKF_update(z, idf, step,gamma_M,Y_M)

global XX PX PARAMS hlm Hlm DATA

% Define some variables
n_L= DATA.numAssoc(step);
n= n_L * PARAMS.m_F;

% If no lms are associated --> return!
if n == 0
    q_D= 0;
    T_D= 0;
    return;
end

% Update detector threshold
T_D= chi2inv(1 - PARAMS.C_REQ, n*(PARAMS.M+1));

% Remove non-associated msmts
z(:, idf == 0)= [];
idf( idf== 0) = [];

% make R a block diagonal with higher dimension
R= kron(eye(n_L), PARAMS.R);

% create the models for the association
h= zeros(n,1);
H= zeros(n,3);
for i= 1:n_L
    ind= i*PARAMS.m_F - 1;
    h(ind:ind+1)= hlm{idf(i)};
    H(ind:ind+1,:)= Hlm{idf(i)};
end

% Compute innovations
gamma= z(:) - h;
gamma(2:2:end)= pi_to_pi(gamma(2:2:end));
Y_k= H*PX*H' + R;

% Update the estimate
K= PX*H'/Y_k;
PX= PX - K*H*PX;
XX= XX + K*gamma;

gamma_M( n+1:end )= gamma_M( 1:end-n );
gamma_M(1:n)= gamma;

% Detector
q_D= gamma_M' / Y_M * gamma_M;
