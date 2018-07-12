function [q_D,T_D,gamma] = EKF_update_shrinking_horizon(z, idf, step)

global XX PX PARAMS hlm Hlm DATA

% Define some variables
n_L= DATA.numAssoc(step);
n= n_L * PARAMS.m_F;

% Update detector threshold
T_D= chi2inv(1 - PARAMS.C_REQ, n*(PARAMS.M+1));

% If no lms are associated --> return!
if n == 0
    q_D= 0;
    T_D= 0;
    return;
end

% Remove non-associated msmts
z(:, idf == 0)= [];
idf( idf== 0) = [];

% make R a block diagonal with higher dimension
R= kron(eye(n_L), PARAMS.R);

% create the models for the association
h= zeros(PARAMS.m_F * n_L, 1);
H= zeros(PARAMS.m_F * n_L, 3);
for i= 1:n_L
    ind= i*PARAMS.m_F - 1;
    h(ind:ind+1)= hlm{idf(i)};
    H(ind:ind+1,:)= Hlm{idf(i)};
end

% Compute innovations
gamma= z(:) - h;
gamma(2:2:end)= pi_to_pi(gamma(2:2:end));
Y= H*PX*H' + R;

% Update the estimate
L= PX*H'/Y;
PX= PX - L*H*PX;
XX= XX + L*gamma;

% Detector
q_D= gamma'/Y*gamma;