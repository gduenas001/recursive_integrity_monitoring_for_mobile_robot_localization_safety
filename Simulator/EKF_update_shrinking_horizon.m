function [T_RB,gamma] = EKF_update_shrinking_horizon(z, idf, step)

global XX PX PARAMS hlm Hlm DATA

% Define some variables
n_L= DATA.numAssoc(step);
n= n_L * PARAMS.m_F;

% Update detector threshold
T_RB= chi2inv(1 - PARAMS.C_REQ, n*(PARAMS.Preceding_Horizon+1));

% If no lms are associated --> return!
if n == 0
    q_RB= 0;
    T_RB= 0;
    return;
end

% Remove non-associated msmts
z(:, idf == 0)= [];
idf( idf== 0) = [];

% make R a block diagonal with higher dimension
R= kron(eye(n_L), PARAMS.R);

% create the models for the association
h= zeros(PARAMS.m_F * n_L,1);
H= zeros(PARAMS.m_F * n_L,3);
for i= 1:n_L
    ind= i*PARAMS.m_F - 1;
    h(ind:ind+1)= hlm{idf(i)};
    H(ind:ind+1,:)= Hlm{idf(i)};
end

% Compute innovations
gamma= z(:) - h;
gamma(2:2:end)= pi_to_pi(gamma(2:2:end));
Y= H*PX*H' + R;

% Save previous estimate
XX_bar= XX;
PX_bar= PX;

% Update the estimate
K= PX*H'/Y;
PX= PX - K*H*PX;
XX= XX + K*gamma;