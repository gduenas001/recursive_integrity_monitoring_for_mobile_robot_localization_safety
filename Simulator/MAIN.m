dbstop if error
dbclear if error

clear; close all; configfile;

n_L= size(LM,2); % Number of landmarks
M= PARAMS.Preceding_Horizon;
P_CA= 1; % Prior probability of correct association

h= setup_animations();

% Innovation vector over the horizon
y=ones((PARAMS.Preceding_Horizon+1)*n_L*PARAMS.m_F,1)*inf;
% Innovation vector covarience matrix over the horizon
Y=zeros((PARAMS.Preceding_Horizon+1)*n_L*PARAMS.m_F,(PARAMS.Preceding_Horizon+1)*n_L*PARAMS.m_F);
% State transition matrix over the horizon
Phi_M=ones((PARAMS.Preceding_Horizon+1)*PARAMS.m,PARAMS.m)*inf;
% Measurement Jacobian matrix over the horizon
H_M=ones((PARAMS.Preceding_Horizon)*n_L*PARAMS.m_F,PARAMS.m)*inf;
% Kalman gain matrix over the horizon
L_M=ones((PARAMS.Preceding_Horizon)*PARAMS.m,n_L*PARAMS.m_F)*inf;
L_p_M=ones((PARAMS.Preceding_Horizon)*PARAMS.m,PARAMS.m)*inf;
L_pp_M=ones((PARAMS.Preceding_Horizon)*PARAMS.m,PARAMS.m)*inf;

% *****************    MAIN LOOP    *****************
for epoch= 2:PARAMS.numEpochs
    disp(['Step: ',num2str(epoch)]);
    
    if epoch <= PARAMS.Preceding_Horizon+1 % avoid the use of shrinking horizon
        % Compute controls
        % [G,iwp]= compute_steering(xtrue, iwp, G);
        % EKF predict step
        [xtrue,XX,PX,Gx,alpha]= predict (xtrue,XX,PX,-deg2rad(3));%G); % Gx is the current state transition matrix
        
        H_k=zeros(n_L*PARAMS.m_F, PARAMS.m); % Observation matrix at the current time step
        h_k=zeros(n_L*PARAMS.m_F, 1); % Expected measurement
        
        for t= 1:n_L
            [h_t,H_t]= compute_lm_model(LM(:,t));
            H_k(((t-1)*PARAMS.m_F)+1:t*PARAMS.m_F,:)= H_t;
            h_k(((t-1)*PARAMS.m_F)+1:t*PARAMS.m_F,:)= h_t;
        end
        
        V= kron(eye(n_L),PARAMS.R); % current measurements covarience matrix
        Y_k= H_k * PX * H_k' + V; % Current IV covarience matrix
        Lk= PX * H_k' / Y_k; % Kalman gain
        Lk_p= eye(PARAMS.m) - Lk*H_k; % current ()
        Lk_pp= Lk_p * Gx; % current ()
        
        % Get measurements
        [z,idft]= get_observations(xtrue);
        z= add_observation_noise(z);
        
        % DA
        [idf,DATA.numAssoc(epoch)]= data_associate_LNN_LS(z);
        
        % Store associations data
        DATA.IA(epoch)= any( (idft - idf).*idf );
        
        [T_D,gamma] = EKF_update_shrinking_horizon(z, idf, epoch);
        
        % Store DATA
        store_data(epoch, 1, 1, xtrue, 0, T_D);
        
        % Plots
        do_plots(xtrue, z, h, epoch);
        
        % Update horizon matrices
        Y((PARAMS.Preceding_Horizon-epoch+1)*n_L*PARAMS.m_F+1:(PARAMS.Preceding_Horizon-epoch+2)*n_L*PARAMS.m_F,(PARAMS.Preceding_Horizon-epoch+1)*n_L*PARAMS.m_F+1:(PARAMS.Preceding_Horizon-epoch+2)*n_L*PARAMS.m_F)=Y_k;
        y((PARAMS.Preceding_Horizon-epoch+1)*n_L*PARAMS.m_F+1:(PARAMS.Preceding_Horizon-epoch+2)*n_L*PARAMS.m_F,:)=gamma;
        Phi_M((PARAMS.Preceding_Horizon-epoch+1)*PARAMS.m+1:(PARAMS.Preceding_Horizon-epoch+2)*PARAMS.m,:)=Gx;
        H_M((PARAMS.Preceding_Horizon-epoch+1)*n_L*PARAMS.m_F+1:(PARAMS.Preceding_Horizon-epoch+2)*n_L*PARAMS.m_F,:)=H_k;
        L_M((PARAMS.Preceding_Horizon-epoch+1)*PARAMS.m+1:(PARAMS.Preceding_Horizon-epoch+2)*PARAMS.m,:)=Lk;
        L_p_M((PARAMS.Preceding_Horizon-epoch+1)*PARAMS.m+1:(PARAMS.Preceding_Horizon-epoch+2)*PARAMS.m,:)=Lk_p;
        L_pp_M((PARAMS.Preceding_Horizon-epoch+1)*PARAMS.m+1:(PARAMS.Preceding_Horizon-epoch+2)*PARAMS.m,:)=Lk_pp;
    else
        % Compute controls
        [G,iwp]= compute_steering(xtrue, iwp, G); G= -deg2rad(3);
        
        % EKF predict step
        [xtrue,XX,PX,Gx,alpha]= predict (xtrue,XX,PX,G); % Gx is the current state transition matrix
        
        % update the states transition matrix over the horizon
        Phi_M(PARAMS.m+1:end,:)=Phi_M(1:end-PARAMS.m,:);
        Phi_M(1:PARAMS.m,:)=Gx;
        
        % Integrity Monitoring
        [P_HMI_H,P_CA,P_HMI_worst,F_mag,Fault_slope,H_M,L_M,L_p_M,L_pp_M,Y]= integrity_monitoring_fault...
            (PX, Phi_M ,H_M ,L_M ,L_p_M ,L_pp_M ,[] ,Y , P_CA, alpha);
        
        % Get measurements
        [z,idft]= get_observations(xtrue);
        z= add_observation_noise(z);
        
        % DA
        if  ~isempty(z)
            [idf,DATA.numAssoc(epoch)]= data_associate_LNN_LS(z);
            
            % Store associations data
            DATA.IA(epoch)= any( (idft - idf).*idf );
        else
            idf= [];
            DATA.numAssoc(epoch)= 0;
            DATA.IA(epoch)= 0;
        end
        
        [q_RB, T_D, gamma,y]= EKF_update(z, idf, epoch,y,Y);
        
        % Store DATA
        store_data(epoch, P_HMI_worst, P_CA, xtrue, q_RB, T_D);
        
        % Plots
        do_plots(xtrue, z, h, epoch);
        
    end
end
% *****************  END OF MAIN LOOP    *****************

% post_processing_and_plots(step)