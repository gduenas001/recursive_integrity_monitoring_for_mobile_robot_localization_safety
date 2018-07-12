dbstop if error
dbclear if error

clear; close all; configfile;

h= setup_animations();


% *****************    MAIN LOOP    *****************
for epoch= 2:PARAMS.numEpochs
    disp(['Step: ',num2str(epoch)]);
    
    %% Growing PH --> store matrices until full-size
    if epoch <= PARAMS.M + 1
        % Compute controls
        [G,iwp]= compute_steering(xtrue, iwp, G); G= -deg2rad(3);
        
        % EKF predict step
        [xtrue,XX,PX,Phi_k,alpha]= predict (xtrue,XX,PX,G); % Gx is the current state transition matrix
        
        h_k= zeros(n_L*PARAMS.m_F, 1); % Expected measurement
        H_k= zeros(n_L*PARAMS.m_F, PARAMS.m); % Observation matrix at the current time step
        
        for t= 1:n_L
            [h_t,H_t]= compute_lm_model( LM(:,t) );
            idx= ( (t-1)*PARAMS.m_F ) + 1 : t*PARAMS.m_F;
            H_k( idx , : )= H_t;
            h_k( idx , : )= h_t;
        end
        
        V= kron(eye(n_L),PARAMS.R);  % current measurements covarience matrix
        Y_k= H_k * PX * H_k' + V;    % Current innovation covarience matrix
        Lk= PX * H_k' / Y_k;         % Kalman gain
        Lk_pp= Phi_k - Lk*H_k*Phi_k; % Kalman gain prime-prime
        
        
        % Get measurements
        [z,idft]= get_observations(xtrue);
        z= add_observation_noise(z);
        
        % DA
        [idf,DATA.numAssoc(epoch)]= data_associate_LNN_LS(z);
        
        % Store associations data
        DATA.IA(epoch)= any( (idft - idf).*idf );
        
        
        [q_D,T_D,gamma] = EKF_update_shrinking_horizon(z, idf, epoch);
        
        % Store DATA
        store_data(epoch, 0, xtrue, q_D, T_D);
        
        % Plots
%         do_plots(xtrue, z, h, epoch);
        
        % Update horizon matrices
        idx= (PARAMS.M+1-epoch)*n_L*PARAMS.m_F+1 : (PARAMS.M+2-epoch)*n_L*PARAMS.m_F;
        Y_M( idx, idx)= Y_k;
        gamma_M( idx, :)= gamma;
        H_M( idx, :)= H_k;
        
        idx= (PARAMS.M-epoch+1)*PARAMS.m+1:(PARAMS.M-epoch+2)*PARAMS.m;
        Phi_M(  idx, :)= Phi_k;
        L_M(    idx, :)= Lk;
        L_pp_M( idx, :)= Lk_pp;
    
        
    %% Full-size PH --> start integrity monitoring
    else 
        % Compute controls
        [G,iwp]= compute_steering(xtrue, iwp, G); G= -deg2rad(3);
        
        % EKF predict step
        [xtrue,XX,PX,Phi_k,alpha]= predict (xtrue,XX,PX,G);
        
        % update the states transition matrix over the horizon
        Phi_M(PARAMS.m+1:end,:)= Phi_M(1:end-PARAMS.m,:);
        Phi_M(1:PARAMS.m,:)= Phi_k;
        
        % Integrity Monitoring
        [P_HMI_worst,H_M,L_M,L_pp_M,Y_M,A_k]= IM (Phi_M ,H_M ,L_M ,L_pp_M , A_k ,Y_M , alpha);
        
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
        
        [q_D, T_D, gamma,gamma_M]= EKF_update(z, idf, epoch,gamma_M,Y_M);
        
        % Store DATA
        store_data(epoch, P_HMI_worst, xtrue, q_D, T_D);
        
        % Plots
        do_plots(xtrue, z, h, epoch);
        
    end
end
%% *****************  END OF MAIN LOOP    *****************

% post_processing_and_plots(step)