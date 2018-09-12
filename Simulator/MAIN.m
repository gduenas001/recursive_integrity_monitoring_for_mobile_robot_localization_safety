dbstop if error
% dbclear if error

clear; close all; configfile;

h= setup_animations();


% *****************    MAIN LOOP    *****************
for epoch= 1:PARAMS.numEpochs
    disp(['Step: ',num2str(epoch)]);
    
    %% Growing PH --> store matrices until full-size
    if epoch <= PARAMS.M + 1
        % Compute controls
        [G,iwp]= compute_steering(xtrue, iwp, G); G= -deg2rad(3);
        
        % EKF predict step
        [xtrue,XX,PX,Phi_k,alpha]= predict (xtrue,XX,PX,G); 
        
        % Get measurements
        [z,idfTrue]= get_observations(xtrue);
        z= add_observation_noise(z);
        
        % DA
        [idf,DATA.numAssoc(epoch)]= data_associate_LNN_LS(z);
        nk= DATA.numAssoc(epoch)* PARAMS.m_F;
        
        % Store associations data
        DATA.IA(epoch)= any( (idfTrue - idf).*idf );
        
        % EKF update while the PH is increasing
        [q_D,T_D,Hk,Lk,Y_k,gamma_k] = EKF_update_shrinking_horizon(z, idf, epoch);
        Lk_pp= Phi_k - Lk*Hk*Phi_k; % Kalman gain prime-prime
        
        % Store DATA
        store_data(epoch, 0, xtrue, q_D, T_D);
        
        % Increase preceding horizon
        idx= (PARAMS.M-epoch+1)*PARAMS.m+1:(PARAMS.M-epoch+2)*PARAMS.m;
        Phi_M(  idx, :)= Phi_k;
        
        if length(n_M_array) == PARAMS.M
            n_M_array= [nk, n_M_array];
            n_M_array(end)= [];
        end
        
        if ~isempty(gamma_M)
            gamma_M= [gamma_k; gamma_M];
            Y_M= [ Y_k, zeros(nk,size(Y_M,1))      ;
                       zeros(size(Y_M,1),nk), Y_M ];
        else
            gamma_M= gamma_k;
            Y_M= Y_k;
        end
        
        % Increase preceding horizon -- Using Cells
        if length(Y_M_cell) == PARAMS.M + 1
            Y_M_cell(end)= [];
        end
        Y_M_cell= [ {Y_k} ,Y_M_cell];
        
        if length(L_M_cell) == PARAMS.M + 1
            L_M_cell(end)= [];
        end
        L_M_cell= [ {Lk} ,L_M_cell];
        
        if length(Lpp_M_cell) == PARAMS.M + 1
            Lpp_M_cell(end)= [];
        end
        Lpp_M_cell= [ {Lk_pp} ,Lpp_M_cell];
        
        if length(H_M_cell) == PARAMS.M + 1
            H_M_cell(end)= [];
        end
        H_M_cell= [ {Hk}, H_M_cell];
        
        % Plots
        do_plots(xtrue, z, h, epoch);

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
        [P_HMI_worst, H_M_cell, Y_M, Y_M_cell, A_k, L_M_cell, Lpp_M_cell, n_M_array]=...
            IM (Phi_M , H_M_cell, A_k , Y_M, Y_M_cell, alpha, L_M_cell, Lpp_M_cell,n_M_array,epoch);
        
        % Get measurements
        [z,idfTrue]= get_observations(xtrue);
        z= add_observation_noise(z);
        
        % Data Association
        if  ~isempty(z)
            [idf,DATA.numAssoc(epoch)]= data_associate_LNN_LS(z);
            
            % Store associations data
            DATA.IA(epoch)= any( (idfTrue - idf).*idf );
        else
            idf= [];
            DATA.numAssoc(epoch)= 0;
            DATA.IA(epoch)= 0;
        end
        
        % EKF update
        [q_D, T_D, gamma_M]= EKF_update(z, idf, epoch,gamma_M,Y_M,n_M_array);
        
        % Store DATA
        store_data(epoch, P_HMI_worst, xtrue, q_D, T_D);
        
        % Plots
        do_plots(xtrue, z, h, epoch);
        
    end
end
%% *****************  END OF MAIN LOOP    *****************

% post_processing_and_plots(step)