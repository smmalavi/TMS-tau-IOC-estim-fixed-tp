%==========================================================================
% Stopping rule for uniform sampling 
%==========================================================================
%
% Seyed Mohammad Mahdi Alavi+, Stellantis (Chrysler), Canada 
% Fidel Vila-Rodriguez, Unitverisyt of British Columbia, Canada 
% Adam Mahdi, University of Oxford, UK
% Stefan M. Goetz, University of Cambridge (UK), Duke University (USA)
% +: code written by
% e-mail: mahdi.alavi.work@gmail.com
%
% April 2022
%==========================================================================



if n>no_ini_pulses+4
    rel_er_1= abs((t_est_u(n,:)-t_est_u(n-1,:))./t_est_u(n-1,:));
    rel_er_2= abs((t_est_u(n-1,:)-t_est_u(n-2,:))./t_est_u(n-2,:));
    rel_er_3= abs((t_est_u(n-2,:)-t_est_u(n-3,:))./t_est_u(n-3,:));
    rel_er_4= abs((t_est_u(n-4,:)-t_est_u(n-5,:))./t_est_u(n-5,:));
    rel_er_5= abs((t_est_u(n,:)-t_est_u(n-1,:))./t_est_u(n-1,:));
    
    r_e_tm1=abs(taum_est_u(n)-taum_est_u(n-1))/taum_est_u(n-1); 
    r_e_tm2=abs(taum_est_u(n-1)-taum_est_u(n-2))/taum_est_u(n-2); 
    r_e_tm3=abs(taum_est_u(n-2)-taum_est_u(n-3))/taum_est_u(n-3); 
    r_e_tm4=abs(taum_est_u(n-3)-taum_est_u(n-4))/taum_est_u(n-4); 
    r_e_tm5=abs(taum_est_u(n-4)-taum_est_u(n-5))/taum_est_u(n-5); 
    r_e_tm=[r_e_tm1 r_e_tm2 r_e_tm3 r_e_tm4 r_e_tm5];
    
    if rel_er_1<er_tol & rel_er_2<er_tol & rel_er_3<er_tol & rel_er_4<er_tol & rel_er_5<er_tol
        if r_e_tm < er_tol_taum.*ones(1,5)
            n_conv_u=[n_conv_u n];
        end
    end
    
end
    