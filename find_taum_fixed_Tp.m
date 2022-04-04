%==========================================================================
% Solve theta3 equation and find tau 
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




function taum_sol = find_taum_fixed_Tp(Tp,taum, theta_est3,g_normal, k1, mu, sigma, w)

rp_=k1/(mu*taum^2-2*sigma*taum+1)*...
    (((mu*taum-sigma)*sin(w*Tp)+w*cos(w*Tp))*exp(-sigma*Tp)-w*exp(-Tp/taum));

taum_sol=abs(rp_ -g_normal/theta_est3);


% taum_sol=double( vpa(k1)/(vpa(mu)*vpa(taum)^2-2*vpa(sigma)*vpa(taum)+1)*...
%     (((vpa(mu)*vpa(taum)-vpa(sigma))*sin(vpa(w)*vpa(Tp))+vpa(w)*cos(vpa(w)*vpa(Tp)))*exp(-vpa(sigma)*vpa(Tp))-...
%     vpa(w)*exp(-vpa(Tp)/vpa(taum)))-vpa(theta_est5));
end
