%==========================================================================
% FIM optimization to find Vc(n+1)
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


function [y,FIM] = find_next_Vc(x,Vc_f_1tonm1,theta)

    den=1+(Vc_f_1tonm1/theta(3)).^theta(4);
    Si2nm1=1-1./den;
    Si3nm1=(theta(1)-theta(2))*theta(4)*(Vc_f_1tonm1.^theta(4))*theta(3).^(-theta(4)-1)./den.^2;
    Si4nm1=-(theta(1)-theta(2))*log(Vc_f_1tonm1./theta(3)).*(Vc_f_1tonm1./theta(3)).^theta(4)./den.^2;
   
    FIMnm1=[
              dot(Si2nm1,Si2nm1) dot(Si2nm1,Si3nm1) dot(Si2nm1,Si4nm1);
              dot(Si3nm1,Si2nm1) dot(Si3nm1,Si3nm1) dot(Si3nm1,Si4nm1);
              dot(Si4nm1,Si2nm1) dot(Si4nm1,Si3nm1) dot(Si4nm1,Si4nm1)];
    
    den=1+(x/theta(3))^theta(4);
    Sn2=1-1/den;
    Sn3=(theta(1)-theta(2))*theta(4)*(x^theta(4))*theta(3)^(-theta(4)-1)/den^2;
    Sn4=-(theta(1)-theta(2))*log(x/theta(3))*(x/theta(3))^theta(4)/den^2;
       
    FIMn=[dot(Sn2,Sn2) dot(Sn2,Sn3) dot(Sn2,Sn4);
          dot(Sn3,Sn2) dot(Sn3,Sn3) dot(Sn3,Sn4);
          dot(Sn4,Sn2) dot(Sn4,Sn3) dot(Sn4,Sn4)];    
      
    FIM=FIMnm1+FIMn;
    y=real(-det(FIM));

end
    
        
 
