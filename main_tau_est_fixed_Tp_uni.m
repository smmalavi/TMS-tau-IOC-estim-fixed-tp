%==========================================================================
% uniform sampling method
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
% uniform 


Vc_u=linspace(min_Vc,max_Vc,n);

y_u=real(true_theta(2)+(true_theta(1)-true_theta(2))./...
    (1+((Vc_u+sigma_x*randn(1,length(Vc_u)))/true_theta(3)).^true_theta(4))+...
    sigma_y*randn(1,length(Vc_u)));


[xData, yData] = prepareCurveData( [Vc_base Vc_u], [y_base y_u]);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Algorithm = 'Trust-Region';
opts.Robust ='LAR';%'LAR';% 'Bisquare';
opts.Lower = paramLB1;
opts.Upper = paramUB1;
opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
[fitresult, gof] = fit( xData, yData, ft1, opts );
t_est_u(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
    
    
taum_func = @(taum) find_taum_fixed_Tp(true_Tp,taum, t_est_u(n,3), g_normal, k1, mu, sigma, w);% 
opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
        'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);

    
    
pts = min_taum+ (max_taum-min_taum).*rand(200,1);
tpoints = CustomStartPointSet(pts);
allpts = {tpoints};
                
%[taum_est_u(n),fval_u(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
  
[taum_est_u(n),fval_u(n),flagm_ms,outptm_ms,manyminsm_ms] = run(MultiStart,problem,allpts);




%% BAd-fit

if n>no_ini_pulses+4
        rel_er_1= abs((t_est_u(n,:)-t_est_u(n-1,:))./t_est_u(n-1,:));
        r_e_tm1=abs(taum_est_u(n)-taum_est_u(n-1))/taum_est_u(n-1); 
           
        if ~isempty(find(rel_er_1 > 20*er_tol)) || r_e_tm1 > 20*er_tol_taum(1)           
            opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
            [fitresult, gof] = fit( xData, yData, ft1, opts );
            t_est_u(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
    
            taum_func = @(taum) find_taum_fixed_Tp(true_Tp,taum, t_est_u(n,3), g_normal, k1, mu, sigma, w);% 
            opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
            problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
                'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
                
            %[taum_est_u(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
  
            pts = min_taum+ (max_taum-min_taum).*rand(200,1);
            tpoints = CustomStartPointSet(pts);
            allpts = {tpoints};
                
            [taum_est_u(n),fval_f(n),flagm_ms,outptm_ms,manyminsm_ms] = run(MultiStart,problem,allpts);    
        end
        
        
        % check bad-fit for the 2nd time
        
        rel_er_1= abs((t_est_u(n,:)-t_est_u(n-1,:))./t_est_u(n-1,:));
        r_e_tm1=abs(taum_est_u(n)-taum_est_u(n-1))/taum_est_u(n-1); 
           
        if ~isempty(find(rel_er_1 > 20*er_tol)) || r_e_tm1 > 20*er_tol_taum(1)
            
            opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
            [fitresult, gof] = fit( xData, yData, ft1, opts );
            t_est_u(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
    
    
            taum_func = @(taum) find_taum_fixed_Tp(true_Tp,taum, t_est_u(n,3), g_normal, k1, mu, sigma, w);% 
            opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
            problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
                'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
                
            %[taum_est_u(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
  
            pts = min_taum+ (max_taum-min_taum).*rand(200,1);
            tpoints = CustomStartPointSet(pts);
            allpts = {tpoints};
                
            [taum_est_u(n),fval_f(n),flagm_ms,outptm_ms,manyminsm_ms] = run(MultiStart,problem,allpts);
        end
        
        
        
    end

%%

check_stopping_u;    


Vc_u_matrix(n,1:n)=Vc_u;
y_u_matrix(n,1:n)=y_u;