%==========================================================================
% Main Matlab file for 
% Estimation of neural membrane time-constant and IO curve in TMS
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

clear all
close all
clc

rng('shuffle'); 



n_max=150;



% cTMS pulse parameters 
R=0.1;%Ohm
r=20e-3;%mOhm
L=16e-6;%H
C=716e-6;%F
delta=3.2e-6;%(V/m)(A/s)
sigma=r/(2*L);
w=sqrt(1/(L*C)-sigma^2);
mu=(w^2+sigma^2);
k1=delta/(L*w);
min_Tp=10e-6;
max_Tp=200e-6;
min_Vc=0.01;
max_Vc=1;
Vc_val = linspace(0, 1, 100);
size_Vc_val=length(Vc_val);

min_taum=90e-6; %1e-6;
max_taum=250e-6;
x_dg=linspace(0,1,1000);


sigma_y=0.1;
sigma_x=0.05;
true_yl=-6.5+ (-5.5-(-6.5))*rand(1,1);
true_yh=-3+ (-2-(-3))*rand(1,1);
true_s=1+ (50-1)*rand(1,1);
true_taum=min_taum+(max_taum-min_taum)*rand(1,1);
true_Tp=10e-6+ (200e-6-10e-6).*rand(1,1);
g_normal=0.025;

true_theta(1)=true_yl;
true_theta(2)=true_yh;
true_theta(4)=true_s;

true_rp=k1/(mu*true_taum^2-2*sigma*true_taum+1)*...
    (((mu*true_taum-sigma)*sin(w*true_Tp)+w*cos(w*true_Tp))*exp(-sigma*true_Tp)-...
    w*exp(-true_Tp/true_taum));


true_theta(3)=g_normal/true_rp;

myy=real(true_theta(2)+(true_theta(1)-true_theta(2))./...
    (1+((x_dg+sigma_x*randn(1,length(x_dg)))/true_theta(3)).^true_theta(4))+...
    sigma_y*randn(1,length(x_dg)));



ycurve_true=real(true_theta(2)+(true_theta(1)-true_theta(2))./...
    (1+(Vc_val/true_theta(3)).^true_theta(4)));

%end
%
hold on
plot(x_dg, 10.^myy,'xk')
plot(Vc_val,10.^ycurve_true,'k')
xlabel('$V_C$ (normalized)', 'interpreter','latex')
ylabel('$y (V_{p-p})$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'log')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = 18;
box on
grid on

%%
% Get basline 

n_base=50;% number of baseline data
Vc_base=zeros(1,n_base);

y_base=real(true_theta(2)+(true_theta(1)-true_theta(2))./...
    (1+((Vc_base+sigma_x*randn(1,length(Vc_base)))/true_theta(3)).^true_theta(4))+...
    sigma_y*randn(1,length(Vc_base)));


hold on
plot(Vc_base, 10.^y_base, 'dc')


% intialization
no_ini_pulses=3;
%Tp_f_ini=repmat(Tp_min+(Tp_max-Tp_min)*rand(1,1),1,3);
ini_Vc=min_Vc+(max_Vc-min_Vc)*rand(1,no_ini_pulses);
ini_taum=min_taum+(max_taum-min_taum)*rand(1,1);

y_ini=real(true_theta(2)+(true_theta(1)-true_theta(2))./...
    (1+((ini_Vc+sigma_x*randn(1,length(ini_Vc)))/true_theta(3)).^true_theta(4))+...
    sigma_y*randn(1,length(ini_Vc)));
figure(1)
plot(ini_Vc, 10.^y_ini, 'sg')


n=no_ini_pulses;
Vc_f=ini_Vc;
y_f=y_ini;


% curve fitting 
[xData, yData] = prepareCurveData( [Vc_base Vc_f], [y_base y_f]);
cf_func1='b+(a-b)./(1+(x/c)^d)';
ft1 = fittype( cf_func1, 'independent', 'x', 'dependent', 'y' );
paramLB1=[-7   -3 0  1];%yl, yh, m, s
paramUB1=[-5 -2 1  100];%yl, yh, m, s
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Algorithm = 'Trust-Region';
opts.Lower = paramLB1;
opts.Upper = paramUB1;
opts.Robust ='LAR';%'LAR';% 'Bisquare';
ini_guess=paramLB1+ (paramUB1-paramLB1).*rand(1,4) ; 
opts.StartPoint =ini_guess;
         
[fitresult, gof] = fit( xData, yData, ft1, opts );

t_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d];

ycurve_est_fim_nmax=real(t_est_f(n,2)+(t_est_f(n,1)-t_est_f(n,2))./(1+(Vc_val/t_est_f(n,3)).^t_est_f(n,4)));

% figure(1)
% plot(Vc_val, 10.^(ycurve_est_f), 'g')


%
taum_func = @(taum) find_taum_fixed_Tp(true_Tp,taum, t_est_f(n,3), g_normal, k1, mu, sigma, w);% 
    opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
        'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
                
    [taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
   

n_conv_f=[];
n_conv_u=[];
er_tol=.01*ones(1,4);
er_tol_taum=.01*ones(1,5);

Vc_u_matrix=NaN(n_max,n_max);
y_u_matrix=NaN(n_max,n_max);

%bad_fit_flag=0;

for n=no_ini_pulses+1:n_max
    n
    ObjFunc_FIM = @(x) find_next_Vc(x,Vc_f,t_est_f(n-1,:));% 
    optsFIM = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','x0',max_Vc*rand,...
                    'objective',ObjFunc_FIM,'lb',min_Vc,'ub',max_Vc,'options',optsFIM);
    [Vc_f(n),fval_fimopt,flagm_fimopt,outptm_fimopt,manyminsm_fimopt] = run(GlobalSearch,problem);
        
    %y_f(n) = virtualsubjectEIVStimulate_01(Vc_f(n), myp);

    y_f(n)=real(true_theta(2)+(true_theta(1)-true_theta(2))./...
        (1+((Vc_f(n)+sigma_x*randn(1,1))/true_theta(3)).^true_theta(4))+...
        sigma_y*randn(1,1));


    
    [xData, yData] = prepareCurveData( [Vc_base Vc_f], [y_base y_f]);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Algorithm = 'Trust-Region';
    opts.Robust ='LAR';%'LAR';% 'Bisquare';
    opts.Lower = paramLB1;
    opts.Upper = paramUB1;
    opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
    [fitresult, gof] = fit( xData, yData, ft1, opts );
    t_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
    
    
    taum_func = @(taum) find_taum_fixed_Tp(true_Tp,taum, t_est_f(n,3), g_normal, k1, mu, sigma, w);% 
    opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
        'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
                
    %[taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
  
    pts = min_taum+ (max_taum-min_taum).*rand(200,1);
    tpoints = CustomStartPointSet(pts);
    allpts = {tpoints};
                
    [taum_est_f(n),fval_f(n),flagm_ms,outptm_ms,manyminsm_ms] = run(MultiStart,problem,allpts);
    
    
    % check bad-fit first time
    
    if n>no_ini_pulses+4
        rel_er_1= abs((t_est_f(n,:)-t_est_f(n-1,:))./t_est_f(n-1,:));
        r_e_tm1=abs(taum_est_f(n)-taum_est_f(n-1))/taum_est_f(n-1); 
           
        if ~isempty(find(rel_er_1 > 20*er_tol)) || r_e_tm1 > 20*er_tol_taum(1)           
            opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
            [fitresult, gof] = fit( xData, yData, ft1, opts );
            t_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
    
            taum_func = @(taum) find_taum_fixed_Tp(true_Tp,taum, t_est_f(n,3), g_normal, k1, mu, sigma, w);% 
            opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
            problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
                'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
                
            %[taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
  
            pts = min_taum+ (max_taum-min_taum).*rand(200,1);
            tpoints = CustomStartPointSet(pts);
            allpts = {tpoints};
                
            [taum_est_f(n),fval_f(n),flagm_ms,outptm_ms,manyminsm_ms] = run(MultiStart,problem,allpts);    
        end
        
        
        % check bad-fit for the 2nd time
        
        rel_er_1= abs((t_est_f(n,:)-t_est_f(n-1,:))./t_est_f(n-1,:));
        r_e_tm1=abs(taum_est_f(n)-taum_est_f(n-1))/taum_est_f(n-1); 
           
        if ~isempty(find(rel_er_1 > 20*er_tol)) || r_e_tm1 > 20*er_tol_taum(1)
            
            opts.StartPoint =paramLB1+ (paramUB1-paramLB1).*rand(1,4);
            [fitresult, gof] = fit( xData, yData, ft1, opts );
            t_est_f(n,:)=[fitresult.a fitresult.b fitresult.c fitresult.d ];
    
    
            taum_func = @(taum) find_taum_fixed_Tp(true_Tp,taum, t_est_f(n,3), g_normal, k1, mu, sigma, w);% 
            opts_taum = optimoptions(@fmincon,'Algorithm','interior-point');
            problem = createOptimProblem('fmincon','x0',min_taum+ (max_taum-min_taum)*rand,...
                'objective',taum_func,'lb',min_taum,'ub',max_taum,'options',opts_taum);
                
            %[taum_est_f(n),fval(n),flagm,outptm,manyminsm] = run(GlobalSearch,problem);
  
            pts = min_taum+ (max_taum-min_taum).*rand(200,1);
            tpoints = CustomStartPointSet(pts);
            allpts = {tpoints};
                
            [taum_est_f(n),fval_f(n),flagm_ms,outptm_ms,manyminsm_ms] = run(MultiStart,problem,allpts);
        end
        
        
        
    end
    
       
    check_stopping_fim;
    
    
    main_tau_est_fixed_Tp_uni;
    
end

main_plots

