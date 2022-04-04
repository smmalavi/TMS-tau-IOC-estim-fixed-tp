%==========================================================================
% Plot io data, estimated parameters and time-constant 
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



fig_font_size=16;
leg_font_size=16;


% ycurve_true=real(true_theta(2)+(true_theta(1)-true_theta(2))./...
%     (1+(Vc_val/true_theta(3)).^true_theta(4)));
% 

%% io
ylim=[2*10^-7 10^-2];
xlim=[0 1.2];
Plot_io_curve_data;
ax=gca;
ax.XLim=xlim;
ax.YLim=ylim;
file_path_name=['fig-io-curve-data'];
saveas(gcf,file_path_name,'fig')
saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')


Plot_io_curve_fim_uni_nmax;
ax=gca;
ax.XLim=xlim;
ax.YLim=ylim;

file_path_name=['fig-io-nmax'];
saveas(gcf,file_path_name,'fig')
saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')


Plot_io_curve_fim_uni_nf;
ax=gca;
ax.XLim=xlim;
ax.YLim=ylim;

file_path_name=['fig-io-nf'];
saveas(gcf,file_path_name,'fig')
saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')


%% io parameters 
figure
hold on
plot([0 n],[true_yl true_yl],'k','LineWidth',1)
plot(t_est_u(:,1),'--b','LineWidth',1)
plot(t_est_f(:,1),'--r','LineWidth',1)
legend('True','Uniform Sampling', 'FIM SPE','Location','northwest','Box','off')
xlabel('$n$', 'interpreter','latex')
ylabel('$\theta_1$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

strmax = ['a'];
text(140,-.5,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);


file_path_name=['fig-yl'];
saveas(gcf,file_path_name,'fig')
saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')

%%

figure
hold on
plot([0 n],[true_yh true_yh],'k')
plot(t_est_f(:,2),'--r','LineWidth',1)
plot(t_est_u(:,2),'--b','LineWidth',1)

xlabel('$n$', 'interpreter','latex')
ylabel('$\theta_2$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

strmax = ['b'];
text(140,-.2,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);


file_path_name=['fig-yh'];
saveas(gcf,file_path_name,'fig')
saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')


%%
figure
hold on
plot([0 n],[true_theta(3) true_theta(3)],'k')
plot(t_est_f(:,3),'--r','LineWidth',1)
plot(t_est_u(:,3),'--b','LineWidth',1)

xlabel('$n$', 'interpreter','latex')
ylabel('$\theta_3$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

strmax = ['c'];
text(140,0.56,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);


file_path_name=['fig-m'];
saveas(gcf,file_path_name,'fig')
saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')


%%
figure
hold on
plot([0 n],[true_s true_s],'k','LineWidth',1)
plot(t_est_f(:,4),'--r','LineWidth',1)
plot(t_est_u(:,4),'--b','LineWidth',1)

%legend('True','Estimate')
xlabel('$n$', 'interpreter','latex')
ylabel('$\theta_4$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

strmax = ['d'];
text(140,92,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);


file_path_name=['fig-s'];
saveas(gcf,file_path_name,'fig')
saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')

%%
figure
hold on
plot([0 n],[true_taum true_taum],'k','LineWidth',1)
plot(taum_est_u,'--b','LineWidth',1)
plot(taum_est_f,'--r','LineWidth',1)

legend('True','Uniform Sampling', 'FIM SPE','Location','southeast','Box','off')
xlabel('$n$', 'interpreter','latex')

yt = get(gca, 'YTick');                                 % 'XTick' Values
set(gca, 'YTick', yt, 'YTickLabel', yt/1e-6) 
ylabel('$\tau~ (\mu s)$','interpreter','latex')
%ylabel('$\tau_m$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'lin')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on

file_path_name=['fig-taum'];
saveas(gcf,file_path_name,'fig')
saveas(gcf,file_path_name,'pdf')
saveas(gcf,file_path_name,'png')
saveas(gcf,file_path_name,'epsc')
