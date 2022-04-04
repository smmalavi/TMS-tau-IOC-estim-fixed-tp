%==========================================================================
% Plot io data 
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


figure
hold on

x_dg=linspace(0,1,10000);

myy=real(true_theta(2)+(true_theta(1)-true_theta(2))./...
    (1+((x_dg+sigma_x*randn(1,length(x_dg)))/true_theta(3)).^true_theta(4))+...
    sigma_y*randn(1,length(x_dg)));


plot(x_dg, 10.^myy,'x','Color', [0.5 0.5 0.5],'MarkerSize',8)

% ycurve_true=real(true_theta(2)+(true_theta(1)-true_theta(2))./...
%     (1+(Vc_val/true_theta(3)).^true_theta(4)));

plot(Vc_val,10.^ycurve_true,'k','LineWidth',1)

%xlabel('$V_C$ (normalized)', 'interpreter','latex')
ylabel('$y (V_{p-p})$', 'interpreter','latex')
set(gca, 'XScale', 'lin', 'YScale', 'log')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = fig_font_size;
box on
grid on

 
%fig_font_size=16;
% legend
strmax = ['a'];
text(1.1,5e-3,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);




%% % legend
% strmax = ['Reference MEP data'];
% text(.17,.004,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);
% strmax = ['Reference IO curve'];
% text(.17,.0022,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);
% strmax = ['Baseline sample'];
% text(.17,.0012,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);
% strmax = ['Uniform sample'];
% text(.17,.00065,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);
% strmax = ['IO curve: uniform sampling'];
% text(.17,.00038,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);
% strmax = ['FIM SPE sample'];
% text(.17,.00022,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);
% strmax = ['IO curve: FIM SPE'];
% text(.17,.000124,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);
% 
% 
% line([0.05 .15],[0.0022 0.0022],'Color','k','LineStyle','-');
% line([0.05 .15],[0.00038 0.00038],'Color','b','LineStyle','--');
% line([0.05 .15],[.000124 .000124],'Color','r','LineStyle','--');
% 
% plot(.1,0.004, 'Marker','x', 'MarkerSize',8,'MarkerEdgeColor',[0.5 0.5 0.5]);
% plot(.1, .0012, 'Marker','d', 'MarkerSize',8,'MarkerEdgeColor','m');
% plot(.1, .00065, 'Marker','^', 'MarkerSize',8,'MarkerEdgeColor','b');
% plot(.1, .00022, 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor','r');
% 

%% % legend
strmax = ['Reference data'];
text(0.72,8e-5,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);
strmax = ['Reference IO curve'];
text(0.72,3.8e-5,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);
strmax = ['Baseline sample'];
text(0.72,1.8e-5,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);
strmax = ['Uniform sample'];
text(0.72,8.5e-6,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);
strmax = ['uniform sampling'];
text(0.72,4.2e-6,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);
strmax = ['FIM SPE sample'];
text(0.72,2e-6,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);
strmax = ['FIM SPE'];
text(0.72,1e-6,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', leg_font_size);


line([0.6 0.7],[3.8e-5 3.8e-5],'Color','k','LineStyle','-');
line([0.6 0.7],[4.2e-6 4.2e-6],'Color','b','LineStyle','--');
line([0.6 0.7],[1e-6 1e-6],'Color','r','LineStyle','--');

plot(0.65, 8e-5, 'Marker','x', 'MarkerSize',8,'MarkerEdgeColor',[0.5 0.5 0.5]);
plot(0.65, 1.8e-5, 'Marker','d', 'MarkerSize',8,'MarkerEdgeColor','m');
plot(0.65, 8.5e-6, 'Marker','^', 'MarkerSize',8,'MarkerEdgeColor','b');
plot(0.65, 2e-6, 'Marker','o', 'MarkerSize',8,'MarkerEdgeColor','r');

