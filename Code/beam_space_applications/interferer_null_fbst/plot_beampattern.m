clc;
clear;
close all;

load('beam_pattern_fbst.mat');
load('beam_pattern_fds.mat');
load('beam_pattern_subBand_nulling.mat');

theta_beam = linspace(-pi/2,pi/2,4000);

figure(1)
p0 = plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--');
hold on
grid on
p2 = plot(rad2deg(theta_beam),db(beam_pattern_fds_nulled_beamspace),'Color',[0 0 1 0.6],'LineWidth',1);
hold on
grid on
p3 = plot(rad2deg(theta_beam),db(beam_pattern_ds_nulled_beamspace),'Color',[0 1 0 0.2],'LineWidth',1);
hold on
grid on
p4 = plot(rad2deg(theta_beam),db(beam_pattern_subBand),'Color',[1 0 1 0.5]);
hold on
grid on
p1 = plot(rad2deg(theta_beam),db(beam_pattern_fbst_nulled),'Color',[0 0 0 0.5],'LineWidth',1);

xlabel('$\theta$ (degrees)','Interpreter','latex','FontSize',12)
ylabel('Response(dB)','Interpreter','latex','FontSize',12)
legend_p0 = plot(-1, -1, 'r--', 'LineWidth', 2);  % Solid red
legend_p1 = plot(-1, -1, 'k-', 'LineWidth', 2); % Dashed blue
legend_p2 = plot(-1, -1, 'b-', 'LineWidth', 2); % Dashed blue
legend_p3 = plot(-1, -1, 'g-','LineWidth',2);
legend_p4 = plot(-1, -1, '-', 'LineWidth', 2,'Color',[1 0 1 1]); % Dashed blue
% Legend
legend([legend_p0, legend_p1, legend_p2, legend_p3, legend_p4], {'Distortionless response', 'FBST', 'FDS (R=16)', 'DS (R=16)', 'Sub-Band Processing'},'location','southwest','Interpreter','latex','FontSize',12);
xlim([-90,90])
% ylim([-45,.5])
% y1 = ylim;
% ylim([-30, 1.5]);
% exportgraphics(gcf, 'super_res_ula_beam_pattern.pdf', 'ContentType', 'vector')
% legend({'Distortionless response','No Nulling', 'Projection Nulling'},'location','northwest','Interpreter','latex')
% First Inset
% ax1 = axes('Position', [0.55, 0.7, 0.2, 0.2]); % [x, y, width, height]
% box on;
% plot(rad2deg(theta_beam(2455:2470)),zeros(length(theta_beam(2455:2470)),1),'r--');
% hold on;
% plot(rad2deg(theta_beam(2455:2470)), db(beam_pattern_fbst_nulled(:,2455:2470)), 'Color',[0 0 0 1], 'LineWidth', 1);
% hold on
% plot(rad2deg(theta_beam(2455:2470)), db(beam_pattern_fds_nulled_beamspace(:,2455:2470)), 'Color',[0 0 1 0.2],'LineWidth',1);
% % y1 = ylim;
% ylim([-90, -30]);
% set(gca, 'YTick', [-90 -60 -30]);

ax1 = axes('Position', [0.45, 0.7, 0.2, 0.2]); % [x, y, width, height]
box on;
plot(rad2deg(theta_beam(2400:2500)),zeros(length(theta_beam(2400:2500)),1),'r--');
hold on;
plot(rad2deg(theta_beam(2400:2500)), db(beam_pattern_fbst_nulled(1:5,2400:2500)), 'Color',[0 0 0 1], 'LineWidth', 1);
hold on
plot(rad2deg(theta_beam(2400:2500)), db(beam_pattern_fds_nulled_beamspace(1:5,2400:2500)), 'Color',[0 0 1 0.4],'LineWidth',1);
hold on
plot(rad2deg(theta_beam(2400:2500)), db(beam_pattern_subBand(1:5,2400:2500)),'Color',[1 0 1 0.5],'LineWidth',1);
% y1 = ylim;
ylim([-90, -30]);
set(gca, 'YTick', [-90 -60 -30]);
exportgraphics(gcf, 'ula_interference_timeBeam_pattern_subBand.pdf', 'ContentType', 'vector')