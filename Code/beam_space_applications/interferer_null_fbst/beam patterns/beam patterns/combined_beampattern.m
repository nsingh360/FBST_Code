clear
clc
close all



figure(1)

load('mpdr_beam.mat')
beam_mpdr = beam_pattern;
p1 = plot(rad2deg(theta_beam),zeros(length(theta_beam),1),'r--');
hold on
p2 = plot(rad2deg(theta_beam),db(beam_mpdr),'b');
load('slepian_beam.mat')
beam_slepian = beam_pattern;
p3 = plot(rad2deg(theta_beam),db(beam_slepian),'k');
grid on
xlabel('$\theta$ (degrees)','Interpreter','latex')
ylabel('Response(dB)','Interpreter','latex')
legend([p1(1,:); p2(1,:); p3(1,:)],{'Distortionless response','8-bit MPDR','1-bit Slepian MPDR'},'Location','Northwest','Interpreter','latex')
xlim([-90,90])
ylim([-45,1])
set(gca,'Fontsize',14)

saveas(gca,'figs/mpdr_beampattern_combine.png')

