B = 2^5;
%% NOTE: The angle sampling here covers only [0,pi/2] in theta and [-pi/2,pi/2] in phi
%% To cover the remaining (i.e. [0,pi/2] in theta and [pi/2,3pi/2] in phi just add pi to the phi angle)
%% There half of the beams will be in the [-pi/2,pi/2] phi group and the other half will be in [pi/2,3pi/2] phi group
theta_upa = zeros(B*(2*B+1),1);
phi_upa = zeros(B*(2*B+1),1);

theta_phi = zeros(B,B);
[K,L] = meshgrid(-B:B,-B:B);
Theta = acos(sqrt((K.^2+L.^2)/(2*B^2)));
Phi = asin(L./sqrt(L.^2 + K.^2));
count=1;
for ll=1:2*B+1
    l = ll-B-1;
    for k=0:B
        if isreal(asin(sqrt((l^2+k^2)/(B^2))))
            theta_upa(count) = asin(sqrt((l^2+k^2)/(B^2)));
            phi_upa(count) = asin(l/sqrt(l^2 + k^2));
            count = count+1;
        end
    end
end
plot(theta_upa*180/pi,phi_upa*180/pi,'.')
% surf(K,L,Theta.*(180/pi),'FaceColor','r');
% hold on;
% surf(K,L,Phi.*(180/pi),'FaceColor','g');
% 
% xlabel('K','Interpreter','latex')
% ylabel('L','Interpreter','latex')
% legend({'$$\Theta$$','$$\Phi$$'},'Location','northwest','Interpreter','latex')