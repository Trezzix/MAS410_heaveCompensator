% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % simulation - WIP
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% time = 0; %s
% timeend = 10; %s
% timestep = 1e-5; %s
% counter = 1;
% 
% z= 0;
% theta = 0;
% thetadot = 0;
% thetadotdto = 0;
% zPl = 0;
% p1 = 0; %TODO needs to be eq
% p2 = 0;
% p1dot = 0;
% p2dot = 0;
% Q1 = 0;
% Q2 = 0;
% V1 = 1e-6; %m^3, temporary, = 1 L
% V2 = 1e-6; %m^3
% V3 = 1e-6; %m^3
% u = 1; % for now
% 
% tic
% while time<timeend
%     %simulate wave
%     z = Zw * sin(((2*pi)/Tw)*time);
%     zdot = Zw * cos((2*pi/Tw)*time) * ((2*pi)/Tw);
%     zdotdot = ((-Zw*(2*pi)^2) / (Tw^2))* sin(((2*pi)/Tw) * time);
%     %flow
%     Q1 = CdAd*u*sign(ps - p1)*sqrt((2/rho) * abs(ps - p1));
%     Q2 = CdAd*u*sign(p2)*sqrt((2/rho) * abs(p2));
%     Qm = thetadot * Dm; %%%%%%%%%%%%%%?
%     %diff eqs
%     p1dot = (beta/V2) * (Q1 - Qm); %?
%     p2dot = (beta/V3) * (Qm - Q2);
%     PL_sim = p1 - p2;
%     %M_M_sim = (mpl * g) / star;
%     thetadotdot = (((PL_sim * Dm) / 2 * pi) - M_M_max)/Jtot; %TODO: M_M will change over time
%     zPl = z + (theta * star); % this is correct
% 
% %save values
%     z_graph(counter) = z;
%     zpl_graph(counter) = zPl;
%     PL_graph(counter) = PL_sim * 1e-5; % to bar
%     P1_graph(counter) = p1 * 1e-5; % to bar
%     P2_graph(counter) = p2 * 1e-5; % to bar
%     Q1_graph(counter) = Q1;
%     Q2_graph(counter) = Q2;
%     time_graph(counter) = time;
% 
%     p1 = p1 + p1dot * timestep;
%     p2 = p2 + p2dot * timestep;
%     thetadot = thetadot + thetadotdot * timestep;
%     theta = theta + thetadotdot * timestep;
% 
%     counter = counter + 1;
%     time = time + timestep;
% 
% end
% toc
% 
% figure
% subplot(2,2,1)
% plot(time_graph,z_graph)
% hold on
% plot(time_graph,zpl_graph)
% grid on
% legend('z-platform','z-payload')
% 
% subplot(2,2,2)
% plot(time_graph,PL_graph)
% hold on
% plot(time_graph,P1_graph)
% plot(time_graph,P2_graph)
% grid on
% legend('PL','P1','P2')
% 
% subplot(2,2,3)
% plot(time_graph,Q1_graph)
% hold on 
% plot(time_graph,Q2_graph)
% grid on
% legend('Q1','Q2')