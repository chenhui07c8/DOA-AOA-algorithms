close all;
clear all;
clc;
% load data before plotting
%% Fig 1
load fig_8_1_data_3.mat;
% close all;
plot_ind = 7:17;
linewidth = 1.5;
markersize = 6;
figure;
semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(rvec2(plot_ind)),'-.+b','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec8(plot_ind)),'-sg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec6(plot_ind)),'-.>m','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec3(plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec7(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec5(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind), sqrt(crlb(plot_ind)),'--k','LineWidth',2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);

grid on;
xlabel('SNR [dB]');
ylabel('RMSE [Deg]');
% set(gca,'FontSize',12)

% Zoomed plot ****************************************
bounds = [32 33 10^-1.35 10^-1.15]; % area to be zoomed
pos = [0.16 0.15 0.17 0.2];      % projected area [x, y, x_width, y_width]
vertex = [2 3];                 % vertex to be linked
p = gca;
% Calculate x,y points of zoomPlot
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = (pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(p.XLim)+(p.XLim(1));
y1 = 10^((pos(2)-(p.Position(2)))/(p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));
y2 = 10^(((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));


rectangle('Position',[bounds(1) bounds(3) bounds(2)-bounds(1) bounds(4)-bounds(3)]);
hold on;
if any(vertex==1)
    plot([bounds(1) x1], [bounds(4) y2], '-.k'); % Line to vertex 1
end
if any(vertex==2)
    plot([bounds(2) x2], [bounds(4) y2], '-.k'); % Line to vertex 2
end
if any(vertex==3)
    plot([bounds(2) x2], [bounds(3) y1], '-.k'); % Line to vertex 4
end
if any(vertex==4)
    plot([bounds(1) x1], [bounds(3) y1], '-.k'); % Line to vertex 3
end
% legend hide

z = axes('position',pos);
box on; % put box around new pair of axes
semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec8(plot_ind)),'-sg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec6(plot_ind)),'-.>m','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec3(plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec7(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec5(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind), sqrt(crlb(plot_ind)),'--k','LineWidth',2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);

axis(bounds);

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})

% set(gcf,'position',[100,100,420*1.2,400*1.2])
% print -dpng -r600 sim1-3sensors.png

legend('PDP','Two-Step','2Q-Order','MUSIC (0.01^\circ)','MLE (0.01^\circ)','MLE (0.2^\circ)','CRLB', 'Location','northeast')
set(gca,'FontSize',12)


%% Fig 2
% load fig_8_1_data_5.mat;
plot_ind = 3:13;
linewidth = 1.5;
markersize = 6;
figure;
semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec8(plot_ind)),'-sg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec5(plot_ind)),'-.>m','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec3(plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec7(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec4(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind), sqrt(crlb(plot_ind)),'--k','LineWidth',2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);

% axis([0 60 10^-7 10^4])
xlabel('SNR [dB]');
ylabel('RMSE [Deg]');
grid on;
% set(gcf,'position',[100,100,420*1.2,400*1.2])

% Zoomed plot ****************************************
bounds = [24.5 25.5 10^-1.43 10^-1]; % area to be zoomed
pos = [0.17 0.15 0.17 0.2];      % projected area [x, y, x_width, y_width]
vertex = [2 3];                 % vertex to be linked
p = gca;
% Calculate x,y points of zoomPlot
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = (pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(p.XLim)+(p.XLim(1));
y1 = 10^((pos(2)-(p.Position(2)))/(p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));
y2 = 10^(((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));

rectangle('Position',[bounds(1) bounds(3) bounds(2)-bounds(1) bounds(4)-bounds(3)]);
hold on;
if any(vertex==1)
    plot([bounds(1) x1], [bounds(4) y2], '-.k'); % Line to vertex 1
end
if any(vertex==2)
    plot([bounds(2) x2], [bounds(4) y2], '-.k'); % Line to vertex 2
end
if any(vertex==3)
    plot([bounds(2) x2], [bounds(3) y1], '-.k'); % Line to vertex 4
end
if any(vertex==4)
    plot([bounds(1) x1], [bounds(3) y1], '-.k'); % Line to vertex 3
end
% legend hide
z = axes('position',pos);

box on % put box around new pair of axes
semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec8(plot_ind)),'-sg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec5(plot_ind)),'-.>m','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec3(plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec7(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec4(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind), sqrt(crlb(plot_ind)),'--k','LineWidth',2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);

axis(bounds);

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})

legend('PDP','Two-Step','2Q-Order','MUSIC (0.01^\circ)','MLE (0.01^\circ)','MLE (0.2^\circ)','CLRB');
set(gca,'FontSize',12)


% print -dpng -r600 sim1-5sensors.png

%% Fig 3
% close all
% load fig_8_1_data_8.mat;
plot_ind = 1:13;
linewidth = 1.5;
markersize = 6;
figure;
semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec8(plot_ind)),'-sg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec5(plot_ind)),'-.>m','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec3(plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec7(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec4(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind), sqrt(crlb(plot_ind)),'--k','LineWidth',2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);

% axis([0 60 10^-7 10^4])
xlabel('SNR [dB]');
ylabel('RMSE [Deg]');
grid on;
% set(gcf,'position',[100,100,420*1.2,400*1.2])

% Zoomed plot ****************************************
bounds = [24.5 25.5 10^-1.8 10^-1.1]; % area to be zoomed
pos = [0.17 0.15 0.17 0.2];      % projected area [x, y, x_width, y_width]
vertex = [2 3];                 % vertex to be linked
p = gca;
% Calculate x,y points of zoomPlot
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = (pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(p.XLim)+(p.XLim(1));
y1 = 10^((pos(2)-(p.Position(2)))/(p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));
y2 = 10^(((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(log10(p.YLim))+log10(p.YLim(1)));

rectangle('Position',[bounds(1) bounds(3) bounds(2)-bounds(1) bounds(4)-bounds(3)]);
hold on;
if any(vertex==1)
    plot([bounds(1) x1], [bounds(4) y2], '-.k'); % Line to vertex 1
end
if any(vertex==2)
    plot([bounds(2) x2], [bounds(4) y2], '-.k'); % Line to vertex 2
end
if any(vertex==3)
    plot([bounds(2) x2], [bounds(3) y1], '-.k'); % Line to vertex 4
end
if any(vertex==4)
    plot([bounds(1) x1], [bounds(3) y1], '-.k'); % Line to vertex 3
end
% legend hide
z = axes('position',pos);

box on % put box around new pair of axes
semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec8(plot_ind)),'-sg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec5(plot_ind)),'-.>m','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec3(plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec7(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec4(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind), sqrt(crlb(plot_ind)),'--k','LineWidth',2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);

axis(bounds);

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})

legend('PDP','Two-Step','2Q-Order','MUSIC (0.01^\circ)','MLE (0.01^\circ)','MLE (0.2^\circ)','CLRB');
set(gca,'FontSize',12)


% print -dpng -r600 sim1-8sensors.png

%% Fig 4
% close all
% load fig_8_2_data_3.mat;
% plot_ind = 7:17;

% load fig_8_2_data_5.mat;
% plot_ind = 7:17;

% load fig_8_2_data_8.mat;
plot_ind = 5:17;

linewidth = 1.5;
markersize = 6;
figure;
semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec8(plot_ind)),'-dg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec3(plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec5(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec4(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on; semilogy(snrdB(plot_ind), sqrt(crlb(plot_ind)),'--k','LineWidth',2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);


xlabel('SNR [dB]');
ylabel('RMSE [Deg]');
legend('PDP','Two-Step','MUSIC (0.01^\circ)','MLE (0.01^\circ)','MLE (0.2^\circ)','CRLB');
grid on; set(gca,'FontSize',12)


% print -dpng -r600 sim2-3sensors.png
% print -dpng -r600 sim2-5sensors.png
% print -dpng -r600 sim2-8sensors.png

