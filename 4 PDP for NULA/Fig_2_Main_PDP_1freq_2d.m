% plot wrapped phase difference (WPD) pattern, for 2 pair of sensors
% author: hui.chen@kaust.edu.sa

close all;
clear all;
clc;
%% Simulation Model Building for Phase Unwrapping
global v;
v = 343;
freq = 20000;
D = 2.3;
delta = 1.25;
lambda = v/freq;                % wavelength
S = [0 1 1+delta]*D*lambda/2;     % sensor location, in [m]

comb = combnk(1:length(S),2);
comb = comb([1 3],:);           % pair 1 & 3
d = abs(S(comb(:,1)) - S(comb(:,2)));   % sensor inner-distance of these pairs
thetam = 90;                    % maximum DOA, detection range [-thetam, thetam]
ang_range = (-thetam):(thetam);

angle_train = [];
phase_real = [];
phase_wrapped = [];
phi = zeros(1,length(freq)-1);
p = [];
for x_i = -thetam:1:thetam
    angle0 = x_i;
    phi = 2*pi*freq*sin(angle0/180*pi)*d/v;
    phase_real = [phase_real; (phi)];
    phase_wrapped = [phase_wrapped; wrapping(phi)];
end

% if(length(d)==3)
%     figure;plot3(phase_wrapped(:,1),phase_wrapped(:,2),phase_wrapped(:,3),'.');
%     character = cellstr(num2str((1:141)'-71));
%     text(phase_wrapped(:,1),phase_wrapped(:,2),phase_wrapped(:,3),character);
%     xlabel('x');ylabel('y');zlabel('z');
% elseif(length(d)==2)
%     figure;plot(phase_wrapped(:,1),phase_wrapped(:,2),'.');
%     character = cellstr(num2str((1:141)'-71));
%     text(phase_wrapped(:,1),phase_wrapped(:,2),character);
% end
%%
figure;plot(phase_wrapped(:,1),phase_wrapped(:,2),'.','MarkerSize',10);
hold on; plot(phase_wrapped((91-20):(91+20),1),phase_wrapped((91-20):(91+20),2),'b.','MarkerSize',10);
hold on; plot(phase_wrapped((91+21):(91+25),1),phase_wrapped((91+21):(91+25),2),'m.','MarkerSize',10);
hold on; plot(phase_wrapped((91-25):(91-21),1),phase_wrapped((91-25):(91-21),2),'m.','MarkerSize',10);
hold on; plot(phase_wrapped((91+26):(91+90),1),phase_wrapped((91+26):(91+90),2),'c.','MarkerSize',10);
hold on; plot(phase_wrapped((91-90):(91-26),1),phase_wrapped((91-90):(91-26),2),'c.','MarkerSize',10);

set(gca,'FontSize',16)
xlabel('{$\psi_1(\theta)$ [Rad]}','Interpreter','Latex');
ylabel('{$\psi_2(\theta)$ [Rad]}','Interpreter','Latex');
set(gcf,'position',[100,100,420*1.2,400*1.2])
axis([-4 4 -3.3 3.3])

% projections
[P, U] = alg_get_projections(max(ang_range), d, v, freq);
hold on;
plot(P(:,1),P(:,2),'-g','LineWidth',1.5);
hold on;
plot(P(:,1),P(:,2),'sr','MarkerSize',8,'MarkerFaceColor','r');

character = cell(1,3);
character{1} = {'{$\mathbf{p}_1$}'};
character{2} = {'{$\mathbf{p}_3$}'};
character{3} = {'{$\mathbf{p}_5$}'};
% character{5} = {'p_5'};
text(P([1 3 5],1)+0.25,P([1 3 5],2)+0.2,character,'FontSize',16,'Color','r','Interpreter','Latex');
text(P(2,1)-0.05,P(2,2)-0.4,'{$\mathbf{p}_2$}','FontSize',16,'Color','r','Interpreter','Latex');
text(P(4,1)+0.15,P(4,2)+0.4,'{$\mathbf{p}_4$}','FontSize',16,'Color','r','Interpreter','Latex');

% boundaryc
hold on;
plot([-pi pi],[-pi -pi],'--b','LineWidth',1);
plot([pi pi],[pi -pi],'--b','LineWidth',1);
plot([-pi pi],[pi pi],'--b','LineWidth',1);
plot([-pi -pi],[pi -pi],'--b','LineWidth',1);

% highlight
highlight = [-90 -25 -26 26 25 90];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
% character = cellstr([num2str(highlight')]);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
% s1 = num2str(highlight');
% s2 = '^\circ';
% s = strcat(s1,s2)
% character = {{'-90^\circ'}, {'-25^\circ'}, {'-26^\circ'}, {'26^\circ'}, {'25^\circ'}, {'90^\circ'}}
text(phase_wrapped(highlight+thetam+1,1)+0.1,phase_wrapped(highlight+thetam+1,2)-0.1,character,'FontSize',16,'Interpreter','Latex');

highlight = [-20 0 21];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
% character = cellstr(num2str(highlight'));
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)-0.20,phase_wrapped(highlight+thetam+1,2)-0.35,character,'FontSize',16,'Interpreter','Latex');

set(gca,'FontSize',16)

highlight = [-21 20];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
% character = cellstr(num2str(highlight'));
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)+0.0,phase_wrapped(highlight+thetam+1,2)+0.35,character,'FontSize',16,'Interpreter','Latex');

set(gca,'FontSize',16)
% an example of phase unwrapping
% ex1 = phase_wrapped(22+thetam+1,:)+[-0.3 0.24];
% ex2 = phase_wrapped(22+thetam+1,:);
% ex3 = phase_real(22+thetam+1,:);

% plot(ex1(1),ex1(2),'^r','MarkerSize',5,'LineWidth',1.2);
% text(ex1(1)-0.94,ex1(2)+0.2,'$\hat{\psi}(22)$','Interpreter','latex','FontSize',12,'Color','r');
% plot(ex2(1),ex2(2),'^r','MarkerSize',5,'LineWidth',1.2);
% text(ex2(1)+0.2,ex2(2)-.1,'$\tilde{\psi}(22)$','Interpreter','latex','FontSize',12,'Color','r');
% plot(ex3(1),ex3(2),'^r','MarkerSize',5,'LineWidth',1.2);
% text(ex3(1)+0.2,ex3(2),'$\hat{\phi}(22)$','Interpreter','latex','FontSize',12,'Color','r');

% plot(-1.5,-2.9,'xb','MarkerSize',7,'LineWidth',1.2);
% text(-1.3,-2.9,'$\hat{\psi}(-15)$','Interpreter','latex','FontSize',12,'Color','b');

text('position',[P(1,1)+1, P(1,2)+1],'Interpreter','latex','String','{\boldmath$\Theta_1$}','FontSize',16,'Color','c');
text('position',[P(5,1)+1, P(5,2)+1],'Interpreter','latex','String','{\boldmath$\Theta_5$}','FontSize',16,'Color','c');
text('position',[P(3,1)+1, P(3,2)+1],'Interpreter','latex','String','{\boldmath$\Theta_3$}','FontSize',16,'Color','b');
text('position',[P(4,1)-0.9, P(4,2)-.1],'Interpreter','latex','String','{\boldmath$\Theta_4$}','FontSize',16,'Color','m');
text('position',[P(2,1)-.7, P(2,2)+0.3],'Interpreter','latex','String','{\boldmath$\Theta_2$}','FontSize',16,'Color','m');

axis([-4 4 -4 4])

set(gcf, 'Position', [10 10 500 450])


% xlabel('{Number of Antennas/SAs at BS $N_\mathrm{B}$}','Interpreter','Latex')
% 

% print -dpng -r600 fig_2.png
