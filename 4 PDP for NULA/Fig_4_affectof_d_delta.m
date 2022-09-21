% Search based phase unwrapping algorithm tests
close all;
clear all;
clc;
%%
% project_3d([a b c], p)
%% Simulation Model Building for Phase Unwrapping
global v;
v = 343;
freq = 20000;
D = 1.8;
delta = 1.25;
lambda = v/freq;
s = [0 1 1+delta]*D*lambda; % sensor location, in [lambda]
comb = combnk(1:length(s),2);
comb = comb([1 3],:);
d = abs(s(comb(:,1)) - s(comb(:,2)));
thetam = 90;
ang_range = (-thetam):(thetam);
Ns = floor(2*freq*sin(thetam/180*pi)*d/v)

angle_train = [];
phase_real = [];
phase_wrapped = [];
phi = zeros(1,length(freq)-1);
p = [];
for x_i = -thetam:1:thetam
    angle0 = x_i;
    % angle = -33.3
    phi = 2*pi*freq*sin(angle0/180*pi)*d/v;
    phase_real = [phase_real; (phi)];
    phase_wrapped = [phase_wrapped; wrapping(phi)];
    % phi = sin(phi);
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


figure;plot(phase_wrapped(:,1),phase_wrapped(:,2),'.','MarkerSize',10);
set(gca,'FontSize',16)
xlabel('{$\psi_1(\theta)$ [Rad]}','Interpreter','Latex');
ylabel('{$\psi_2(\theta)$ [Rad]}','Interpreter','Latex');
set(gcf,'position',[100,100,420*1.2,400*1.2])
axis([-4 4 -4 4])
% projections
[P, U] = alg_get_projections(max(ang_range), d, v, freq);
hold on;
plot(P(:,1),P(:,2),'-g','LineWidth',1.5);
hold on;
plot(P(:,1),P(:,2),'sr','MarkerSize',8,'MarkerFaceColor','r');
character = [];
for i = 1:length(P(:,1))
    character = [character; ['{$' '\mathbf{p}_' num2str(i) '$}']]
end
character = cellstr(character);
text(P(:,1)-0.07,P(:,2)-0.3,character,'FontSize',16,'Color','r', 'Interpreter', 'Latex');

% boundary
hold on;
plot([-pi pi],[-pi -pi],'--b','LineWidth',1);
plot([pi pi],[pi -pi],'--b','LineWidth',1);
plot([-pi pi],[pi pi],'--b','LineWidth',1);
plot([-pi -pi],[pi -pi],'--b','LineWidth',1);
% highlight


highlight = [-41 -12 13 42 -90];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)-0.2,phase_wrapped(highlight+thetam+1,2)-0.4,character,'FontSize',16, 'Interpreter', 'Latex');
set(gca,'FontSize',14)

highlight = -[-41 -12 13 42 -90];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)-0.2,phase_wrapped(highlight+thetam+1,2)+0.4,character,'FontSize',16, 'Interpreter', 'Latex');
set(gca,'FontSize',14)

highlight = [0];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)+0.2,phase_wrapped(highlight+thetam+1,2)+0.1,character,'FontSize',16, 'Interpreter', 'Latex');
set(gca,'FontSize',14)

highlight = [-57 -17 16 56];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)+0.17,phase_wrapped(highlight+thetam+1,2)-0.1,character,'FontSize',16, 'Interpreter', 'Latex');
set(gca,'FontSize',14)

highlight = -[-57 -17 16 56];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)-0.77,phase_wrapped(highlight+thetam+1,2)-0.1,character,'FontSize',16, 'Interpreter', 'Latex');
set(gca,'FontSize',14)

% print -dpng -r600 fig-3-1.png

%% Simulation Model Building for Phase Unwrapping
freq = 20000;
D = 1.15;
delta = 1.45;
lambda = v/freq;
s = [0 1 1+delta]*D*lambda; % sensor location, in [lambda]
comb = combnk(1:length(s),2);
comb = comb([1 3],:);
d = abs(s(comb(:,1)) - s(comb(:,2)));
thetam = 90;
ang_range = (-thetam):(thetam);
Ns = floor(2*freq*sin(thetam/180*pi)*d/v);

% distance = 36000*10^3;
lambda = v/freq;
angle_train = [];
phase_real = [];
phase_wrapped = [];
phi = zeros(1,length(freq)-1);
p = [];
for x_i = -thetam:1:thetam
    angle0 = x_i;
    % angle = -33.3
    phi = 2*pi*freq*sin(angle0/180*pi)*d/v;
    phase_real = [phase_real; (phi)];
    phase_wrapped = [phase_wrapped; wrapping(phi)];
    % phi = sin(phi);
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


figure;plot(phase_wrapped(:,1),phase_wrapped(:,2),'.','MarkerSize',10);
set(gca,'FontSize',16)
xlabel('{$\psi_1(\theta)$ [Rad]}','Interpreter','Latex');
ylabel('{$\psi_2(\theta)$ [Rad]}','Interpreter','Latex');
set(gcf,'position',[100,100,420*1.2,400*1.2])
axis([-4 4 -4 4])
% projections
[P, U] = alg_get_projections(max(ang_range), d, v, freq);
hold on;
plot(P(:,1),P(:,2),'-g','LineWidth',1.5);
hold on;
plot(P(:,1),P(:,2),'sr','MarkerSize',8,'MarkerFaceColor','r');
% character = [];

for i = 3:7
    character = cellstr(['{$' '\mathbf{p}_' num2str(i) '$}']);
    text(P(i,1)-0.15,P(i,2)-0.4,character,'FontSize',16,'Color','r','Interpreter','Latex');
end

for i = 1:2
    character = cellstr(['{$' '\mathbf{p}_' num2str(i) '$}']);
    text(P(i,1)-0.3,P(i,2)+0.42,character,'FontSize',16,'Color','r','Interpreter','Latex');
end


% boundary
hold on;
plot([-pi pi],[-pi -pi],'--b','LineWidth',1);
plot([pi pi],[pi -pi],'--b','LineWidth',1);
plot([-pi pi],[pi pi],'--b','LineWidth',1);
plot([-pi -pi],[pi -pi],'--b','LineWidth',1);
% highlight


highlight = [0];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)+0.2,phase_wrapped(highlight+thetam+1,2)+0.1,character,'FontSize',16,'Interpreter','Latex');
set(gca,'FontSize',14)


highlight = -[-18 17 -65 -64];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)-0.2,phase_wrapped(highlight+thetam+1,2)-0.4,character,'FontSize',16,'Interpreter','Latex');
set(gca,'FontSize',14)
% 
highlight = [-18 17 -65 -64];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)-0.35,phase_wrapped(highlight+thetam+1,2)+0.4,character,'FontSize',16,'Interpreter','Latex');
set(gca,'FontSize',14)


highlight = [25 -26];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)+0.17,phase_wrapped(highlight+thetam+1,2)-0.1,character,'FontSize',16,'Interpreter','Latex');
set(gca,'FontSize',14)

highlight = -[25 -26];
hold on;plot(phase_wrapped(highlight+thetam+1,1),phase_wrapped(highlight+thetam+1,2),'.r','MarkerSize',16);
character = cellstr(strcat([num2str(highlight')], '^\circ'));
for i = 1:length(character)
    character{i} = ['{$' character{i} '$}'];
end
text(phase_wrapped(highlight+thetam+1,1)-0.82,phase_wrapped(highlight+thetam+1,2)+0.15,character,'FontSize',16,'Interpreter','Latex');
set(gca,'FontSize',14)

% print -dpng -r600 fig-3-2.png
