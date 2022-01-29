clear;clc;close all;

addpath(genpath('..\CR3BP\'));
addpath(genpath('..\conversions\'));
addpath(genpath('..\generic libs\'));

% - mass of primary
m1 = 5.97219e24;
% - mass of secondary
m2 = 7.34767e22;

% - mean distance 
r12 = 384400;

% - Tolerance
tol = 1e-6;
% - Max iteration
iterMax = 100;

% - computing characteristic values
[mu,m_star,t_star,l_star,v_star,n_star] = getCharValues_CR3BP(m1,m2,r12);

% - Lagrange Points computation
LP = getLagrangePoints_CR3BP(mu);

%% -- Initialzing Database

DATASET = table([],[],[],[],[],[],[],[],...
    'VariableNames',{'Center','Family','Orient','SV','T','C','nu','k'});

%% -- Initializing Plots

figure();
title('Synodic Frame');
xlabel('x-axis [nd]');
ylabel('y-axis [nd]');
zlabel('z-axis [nd]');
grid on;
hold on;
axis equal;
plot3(-mu,0,0,'.k');
text(-mu,0,0,'  Earth');
plot3(1-mu,0,0,'.k');
text(1-mu,0,0,'  Moon');
plot3(LP(1,4),LP(2,4),LP(3,4),'.g');
text(LP(1,4),LP(2,4),LP(3,4),sprintf('L%g',4));

%% -- Computing Initial Orbit from Guess Solution

% -- Guess State and Time period [nd]
SV0 = [0.080788303911618;0.549925532615439;0.807709030651238;-0.214506093111874;0.529406107790107;-0.249472850183398];
tf = 6.29583971028002;

% - General Info
Family = 'Axial';
Orientation = 'North';
Center = 'L4';

% -- Correction of Initial Guess
[SV0,tf] = newton_assymPeriodic_3DXY_fixedZ_CR3BP(SV0,tf,mu,tol);

% -- Computing Solution
[SV,t,STM] = propagateTrajectory_CR3BP(SV0,tf,mu,0,...
    'ComputeSTM',true,'InitialSTM',eye(6,6));

% - Ploting solution
plot3(SV(1,:),SV(2,:),SV(3,:),'b');

% - Monodromy
M = STM(:,:,end);
% - Stability Analysis
[nu,k] = getStability_CR3BP(M);
% - Jacobi Constant
C = jacobiConstant(SV0,mu);

% -- Adding to family
DATASET = [DATASET;{Center,Family,Orientation,SV0',tf,C,nu,k}];

% -- Adding SV0 and tf to origin values
SV0_origin = SV0;
tf_origin = tf;

%% --- Pseudo Arc-Length Continuaiton (positive Step)

% - Step size
DeltaS = 0.025;
% - Number of Orbits to compute
N = 10;


% - performing continuation
[SV0,tf] = continuation_PAL_CR3BP(SV0_origin,tf_origin,DeltaS,N,mu,6,...
    'DirectionalIncrement',true,'TargetVector',1,'TargetDirection',+1);

% - Computing orbits
for n = 1:N
    % - computing the 3-body problem solution
    [SV_nd,t_nd,STM] = propagateTrajectory_CR3BP(SV0(:,n),tf(n),mu,0,...
        'ComputeSTM',true,'InitialSTM',eye(6,6));
    % - Monodromy
    M = STM(:,:,end);
    % - Stability Analysis
    [nu,k,~,BC,region] = getStability_CR3BP(M);
    % - Jacobi Constant
    C = jacobiConstant(SV0(:,n),mu);
    % - Adding Set to Database
    DATASET = [DATASET;{Center,Family,Orientation,SV0(:,n)',tf(n),C,nu,k}];
    % - Plotting orbit
    plot3(SV_nd(1,:),SV_nd(2,:),SV_nd(3,:),'m');
    drawnow;
end

%% --- Pseudo Arc-Length Continuaiton (negetive Step)

% - Step size
DeltaS = -0.025;
% - Number of Orbits to compute
N = 53;

% - performing continuation
[SV0,tf] = continuation_PAL_CR3BP(SV0_origin,tf_origin,DeltaS,N,mu,6,...
    'DirectionalIncrement',true,'TargetVector',1,'TargetDirection',-1);

% - Computing orbits
for n = 1:N
    % - computing the 3-body problem solution
    [SV_nd,t_nd,STM] = propagateTrajectory_CR3BP(SV0(:,n),tf(n),mu,0,...
        'ComputeSTM',true,'InitialSTM',eye(6,6));
    % - Monodromy
    M = STM(:,:,end);
    % - Stability Analysis
    [nu,k,~,BC,region] = getStability_CR3BP(M);
    % - Jacobi Constant
    C = jacobiConstant(SV0(:,n),mu);
    % - Adding Set to Database
    DATASET = [{Center,Family,Orientation,SV0(:,n)',tf(n),C,nu,k};DATASET];
    % - Plotting orbit
    plot3(SV_nd(1,:),SV_nd(2,:),SV_nd(3,:),'m');
    drawnow;
end


SV_last = SV0(:,n);
tf_last = tf(n);

% - Step size
DeltaS = -0.025;
% - Number of Orbits to compute
N = 45;

% - performing continuation
[SV0,tf] = continuation_PAL_CR3BP(SV_last,tf_last,DeltaS,N,mu,6,...
    'DirectionalIncrement',true,'TargetVector',3,'TargetDirection',-1);

% - Computing orbits
for n = 1:N
    % - computing the 3-body problem solution
    [SV_nd,t_nd,STM] = propagateTrajectory_CR3BP(SV0(:,n),tf(n),mu,0,...
        'ComputeSTM',true,'InitialSTM',eye(6,6));
    % - Monodromy
    M = STM(:,:,end);
    % - Stability Analysis
    [nu,k,~,BC,region] = getStability_CR3BP(M);
    % - Jacobi Constant
    C = jacobiConstant(SV0(:,n),mu);
    % - Adding Set to Database
    DATASET = [{Center,Family,Orientation,SV0(:,n)',tf(n),C,nu,k};DATASET];
    % - Plotting orbit
    plot3(SV_nd(1,:),SV_nd(2,:),SV_nd(3,:),'m');
    drawnow;
end


%% --- Saving Database

save('..\database\DATASET_RAW_L4Axial_North.mat','DATASET');