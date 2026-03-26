 clc
clear
close all

addpath(genpath('../'))

%% setup
% flags and constants
% ALL_IN_DISC = true : states are mapped from the real world to the
%                      disc world, and nominal and CBF-QP controllers are
%                      evaluated in the disc world (this way we can
%                      literally get the obstacles out of the way)
% ALL_IN_BALL = false : nominal controllers are evaluated in the real
%                       world and mapped to the ball world through the
%                       jacobian of the diffeomorphism, then CBF-QP
%                       controller is evaluated in the ball world and
%                       mapped back to the real world through the
%                       pseudo-inverse Jacobian
ALL_IN_BALL = true;
DT = 0.001;
T_MAX = 1e4;
POSITION_RADIUS = 'equal'; % 'move' or 'radius' or 'equal'
LAMBDA = 1e0;

% build-up world
realWorld.domain.type = 'qc';
realWorld.domain.contour = 2.5*[-1 1 1 -1;
    -1 -1 1 1];
realWorld.domain.goal = [0;-2];
realWorld.obstacles{1}.type = 'qc';
realWorld.obstacles{1}.contour = [-0.5; 0.25] + 1.1*[-0.5 -0.4 -0.4 0.4 0.4 0.5 0.5 -0.5;
    -0.5 -0.5 0 0 -0.5 -0.5 0.1 0.1];
ballWorld.domain.center = [0;0];
ballWorld.domain.radius = 5;
ballWorld.domain.goal = realWorld.domain.goal;
ballWorld.obstacles{1}.center = [-0.95;0.25];
ballWorld.obstacles{1}.centerOriginal = ballWorld.obstacles{1}.center;
ballWorld.obstacles{1}.radius = 0.75;
ballWorld.obstacles{1}.radiusOriginal = ballWorld.obstacles{1}.radius;

%% init
% real world
wm = WorldMapping(realWorld, ballWorld);
wm.evaluateMappings(LAMBDA)
[r2bMap, b2rMap, r2bJac, b2rJac] = wm.getMappings();

% robot, goal, and trajectory variables
x = [2;2];
xBall = r2bMap(x);
xG = [-0.5;-1];
xGBall = r2bMap(xG);
xTraj = zeros(2,T_MAX);
xTraj(:,1) = x;
xTrajBall = zeros(2,T_MAX);
xTrajBall(:,1) = xBall;

% plots
figure('units','normalized','position',[0 0 1 1])

subplot(1,2,1), hold on, axis equal, axis([-3 3 -3 3]), set(gca, 'Visible', 'off')
hS = scatter(x(1), x(2), 1000, '.');
hSTraj = line(x(1), x(2), 'LineWidth', 2);
plot(realWorld.domain.contour(1,[1:end,1]), realWorld.domain.contour(2,[1:end,1]), 'LineWidth', 2)
for i = 1 : numel(realWorld.obstacles)
    plot(realWorld.obstacles{i}.contour(1,[1:end,1]), realWorld.obstacles{i}.contour(2,[1:end,1]), 'LineWidth', 2)
end
scatter(xG(1), xG(2), 1000, '.');

subplot(1,2,2), hold on, axis equal, axis([-6 6 -6 6]), set(gca, 'Visible', 'off')
hSBall = scatter(xBall(1), xBall(2), 1000, '.');
hSTrajBall = line(xBall(1), xBall(2), 'LineWidth', 2);
domainBall = ballWorld.domain.center + ballWorld.domain.radius * [cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))];
plot(domainBall(1,[1:end,1]), domainBall(2,[1:end,1]), 'LineWidth', 2)
hObstBall = cell(1,numel(ballWorld.obstacles));
for i = 1 : numel(ballWorld.obstacles)
    obstBall = ballWorld.obstacles{i}.center + ballWorld.obstacles{i}.radius * [cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))];
    hObstBall{i} = plot(obstBall(1,[1:end,1]), obstBall(2,[1:end,1]), 'LineWidth', 2);
end
scatter(xGBall(1), xGBall(2), 1000, '.');

drawnow












    %%% controller synthesis
    % robot does not care about the obstacles
    uBall = uNomBall;
    % obstacles care about the robot
    safetyScale = 3;
    roimin = 0.5;
    xDotObstNomBall = zeros(2*Nobst,1);
    rDotObstNomBall = zeros(Nobst,1);
    Acbf = [];
    bcbf = [];
    for i = 1 : Nobst
        xoiHat = ballWorld.obstacles{i}.centerOriginal;
        roiHat = ballWorld.obstacles{i}.radiusOriginal;
        xoi = ballWorld.obstacles{i}.center;
        roi = ballWorld.obstacles{i}.radius;
        xd = ballWorld.domain.center;
        rd = ballWorld.domain.radius;
        xDotObstNomBall(2*i-1:2*i) = 100*(xoiHat - xoi);
        rDotObstNomBall(i) = 100*(roiHat - roi);
        % do not collide with robot
        Acbf(end+1,[2*i-1:2*i,2*Nobst+i]) = [-2*(xoi-xBall)', 2*safetyScale*roi];
        bcbf(end+1,1) = -2*(xoi-xBall)'*uBall + 1e2 * (norm(xoi-xBall)^2-(safetyScale*roi)^2);
        % do not go out of environment
        Acbf(end+1,[2*i-1:2*i,2*Nobst+i]) = [2*(xoi-xd)', 2*(rd-roi)];
        bcbf(end+1,1) = 1e2 * ((rd-roi)^2-norm(xoi-xd)^2);
        % do not collide with other obstacles
        for j = i+1 : Nobst
            xoj = ballWorld.obstacles{j}.center;
            roj = ballWorld.obstacles{j}.radius;
            Acbf(end+1,[2*i-1:2*i,2*j-1:2*j,2*Nobst+i,2*Nobst+j]) = [-2*(xoi-xoj)', 2*(xoi-xoj)', 2*(roi+roj), 2*(roi+roj)];
            bcbf(end+1,1) = 1e2 * (norm(xoi-xoj)^2-(roi+roj)^2);
        end
        % keep minimum obstacle radius
        Acbf(end+1,2*Nobst+i) = -2*roi;
        bcbf(end+1,1) = 1e2 * (roi^2-roimin^2);
    end
    if strcmp(POSITION_RADIUS, 'move')
        Wcenter = 1;
        Wradius = 1e3;
    elseif strcmp(POSITION_RADIUS, 'radius')
        Wcenter = 1e3;
        Wradius = 1;
    elseif strcmp(POSITION_RADIUS, 'equal')
        Wcenter = 1;
        Wradius = 1;
    end
    xDotRDotObstBall = quadprog(2*blkdiag(Wcenter*eye(2*Nobst), Wradius*eye(Nobst)), ...
        -2*[Wcenter*xDotObstNomBall; Wradius*rDotObstNomBall]', ...
        Acbf, ...
        bcbf, ...
        [],[],[],[],[],optimoptions(@quadprog,'Display','off'));
    for i = 1 : Nobst
        xDotI = xDotRDotObstBall(2*i-1:2*i);
        rDotI = xDotRDotObstBall(2*Nobst+i);
        ballWorld.obstacles{i}.center = ballWorld.obstacles{i}.center + xDotI*DT;
        ballWorld.obstacles{i}.radius = ballWorld.obstacles{i}.radius + rDotI*DT;
    end
    















































