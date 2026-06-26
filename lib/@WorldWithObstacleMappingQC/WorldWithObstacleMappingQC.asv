classdef WorldWithObstacleMappingQC < handle
    properties (Access=private)
        interior = false;
        vertices
        fx
        fy
        R
        xPrev
    end
    methods
        % constructor
        function this = WorldWithObstacleMappingQC(vertices, ballDomain, ballObstacles)
            % function this = WorldWithObstacleMappingQC(vertices, ballDomain, ballObstacles, ballRobot)
            this.vertices = polyshape(vertices{1}, vertices{2}); % construct a polygon with the given list of coordinates
            tr = triangulation(this.vertices); % %% a default (delunay?) triangulation of the polygon, connectivity list T not pre-determined
            model = createpde; % create empty scalar pde structure (1 equation) 
            geometryFromMesh(model, tr.Points', tr.ConnectivityList'); % define geometry info for model: faces, edges, No. of vertices and 3d coordinates of vertices, this automatically generates a mesh
            generateMesh(model,'Hmax',0.1); % h_max: maximum element size is a tuning parameter (precision VS speed) 
            %% but why only sample h_max from vector [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6]?
            % update model.Mesh with given maximum element size h_max
            v = model.Mesh.Nodes';
            f = model.Mesh.Elements';
            f = f(:,1:3); % Why only consider first three columns? (last three not useful?)
            [f,v] = gpp_clean_mesh(f,v); % Cleaning mesh and left with less necessary repeated vertices and faces
            bd = meshboundaries(f); % check meshboundaries.m
            c = [ballDomain.center'; cell2mat(cellfun(@(x) x.center, ballObstacles, 'UniformOutput', false))'];
            r = [ballDomain.radius; cell2mat(cellfun(@(x) x.radius, ballObstacles, 'UniformOutput', false))'];
            % c_robot = [ballDomain.center'; cell2mat(cellfun(@(x) x.center, ballRobot, 'UniformOutput', false))'];
            % r_robot = [ballDomain.radius; cell2mat(cellfun(@(x) x.radius, ballRobot, 'UniformOutput', false))'];
            if numel(bd) > 2
                idx_ordered = zeros(numel(bd)-1,1);
                for i = 2 : numel(bd)
                    for j = 1 : numel(vertices{1})
                        if inpolygon(v(bd{i},1), v(bd{i},2), vertices{1}{j}, vertices{2}{j})
                            idx_ordered(i-1) = j;
                        end
                    end
                end
                c(2:end,:) = c(idx_ordered,:);
                r(2:end) = r(idx_ordered);
            end
            map = polygonal_to_ball_qc_map_prescribed_holes(v, f, bd, c, r); % polygonal_to_ball_qc_map_prescribed_holes.m
            Nmax = size(map,1);
            if size(map,1) > size(v,1)
                Nmax = size(map,1)-1;
            end
            this.fx = TriScatteredInterp(v(:,1),v(:,2),map(1:Nmax,1));
            this.fy = TriScatteredInterp(v(:,1),v(:,2),map(1:Nmax,2));
        end
        function storeState(this, statePrev)
            this.xPrev = statePrev;
        end
        % getters
        function mapHandle = getReal2BallMapHandle(this)
            mapHandle = @(q) this.real2ball(q);
        end
        function mapHandle = getBall2RealMapHandle(this)
            mapHandle = @(xBall, xRealPrev) this.ball2real(xBall, xRealPrev);
        end
        function mapHandle = getBall2RealMapHandle_new(this)
            mapHandle = @(xBall) this.ball2real_new(xBall);
        end
    end
    methods (Access=private)
        function xBall = real2ball(this, xReal)
            xBall = [this.fx(xReal'); this.fy(xReal')];
        end
        function xReal = ball2real(this, xBall, xRealPrev)
            options = optimset('MaxFunEvals', 100, ...
                'MaxIter', 100, ...
                'TolFun', 1e-6, ...
                'TolX', 1e-6, ...
                'Display', 'off');
            directMapping = @(p) this.real2ball(p)-xBall;
            xReal = fsolve(directMapping, xRealPrev, options);
        end
        function xReal = ball2real_new(this, xBall)
            options = optimset('MaxFunEvals', 100, ...
                'MaxIter', 100, ...
                'TolFun', 1e-6, ...
                'TolX', 1e-6, ...
                'Display', 'off');
            directMapping = @(p) this.real2ball(p)-xBall;
            xReal = fsolve(directMapping, this.xPrev, options);
        end
    end
end