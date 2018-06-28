classdef Particle
    %PARTICLE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ID = -1;
        Lumi = 1;
        Diameter = 10;
        Life = 100000000000000;
        InWindow = 0;
        Spawned = 0;
        RecentlySpawned = 0;
        XPos
        YPos
        XTrack = {}
        YTrack = {}
    end

    methods
        function particle = Particle(diameter, ID, lumi0)
            if nargin == 1
                particle.Diameter = round(diameter);
            end
            if nargin == 2
                particle.Diameter = round(diameter);
                particle.ID = ID;
            end
            if nargin == 3
                particle.Diameter = round(diameter);
                particle.ID = ID;
                particle.Lumi = lumi0;
            end
%             if isnumeric(val)
%                obj.Value = val;
%             else
%                error('Value must be numeric')
%             end
        end

        function obj = set.ID(obj, ID)
            obj.ID = ID;
        end
        function ID = get.ID(obj)
            ID = obj.ID;
        end
        function obj = set.Lumi(obj, lumi)
            obj.Lumi = lumi;
        end
        function lumi = get.Lumi(obj)
            lumi = obj.Lumi;
        end
        function obj = set.Diameter(obj, diameter)
            obj.Diameter = diameter;
        end
        function diameter = get.Diameter(obj)
            diameter = obj.Diameter;
        end
        function obj = set.Life(obj, life)
            obj.Life = life;
        end
        function life = get.Life(obj)
            life = obj.Life;
        end
        function obj = set.InWindow(obj, inWindow)
            obj.InWindow = inWindow;
        end
        function inWindow = get.InWindow(obj)
            inWindow = obj.InWindow;
        end
        function obj = set.Spawned(obj, spawned)
            obj.Spawned = spawned;
        end
        function spawned = get.Spawned(obj)
            spawned = obj.Spawned;
        end
        function obj = set.RecentlySpawned(obj, recentlySpawned)
            obj.RecentlySpawned = recentlySpawned;
        end
        function recentlySpawned = get.RecentlySpawned(obj)
            recentlySpawned = obj.RecentlySpawned;
        end
        function obj = set.XPos(obj, xPos)
            obj.XPos = xPos;
        end
        function xPos = get.XPos(obj)
            xPos = obj.XPos;
        end
        function obj = set.YPos(obj, yPos)
            obj.YPos = yPos;
        end
        function yPos = get.YPos(obj)
            yPos = obj.YPos;
        end
        function obj = set.XTrack(obj, xTrack)
            obj.XTrack = xTrack;
        end
        function xTrack = get.XTrack(obj)
            xTrack = obj.XTrack;
        end
        function obj = set.YTrack(obj, yPos)
            obj.YTrack = yPos;
        end
        function yTrack = get.YTrack(obj)
            yTrack = obj.YTrack;
        end
        
        function pl = luminosity(obj)
            x = linspace(-obj.Diameter, obj.Diameter, obj.Diameter);
            y = linspace(-obj.Diameter, obj.Diameter, obj.Diameter);
            [X, Y] = meshgrid(x,y);
            pl = obj.Lumi*exp((-(X).^2-(Y).^2)/((1/8.0)*obj.Diameter.^2));
        end
    end
end

