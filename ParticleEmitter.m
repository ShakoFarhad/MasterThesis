classdef ParticleEmitter < handle
    %PARTICLEEMITTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        AmountOfParticles = 10;
        BGX = 50;
        BGY = 50;
        VX
        VY
        EmitXArray = [1,1];
        EmitYArray = [1,1];
        ParticleArray = -1;
        ParticleSpawnRate = 1;
        LocalXVelocity = 0;
        LocalYVelocity = 0;
    end
    
    properties(Access = private)
        Window = -1;
        ParticleCount
        ParticleSpawnNow = 1;
        PadWindow = 0;
        ParticleSpawnCount = 0;
    end
    

    methods
        function particleE = ParticleEmitter(AmountOfParticles, particleSpawnRate, BGX, BGY, emitXArray, emitYArray, VX, VY, localXVelocity, localYVelocity)
            if nargin == 0
                mu = [0 0]; Sigma = [2 1; 1 2];
                x1 = linspace(-2,2,particleE.BGX); x2 = linspace(-2,2,particleE.BGY);
                [X1,X2] = meshgrid(x1,x2);
                F = mvnpdf([X1(:) X2(:)],mu,Sigma)*(particleE.BGX+particleE.BGY)/2;
                particleE.VX = reshape(F,length(x2),length(x1));
                particleE.VY = reshape(F,length(x2),length(x1));
            end
            if nargin == 1
                particleE.AmountOfParticles = AmountOfParticles;
                mu = [0 0]; Sigma = [2 1; 1 2];
                x1 = linspace(-2,2,particleE.BGX); x2 = linspace(-2,2,particleE.BGY);
                [X1,X2] = meshgrid(x1,x2);
                F = mvnpdf([X1(:) X2(:)],mu,Sigma)*(particleE.BGX+particleE.BGY)/2;
                particleE.VX = reshape(F,length(x2),length(x1));
                particleE.VY = reshape(F,length(x2),length(x1));
            end
            if nargin == 2
                particleE.AmountOfParticles = AmountOfParticles;
                particleE.ParticleSpawnRate = particleSpawnRate;
                mu = [0 0]; Sigma = [2 1; 1 2];
                x1 = linspace(-2,2,particleE.BGX); x2 = linspace(-2,2,particleE.BGY);
                [X1,X2] = meshgrid(x1,x2);
                F = mvnpdf([X1(:) X2(:)],mu,Sigma)*(particleE.BGX+particleE.BGY)/2;
                particleE.VX = reshape(F,length(x2),length(x1));
                particleE.VY = reshape(F,length(x2),length(x1));
            end
            if nargin == 3
                particleE.AmountOfParticles = AmountOfParticles;
                particleE.ParticleSpawnRate = particleSpawnRate;
                particleE.BGX = BGX;
                particleE.BGY = BGX;
                mu = [0 0]; Sigma = [2 1; 1 2];
                x1 = linspace(-2,2,BGX); x2 = linspace(-2,2,BGY);
                [X1,X2] = meshgrid(x1,x2);
                F = mvnpdf([X1(:) X2(:)],mu,Sigma)*(BGX+BGY)/2;
                particleE.VX = reshape(F,length(x2),length(x1));
                particleE.VY = reshape(F,length(x2),length(x1));
            end
            if nargin == 4
                particleE.AmountOfParticles = AmountOfParticles;
                particleE.ParticleSpawnRate = particleSpawnRate;
                particleE.BGX = BGX;
                particleE.BGY = BGY;
                mu = [0 0]; Sigma = [2 1; 1 2];
                x1 = linspace(-2,2,BGX); x2 = linspace(-2,2,BGY);
                [X1,X2] = meshgrid(x1,x2);
                F = mvnpdf([X1(:) X2(:)],mu,Sigma)*(BGX+BGY)/2;
                particleE.VX = reshape(F,length(x2),length(x1));
                particleE.VY = reshape(F,length(x2),length(x1));
            end
            if nargin == 5
                particleE.AmountOfParticles = AmountOfParticles;
                particleE.ParticleSpawnRate = particleSpawnRate;
                particleE.BGX = BGX;
                particleE.BGY = BGY;
                
                emitYArray = emitXArray;
                if(emitXArray(1) < 1)
                    emitXArray(1) = 1;
                end
                if(emitXArray(2) > particleE.BGX)
                    emitXArray(2) = particleE.BGX;
                end
                if(emitYArray(1) < 1)
                    emitYArray(1) = 1;
                end
                if(emitYArray(2) > particleE.BGY)
                    emitYArray(2) = particleE.BGY;
                end
                particleE.EmitXArray = emitXArray;
                particleE.EmitYArray = emitYArray;
                
                mu = [0 0]; Sigma = [2 1; 1 2];
                x1 = linspace(-2,2,BGX); x2 = linspace(-2,2,BGY);
                [X1,X2] = meshgrid(x1,x2);
                F = mvnpdf([X1(:) X2(:)],mu,Sigma)*(BGX+BGY)/2;
                particleE.VX = reshape(F,length(x2),length(x1));
                particleE.VY = reshape(F,length(x2),length(x1));
            end
            if nargin == 6
                particleE.AmountOfParticles = AmountOfParticles;
                particleE.ParticleSpawnRate = particleSpawnRate;
                particleE.BGX = BGX;
                particleE.BGY = BGY;
                
                if(emitXArray(1) < 1)
                    emitXArray(1) = 1;
                end
                if(emitXArray(2) > particleE.BGX)
                    emitXArray(2) = particleE.BGX;
                end
                if(emitYArray(1) < 1)
                    emitYArray(1) = 1;
                end
                if(emitYArray(2) > particleE.BGY)
                    emitYArray(2) = particleE.BGY;
                end
                particleE.EmitXArray = emitXArray;
                particleE.EmitYArray = emitYArray;
                
                mu = [0 0]; Sigma = [2 1; 1 2];
                x1 = linspace(-2,2,BGX); x2 = linspace(-2,2,BGY);
                [X1,X2] = meshgrid(x1,x2);
                F = mvnpdf([X1(:) X2(:)],mu,Sigma)*(BGX+BGY)/2;
                particleE.VX = reshape(F,length(x2),length(x1));
                particleE.VY = reshape(F,length(x2),length(x1));
            end
            if nargin == 7
                particleE.AmountOfParticles = AmountOfParticles;
                particleE.ParticleSpawnRate = particleSpawnRate;
                particleE.BGX = BGX;
                particleE.BGY = BGY;
                if(emitXArray(1) < 1)
                    emitXArray(1) = 1;
                end
                if(emitXArray(2) > particleE.BGX)
                    emitXArray(2) = particleE.BGX;
                end
                if(emitYArray(1) < 1)
                    emitYArray(1) = 1;
                end
                if(emitYArray(2) > particleE.BGY)
                    emitYArray(2) = particleE.BGY;
                end
                particleE.EmitXArray = emitXArray;
                particleE.EmitYArray = emitYArray;
                particleE.VX = VX;
                particleE.VY = VX;
            end
            if nargin == 8
                particleE.AmountOfParticles = AmountOfParticles;
                particleE.ParticleSpawnRate = particleSpawnRate;
                particleE.BGX = BGX;
                particleE.BGY = BGY;
                if(emitXArray(1) < 1)
                    emitXArray(1) = 1;
                end
                if(emitXArray(2) > particleE.BGX)
                    emitXArray(2) = particleE.BGX;
                end
                if(emitYArray(1) < 1)
                    emitYArray(1) = 1;
                end
                if(emitYArray(2) > particleE.BGY)
                    emitYArray(2) = particleE.BGY;
                end
                particleE.EmitXArray = emitXArray;
                particleE.EmitYArray = emitYArray;
                particleE.VX = VX;
                particleE.VY = VY;
            end
            if nargin == 9
                particleE.AmountOfParticles = AmountOfParticles;
                particleE.ParticleSpawnRate = particleSpawnRate;
                particleE.BGX = BGX;
                particleE.BGY = BGY;
                if(emitXArray(1) < 1)
                    emitXArray(1) = 1;
                end
                if(emitXArray(2) > particleE.BGX)
                    emitXArray(2) = particleE.BGX;
                end
                if(emitYArray(1) < 1)
                    emitYArray(1) = 1;
                end
                if(emitYArray(2) > particleE.BGY)
                    emitYArray(2) = particleE.BGY;
                end
                particleE.EmitXArray = emitXArray;
                particleE.EmitYArray = emitYArray;
                particleE.VX = VX;
                particleE.VY = VY;
                particleE.LocalXVelocity = localXVelocity;
                particleE.LocalYVelocity = localXVelocity;
            end
            if nargin == 10
                particleE.AmountOfParticles = AmountOfParticles;
                particleE.ParticleSpawnRate = particleSpawnRate;
                particleE.BGX = BGX;
                particleE.BGY = BGY;
                if(emitXArray(1) < 1)
                    emitXArray(1) = 1;
                end
                if(emitXArray(2) > particleE.BGX)
                    emitXArray(2) = particleE.BGX;
                end
                if(emitYArray(1) < 1)
                    emitYArray(1) = 1;
                end
                if(emitYArray(2) > particleE.BGY)
                    emitYArray(2) = particleE.BGY;
                end
                particleE.EmitXArray = emitXArray;
                particleE.EmitYArray = emitYArray;
                particleE.VX = VX;
                particleE.VY = VY;
                particleE.LocalXVelocity = localXVelocity;
                particleE.LocalYVelocity = localYVelocity;
            end
        end
        
        function initializeParticles(obj, minmaxDia, minmaxLumi0, life)
            particleArray(1:obj.AmountOfParticles) = Particle();
            obj.ParticleArray = particleArray;
            tmp = size(minmaxDia);
            if(tmp(2) == 2)
                dia = num2cell(round((minmaxDia(2)-minmaxDia(1)).*rand(1,obj.AmountOfParticles) + minmaxDia(1)));
            else
                dia = num2cell(round(zeros(1,numel(particleArray))+minmaxDia));
            end
            tmp = size(minmaxLumi0);
            if(tmp(2) == 2)
                lumi = num2cell((minmaxLumi0(2)-minmaxLumi0(1)).*rand(1,obj.AmountOfParticles) + minmaxLumi0(1));
            else
                lumi = num2cell(zeros(1,numel(particleArray))+minmaxLumi0);
            end
            if exist('life', 'var')
                tmp = size(life);
                if(tmp(2) == 2)
                    l = num2cell(round((life(2)-life(1)).*rand(1,obj.AmountOfParticles) + life(1)));
                else
                    l = num2cell(round(zeros(1,numel(particleArray))+life));
                end
                [particleArray.Life] = l{:};
            end
            id = num2cell(randperm(obj.AmountOfParticles));
            [particleArray.Diameter] = dia{:};
            [particleArray.Lumi] = lumi{:};
            [particleArray.ID] = id{:};

            obj.ParticleArray = particleArray;
            
            obj.ParticleCount = obj.AmountOfParticles;
            
            obj.PadWindow = max([obj.ParticleArray.Diameter]);
            obj.Window = zeros(obj.BGX + obj.PadWindow*2, obj.BGY + obj.PadWindow*2);
        end
        function [window, particleArray] = update(obj)
            if(size(obj.ParticleArray,2) > 1)
                if(sum([obj.ParticleArray.Spawned]) < obj.AmountOfParticles)
                    if(obj.ParticleSpawnNow >= 1)
                        for i=1:obj.ParticleCount
                            if(obj.ParticleArray(i).Spawned == 0)
                                pl = obj.ParticleArray(i).luminosity();
                                pl = pl/max(max(pl))*255;
                                plRowSize = size(pl,1);
                                plColSize = size(pl,2);
                                obj.ParticleArray(i).XPos = round((obj.EmitXArray(2)-obj.EmitXArray(1)).*rand(1) + obj.EmitXArray(1));
                                obj.ParticleArray(i).YPos = round((obj.EmitYArray(2)-obj.EmitYArray(1)).*rand(1) + obj.EmitYArray(1));
                                obj.ParticleArray(i).XTrack{end+1} = obj.ParticleArray(i).XPos;
                                obj.ParticleArray(i).YTrack{end+1} = obj.ParticleArray(i).YPos;
                                obj.Window(obj.ParticleArray(i).XPos:obj.ParticleArray(i).XPos + plRowSize -1, obj.ParticleArray(i).YPos:obj.ParticleArray(i).YPos + plColSize-1) = obj.Window(obj.ParticleArray(i).XPos:obj.ParticleArray(i).XPos + plRowSize -1, obj.ParticleArray(i).YPos:obj.ParticleArray(i).YPos + plColSize-1) + pl;

                                obj.ParticleArray(i).Spawned = 1;
                                obj.ParticleArray(i).RecentlySpawned = 1;
                                obj.ParticleArray(i).InWindow = 1;

                                obj.ParticleSpawnNow = obj.ParticleSpawnRate;
                                obj.ParticleSpawnCount = obj.ParticleSpawnCount + 1;

                                if(obj.ParticleSpawnCount >= obj.ParticleSpawnNow)
                                    obj.ParticleSpawnCount = 0;
                                    %Break out of loop
                                    break
                                end
                            end
                        end
                    else
                        obj.ParticleSpawnNow = obj.ParticleSpawnNow + obj.ParticleSpawnRate;
                    end
                end
                i = 1;
                while(i <= obj.ParticleCount)
                    if(obj.ParticleArray(i).RecentlySpawned == 1)
                        obj.ParticleArray(i).RecentlySpawned = 0;
                        i = i + 1;
                    elseif(obj.ParticleArray(i).Spawned == 1 && obj.ParticleArray(i).RecentlySpawned == 0)
                        pl = obj.ParticleArray(i).luminosity();
                        pl = pl/max(max(pl))*255;
                        plRowSize = size(pl,1);
                        plColSize = size(pl,2);
                        
                        obj.ParticleArray(i).XPos = obj.ParticleArray(i).XPos + round(obj.ParticleArray(i).XPos*(0.01*rand(1)-0.005)*obj.LocalXVelocity) + round(obj.VX(obj.ParticleArray(i).XPos,obj.ParticleArray(i).YPos));
                        if((obj.ParticleArray(i).XPos < 1) || (obj.ParticleArray(i).XPos > obj.BGX))
                            obj.ParticleArray(i).YTrack{end+1} = obj.ParticleArray(i).YPos;
                            obj.ParticleArray(i).InWindow = 0;
                            %Saving the last coordinate of the out of bound
                            %particle.
                            tmpXPos = obj.ParticleArray(i).XPos;
                            if(tmpXPos < 1)
                                tmpXPos = 1;
                                obj.ParticleArray(i).XTrack{end+1} = tmpXPos;
                            elseif(tmpXPos > obj.BGX)
                                tmpXPos = obj.BGX;
                                obj.ParticleArray(i).XTrack{end+1} = tmpXPos;
                            else
                                obj.ParticleArray(i).XTrack{end+1} = tmpXPos;
                            end
                        else
                            obj.ParticleArray(i).YPos = obj.ParticleArray(i).YPos + round(obj.ParticleArray(i).YPos*(0.01*rand(1)-0.005)*obj.LocalYVelocity) + round(obj.VY(obj.ParticleArray(i).XPos,obj.ParticleArray(i).YPos));
                            if((obj.ParticleArray(i).YPos < 1) || (obj.ParticleArray(i).YPos > obj.BGY))
                                obj.ParticleArray(i).XTrack{end+1} = obj.ParticleArray(i).XPos;
                                obj.ParticleArray(i).InWindow = 0;
                                %Saving the last coordinate of the out of bound
                                %particle.
                                tmpYPos = obj.ParticleArray(i).YPos;
                                if(tmpYPos < 1)
                                    tmpYPos = 1;
                                    obj.ParticleArray(i).YTrack{end+1} = tmpYPos;
                                elseif(tmpYPos > obj.BGY)
                                    tmpYPos = obj.BGY;
                                    obj.ParticleArray(i).YTrack{end+1} = tmpYPos;
                                else
                                    obj.ParticleArray(i).YTrack{end+1} = tmpYPos;
                                end
                            else
                                obj.Window(obj.ParticleArray(i).XPos:obj.ParticleArray(i).XPos + plRowSize -1, obj.ParticleArray(i).YPos:obj.ParticleArray(i).YPos + plColSize-1) = obj.Window(obj.ParticleArray(i).XPos:obj.ParticleArray(i).XPos + plRowSize -1, obj.ParticleArray(i).YPos:obj.ParticleArray(i).YPos + plColSize-1) + pl;

                                obj.ParticleArray(i).XTrack{end+1} = obj.ParticleArray(i).XPos;
                                obj.ParticleArray(i).YTrack{end+1} = obj.ParticleArray(i).YPos;

                                obj.ParticleArray(i).Life = obj.ParticleArray(i).Life - 1;
                            end
                        end

                        if(obj.ParticleArray(i).Life <= 0 || obj.ParticleArray(i).InWindow ~= 1)
                            deadParticle = obj.ParticleArray(i);
                            lastAliveParticle = obj.ParticleArray(obj.ParticleCount);
                            obj.ParticleArray(obj.ParticleCount) = deadParticle;
                            obj.ParticleArray(i) = lastAliveParticle;
                            obj.ParticleCount = obj.ParticleCount - 1;
                        else
                            i = i + 1;
                        end
                    elseif(obj.ParticleArray(i).Spawned == 0)
                        i = i + 1;
                    end
                end
                particleArray = obj.ParticleArray;
                window = obj.Window;
                
                %Resetting the window
                obj.Window = zeros(obj.BGX + obj.PadWindow*2, obj.BGY + obj.PadWindow*2);
                %Removing the padding
                window = window(round(obj.PadWindow/2) + 1:end - round(obj.PadWindow/2), round(obj.PadWindow/2) + 1:end - round(obj.PadWindow/2));
            else
                disp('Need to initialize particles first.')
            end
        end
        
        function set.BGX(obj, BGX)
            obj.BGX = BGX;
        end
        function set.BGY(obj, BGY)
            obj.BGY = BGY;
        end
        function bgx = get.BGX(obj)
            bgx = obj.BGX;
        end
        function bgy = get.BGY(obj)
            bgy = obj.BGY;
        end
        function set.VX(obj, VX)
            obj.VX = VX;
        end
        function set.VY(obj, VY)
            obj.VY = VY;
        end
        function vx = get.VX(obj)
            vx = obj.VX;
        end
        function vy = get.VY(obj)
            vy = obj.VY;
        end
        function set.EmitXArray(obj, EmitXArray)
            obj.EmitXArray = EmitXArray;
        end
        function set.EmitYArray(obj, EmitYArray)
            obj.EmitYArray = EmitYArray;
        end
        function emitXArray = get.EmitXArray(obj)
            emitXArray = obj.EmitXArray;
        end
        function emitYArray = get.EmitYArray(obj)
            emitYArray = obj.EmitYArray;
        end
        function set.ParticleArray(obj, particleArray)
            obj.ParticleArray = particleArray;
        end
        function particleArray = get.ParticleArray(obj)
            particleArray = obj.ParticleArray;
        end
    end
end
