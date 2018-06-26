function [particle] = particleLuminosity(x, y, luminosity, particleDiameter, x0, y0)
    particle = luminosity*exp((-(x-x0).^2-(y-y0).^2)/((1/8)*particleDiameter.^2));
end
