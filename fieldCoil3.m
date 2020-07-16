function [Bx,By,Bz]=fieldCoil3(x,y,z,C)
%FIELDCOIL Summary of this function goes here
%   Detailed explanation goes here

if nargin~=4
    x=0;    % x position to evaluaate at
    y=0;    % y position to evaluate at
    z=0;    % z position to evaluate at
    
    % Coils are assumed to be symmetric along the x direction.
    C(1)=0; % coil x position
    C(2)=0; % coil z position
    C(3)=1; % coil radius
end
% Redefine the coil vector array so it's easier to read

    x0=C(1);
    z0=C(2);
    R=C(3);

    % permeability of free space
    mu0=4*pi*1E-7;

    % center magnetic field
    B0=mu0./(2*R);

    % Define relative coordinates
    dZ=z-z0;        % differential z distance from center
    dX=x-x0;        % differential x distance from center
    dY=y;           % differential y distance from center

    r=sqrt(dX.^2+dY.^2);    % Radial distance 

    % Define normalized coordinates
    alpha=r./R;
    beta=dZ./R;
    gamma=dZ./r;

    Q=(1+alpha).^2+beta.^2;
    k2=4*alpha./Q;

    % Elliptic integrals
    [K,E]=ellipke(k2);

    % Output the field
    Bz=B0.*(1./(pi.*sqrt(Q))).*(E.*(1-alpha.^2-beta.^2)./(Q-4.*alpha)+K);

    Br=B0.*(gamma./(pi.*sqrt(Q))).*(E.*(1+alpha.^2+beta.^2)./(Q-4.*alpha)-K);

    Bx=Br.*(dX./r);
    By=Br.*(dY./r);

    if r==0
        Bx=0;
        By=0;
    end




end

