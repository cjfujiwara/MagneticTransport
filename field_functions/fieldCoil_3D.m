function [Bx,By,Bz]=fieldCoil_3D(x,y,z,C)
% The functions generates the magnetic field at position (x,y,z) with a
% coil descriptor of C.
%
% C = [X_c,Z_c,R,A]
%
% Without loss of generality, the center point of all coils lie in the x-z
% plane with axis of clyndrical symmetry pointed in the Z direction.
%
% It is assumed that all units are in SI, so that the positins are in
% meters. The output magnetic field is in Tesla.
%
% Because it is often useful to specify multiple current loops at one time,
% (such as a wound inductor), the descriptor C may be a matrix.  Here, each
% row of the matrix defines a different current loop.
%
% The magnetic fields are calculated in general terms which allow for off
% axis values.  See
% 
% https://tiggerntatie.github.io/emagnet/offaxis/iloopoffaxis.htm
% 
% For more details about how these functions are generated.

if nargin~=4
    x=0;    % x position to evaluaate at
    y=0;    % y position to evaluate at
    z=0;    % z position to evaluate at
    
    % Coils are assumed to be symmetric along the x direction.
    C(1)=0; % coil x position
    C(2)=0; % coil z position
    C(3)=1; % coil radius
    C(4)=1; % coil current (amps)
end

Bx_total = {};
By_total = {};
Bz_total = {};

for nn = 1:size(C,1)
    % Redefine the coil vector array so it's easier to read
    x0 = C(nn,1);
    z0 = C(nn,2);
    R  = C(nn,3);
    A  = C(nn,4);

    % permeability of free space (in SI)
    mu0=4*pi*1E-7;

    % center magnetic field (normalized to 1A)
    B0=A*mu0./(2*R);

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

    Bx(r==0) = 0;
    By(r==0) = 0;
 
    Bx_total{nn} = Bx;
    By_total{nn} = By;
    Bz_total{nn} = Bz;
    
end

Bx = Bx_total{1};
By = By_total{1};
Bz = Bz_total{1};
for ii=2:length(Bx_total)
    Bx = Bx + Bx_total{ii};
    By = By + By_total{ii};
    Bz = Bz + Bz_total{ii};
end

end

