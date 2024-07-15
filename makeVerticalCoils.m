function coils = makeVerticalCoils
% Coils from from Yee thesis
%
% This functin creates a structure array where each element is a descriptor
% for each coil.  These coils are called
%
% Coil 12a, Coil 12b, Coil 13, Coil 14, Coil 15, and Coil 16.
%
% Coil 12a and Coil 12b are named as such because during horizontal
% transport they act as a single pair of coils for the final horizontal
% transport coil pair.

%
% Coil 12b, Coil 13, and Coil 14 are used with H-bridges to change their
% polarity.
%
% Coil 15 and Coil 16 are the final trapping coils which are used in the
% quadrupolar trap.  Coil 15 changes polarity but is controlled with a
% custom set of FETs in order to maintain that Coil 15 and Coil 16 are in
% series at the end of transport.  This is important from a noise
% perspective for the quadrupolar trap.
%
% Each coil is specified by the number of radial and vertical layers.  
%
% ID : The inner diameter
% OD : The outer diameter
%
% Nrad : The number of radial layers
% Nax : The number of axial layers
% WireDim : The dimensions of the square wire in meters (width x height)
% Height : The total thickness of hte coil
% zbot : THe bottom position of the coil.  zbot = 0 is defined as the
% atomic location at the end of horizontal transport.

%% Coil 12a ("Vert1/H10")

c12a = struct;
c12a.Name = 'Coil 12a';
c12a.Name2 = 'Vert 1';
c12a.ID = 50e-3;
c12a.OD = 104.3e-3;
c12a.Nrad = 19;
c12a.Nax = 4;
c12a.WireDim = [1.43 2.33]*1e-3;
c12a.Height = c12a.Nax*c12a.WireDim(2);
c12a.zbot = -20.5e-3-c12a.Height;

%% Coil 12b ("Vert2/H10")

c12b = struct;
c12b.Name = 'Coil 12a';
c12b.Name2 = 'Vert 2';
c12b.ID = 50e-3;
c12b.OD = 104.3e-3;
c12b.Nrad = 19;
c12b.Nax = 4;
c12b.WireDim = [1.43 2.33]*1e-3;
c12b.Height = c12b.Nax*c12b.WireDim(2);
c12b.zbot = 20.5e-3;

%% Coil 13 ("Vert3")

c13 = struct;
c13.Name = 'Coil 13';
c13.Name2 = 'Vert 3';
c13.ID = 50e-3;
c13.OD = 104.3e-3;
c13.Nrad = 19;
c13.Nax = 4;
c13.WireDim = [1.43 2.33]*1e-3;
c13.Height = c13.Nax*c13.WireDim(2);
c13.zbot = c12b.zbot + c12b.Height + 34e-3;

%% Coil 14 ("Vert4")

c14 = struct;
c14.Name = 'Coil 14';
c14.Name2 = 'Vert 4';
c14.ID = 50e-3;
c14.OD = 104.3e-3;
c14.Nrad = 19;
c14.Nax = 4;
c14.WireDim = [1.43 2.33]*1e-3;
c14.Height = c14.Nax*c14.WireDim(2);
c14.zbot = c13.zbot + c13.Height + 30e-3;

%% Coil 15 ("Vert5")

c15 = struct;
c15.Name = 'Coil 15';
c15.Name2 = 'Vert 5';
c15.ID = 56.2e-3;
c15.OD = 96.8e-3;
c15.Nrad = 18;
c15.Nax = 8;
c15.WireDim = [1.13 2.63]*1e-3;
c15.Height = c15.Nax*c15.WireDim(2);
c15.zbot = c14.zbot + c14.Height + 21.3e-3;

%% Coil 16 ("Vert6")

c16 = struct;
c16.Name = 'Coil 16';
c16.Name2 = 'Vert 5';
c16.ID = 56.2e-3;
c16.OD = 96.8e-3;
c16.Nrad = 18;
c16.Nax = 8;
c16.WireDim = [1.13 2.63]*1e-3;
c16.Height = c16.Nax*c16.WireDim(2);
c16.zbot = c15.zbot + c15.Height + 34.6e-3; % WRONG
c16.zbot = 191.5e-3; % Yee Thesis. OUtdated
c16.zbot = 191.5e-3 + 14e-3; % "raise top QP" in ?? year
c16.zbot = c15.zbot + .0701; % Using levitation gradient
%% Collect all Coils
coils = [c12a c12b c13 c14 c15 c16];

for kk=1:length(coils)
    coils(kk).Coil = makeVerticalCoil(coils(kk));
end

end

