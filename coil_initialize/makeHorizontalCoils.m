function coils_out = makeHorizontalCoils
% Coils from from Yee thesis
%
% This functin creates a structure array where each element is a descriptor
% for each coil.  These coils are called
%
% PUSH, MOT, Coil 3, Coil 4, Coil 5, Coil 6, Coil 7, Coil 8, Coil 9,
% Coil 10, and Coil 11
%
% Coil 12a and Coil 12b are named as such because during horizontal
% transport they act as a single pair of coils for the final horizontal
% transport coil pair.
%
% Coil 11 Extra is an additional coil that was installed by DM/DJ before I
% joined the group.  The properties of the coil were determined by a
% mathematica notebook in Dave's files.

%% Coil 1 ("Push")

cpush = makePushCoil;

%% Coil 2 ("MOT")
c2 = struct;
c2.Name = 'MOT';
c2.Name2 = 'MOT';
c2.ID = 45.7*1e-3;
c2.OD = 96*1e-3;
c2.Nrad = 33;
c2.Nax = 2;

c2.WireDim = [0.73 2.63]*1e-3;
c2.InnerSeparation = 50*1e-3;
c2.Position = 0*1e-3;

%% Coil 3 ("Horiz 1")
c3 = struct;
c3.Name = 'Coil 3';
c3.Name2 = 'Horiz 1';
c3.ID = 24.8*1e-3;
c3.OD = 61*1e-3;
c3.Nrad = 16;
c3.Nax = 2;
% c3.Nrad = 1;
% c3.Nax = 1;
c3.WireDim = [1.13 2.63]*1e-3;
c3.InnerSeparation = 63.48*1e-3;
c3.Position = 53.3*1e-3;

%% Coil 4 ("Horiz 2")
c4 = struct;
c4.Name = 'Coil 4';
c4.Name2 = 'Horiz 2';
c4.ID = 24.8*1e-3;
c4.OD = 61*1e-3;
c4.Nrad = 16;
c4.Nax = 2;
c4.WireDim = [1.13 2.63]*1e-3;
c4.InnerSeparation = 51.48*1e-3;
c4.Position = 84.8*1e-3;

%% Coil 5 ("Horiz 3")
c5 = struct;
c5.Name = 'Coil 3';
c5.Name2 = 'Horiz 3';
c5.ID = 24.8*1e-3;
c5.OD = 61*1e-3;
c5.Nrad = 16;
c5.Nax = 2;
c5.WireDim = [1.13 2.63]*1e-3;
c5.InnerSeparation = 63.48*1e-3;
c5.Position = 116.3*1e-3;

%% Coil 6 ("Horiz 4")
c6 = struct;
c6.Name = 'Coil 6';
c6.Name2 = 'Horiz 4';
c6.ID = 24.8*1e-3;
c6.OD = 61*1e-3;
c6.Nrad = 16;
c6.Nax = 2;
c6.WireDim = [1.13 2.63]*1e-3;
c6.InnerSeparation = 51.48*1e-3;
c6.Position = 147.8*1e-3;

%% Coil 7 ("Horiz 5")
c7 = struct;
c7.Name = 'Coil 7';
c7.Name2 = 'Horiz 5';
c7.ID = 24.8*1e-3;
c7.OD = 61*1e-3;
c7.Nrad = 16;
c7.Nax = 2;
c7.WireDim = [1.13 2.63]*1e-3;
c7.InnerSeparation = 63.48*1e-3;
c7.Position = 179.3*1e-3;

%% Coil 8 ("Horiz 6")
c8 = struct;
c8.Name = 'Coil 8';
c8.Name2 = 'Horiz 6';
c8.ID = 24.8*1e-3;
c8.OD = 61*1e-3;
c8.Nrad = 16;
c8.Nax = 2;
c8.WireDim = [1.13 2.63]*1e-3;
c8.InnerSeparation = 51.48*1e-3;
c8.Position = 210.8*1e-3;

%% Coil 9 ("Horiz 7")
c9 = struct;
c9.Name = 'Coil 9';
c9.Name2 = 'Horiz 7';
c9.ID = 24.8*1e-3;
c9.OD = 61*1e-3;
c9.Nrad = 16;
c9.Nax = 2;
c9.WireDim = [1.13 2.63]*1e-3;
c9.InnerSeparation = 63.48*1e-3;
c9.Position = 242.3*1e-3;

%% Coil 10 ("Horiz 8")
c10 = struct;
c10.Name = 'Coil 10';
c10.Name2 = 'Horiz 8';
c10.ID = 24.8*1e-3;
c10.OD = 61*1e-3;
c10.Nrad = 16;
c10.Nax = 2;
c10.WireDim = [1.13 2.63]*1e-3;
c10.InnerSeparation = 51.48*1e-3;
c10.Position = 273.8*1e-3;

%% Coil 11 ("Horiz 9")
c11 = struct;
c11.Name = 'Coil 11';
c11.Name2 = 'Horiz 9';
c11.ID = 24.8*1e-3;
c11.OD = 61*1e-3;
c11.Nrad = 16;
c11.Nax = 2;
c11.WireDim = [1.13 2.63]*1e-3;
c11.InnerSeparation = 63.48*1e-3;
c11.Position = 305.3*1e-3;


%% Coil 11 Extra
% This is an extra coil pair which I think D. McKay wound during his PhD, there
% is little information. The idea of this coil is that the distane
% from Coil 10 to the Center of 12ab is too large because of the fact that
% the coils cannot intersect with the vacuum chamber.  Michael's intial
% calculation suggested this was okay,
%
% "4.4.6 Horizontal to vertical transport stage" Pp.22

c11_extra = struct;
c11_extra.Name = 'Coil 11 Extra';
c11_extra.Name2 = 'NEW COIL';

% Estimated as of 2024/08/13
c11_extra.ID = 2*7.6*1e-3;
c11_extra.OD = 61*1e-3;
c11_extra.Nrad = 22;
c11_extra.Nax = 2;
c11_extra.WireDim = [1.02 2.29]*1e-3;
c11_extra.InnerSeparation = 87*1e-3;
c11_extra.Position = 318*1e-3; 

%% Coil 12

v_coils = makeVerticalCoils;
c12a_vert = v_coils(1);
c12b_vert = v_coils(2);

c12 = struct;
c12.Name = 'Coil 12';
c12.Name2 = 'Vert 1/2';

% Estimated as of 2024/08/13
c12.ID = c12a_vert.ID;
c12.OD = c12a_vert.OD;
c12.Nrad = c12a_vert.Nrad;
c12.Nax = c12a_vert.Nax;
c12.WireDim = c12a_vert.WireDim;
c12.InnerSeparation = c12b_vert.zbot - (c12a_vert.zbot+c12a_vert.Height);
c12.Position = 360*1e-3; 

% Calculate coil 11 position
D = 19.2e-3; % diameter of vertical tube
dx = 2e-3;
c11_extra.Position = c12.Position(1) - (D/2+dx+c11_extra.OD/2);


%% Assemble Coils

coils = [c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c11_extra c12];

for kk=1:length(coils)
    coils(kk).Coil = makeCoilPair(coils(kk));
end

%% Combine with Push
coils_out = [cpush coils];

end


