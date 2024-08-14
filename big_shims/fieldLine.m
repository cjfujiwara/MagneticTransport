function [Bx,By,Bz] = fieldLine(x,y,z,L,R,A)

% permeability of free space
mu0=4*pi*1E-7;

Bx_total = {};
By_total = {};
Bz_total = {};

for nn = 1:size(L,1)
    Lthis = L(nn,:);
    Rthis = R(nn,:);
    
    
    Lhat = Lthis/norm(Lthis);
    
    % Redefine the coil vector array so it's easier to read
    x0 = Rthis(1);
    y0 = Rthis(2);
    z0 = Rthis(3);

    % Vector from wire base to the point of inquiry
    x_prime = x - x0;
    y_prime = y - y0;
    z_prime = z - z0;  
        
    % Project r_prime onto the Lhat direction
    P = x_prime*Lhat(1) + y_prime*Lhat(2) + z_prime*Lhat(3); 
    
    % Calculate the cylindrical distance from the wire endpoints to the
    % point of inquiry
    zL1 = -P;
    zL2 = zL1+norm(Lthis);
    
    % Calculate the nearest vector distance.
    rho0_x = x_prime - P*Lhat(1);
    rho0_y = y_prime - P*Lhat(2);
    rho0_z = z_prime - P*Lhat(3);   
    
    % Norm of closest approach
    rho0 = sqrt(rho0_x.^2+rho0_y.^2+rho0_z.^2);
    
    % Calculate direction of magnetic field
    phihat_x = (Lhat(2)*rho0_z - Lhat(3)*rho0_y)./rho0;
    phihat_y = (Lhat(3)*rho0_x - Lhat(1)*rho0_z)./rho0;
    phihat_z = (Lhat(1)*rho0_y - Lhat(2)*rho0_x)./rho0;
    
    % Magnetic field prefactor
    B0 = mu0*A./(4*pi*rho0);
    B0 = B0*1e4; % Convert T to Gauss
    
    % Calcluate sin(theta) which are the tips of the wire
    f1 = zL1./sqrt(rho0.^2+zL1.^2);
    f2 = zL2./sqrt(rho0.^2+zL2.^2);

    % Value of magnetic field in phihat direction
    B = B0.*abs(f1-f2);   

    % Convert magnetic field into cartersian
    Bx = B.*phihat_x;
    By = B.*phihat_y;
    Bz = B.*phihat_z;
    
     
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

