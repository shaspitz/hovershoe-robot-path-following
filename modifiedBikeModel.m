function zp = modifiedBikeModel(z, u, TS)

x = z(1);
y = z(2); 
psi = z(3);

v = u(1);
psidot = u(2);

xdot = v*cos(psi);
ydot = v*sin(psi);

xp = x + TS*xdot;
yp = y + TS*ydot;
psip = psi + TS*psidot;

zp = [xp; yp; psip];

end

