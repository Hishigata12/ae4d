%% Flattop
t = 'Flattop'
N=101;
[X,Y] = meshgrid(-1:2/(N-1):1);
Z = sqrt(X.^2 + Y.^2)<=1;
Z = Z*1;

Nmid = ceil((N-1)/2+1);
X = X(1,:);
Y = Y(:,1);
Zx = Z(Nmid,Nmid:end);
Zy = Z(Nmid:end,Nmid)';
M=floor((N-1)/2+1);
[Zr,r] = radialavg(Z,M);
err = Zr-(Zx+Zy)/2;

figure;
subplot(2,1,1)
	imagesc(X,Y,Z);colormap(jet);colorbar;title(t);
	xlabel('X');ylabel('Y');
subplot(2,1,2)
	plot(r,Zr,'o', X(Nmid:end),Zx,'.-', Y(Nmid:end),Zy,'x');
	title(t)
	xlabel('r'); ylabel('Z');
	legend('Zr', 'Zx', 'Zy');
	axis([0 1 0 1])

figure;plot(r,err);
title(sprintf('%s Error=%0.1g', t, std(err)))
xlabel('rho'); ylabel('Z');
%% Cone
t = 'Cone'
N=201;
[X,Y] = meshgrid(-1:2/(N-1):1);
Nmid = (N-1)/2+1;

Z = 1-sqrt(X.^2 + Y.^2);

X = X(1,:);
Y = Y(:,1);
Zx = Z(Nmid,Nmid:end);
Zy = Z(Nmid:end,Nmid)';
M=(N-1)/2+1;
[Zr,r] = radialavg(Z,M);
err = Zr-(Zx+Zy)/2;

figure;
subplot(2,1,1)
	imagesc(X,Y,Z);colormap(jet);colorbar;title(t);
	xlabel('X');ylabel('Y');
subplot(2,1,2)
	plot(r,Zr,'o', X(Nmid:end),Zx,'.-', Y(Nmid:end),Zy,'x');
	title(t)
	xlabel('r'); ylabel('Z');
	legend('Zr', 'Zx', 'Zy');
	axis([0 1 0 1])

figure;plot(r,err);
title(sprintf('%s Error=%0.1g', t, std(err)))
xlabel('rho'); ylabel('Z');
%% Sine
t = 'Radially Sinusoidal'
N=201;
[X,Y] = meshgrid(-1:2/(N-1):1);
[Q,R] = cart2pol(X,Y);
Z = sin(2*pi*R);

Nmid = (N-1)/2+1;
X = X(1,:);
Y = Y(:,1);
Zx = Z(Nmid,Nmid:end);
Zy = Z(Nmid:end,Nmid)';
M=(N-1)/2+1;
[Zr,r] = radialavg(Z,M);
err = Zr-(Zx+Zy)/2;

figure;
subplot(2,1,1)
	imagesc(X,Y,Z);colormap(jet);colorbar;title(t);
	xlabel('X');ylabel('Y');
subplot(2,1,2)
	plot(r,Zr,'o', X(Nmid:end),Zx,'.-', Y(Nmid:end),Zy,'x');
	title(t)
	xlabel('r'); ylabel('Z');
	legend('Zr', 'Zx', 'Zy');
	axis([0 1 -1 1])

figure;plot(r,err);
title(sprintf('%s Error=%0.1g', t, std(err)))
xlabel('rho'); ylabel('Z');
%% Sine
t = 'Azimuthally Sinusoidal'
N=201;
[X,Y] = meshgrid(-1:2/(N-1):1);
[Q,R] = cart2pol(X,Y);
Z = sin(Q);

Nmid = (N-1)/2+1;
X = X(1,:);
Y = Y(:,1);
Zx = Z(Nmid,:);
Zy = Z(:,Nmid)';
M=(N-1)/2+1;
[Zr,r] = radialavg(Z,M);
err = Zr - 0;

figure;
subplot(2,1,1)
	imagesc(X,Y,Z);colormap(jet);colorbar;title(t);
	xlabel('X');ylabel('Y');
subplot(2,1,2)
	plot(r,Zr,'o', X,Zx,'.-', Y,Zy,'x');
	title(t)
	xlabel('r'); ylabel('Z');
	legend('Zr', 'Zx', 'Zy');
	axis([-1 1 -1 1])

figure;plot(r,err);
title(sprintf('%s Error=%0.1g', t, std(err)))
xlabel('rho'); ylabel('Z');
%% Flattop with NaNs
t = 'Flattop with NaNs'
N=201;
[X,Y] = meshgrid(-1:2/(N-1):1);
Nmid = ceil((N-1)/2+1);

Z = sqrt(X.^2 + Y.^2)<=1;
Z = Z*1;
Z(Nmid-10:Nmid+10,:) = NaN;
Z(sqrt(X.^2 + Y.^2)>=0.7 & sqrt(X.^2 + Y.^2)<=0.8) = NaN;

X = X(1,:);
Y = Y(:,1);
Zx = Z(Nmid,Nmid:end);
Zy = Z(Nmid:end,Nmid)';
M=floor((N-1)/2+1);
[Zr,r] = radialavg(Z,M);
% err = Zr-(Zx+Zy)/2;
err = Zr - 1;

figure;
subplot(2,1,1)
	imagesc(X,Y,Z);colormap(jet);colorbar;title(t);
	xlabel('X');ylabel('Y');
subplot(2,1,2)
	plot(r,Zr,'o', X(Nmid:end),Zx,'.-', Y(Nmid:end),Zy,'x');
	title(t)
	xlabel('r'); ylabel('Z');
	legend('Zr', 'Zx', 'Zy');
	axis([0 1 0 1])

figure;plot(r,err);
title(sprintf('%s Error=%0.1g', t, std(err)))
xlabel('rho'); ylabel('Z');
%% Cone with Offsets
t = 'Cone with Offsets'
N=201;
[X,Y] = meshgrid(-1:2/(N-1):1);
xo = +0.25;
yo = -0.25;

X = X-xo;
Y = Y-yo;
Z = 1-sqrt(X.^2 + Y.^2);

X = X(1,:);
Y = Y(:,1);
Nmid = ceil((N-1)/2+1);
Nmidx = find(X==0);
Nmidy = find(Y==0);
Zx = Z(Nmidy,Nmidx:end);
Zy = Z(Nmidy:end,Nmidx)';
M=(N-1)/2+1;
[Zr,r] = radialavg(Z,M,xo,yo);

figure;
subplot(2,1,1)
	imagesc(X,Y,Z);colormap(jet);colorbar;title(t);
	xlabel('X');ylabel('Y');
subplot(2,1,2)
	plot(r,Zr,'o', X(Nmidx:end),Zx,'.-', Y(Nmidy:end),Zy,'x');
	title(t)
	xlabel('r'); ylabel('Z');
	legend('Zr', 'Zx', 'Zy');
	axis([0 1 0 1])
%% Radially Sinusoidal with Offsets
t = 'Radially Sinusoidal with Offsets'
N=201;
[X,Y] = meshgrid(-1:2/(N-1):1);
xo = +0.25;
yo = -0.25;
X = X-xo;
Y = Y-yo;
[Q,R] = cart2pol(X,Y);
Z = sin(2*pi*R);

X = X(1,:);
Y = Y(:,1);
Nmid = ceil((N-1)/2+1);
Nmidx = find(X==0);
Nmidy = find(Y==0);
Zx = Z(Nmidy,Nmidx:end);
Zy = Z(Nmidy:end,Nmidx)';
M=(N-1)/2+1;
[Zr,r] = radialavg(Z,M,xo,yo);

figure;
subplot(2,1,1)
	imagesc(X,Y,Z);colormap(jet);colorbar;title(t);
	xlabel('X');ylabel('Y');
subplot(2,1,2)
	plot(r,Zr,'o', X(Nmidx:end),Zx,'.-', Y(Nmidy:end),Zy,'x');
	title(t)
	xlabel('r'); ylabel('Z');
	legend('Zr', 'Zx', 'Zy');
	axis([0 1 -1 1])