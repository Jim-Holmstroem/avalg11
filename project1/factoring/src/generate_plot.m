
%only write down the diagonal

B = A + A' - diag(diag(A))
[X,Y] = meshgrid(1:12);
[Xp,Yp] = meshgrid(1:1/4:12);

Bp = interp2(X,Y,B,Xp,Yp,'linear');

Xp4 = 4*Xp;
Yp4 = 4*Yp;

figure();
    mesh(Xp4,Yp4,log(0.1+Bp));
figure();
    mesh(Xp4,Yp4,Bp);
figure();
    mesh(Xp4,Yp4,min(Xp4,Yp4));
