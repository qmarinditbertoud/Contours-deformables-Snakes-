clear all;
close all;
I=double(imread('plane.png'));
%I2=double(I);
imagesc(I);
colormap gray;

n=50;
R=50;
Cx=400/2; %527 pixels de largeur
Cy=300/2; %451 pixels de hauteur
t=[0:n-1]/n;
t=t';
vx=R*cos(2*pi*t)+Cx;
vy=R*sin(2*pi*t)+Cy;
hold on;
plot([vx;vx(1)],[vy;vy(1)], 'r');

D2=zeros(n);

D2=diag(-2*ones(1,n)) + diag(ones(1,n-1),1) + diag(ones(1,n-1),-1) + diag(1,n-1) + diag(1,-n+1);
D4=diag(6*ones(1,n))  + diag(-4*ones(1,n-1),1) + diag(-4*ones(1,n-1),-1) + diag(ones(1,n-2),2)+ diag(ones(1,n-2),-2) + diag(ones(1,2),n-2) + diag(ones(1,2),-n+2) + (-4*diag(1,n-1)) + (-4*diag(1,-n+1));

D2x=D2*vx;
D2y=D2*vy;
D4x=D4*vx;
D4y=D4*vy;

figure(2);
hold on;
plot(t,D2x);
plot(t,D2y);
legend('courbe D2x en fonction de t','courbe D2y en fonction de t');

figure(3);
hold on;
plot(t,D4x);
plot(t,D4y);
legend('courbe D4x en fonction de t','courbe D4y en fonction de t');

deltat=0.05;
lambda1=7;%influe sur la longueur
lambda2=150;%Courbure
lambda3=0.1;%pour converger sur les contours


A=diag(ones(1,n))+2*deltat*(-lambda1*D2+lambda2*D4);
IA=inv(A);



D=diag(-1*ones(1,n-1),-1) + diag(ones(1,n-1),1)+ diag(1,n-1) + diag(1,-n+1); %matrice dérivée
Einterne=[];

H=[-1,0,1];
G=[-1,0,1]';
Gh=imfilter(I,H);
Gg=imfilter(I,G);
Gradient=(Gg.^2)+(Gh.^2);

%clf;
figure(4);
subplot(1,3,1);
imagesc(Gradient);colormap gray;
title('norme du gradient')

Kx=-imfilter(Gradient,H);
Ky=-imfilter(Gradient,G);

subplot(1,3,2);
imagesc(Kx);
title('Kx');
subplot(1,3,3);
imagesc(Ky);
title('Ky');

figure(5)



for k=[1:50]
 


 %interpolation linéaire calculant P=K(vx,vy)
    Px=interp2(Kx,vx,vy);
    Py=interp2(Ky,vx,vy);
   
    vx=IA*(vx+(deltat*lambda3*Px));
    vy=IA*(vy+(deltat*lambda3*Py));
  
    
    
    
    Dx=D*vx;
    Dy=D*vy;
    D2x=D2*vx;
    D2y=D2*vy;
    Einterne=[Einterne, sum(lambda1*((Dx).^2+(Dy).^2)+lambda2*((D2x).^2+(D2y).^2))]; % expression de l'energie interne

    
    
    if (mod(k,2)==0)
        
        clf;
        imagesc(I);
        colormap gray;
        hold on;       
   
        plot([vx;vx(1)],[vy;vy(1)]);
        pause(0.01);
        
    end
end

figure(6);
plot(Einterne);
title('energie interne en fonction de lindice');