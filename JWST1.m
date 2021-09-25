%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Shun Qin, PhD, University of G?ttingen
% Email: shun.qin@outlook.com
% reference: 
% [1] Shun Qin,et al.,A Tip–Tilt and Piston Detection Approach for
%     Segmented Telescopes,Photonics, 8(1), 3, 2021.
%
% sample and comments: 
% N=512/2; % the sampling pixels
% Dia=1;   % the diameter of the circumcircle of a hexaggon;
% nring=1; % numbers of ring of hexaggons,without counting the core hexaggon
% gap=0.01;% the width of the gap between segments
% times=1; % the factor of zeros-padding outside the area of the wavefront

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W=JWST1(Dia,gap,nring,times,N)

[x,y]=meshgrid(-times*Dia/2:Dia/(N-1):times*Dia/2,times*Dia/2:-Dia/(N-1):-times*Dia/2);%
D=Dia/(nring+sqrt(3)/3)/2;%两子镜之间中心的距离
d=D/sqrt(3);%两蜂窝元间的距离比例系数
co=(D/2-gap/2)/(sqrt(3)/2);
%z=singleface(30,30,x,y);
%surf(x,y,z); 
c=hexaggon(nring);%取蜂窝平面的坐标
%surf(x,y,z);
num=length(c);%制动器的个数,不包括中心那个
% central=polygon(x,y,co,[0 0]);
% for i=1:L;
% central=polygon(x,y,co,c{i}*d) | central;
% end
lamda=632.8e-9;%单位m

% A, B, C store the slopes with respective to 3 axises of all segments
A0=0.1*lamda/D;A=normrnd(0,A0/2,num+1,1);A=tan(A);%A=0;
B0=0.1*lamda/D;B=normrnd(0,B0/2,num+1,1);B=tan(B);%B=0;
C0=100e-9;C=normrnd(0,C0/2,num+1,1);%C=0;

z=polygon(x,y,A(1),B(1),C(1),co,[0 0]);
for i=1:num;
    
    z=polygon(x,y,A(i+1),B(i+1),C(i+1),co,c{i}*d)+z;
    
end
W=z;
PV0=max(max(W))-min(min(W));
nlamda=0.7;
W=W*2*pi*nlamda/PV0;
W=(W-sum(sum(W))/sum(sum(W(:,:)~=0)));
% imagesc(W)
% figure;
% subplot(2,2,1);surf(pupil);shading flat
% Z=ft2(pupil,Dia/(N-1));
% subplot(2,2,3);imagesc(log10(abs(Z)^2));
% maxv0=max(max(abs(Z)^2))
% 
% 
% subplot(2,2,2);surf(z);shading flat
% Z=ft2(pupil.*exp(1i*z*2*pi/lamda),Dia/(N-1));
% subplot(2,2,4);imagesc(log10(abs(Z)^2));
% maxv=max(max(abs(Z)^2))
% 
% SR=maxv/maxv0


%%%%%%%%%%%%%%%%%%%%%%%%%%%% 画一个六边形 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z=polygon(x,y,A,B,C,co,vector)%c比例系数原型是单位一的六边形，乘此系数得实际值 vector为平移向量
k=sqrt(3);
x=x-vector(1);y=y+vector(2);
z1=A*x+B*y+C;
z0=( y>=-k*co/2 & y<=k*co/2 & y<k*x+k*co & y<-k*x+k*co & y>-k*x-k*co & y>k*x-k*co )*1;
z=double(z0).*z1;
%surf(x,y,z);

function  v=hexaggon(n0)%n0为圈数v(k+6*(n-1))=
if n0==0
    v=cell(1,1);
    v{1}=[0 0];
    return;
end
all=0;next=0;
for i=1:n0
    all=i*6+all;
end
v=cell(all,1);

for n=1:n0
    next=next+(n-1)*6;
    for k=1:6*n
        
        v{k+next}=cores(n,k);
    end
end

function  A=cores(n,k)%n圈 第k个六边形中心
i=floor(abs((k+n-1)/n));
A=(n-(k-(i-1)*n-1))*core1(i)+(k-(i-1)*n-1)*core1(i+1);
function central=core1(k)
central=sqrt(3)*[cos(pi/6+(k-1)*pi/3),sin(pi/6+(k-1)*pi/3)];

