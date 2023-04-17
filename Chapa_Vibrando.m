%% Gabriel Prieto Paris
% P1 - PLACAS E CASCAS
close all
clear
clc

%% DADOS DE ENTRADA

% Elemento retangular de chapa
nnel=4; % Número de nós por elemento
nglpn=2;% Número de graus de liberdade por nó

% Dimensões do retângulo
compr=0.44; % [m]
alt=0.11; % [m]

% Divisões
ndx=8;ndy=2; % Elementos em x e elementos em y
nno=(ndx+1)*(ndy+1); % Número total de nós
nel=ndx*ndy; % Número de elementos

% Dimensões do problema
nds=nno*nglpn; % Número de deslocamentos do sistema
ndpel=nnel*nglpn;% Número de deslocamentos por elemento
ns=1;% Número de seções diferentes; se ns=1, todas iguais

% Dados físicos dos elementos
t=0.055; % [m]
EM=1e9; % [N/m^2]
nu=0.25;
ro=1000; % [kg/m^3]
EL=EM/(1-nu*nu);G=EM/2/(1+nu);

%% Inicialização de matrizes e vetores
gcoord=zeros(nno,2);
nodel=zeros(nel,nnel);
LN=zeros(nno,nglpn);
LM=zeros(nel,nglpn*nnel);
k=zeros(ndpel,ndpel);
m=zeros(ndpel,ndpel);
K=zeros(nds,nds);
M=zeros(nds,nds);

%% NÓS

% Geracao das coordenadas
dx=compr/ndx;dy=alt/ndy;
x=0;y=0;k=0;
for i=1:ndx+1
    for j=1:ndy+1
        k=k+1;
        gcoord(k,1)=x;gcoord(k,2)=y;
        y=y+dy;
    end
    y=0;
    x=x+dx;
end

% Geracao da conectividade dos elementos
kel=0;kaux=0;
for i=1:ndx
    for j=1:ndy
        kaux=kaux+1;
        kel=kel+1;
        nodel(kel,1)=kaux+ndy+2;nodel(kel,2)=kaux+1;
        nodel(kel,3)=kaux;nodel(kel,4)=kaux+ndy+1;
    end
    kaux=kaux+1;
end

%condições de contorno (n. de nó restrito e direções restritas ou livres)
%se fixa = -1; se livre = 0
%engaste de viga em balanco
LN(2,:)=[-1 -1];
LN(26,:)=[0 -1];

% Determinação das matrizes

% Matriz LN
ngl=0;
for i=1:nno
    for j=1:nglpn
        if  LN(i,j)==0
            ngl=ngl+1;
            LN(i,j)=ngl;
        end
    end
end
ngr=ngl;
for i=1:nno
    for j=1:nglpn
        if LN(i,j)<0
            ngr=ngr+1;
            LN(i,j)=ngr;
        end
    end
end  

% Matriz LM
for iel=1:nel
    kl=0;
    for i=1:nnel
        ino=nodel(iel,i);
        for j=1:nglpn
            kl=kl+1;
            LM(iel,kl)=LN(ino,j);
        end
    end
end


% Matriz de rigidez e massa
for iel=1:nel
    xa=gcoord(nodel(iel,1),1);
    xb=gcoord(nodel(iel,2),1);
    yb=gcoord(nodel(iel,2),2);
    yc=gcoord(nodel(iel,3),2);

    % Dimensoes do retangulo
    a=(xa-xb)/2;b=(yb-yc)/2;

    % Constantes
    c1=EL*t*b/3/a;c2=c1/2;c3=EL*t*nu/4;
    c4=G*t*a/3/b;c5=c4/2;c6=G*t/4;

    kd(1,1)=c1;kd(1,2)=c3;kd(1,3)=-c1;kd(1,4)=c3;kd(1,5)=-c2;kd(1,6)=-c3;kd(1,7)=c2;kd(1,8)=-c3;
    kd(2,2)=c1;kd(2,3)=-c3;kd(2,4)=c2;kd(2,5)=-c3;kd(2,6)=-c2;kd(2,7)=c3;kd(2,8)=-c1;
    kd(3,3)=c1;kd(3,4)=-c3;kd(3,5)=c2;kd(3,6)=c3;kd(3,7)=-c2;kd(3,8)=c3;
    kd(4,4)=c1;kd(4,5)=-c3;kd(4,6)=-c1;kd(4,7)=c3;kd(4,8)=-c2;
    kd(5,5)=c1;kd(5,6)=c3;kd(5,7)=-c1;kd(5,8)=c3;
    kd(6,6)=c1;kd(6,7)=-c3;kd(6,8)=c2;
    kd(7,7)=c1;kd(7,8)=-c3;
    kd(8,8)=c1;

    ks(1,1)=c4;ks(1,2)=c6;ks(1,3)=c5;ks(1,4)=-c6;ks(1,5)=-c5;ks(1,6)=-c6;ks(1,7)=-c4;ks(1,8)=c6;
    ks(2,2)=c4;ks(2,3)=c6;ks(2,4)=-c4;ks(2,5)=-c6;ks(2,6)=-c5;ks(2,7)=-c6;ks(2,8)=c5;
    ks(3,3)=c4;ks(3,4)=-c6;ks(3,5)=-c4;ks(3,6)=-c6;ks(3,7)=-c5;ks(3,8)=c6;
    ks(4,4)=c4;ks(4,5)=c6;ks(4,6)=c5;ks(4,7)=c6;ks(4,8)=-c5;
    ks(5,5)=c4;ks(5,6)=c6;ks(5,7)=c5;ks(5,8)=-c6;
    ks(6,6)=c4;ks(6,7)=c6;ks(6,8)=-c4;
    ks(7,7)=c4;ks(7,8)=-c6;
    ks(8,8)=c4;

    k=kd+ks;

    % Simetria
    for i=2:8
        for j=1:i-1
            k(i,j)=k(j,i);
        end
    end

    m=a*b*t*ro*eye(8);
    
    % Matriz global
    for i=1:ndpel
        ii=LM(iel,i);
        for j=1:ndpel
            jj=LM(iel,j);
            K(ii,jj)=K(ii,jj)+k(i,j);
            M(ii,jj)=M(ii,jj)+m(i,j);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Massas pontuais aplicadas nos nós
massa1 = 5500; % [kg]
M(LN(14,1),LN(14,1))=M(LN(14,1),LN(14,1))+massa1;
M(LN(14,2),LN(14,2))=M(LN(14,2),LN(14,2))+massa1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SOLUÇÃO

% Problema de autovalores e autovetores
[mod,fre]=eig(K(1:ngl,1:ngl),M(1:ngl,1:ngl));
disp('Primeira Frequencia em rad/s')
for i=1:ngl
    om(i)=sqrt(fre(i,i));
end
disp(om(1))

% Grafico do primeiro modo (eixo)
X=zeros(ndx+1,1);mod1=zeros(ndx+1,1);
xe=0;
for i=1:ndx+1
    noeixo=(i-1)*(ndy+1)+ndy/2+1;
    ge=LN(noeixo,2);
    X(i,1)=xe;
    xe=xe+dx;
    if ge<=ngl
    mod1(i,1)=mod(ge,1);
    else
        mod1(i,1)=0;
    end
end

figure(1)
plot(X,mod1,'r',X,-mod1,'g')
axis equal
legend('Máximo','Mínimo')
    
% Grafico do primeiro modo (todos os pontos)
for i=1:nno
    for j = 1:2
        ge=LN(i,j);
        if ge<=ngl
            mod_vib1(i,j)=mod(ge,1);
        else
            mod_vib1(i,j)=0;
        end
    end
end    

figure(2)
hold on
patch('faces', nodel, 'vertices', gcoord+mod_vib1, 'FaceColor', 'none', 'LineStyle', '--', 'EdgeColor','red');
patch('faces', nodel, 'vertices', gcoord-mod_vib1, 'FaceColor', 'none', 'LineStyle', '--', 'EdgeColor','green');
axis equal
legend('Máximo','Mínimo')
hold off

%% ANIMAÇÃO DO PRIMEIRO MODO DE VIBRAÇÃO

n = 20;
dmodo = (mod_vib1+mod_vib1)/n;
altura = gcoord+mod_vib1;

for ii = 1:10
    for i = 1:n
        fig = figure(3);
        clf
        patch('faces', nodel, 'vertices', altura, 'FaceColor', 'blue', 'EdgeColor','blue');
        axis([0 compr -0.2 0.2])
        axis equal
        axis off
        altura = altura-dmodo;
        frame = getframe(fig);
        im{i} = frame2im(frame);
        pause(0.01)
    end
    for i = 1:n
        fig = figure(3);
        clf
        patch('faces', nodel, 'vertices', altura, 'FaceColor', 'blue', 'EdgeColor','blue');
        axis([0 compr -0.2 0.2])
        axis equal
        axis off
        altura = altura+dmodo;
        frame = getframe(fig);
        im{i+n} = frame2im(frame);
        pause(0.01)
    end
end

%% Criando GIF
% filename = 'Chapa Vibrando.gif'; % Specify the output file name
% for idx = 1:40
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.02);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
%     end
% end