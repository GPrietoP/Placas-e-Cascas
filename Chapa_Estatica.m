%% Gabriel Prieto Paris
% P1 - PLACAS E CASCAS
close all
clear
clc

%% DADOS DE ENTRADA

% Elemento retangular de chapa
nnel=4; % Número de nós por elemento
nglpn=2;% Número de graus de liberdade por nó

% Dimensoes do retangulo
compr=44;
alt=11;

% Divisoes
ndx=8;ndy=2; % Elementos em x e elementos em y
nno=(ndx+1)*(ndy+1);% Numero de nos
nel=ndx*ndy;% Numero de elementos

% Dimensões do problema
nds=nno*nglpn; % Número total de deslocamentos do sistema
ndpel=nnel*nglpn;% Número de deslocamentos por elemento
ns=1;% Numero de seções diferentes; se ns=1, todas iguais

% Dados físicos dos elementos, mód elasticidade, area da seção
t=5.5; EM=100000; nu=0.25;
EL=EM/(1-nu*nu);G=EM/2/(1+nu);
disp('# Módulo de elasticidade:')
E=EL*[1 nu 0;nu 1 0;0 0 (1-nu)/2]
PP=0;%peso proprio

%% Inicialização de matrizes e vetores
gcoord=zeros(nno,2);
nodel=zeros(nel,nnel);
LN=zeros(nno,nglpn);
LM=zeros(nel,nglpn*nnel);
k=zeros(ndpel,ndpel);
B=zeros(nnel,ndpel,nel);
K=zeros(nds,nds);
p=zeros(nds,1);
P=zeros(nds,1);
q=zeros(ndpel,1);
P0=zeros(nds,1);
Q0=zeros(nel,ndpel);
AS=zeros(nel);
desloc=zeros(nno,nglpn+1);
forcas=zeros(nno,nglpn+1);

%% NÓS

% Geracao das coordenadas
dx=compr/ndx;dy=alt/ndy;
gcoord=zeros(nno,2);
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
LN(ceil((ndy+1)/2),:)=[-1 -1];
LN(ceil(nno-(ndy+1)/2),:)=[0 -1];

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

% Matriz de rigidez e vetor peso proprio
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

%   Simetria
    for i=2:8
        for j=1:i-1
            k(i,j)=k(j,i);
        end
    end
% fprintf('K%i = \n',iel)
% disp(k)
    
%   Matriz global
    for i=1:ndpel
        ii=LM(iel,i);
        for j=1:ndpel
            jj=LM(iel,j);
            K(ii,jj)=K(ii,jj)+k(i,j);
        end
    end
    
%   Peso próprio
    for i=2:2:8
         Q0(iel,i)=PP;
    end
 
    for i=1:ndpel
        P0(LM(iel,i))=P0(LM(iel,i))-Q0(iel,i);
    end
end
% K(1:ngl,1:ngl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alterar vetor de carregamento P
% Carga horizontal nos nós de cima
M=5500*9.81;
P(LN(ceil(nno/2),2))=-M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SOLUÇÃO

% Vetor de forças nodais
P=P+P0;

% Solução do sistema
AK=K(1:ngl,1:ngl);
disp('# Deslocamentos: ')
b=(P(1:ngl)-K(1:ngl,ngl+1:nds)*p(ngl+1:nds));
p(1:ngl)=AK\b;
for i=1:nno
    desloc(i,:)=[i,p(LN(i,1)),p(LN(i,2))];
end
disp('    Nó        DX        DY')
disp(desloc)

% Tensões
for iel=1:nel
    xa=gcoord(nodel(iel,1),1);
    xb=gcoord(nodel(iel,2),1);
    yb=gcoord(nodel(iel,2),2);
    yc=gcoord(nodel(iel,3),2);

%  DimensÃµes do retangulo
    a=(xa-xb)/2;b=(yb-yc)/2;

%  Constantes
    ca=1/4/a;cb=1/4/b;
    
%  matriz B=L*N calculada no centro do elemento x=y=0
    B=[ca 0 -ca 0 -ca 0 ca 0;0 cb 0 cb 0 -cb 0 -cb;cb ca cb -ca -cb -ca -cb ca];

    for i=1:ndpel
        q(i)=p(LM(iel,i));
    end
    tau(iel,:)=[iel (E*B*q)'];
end
disp('Tensões nos elementos')
disp('    Elemento  sigma_x   sigma_y   tau_xy')
disp(tau)

% Plotagem da figura
figure(1)
hold on
patch('faces', nodel, 'vertices', gcoord+desloc(:,2:3),'FaceVertexCData',vecnorm(desloc(:,2:3),2,2),'FaceColor','interp');
patch('faces', nodel, 'vertices', gcoord, 'FaceColor', 'none', 'LineStyle', '--', 'EdgeColor',[0.5 0.5 0.5]);
colormap(jet)
colorb = colorbar;
colorb.Label.String = 'Deslocamento [cm]';
axis equal
axis off

% figure(2)
% hold on
% patch('faces', nodel, 'vertices', gcoord+desloc(:,2:3),'FaceVertexCData',abs(tau(:,3)),'FaceColor', 'flat');
% patch('faces', nodel, 'vertices', gcoord, 'FaceColor', 'none', 'LineStyle', '--', 'EdgeColor',[0.5 0.5 0.5]);
% colormap(jet)
% colorb = colorbar;
% colorb.Label.String = 'Deslocamento [cm]';
% axis equal
% axis off