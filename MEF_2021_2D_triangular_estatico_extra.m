%nome do arquivo MEF_2021_2D_triangular_estatico_extra.m
%programa de EF para calculo estático
%chapas 2D elemento triangular
% Prof Reyolando Brasil - fevereiro 2021
%
clc
clear
%
%treliça
%nnel=2; %número de nós por elemento
%nglpn=2;%número de graus de liberdade por nó
%
%pórtico
%nnel=2; %número de nós por elemento
%nglpn=3;%número de graus de liberdade por nó
%
%elemento triangular de chapa
nnel=3; %número de nós por elemento
nglpn=2;%número de graus de liberdade por nó
%
%elemento retangular de chapa
%nnel=4; %número de nós por elemento
%nglpn=2;%número de graus de liberdade por nó
%entrada de dados
%
%dimensões do problema
%
nno=40;%número de nós
nds=nno*nglpn; %número de deslocamentos do sistema
nel=36;%número de elementos
ndpel=nnel*nglpn;%número de deslocamentos por elemento
ns=1;%numero de seções diferentes; se ns=1, todas iguais
%
%inicialização de matrizes e vetores
%
gcoord=zeros(nno,2);
nodel=zeros(nel,nnel);
LN=zeros(nno,nglpn);
LM=zeros(nel,nglpn*nnel);
k=zeros(ndpel,ndpel,nel);
B=zeros(nnel,ndpel,nel);
K=zeros(nds,nds);
p=zeros(nds,1);
P=zeros(nds,1);
q=zeros(ndpel,1);
P0=zeros(nds,1);
Q0=zeros(ndpel,nel);
AS=zeros(nel);
desloc=zeros(nno,nglpn+1);
forcas=zeros(nno,nglpn+1);
%
%
%coordenadas X Y dos nós
%
gcoord=[0 0;
0.5 0;
0 0.5;
0.5 0.5;
0 1;
0.5 1;
0 1.5;
0.5 1.5;
0 2;
0.5 2;
0 2.5;
0.5 2.5;
0.5 3;
4.5 0;
5 0;
4.5 0.5;
5 0.5;
4.5 1;
5 1;
4.5 1.5;
5 1.5;
4.5 2;
5 2;
4.5 2.5;
5 2.5;
4.5 3;
1 3;
1 3.5;
1.5 3;
1.5 3.5;
2 3;
2 3.5;
2.5 3;
2.5 3.5;
3 3;
3 3.5;
3.5 3;
3.5 3.5;
4 3
4 3.5];
%
%conectividade por elemento (numeração dos nós de cada elemento)
%
nodel=[1 2 4;
1 4 3;
3 4 6;
3 6 5;
5 6 8;
5 8 7;
7 8 10;
7 10 9;
9 10 12;
9 12 11;
11 12 13;
14 15 16;
16 15 17;
16 17 18;
18 17 19;
18 19 20;
20 19 21;
20 21 22; 
22 21 23;
22 23 24;
24 23 25;
24 25 26;
13 27 28;
27 29 28;
29 30 28;
29 31 30;
31 32 30;
31 33 32;
33 34 32;
33 36 34;
33 35 36;
35 38 36;
35 37 38;
37 40 38;
37 39 40;
39 26 40];
%
% matriz de número de graus de liberdade por nó
%
%condições de contorno (n. de nó restrito e direções restritas ou livres)
%se fixa = -1; se livre = 0
%
LN(1,:)=[-1 -1];
LN(2,:)=[-1 -1];
LN(14,:)=[-1 -1];
LN(15,:)=[-1 -1];
%
%determinação das matrizes
%
%matriz LN
%
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
%
%matriz LM
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
%
%dados físicos dos elementos, mód elasticidade, area da seção
%
t=0.25;%m
EM=10e9;%N/m² 
nu=0.25;
ro=1000;%kg/m³
%
%matriz constitutiva Estado Plano de Tensões
%
E=EM/(1-nu^2)*[1 nu 0;nu 1 0;0 0 (1-nu)/2];
%
%peso próprio dos elementos por vértice
PP=t*0.5^2*ro*10/6;%N, g=10 m/s²
%
%matriz de rigidez e vetor peso proprio
%
for iel=1:nel
    xa=gcoord(nodel(iel,1),1);ya=gcoord(nodel(iel,1),2);
    xb=gcoord(nodel(iel,2),1);yb=gcoord(nodel(iel,2),2);
    xc=gcoord(nodel(iel,3),1);yc=gcoord(nodel(iel,3),2);
%
%área dotriângulo
%
    am=[1 xa ya;1 xb yb;1 xc yc];
    A=0.5*det(am);
%
% matriz B=LN
%
    B(:,:,iel)=1/2/A*[yb-yc 0 yc-ya 0 ya-yb 0;0 xc-xb 0 xa-xc 0 xb-xa;...
        xc-xb yb-yc xa-xc yc-ya xb-xa ya-yb];
%
% matriz de rigidez do elemento k
%
    k=A*t*B(:,:,iel)'*E*B(:,:,iel);
%
%matriz global_
%
    for i=1:ndpel
        ii=LM(iel,i);
        for j=1:ndpel
            jj=LM(iel,j);
            K(ii,jj)=K(ii,jj)+k(i,j);
        end
    end
    %peso próprio
    for i=2:2:6
         Q0(i,iel)=PP;
    end
 %
    for i=1:ndpel
        P0(LM(iel,i))=P0(LM(iel,i))-Q0(i,iel);
    end
end
%
%cargas pontuais aplicadas nos nós
% V=500;%KN
% P(LN(4,2))=-V;
%
%vetor de forças nodais
P=P+P0;
%
% solução do sistema
%
AK=K(1:ngl,1:ngl);
disp('deslocamentos * 10^4')
b=(P(1:ngl)-K(1:ngl,ngl+1:nds)*p(ngl+1:nds));
p(1:ngl)=AK\b;
for i=1:nno
    desloc(i,:)=[i,10000*p(LN(i,1)),10000*p(LN(i,2))];
end
disp('    Nó        DX        DY')
disp(desloc)
%
disp('Esforços Nodais inclusive reações de apoio')
P(ngl+1:nds)=K(ngl+1:nds,1:ngl)*p(1:ngl)+K(ngl+1:nds,ngl+1:nds)...
    *p(ngl+1:nds);
P=P-P0;
for i=1:nno
    forcas(i,:)=[i,P(LN(i,1)),P(LN(i,2))];
end
disp('    Nó        FX        FY')
disp(forcas)
%
%tensões
%
disp('tensões nos elementos')
%
for iel=1:nel    
    
    for i=1:ndpel
        q(i)=p(LM(iel,i));
    end
    tau(iel,:)=(E*B(:,:,iel)*q)';
 
end
disp('sigma_x,    sigma_y,    tau_xy')
disp(tau);
    
% Plotagem da figura
figure(1)
hold on
patch('faces', nodel, 'vertices', gcoord+desloc(:,2:3),'FaceVertexCData',vecnorm(desloc(:,2:3),2,2),'FaceColor','interp');
patch('faces', nodel, 'vertices', gcoord, 'FaceColor', 'none', 'LineStyle', '--');
colormap(jet)
colorb = colorbar;
colorb.Label.String = 'Deslocamento [cm]';
axis equal
axis off

figure(2)
hold on
patch('faces', nodel, 'vertices', gcoord+desloc(:,2:3),'FaceVertexCData',vecnorm(desloc(:,2:3),2,2),'FaceColor','interp');
patch('faces', nodel, 'vertices', gcoord, 'FaceColor', 'none', 'LineStyle', '--');
colormap(jet)
colorb = colorbar;
colorb.Label.String = 'Deslocamento [cm]';
axis equal
axis off
