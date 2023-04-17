%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                  Prova final: Placas e Cascas                        %%
%%                         QUESTÃO 1                                    %%
%%                 Gabriel Prieto Paris RA: 11052516                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Código de difereças finitas para resolver numericamente a deformação   %
% de uma placa. Utilizando a equação de Sophie Germain/Lagrange.         %
% # Resultados:                                                          %
%   - Gráfico de deslocamento vertical máximo e frequência em função     % 
%     da espessura.                                                      %
%   - Espessura mínima para utilizar no projeto.                         %
%   - Distribuição de deslocamento vertical dada a espessura mínima      %
%   - Animação da vibração do modo escolhido.                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

clc
clear
close all

%% DADOS DE ENTRADA

% Propriedade físicas
c = 5.5; % [m] - Média dos último e antepenúltimo algarismos do RA
a = 0.3*c; b = 0.2*c; % [m] - Dimensões do retângulo
EM = 70e6; % [KPa] - Módulo de elasticidade 
nu = 0.25; % Coeficiente de Poisson
rho = 2700; % Massa específica
p = -100; % [KPa] - Carga Uniformemente distribuída

% Divisão da malha
nd=16; % Número de divisões da malha

%Condições de contorno: (paralela a x) isup, iinf - (paralela a y) iesq, idir
% código: -1=apoiado; 1=engastado
isup=-1;iinf=1;iesq=1;idir=-1;

% Modo de vibração
modo_de_vibrar = 1;

% o maior deslocamento e a menor frequência
deslocamento_max  = 27.5; % [mm]
frequencia_min  = 2; % [Hz]

espessura = 0.01:0.00001:0.07; % [m] - Espessura;

%% Cálculos da espessura

% Aplicação do método das diferenças finitas em diversas espessuras
[~, wmax, ~, omega_modo] = Placas(a, b, espessura, EM, nu, rho, p, nd, isup, iinf, iesq, idir, modo_de_vibrar);
frequencia = omega_modo/(2*pi); % Convertendo a frequência de rad/s para Hz

% Primeiramente achamos quais as localizaçoes onde a diferênça entre todos
% os deslocamentos maximos e o deslocamento especificado são menores do que zero. Depois
% encontramos o maior valor de wmax dentro das localizações encontradas.
% Por fim, em que posição dentro de wmax se encontra esse valor.
loc = find(wmax-deslocamento_max<0);
maior = max(wmax(loc));
maior_loc = find(wmax == maior);
% maior_loc = find(abs(wmax - maior) < 0.01);

% Primeiramente achamos quais as localizaçoes onde a diferênça entre todos
% as frequencias do primeiro modo e a frequencia especificada são maiores do que zero. Depois
% encontramos o menor valor de frequencia dentro das localizações encontradas.
% Por fim, em que posição dentro de frequencia se encontra esse valor.
loc2 = find(frequencia-frequencia_min>0);
menor = min(frequencia(loc2));
menor_loc = find(frequencia == menor);
% menor_loc = find(abs(frequencia - menor) < 0.1);

% A espessura utilizada sera a maior espessura entre as duas encontradas
espessura_utilizada = max([espessura(maior_loc) espessura(menor_loc)]);
fprintf('Espessura mínima que deve ser adotada no projeto: %.3f mm \n\n', 1000*espessura_utilizada)

%% PLOTAGEM DAS FIGURAS

% Gráficos de deslocamento vertical máximo e frequência para cada espessura
figure
sgtitle('Deslocamento vertical máximo e frequência para cada espessura')
p1 = plot(espessura,wmax,'b','linewidth',2);
datatip(p1,espessura(maior_loc),wmax(maior_loc));
p1.DataTipTemplate.DataTipRows(1).Label = "Espessura [m]: ";
p1.DataTipTemplate.DataTipRows(2).Label = "Deslocamento vertical máximo [mm]: ";
ylabel('Deslocamento vertical máximo [mm]')
hold on
yyaxis right
p2 = plot(espessura,frequencia,'r','linewidth',2);
p2.DataTipTemplate.DataTipRows(1).Label = "Espessura [m]: ";
p2.DataTipTemplate.DataTipRows(2).Label = "Primeira frequência de vibração livre [Hz]: ";
datatip(p2,espessura(menor_loc),frequencia(menor_loc));
ylabel('Primeira frequência de vibração livre [Hz]')
grid on
grid minor
set(gcf,'Position',[150 150 700 550])
saveas(gcf,'Gráficos de deslocamento vertical máximo e frequência para cada espessura.png')

% Aplicando o método das diferenças finitas para a espessura que será
% utilizada no projeto.
[w, ~, modo, ~] = Placas(a, b, espessura_utilizada, EM, nu, rho, p, nd, isup, iinf, iesq, idir, modo_de_vibrar);

answer = questdlg('Gostaria de plotar os próximos gráficos?', ...
	'Menu', ...
	'Sim','Não','Não');

switch answer
    case 'Sim'
        % Deslocamento
        figure
        [X,Y] = meshgrid(-a/2:a/nd:a/2,-b/2:b/nd:b/2);
        Z = zeros(nd+1);
        c = 0;
        for i = 2:nd+1-1
            for j = 2:nd+1-1
                c = c+1;
                Z(i,j) = w(c);
            end
        end
        surf(X,Y,Z, 'FaceColor', 'interp') % 'EdgeColor', 'none'
        colormap(flipud(jet))
        colorb = colorbar;
        colorb.Label.String = 'Deslocamento [m]';
        title('Distribuição de deslocamento vertical')
        axis equal
        axis off
        saveas(gcf,'Distribuição de deslocamento vertical.png')        
    
        % ANIMAÇÃO DO MODO DE VIBRAÇÃO ESCOLHIDO
        z = reshape(modo(:,1)',nd-1,nd-1)';
        Z = zeros(size(z)+2);
        Z(2:end-1,2:end-1)=z;
        z=Z;
        n = 20;
        dmodo = (z+z)/n;
        altura = z;
        
        for ii = 1:5
            for i = 1:n
                fig = figure(3);
                clf
                surf(X,Y,altura, 'FaceColor', 'interp');
                axis equal
                axis([-a/2 a/2 -b/2 b/2 -max(abs(z),[],'all') max(abs(z),[],'all')])
                axis off
                colormap(jet);
                caxis([-max(abs(z),[],'all') max(abs(z),[],'all')]);
                %         view(0,0)
                altura = altura-dmodo;
                frame = getframe(fig);
                im{i} = frame2im(frame);
                pause(0.01)
            end
            for i = 1:n
                fig = figure(3);
                clf
                surf(X,Y,altura, 'FaceColor', 'interp');
                axis equal
                axis([-a/2 a/2 -b/2 b/2 -max(abs(z),[],'all') max(abs(z),[],'all')])
                axis off
                colormap(jet);
                caxis([-max(abs(z),[],'all') max(abs(z),[],'all')]);
                %         view(0,0)
                altura = altura+dmodo;
                frame = getframe(fig);
                im{i+n} = frame2im(frame);
                pause(0.01)
            end
        end
        %% Criando GIF
%         filename = 'Placa Vibrando.gif'; % Specify the output file name
%         for idx = 1:40
%             [A,map] = rgb2ind(im{idx},256);
%             if idx == 1
%                 imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.02);
%             else
%                 imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
%             end
%         end

end


fprintf('\n * PROGRAMA FINALIZADO * \n')

%% FUNÇÕES

function [w, wmax, modo, omega_modo] = Placas(a, b, espessura, EM, nu, rho, p, nd, isup, iinf, iesq, idir, modo_de_vibrar)

np1=nd+1; % Número de nós em cada direção
np3=nd+3; % Número de nós da malha expandida em cada direção (Nós fictícios)
neq=(nd-1)^2; % Número de equações a resolver
hx=a/nd;hy=b/nd; % Dimensões do elemento
hx4=hx*hx*hx*hx;hy4=hy*hy*hy*hy;hx2hy2=hx*hx*hy*hy;
hx2=hx*hx;hy2=hy*hy;

% Criar Barra de carregamento
steps = length(espessura); % Número total de iterações
f = waitbar(0,'1','Name','Calculando espessura...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

it = 0;
for t = espessura
    it=it+1;
    
    % Checando se o botão Cancel foi apertado
    if getappdata(f,'canceling')
        break
    end
    % Barra de carregamento
    waitbar(it/steps,f,sprintf('%2.1f %%',(it/steps)*100))    
    
    % Calculo de p/D
    D = EM*t^3/(12*(1-nu^2));
    pdD=p/D;
    
    %% GERAÇÃO DA MALHA
    
    % Inicialização
    A=zeros(np3*np3,np3*np3); % Número de nós da malha = np3*np3 (Inclusive nós fictícios)
    B=zeros(np3*np3,1);
    kk=zeros(neq,1);
    aa=zeros(neq,neq);
    bb=zeros(neq,1);
    
    % Geração das equações de diferenças finitas para os nós internos
    % (Os nós são contados da esquerda para direita de cima para baixo)
    % As linhas da matriz A representam os nós e as colunas cada coeficiente na
    % frente de cada icognita.
    m=0;
    for i=3:np1 % Número da linha
        k=(i-1)*np3+2; % número do nó, apenas internos
        for j=3:np1 % Número da coluna
            k=k+1;m=m+1;kk(m,1)=k; % kk é o vetor do número dos nós internos
            % Coeficientes na frente de w
            A(k,k-2)=1/hx4; % w_j,k-2
            A(k,k-1)=-4/hx4-4/hx2hy2; % w_j,k-1
            A(k,k+1)=-4/hx4-4/hx2hy2; % w_j,k+1
            A(k,k+2)=1/hx4; % w_j,k+2
            A(k,k)=6/hx4+6/hy4+8/hx2hy2; % w_j,k
            A(k+2*np3,k)=1/hy4; % w_j-2,k
            A(k+np3,k)=-4/hy4-4/hx2hy2; % w_j-1,k
            A(k-np3,k)=-4/hy4-4/hx2hy2; % w_j+1,k
            A(k-2*np3,k)=1/hy4; % w_j+1,k
            A(k+np3,k+1)=2/hx2hy2; % w_j+1,k+1
            A(k+np3,k-1)=2/hx2hy2; % w_j+1,k-1
            A(k-np3,k+1)=2/hx2hy2; % w_j-1,k+1
            A(k-np3,k-1)=2/hx2hy2; % w_j-1,k-1
            % Depois do igual no sistema
            B(k,1)=pdD;
        end
    end
    
    % Bordo superior
    k=2*np3+2;
    for j=3:np1
        k=k+1;
        A(k,k)=A(k,k)+isup*A(k-2*np3,k); % w_j,k = +- w_j-2,k o ponto fora da malha é igual a +- ao ponto adjacente da fronteira
    end
    % Bordo inferior
    k=nd*np3+2;
    for j=3:np1
        k=k+1;
        A(k,k)=A(k,k)+iinf*A(k+2*np3,k); % w_j,k = +- w_j+2,k
    end
    
    % Bordo esquerdo
    k=2*np3+3;
    for i=3:np1
        A(k,k)=A(k,k)+iesq*A(k,k-2); % w_j,k = +- w_j,k-2
        k=k+np3;
    end
    % Bordo direito
    k=3*np3-2;
    for i=3:np1
        A(k,k)=A(k,k)+idir*A(k,k+2); % w_j,k = +- w_j,k+2
        k=k+np3;
    end
    
    %% SOLUÇÃO
    
    % Montagem sistema de equações
    for i=1:neq
        for j=1:neq
            aa(i,j)=A(kk(i),kk(j));
        end
        bb(i)=B(kk(i));
    end
    
    % Solução do sitema
    w = aa\bb;
    wmax(it) = 1000*max(abs(w)); % Deslocamento máximo em mm
    % fprintf('Deslocamento Máximo:\n %f mm \n\n', wmax)
    
    % Momentos fletores no nó central da malha
    nnc=(neq+1)/2;
    if nnc>1
        Wdir=w(nnc+1);Wcen=w(nnc);Wesq=w(nnc-1);
        Wsup=w(nnc-nd+1);Winf=w(nnc+nd-1);
    else
        Wdir=0;Wcen=w(nnc);Wesq=0;
        Wsup=0;Winf=0;
    end
    Mx=-D*((Wesq-2*Wcen+Wdir)/hx2+nu*(Wsup-2*Wcen+Winf)/hy2);
    My=-D*(nu*(Wesq-2*Wcen+Wdir)/hx2+(Wsup-2*Wcen+Winf)/hy2);
    
    % Solução do problema de autovalores e autovetores
    [modo,lambda] = eig(aa);
    omega = sqrt(lambda*D/rho);
    omega_modo(it) = omega(modo_de_vibrar);
    % fprintf('Frequência do %iº modo de vibração:\n %f rad/s \n\n', modo_de_vibrar, omega_modo)
    K=1;N=0.25;
    omega_t = (pi/a)^2 * sqrt(D*K/(rho*N));
end
delete(f)
end