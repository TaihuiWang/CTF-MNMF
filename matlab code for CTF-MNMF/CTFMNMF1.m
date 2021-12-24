function [Q,  T, V, piaoY] = DerbMNMF(X, N, It, nb)
%% 
% X is a I x J x M tensor;
% Q is a M x N x I tensor;
% T is a I x K x N tensor;
% V is a K x J x N tensor;
%% initialization of parameters
[I,J,M] = size(X);
x = permute(X,[3,1,2]); % M x I x J;
xp = permute(X,[3,2,1]); % M x J x I;
Tao = M/N;
if M > I 
    error('The size of X might be wrong: frequency x frame x channel');
elseif M > J
    error('The size of X might be wrong: frequency x frame x channel');
end
K = nb;
if nargin < 5 
    T = max(eps, rand(I,K,N));
    V = max(eps, rand(K,J,N));
end
if nargin < 7
    tempQ = eye(M) +0.9*(ones(M)-eye(M));
    for i = 1:I
        Q(:,:,i) = tempQ  ;  % M x NTao x I
    end  
end
for i = 1:I
    for j = 1:J
        for n = 1:N
            for tao = 1:Tao
              piaoX(i,j,n,tao) =  abs(  Q((n-1)*Tao+tao ,:,i)*x(:,i,j) ).^2;
            end
        end
    end
end
for n = 1:N
    for tao = 1:Tao
        piaoY(:,:,n,tao) =T(:,:,n)*local_Vtao(V(:,:,n),tao); % I x J x N x Tao
    end
end
UN1 = ones(M,1);
U1J = ones(1,J);
E = eye(M);
fprintf('Iteration:    ');
for it = 1:It
    fprintf('\b\b\b\b%4d', it);
    %% update T  
    ZfenziT = piaoX.*(piaoY.^(-2));  % I x J x N x T
    ZfenmuT = piaoY.^(-1);
    for n = 1:N  
        Tfenzi = zeros(I,K);
        Tfenmu = zeros(I,K);
        for tao = 1:Tao
            Tfenzi = Tfenzi + ZfenziT(:,:,n,tao)*local_Vtao(V(:,:,n),tao)';
            Tfenmu = Tfenmu + ZfenmuT(:,:,n,tao)*local_Vtao(V(:,:,n),tao)';
        end
        T(:,:,n) = max(eps, T(:,:,n).*sqrt(Tfenzi./Tfenmu));
    end
    T = abs(T);
    
    %% update piaoY
    for n = 1:N
        for tao = 1:Tao
            piaoY(:,:,n,tao) =T(:,:,n)*local_Vtao(V(:,:,n),tao); 
        end
    end
    
    %% update V
    ZfenziV = piaoX.*(piaoY.^(-2)); % I x J x N x Tao
    ZfenmuV = piaoY.^(-1);
    tao = 1;
	for n = 1:N
        Vfenzi = T(:,:,n)'*ZfenziV(:,:,n,tao);
        Vfenmu = T(:,:,n)'*ZfenmuV(:,:,n,tao);
        V(:,:,n) = max(eps, V(:,:,n).*sqrt(Vfenzi./Vfenmu));
    end
    
    %% update piaoY
    for n = 1:N
        for tao = 1:Tao
            piaoY(:,:,n,tao) =T(:,:,n)*local_Vtao(V(:,:,n),tao); 
        end
    end
    
    %% update Q using IP
    for n = 1:N
        for i = 1:I
            for tao = 1:Tao
                D = ((xp(:,:,i).*(UN1*(U1J./piaoY(i,:,n,tao)))  )*xp(:,:,i)')/J;                
                q = pinv(Q(:,:,i)*D+10000*eps*eye(M))*E(:,(n-1)*Tao+tao);
                q = q/sqrt((q')*D*q);
                Q((n-1)*Tao+tao,:,i) = q';
                piaoX(i,:,n,tao) = abs((q')*xp(:,:,i)).^2;
            end
        end
    end
end
end

function [ HNew ] = local_Vtao( H,tao )
tao = tao - 1;
[K,J] = size(H);
    for j = 1:J
       if j <= tao
           HNew(:,j) = H(:,j);
       elseif j > tao
           HNew(:,j) = H(:,j-tao);
       else
           error('error in the local_Htao')
       end
    end
end
