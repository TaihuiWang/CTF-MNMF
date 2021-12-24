function [sep,  HH] = bss_CTFMNMF1(mix, ns, nb, fftSize, shiftSize, it, refMic)
[X, window] = STFT( mix, fftSize, shiftSize, 'hamming' );

[Qori, T, V, piaoY] = CTFMNMF1( X, ns, it, nb);

[I,J,M] = size(X);

for i = 1:I
   HH(:,:,i) = inv(Qori(:,:,i));
end
N= ns;
Tao = M/N;
for i = 1:I
	for n = 1:ns
		SCM_MMNF(:,:,n,i) = HH(:,(n-1)*Tao+1:n*Tao,i) * HH(:,(n-1)*Tao+1:n*Tao,i)';
	end
end
piaoY = permute(piaoY,[3,4,1,2]); % N x Tao x I x J;
Xp = permute(X,[3,1,2]); % M x I x J;
for i = 1:I
    for j = 1:J
        lambda_ntao = [];
         for n = 1:N
             for tao = 1:Tao
                 lambda_ntao = [lambda_ntao T(i,:,n)*local_Vtao( V(:,:,n),j,tao )];
             end
         end
         fenmu = HH(:,:,i)*diag(lambda_ntao)*HH(:,:,i)';
         
         for n = 1:N
             lambda_tao = [];
             for tao = 1:Tao
                 lambda_tao = [lambda_tao T(i,:,n)*local_Vtao( V(:,:,n),j,tao )];
             end
             fenzi = HH(:,(n-1)*Tao+1:n*Tao,i) * diag(lambda_tao) * HH(:,(n-1)*Tao+1:n*Tao,i)';
             
             tempZ = fenzi*pinv(fenmu)*Xp(:,i,j);

             Z(i,j,n) = tempZ(refMic);    
         end
    end
end
sep = ISTFT( Z, shiftSize, window, size(mix,1) );
end

function [ HNew ] = local_Vtao( H,j,tao )
tao = tao - 1;
[K,J] = size(H);
if j <= tao
   HNew = H(:,j);
elseif j > tao
   HNew = H(:,j-tao);
else
   error('error in the local_Htao')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%