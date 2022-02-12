x = [0 0 0 0 0 0 1 0 1 0 0 1];
K = 3;
N = 4;
G = 2;
M = 4; % 16-QAM, muc dieu che
%checkPilot = true;
checkPilot = false;
numSymbol = length(x)/K;
A = reshape(x,[K,numSymbol])
B = [];
iter = numSymbol/log2(M);
for i = 1:iter
    temp = [bi2de(A(:,i:i+log2(M)-1),'left-msb')];
    B = [B qammod(temp,M);]
end
if checkPilot == true
    df = 1;
    dt = 2;
    pilot = [1 2 3 4];  % pilot must be not zero
    temp = zeros(1, numel(B) + numel(pilot));
    temp = reshape(temp,K,[]);
    l = size(temp,2);
    row = [1 (1:(K-1)/df)*df+1];
    col = [1 (1:(l-1)/dt)*dt+1];
    temp(row,col) = 1;
    temp(temp > 0 ) = pilot;
    temp(temp == 0) = B;
    B = temp;
end
numSymbol = size(B,2);
C = ones(N,numSymbol);
if (N-K)~= 0   
    col = [1:numSymbol];
    row = [1 (K/2+2:K/2+N-K)];
    C(row,col) = 0;
    C(C ~= 0) = B;
else
    C = B;
end
C
D = ifft(C)
E = [];
Guard = D(end-G+1:end,:);
E = [Guard;D]
F = reshape(E,1,[]);
F1 = real(F); %real part
F2 = imag(F); %imaginary part
G1 = kron(ones(20,1),F1);
G1 = G1(:);
subplot(5,1,1);
plot(G1)
G2 = kron(ones(20,1),F2);
G2 = G2(:);
subplot(5,1,2);
plot(G2)
t = 0:2*pi/20:100*pi;
xs = sin(t);
xc = cos(t);
xs = xs(1:size(G1));
xc = xc(1:size(G2));
H1 = G1.*xs';
subplot(5,1,3);
plot(H1);
H2 = G2.*xc';
subplot(5,1,4);
plot(H2);
I = H1 + H2;
subplot(5,1,5);
plot(I);