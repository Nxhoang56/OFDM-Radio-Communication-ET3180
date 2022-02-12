x = [0 0 0 1 0 0 1 0 0 0 1 1 0 1 0 0];
K = 4;
N = 4;
G = 2;
M = 16; % M-QAM
checkPilot = true;
%checkPilot = false;
df = 1;
dt = 2;
pilot = [1 2 3 4];  % pilot must be not zero
numSymbol = length(x)/K;
A = reshape(x,[K,numSymbol])
B = [];
iter = numSymbol/log2(M);

%const = [1 1i -1i -1]; %star constellation with order 0 1 2 3
const = [-3+3i -3+i -3-3i -3-i -1+3i -1+i -1-3i -1-i 3+3i 3+i 3-3i 3-i 1+3i 1+i 1-3i 1-i]; %chom sao theo 16QAM theo ma gray
for i = 1:iter
    temp = [bi2de(A(:,(i-1)*log2(M)+1:i*log2(M)),'left-msb')];
    B = [B qammod(temp,M,'PlotConstellation',true);];
    %B = [B genqammod(temp,const)];
end
if checkPilot == true
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
B
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
F1 = real(F) %real part
F2 = imag(F) %imaginary part
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