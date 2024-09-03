function  [Tsigma, c_10, c, c_90, m] = FatigueScatter(CP, Nf, d_10, d_90)

% Fatigue life prediction

LogNf(:,2) = log(Nf);
LogNf(:,1) = 1;
SOL = pinv(LogNf)*log(CP);
c = exp(SOL(1,:));
m = SOL(2,:);

Nexpected = (CP./c).^(1./m);
Residual(:,1) = log(Nf) - log(Nexpected);
S = std(Residual);

c_10 = c./exp(S*d_10.*m);
c_90 = c./exp(S*d_90.*m);


CP10 = c_10*1000000^m;

Cp90 = c_90*1000000^m;

Tsigma = CP10/Cp90;

