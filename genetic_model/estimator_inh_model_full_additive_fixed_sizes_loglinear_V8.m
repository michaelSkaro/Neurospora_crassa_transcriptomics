% enter data from Emily
% imposes constraint on Lambda but leaves in the scoring
%used Genetics paper form to IRLS in Avise et al (1987)
% imposes sample size constraint for each cross!
clear all
N(1,1) = 103
N(1,2) = 167
N(1,3) = 126
N(2,1) = 142
N(2,2) =  86
N(2,3) =  58
N(3,1) = 170
N(3,2) = 108
N(3,3) =  88
NN = 0
k = 0
m = 6
for i=1:3
    for j=1:3
        k = k + 1
         O(k) = log(N(i,j));
        NN = NN + N(i,j);
    end
end
for i=1:3
    NR(i) = 0
    for j=1:3
        NR(i) = NR(i) + N(i,j)
    end
end
% use regression to initialize IRLS
X = [ 1 0 0 1 -1 0 1 0 0; 1 0 0 -1 1 0 1 0 0; 1 0 0 -1 -1 0 -1 0 0; 0 1 0 0 -1 -1 0 -1 1; 0 1 0 0 1 -1 0 1 0; 0 1 0 0 -1 1 0 1 0; 0 0 1 1 0 -1 0 0 1; 0 0 1 -1 0 -1 0 0 -1; 0 0 1 1 0 -1 0 0 1]
%X = [ 1 1 -1 0 1 0 0; 1 -1 1 0 1 0 0; 1 -1 -1 0 -1 0 0; 1 0 -1 -1 0 -1 1; 1 0 1 -1; 1 0 -1 1 0 1 0; 1 1 0 -1 0 0 1; 1 -1 0 -1 0 0 -1; 1 1 0 -1 0 0 1]
%X = [ 1 1 -1 0; 1 -1 1 0; 1 -1 -1 0; 1 0 -1 -1; 1 0 1 -1; 1 0 -1 1; 1 1 0 -1; 1 -1 0 -1; 1 -1 0 1]
X = [1 0 0 1 -1 0; 1 0 0 -1 1 0; 1 0 0 -1 -1 0; 0 1 0 0 -1 -1; 0 1 0 0 1 -1; 0 1 0 0 -1 1; 0 0 1 1 0 -1; 0 0 1 -1 0 -1; 0 0 1 -1 0 1]
H = X'*X
O = O'
theta0 = inv(H)*X'*O
Pred = (X*theta0)
Exp=exp(Pred)
Obs=exp(O)
X2 = 0;
for i=1:9
    X2 = X2 + ((Obs(i))-Exp(i))^2/Exp(i)
end
for i=1:3
      k = 3*(i-1)
    ER(i) = 0
    for j=1:3
        ER(i) = ER(i) + Exp(k+j)
    end
end
'initial',X2,theta0
%IRLS
beta = theta0
iter = 100
%start IRLS
error = 1000
while error > 10^-8
   % calculate Ks from model and parameters
K = exp(X*beta)
for i=1:3
   k = 3*(i-1)
   for j=1:3
       K(k+j) = K(k+j)/NR(i)
   end
end
    A = diag(K)
    A = inv(A)
   
         info = X'*A*X
         info = NN*info
         VAR = inv(info)
         S = A*O
         Y = S + NN*A*X*beta
         betastar=VAR*X'*Y
        Pred = (X*betastar)
        Exp = exp(Pred)
        for i=1:3
      k = 3*(i-1)
      % calculate expected marginal counts for each cross
    ER(i) = 0;
    for j=1:3
        ER(i) = ER(i) + Exp(k+j);
    end
        end
% After getting expected margins, update the betastar to enforce these
% margins
         for i=1:3
    betastar(i)=betastar(i) - log(ER(i)/NR(i));
         end
    % update predictions once the margins are enforced for each cross
      Pred = (X*betastar)
        Exp = exp(Pred)   
%calculate error in marginal counts of crosses
        diff = ER(1) -NR(1) + ER(2) - NR(2) + ER(3) - NR(3)
        Obs=exp(O)
        X2 = 0;
        for i=1:9
            X2 = X2 +(Obs(i)-Exp(i))^2/Exp(i);
        end
        k
        X2
        chisq(k) = X2;
        error = (betastar-beta)'*(betastar-beta)
        rel = (beta)'*(beta)
        error = error/rel
        beta = betastar;
        for i=1:m
            sd(i) = sqrt(VAR(i,i));
        end
        sdd = sd'
        parameters = [beta,sdd]
        Score(k) = sum(X'*Y);
        lnL=0;
        for i=1:9
            lnL = lnL + O(i)*log(K(i));
        end
        Like(k) = lnL;
        beta1(k) = beta(1);
        beta2(k) = beta(2);
        beta4(k) = beta(4);
        
end
subplot(2,2,1)
plot(1:k,Score)
ylabel('Score')
subplot(2,2,2)
plot(1:k,chisq)
ylabel('Chi-Squared')
subplot(2,2,3)
plot(1:k,Like)
ylabel('Loglikelihood')
subplot(2,2,4)
xx=1:k;
plot(xx,beta1,xx,beta2,xx,beta4)
legend('Lambda','alpha','alpha-beta')
ylabel('thetas')