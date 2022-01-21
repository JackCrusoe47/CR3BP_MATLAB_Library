function [SVp,tp] = newton_multipleshootingDiffCorrector_endFixed_forward(SVp,tp,mu,Tol,IterMax)
% --- Preprocessing
% -- Number of segments
N = length(tp)-1;
% -- Number of variables
nvar = 7*(N-1)+4; 
% -- Number of constraints
ncon = 6*(N-1)+3;
% -- Initializing constraint vector
F = zeros(ncon,1);
% -- Initializing derivative matrix
dF = zeros(ncon,nvar);
% -- time of flight of each segment
tf = tp(2:end)-tp(1:end-1);
% -- definiting variable vector
X = [SVp(4:6,1); reshape(SVp(:,2:end-1),6*(N-1),1); tf'];
% -- Initializing iteration
Iter = 0;

% --- Newton-Raphson loop
while Iter<IterMax
   % -- Looping through all segments
   for n = 1:N
       % - recovering initial state vector components
       if n==1
           SVp(4:6,n) = X(1:3);
       else
           SVp(:,n) = X(6*(n-2)+4:6*(n-2)+9);
       end
       % - recovering integration time span
       tf(n) = X(6*(N-1)+3+n);
       % - integrating forward with STM
       [SV,~,STM] = propagateTrajectory_CR3BP(SVp(:,n),tf(n),mu,0,...
           'ComputeSTM',true,'InitialSTM',eye(6,6));
       % - State vector at the end of computation
       SV = SV(:,end);
       % - State trasition matrix at end of computation
       STM = STM(:,:,end);
       % - Derivative state vector
       dSV = getSV_dot_CR3BP(SV,mu);
       
       % - adding elements to constraint vector
       if n<N
           F(6*(n-1)+1:6*(n-1)+6) = SV-SVp(:,n+1);
       else
           F(6*(n-1)+1:6*(n-1)+3) = SV(1:3)-SVp(4:6,n+1);
       end
       
       % - adding elements to derivative matrix
       if n==1
          dF(1:6,1:3) = STM(1:6,4:6);
          dF(1:6,4:9) = -eye(6,6);
          dF(1:6,6*(N-1)+3+n) = dSV;
       elseif n<N
           dF(6*(n-1)+1:6*(n-1)+6,6*(n-2)+4:6*(n-2)+9) = STM;
           dF(6*(n-1)+1:6*(n-1)+6,6*(n-2)+10:6*(n-2)+15) = -eye(6,6);
           dF(6*(n-1)+1:6*(n-1)+6,6*(N-1)+3+n) = dSV;
       else
           dF(end-2:end,6*(n-2)+4:6*(n-2)+9) = STM(1:3,1:6);
           dF(end-2:end,6*(N-1)+3+n) = dSV(1:3);
       end
   end
    
   % -- Newton-Raphson Step
   Xnew = X - dF'/(dF*dF')*F;
   
   % -- Breaking Conditon
   if norm(F)<Tol
       break
   end
   
   X = Xnew;
   
   Iter = Iter+1;
end

X = Xnew;

end