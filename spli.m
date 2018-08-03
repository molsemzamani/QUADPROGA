function x= spli(B,b,lb,ub)%obtaing Chebyshev center
%cvx_solver  sedumi
  n=size(B, 2);
  sl=size(lb,1);
  if sl>0 
    B=[B;-eye(n);eye(n)];
    b=[b;-lb;ub];
  end
  m=size(B, 1);
  D=zeros(m,1);
  for  i= 1 : m,
         D(i)=norm(B(i,:));
  end
  c=zeros(n+1,1);
  c(n+1)=-1;
  b1=[b;b];
  B1=[B D;B zeros(m,1)];
  x=cplexlp(c,B1,b1);
  x=x(1:n);
end
%{    
 cvx_begin quiet
variables d(n); 
 variable S(n,n) semidefinite
    maximize(log_det(S));
 
     subject to
      for  i= 1 : m,
        norm(S*B(i,:)')+B(i,:)*d<=b(i);
      end
      
    cvx_end
    [V,D] = eig(S);
    n=V(:,1)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cvx_begin quiet
    variables x(n) alpa; 
     maximize(alpa);
  
      subject to
      for  i= 1 : m,
         B(i,:)*x<=b(i,:)-alpa*norm(B(i,:));
      end
      B*x-b<=0;
    cvx_end
%}