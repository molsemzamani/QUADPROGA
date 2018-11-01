function [alpa,G] = semidf(A,c,B,b,lb,ub)
  n=size(A, 1);
  sl=size(lb,1);
  if sl>0 
    B=[B;-eye(n);eye(n)];
    b=[b;-lb;ub];
  end
  m=size(B, 1);
  M=[-b B]; 
  
  cvx_begin quiet
    variables x(n,m) l(m) alpa;
    variable lambda(m,m) nonnegative;
    
    maximize(alpa);
  
      subject to
       for  i= 1 : m,
         x(:,i)==-B'*lambda(:,i);
       end
      lambda'*b<=l;
      .5*[l';x]*M+.5*M'*[l';x]'+[-alpa c';c  A]== semidefinite(n+1);
    cvx_end
  G=.5*[l';x]*M+.5*M'*[l';x]';
end


