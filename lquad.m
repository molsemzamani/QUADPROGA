function xq=lquad(G,A,c,B, b,lb,ub,x0)% solvig convex QP

  n=size(G, 1)-1;
 
 G2=G(2:n+1,2:n+1)+A;
 f11=G(1,:);
 f1=2*f11(2:n+1)'+2*c;
 G2=G2+G2';
 opts1=optimset('display','off');
 xq=quadprog(G2,f1,B,b,[],[],lb,ub,x0,opts1);

end

