function [xq, adj,D]=lov(A,c,B,b,lb,ub,xq)% obtaining a local vertex optimal
   n=size(A, 1); 
   xq=cplexlp(A*xq+c,B,b,[],[],lb,ub);
   sl=size(lb,1);
  if sl>0 
    B=[B;-eye(n);eye(n)];
    b=[b;-lb;ub];
  end
  m=size(B, 1);
  br=1;
 while(br)
  JJ=[];
  JJb=[];
  B1=B;
  b1=b;
  con=0;
  for o=1:m
     if abs(B(o,:)*xq-b(o))<=10^(-5)
        JJ=[JJ;B(o,:)];
        JJb=[JJb;b(o)];
        B1(o-con,:)=[];
        b1(o-con,:)=[];
        con=con+1;
     end
  end

 D=[];
 adj=[];
 valadj=[];
 if rank(JJ)<n
     D=zeros(n);
     return;
 end
  for i=1:size(JJ,1)
      GH=JJ;
      GH(i,:)=[];
      [Uy Sy Vy] = svd(GH);
       x = Vy(:,end);
     if JJ(i,:)*x>0
         x=-x;
     end
     x=norm(x)^(-1)*x;
     D=[D x]; 
     tt=[];
     for j=1:m
         if B(j,:)*x>=10^(-7)
             tt=[tt (b(j)-B(j,:)*xq)/(B(j,:)*x)];
         end
     end
     ax=xq-max(-tt)*x;
      adj=[adj ax];
      valadj=[valadj ax'*A*ax+2*ax'*c];
 end
  
adj=[xq adj];
valadj=[xq'*A*xq+2*xq'*c valadj];
[fffe idx]=max(-valadj);
  xq=adj(:,idx);

   if idx==1
      br=0;
   end
 end
adj(:,1)=[];
valadj(:,1)=[];
end

