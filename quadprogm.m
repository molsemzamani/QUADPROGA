function [x,fval,time,logq,lb] = quadprogm(A,c,B,b,lb,ub)
%% [x,fval,time,logq] = quadprogbb(H,f,A,b,Aeq,beq,LB,UB,options)
%
% QUADPROGM globally solves the following nonconvex quadratic
% programming problem:
%
%    min      x'*A*x +2c'*x
%    s.t.       B * x <= b
% where the interior of feasible set is not empty and bounded. 
% To apply this code you need install CVX and CPLEX pachage
% x reurns optimal solution
% fval reurns optimal value
% time reurns implementation time
% lb reurns lower bound
 % x reurns nodes details
tic;
m=size(B,1);
n=size(A,1);
i=1;
L=64;% the number of alloable branching
br(i).B=B;
br(i).b=b;
br(i).inx=i;
br(i).lo=inf;
br(i).up=-inf;
br(i).f=0;% father index
eps=10^(-5);
val_p=inf;
val_d=-inf;
VAL_D=inf(1,L);
VAL_P=[];
ly=1;% index saving log
x_m=zeros(n,1);
tic
 while(1)
     tf = arrayfun(@(k) ~isempty(br(k).inx), 1:numel(br));
     sq=size(tf,2);
     VAL_P1=[];
     for kk=1:sq
         k=1;
            [alpa,G] = semidf(A,c, br(k).B, br(k).b,lb,ub);
             if br(k).f>=1
                  VAL_D(br(k).f)=inf;
              end
            br(k).up=alpa;
            VAL_D(br(k).inx)=alpa;
            if alpa<=val_p-eps  %2
               xq=lquad(G, A, c,B,b,lb,ub,x_m);% solving a convex QP
               xq=lov(A,c,B,b,lb,ub,xq);% obtaining a local vertex optimal
                 val_p=min(val_p,xq'*A*xq+2*xq'*c);
                 br(k).lo=xq'*A*xq+2*xq'*c;
                   if val_p>=xq'*A*xq+2*xq'*c-5*eps
                      x_m=xq;
                      VAL_P1=[VAL_P1 xq'*A*xq+2*xq'*c];
                   end
               end                           
            if alpa<=val_p-eps   %branching
                c1=randn(1,n);
                x= spli(br(k).B,br(k).b,lb,ub);%obtaing Chebyshev center
                c1=(-1/norm(c1))*c1;
                br(i+1).B=[br(k).B;-c1];
                br(i+1).b=[br(k).b;-c1*x];
                br(i+1).inx=i+1;
                br(i+1).lo=inf;
                br(i+1).up=-inf;
                br(i+1).f=br(k).inx;% father index
                % second branch
                br(i+2).B=[br(k).B;c1];
                br(i+2).b=[br(k).b;c1*x];
                br(i+2).inx=i+2;
                br(i+2).lo=inf;
                br(i+2).up=-inf;
                 br(i+2).f=br(k).inx;% father index
                i=i+2;
            end
            logq(ly)=br(k);
            ly=ly+1;
            br(k)=[];
            br= br(~cellfun(@isempty,{br.inx}));%delete empty row 
     end
     VAL_P=[VAL_P min(VAL_P1)];
    val_d=max(min(VAL_D));   
    if sq==0 | val_p-val_d<=eps | i>=L  
         break
    end
    m_VAL_P=size(VAL_P,2);
    if m_VAL_P>=3 & VAL_P(m_VAL_P-2)-VAL_P(m_VAL_P)<=5*eps
        break
    end
 end
time=toc;
x=x_m;
fval=val_p;
lb=val_d;

end