function [e,ce,ci,g,ae,ai,hl,indic] = chs(indic,xy,lme,lmi)
  global A B L R S
  
  e=0;
  ce=[];
  ci=[];
  g=[];
  ae=[];
  ai=[];
  hl=[];
  
  if (mod(length(xy),2) ~= 0)
    indic=-1;
    return;
  end  
  
  p=length(R);
  nn=length(xy)/2;
  nb=nn+1;
  x=xy(1:nn);
  y=xy(nn+1:end);
  
  
  
  if indic==1   %Tracé de la chaîne et plancher
   if lmi~=0
    if (p==1)
      X=[-1,A+1];
      Y=S*X+R;
      plot(X,Y,'-ro',[0;x;A],[0;y;B],'-bo');
    else
      r1=[0;R];
      r2=[R;0];
      s1=[S;0];
      s2=[0;S];
      num=r1-r2;
      den=s1-s2;
      inter=num./den     
      X=inter(2:length(inter)-1,1)' %les abscisses des points d'intersection des droites
      X=[-1 X A+1]
      coef=[S(1,1);S];
      ord=[R(1,1);R];
      Y=coef'.*X+ord';
      plot(X,Y,'-ro',[0;x;A],[0;y;B],'-bo');
    end  
   else
      plot([0;x;A],[0;y;B],'-bo');
   end
    %plot([0;x;A],[0;y;B],'-b');
    indic=1;
  end
  if (indic==2)||(indic==4)   %Calcul de e, ce et ci 
    Y1=[0 y']';
    X1=[0 x']';
    Y2=[y' B]';
    X2=[x' A]';
    V=Y2+Y1;
    ly=Y2-Y1;
    lx=X2-X1;
    li2=lx.^2+ly.^2;
    e=sum(L.*(V/2));
    ce=li2-L.^2; 
    
    
    rep1=ones(nn,1);
    rep2=ones(p,1);
    ci=kron(R,rep1)+kron(S,x)-kron(rep2,y);


    
    if indic==2
      indic=1;
      return;
    end
  end
  
  if (indic==4)   %Calcul e,ce,ci,g,ae et ai
    L1=L(1:nn)/2;
    L2=L(2:end)/2;
    gy=L1+L2;
    g=[zeros(1,nn) gy']';
    
    ae=zeros(nb,2*nn);
    for i=1:nb
      v=zeros(1,2*nn);
      if i==1
        v(1,1)=2*x(1,1);
        v(1,nn+1)=2*y(1,1);
        ae(1,:)=v;
      elseif i==nb
        v(1,nn)=-2*(A-x(nn,1));
        v(1,2*nn)=-2*(B-y(nn,1));
        ae(nb,:)=v;
      else
        v(1,i)=2*(x(i,1)-x(i-1,1));
        v(1,i-1)=-v(1,i);
        v(1,i+nn)=2*(y(i,1)-y(i-1,1));
        v(1,i+nn-1)=-v(1,i+nn);
        ae(i,:)=v;
      end  
    end
    
    m=[S -ones(p,1)];
    ai=kron(m,eye(nn));  
    indic=1;
  end


  if indic==5   %Calcul hl
    hl=zeros(2*nn,2*nn);
    for i=1:nb
      h=zeros(2*nn,2*nn);
      if i==1
        h(1,1)=2;
        h(nn+1,nn+1)=2;
      elseif i==nb
        h(nn,nn)=2;
        h(2*nn,2*nn)=2;
      else
        h(i,i)=2;
        h(i-1,i-1)=2;
        h(nn+i,nn+i)=2;
        h(nn+i-1,nn+i-1)=2;
        h(i,i-1)=2;
        h(i-1,i)=2;
        h(i+nn,nn+i-1)=2;
        h(nn+i-1,nn+i)=-2;
      end
      hl=hl+lme(i,1)*h;
    end
    indic=1;
  end


return

