function [xy,lme,lmi,info] = sq (xy,lme,lmi,options)
 
  small = 1.e-5;
  big   = 1.e+5;
  alpha=1;
  
  n=length(xy);
  maxiter=options.maxit;
  eps1=options.tol(1);    
  eps2=options.tol(2);
  eps3=options.tol(3);
  nbiter=0;
  
  
  s=size(lmi,1)
  if s==0
      ci=0;
      lmi=0;
  end
  
  
  
if options.deriv==2  %Méthode Newton 
 
  while 1
    if nbiter==0
     if s~=0
      [e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi);
      [e,~,~,~,~,~,hl,indic]=chs(5,xy,lme,lmi);
     else
       [e,ce,~,g,ae,~,~,indic]=chs(4,xy,lme,lmi);
       [e,~,~,~,~,~,hl,indic]=chs(5,xy,lme,lmi);
     end
    nbiter=nbiter+1;
    else
        
      
      supce=max(ce.*sign(ce));
      if s~=0
        gradl=g+ae'*lme+ai'*lmi;
      else
        gradl=g+ae'*lme;
      end
      
      supl=max(gradl.*sign(gradl));
    
      supci_lm=max(min(lmi,-ci));
    
     if nbiter == maxiter
        info.status=2;
        info.niter=maxiter;
        break;
     end
      if s~=0
          if (supl<=eps1) && (supce<=eps2) && (supci_lm<=eps3)
            info.status=0;
            info.niter=nbiter;
            break;
          end
      else
          if (supl<=eps1) && (supce<=eps2) 
             info.status=0;
             info.niter=nbiter;
             break;
          end
      end
     
      [U,D,flag]=cholmod(hl,small,big);
      M=U*diag(D)*U';
      %M=eye(size(hl));
      if s~=0
          [d, obj,~,~,lambda] = quadprog (M,g,ai,-ci,ae,-ce);
      else
          [d, obj,~,~,lambda] = quadprog (M,g,[],[],ae,-ce);
      end
      xy=xy+d;
      
      lme=lambda.eqlin;
      
      if s~=0
          lme=lambda.eqlin;
          lmi=lambda.ineqlin;
      else
          lme=lambda.eqlin;
      end
    
      if s~=0 
          [e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi);
          [e,~,~,~,~,~,hl,indic]=chs(5,xy,lme,lmi);
      else
          [e,ce,~,g,ae,~,~,indic]=chs(4,xy,lme,lmi);
          [e,~,~,~,~,~,hl,indic]=chs(5,xy,lme,lmi);
      end
      nbiter=nbiter+1;
      
    end
    
     end
    
end

  







if options.deriv==1 %Méthode quasi-Newton
  %Calcul des variables nécessaires 

  
  %M1:
  M= eye(n,n);
  
  while 1
      
      
      
   
   if nbiter==0
      [e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi);
      [e,~,~,~,~,~,hl,indic]=chs(5,xy,lme,lmi);
      nbiter=nbiter+1;
   else
      
      supce=max(ce.*sign(ce));
    if s~=0
      gradl=g+ae'*lme+ai'*lmi;

    else
      gradl=g+ae'*lme;

    end
    
      supl=max(gradl.*sign(gradl));
    
      supci_lm=max(min(lmi,-ci)); 
     %Si le nombre d'itérations maximal est atteint
     if nbiter == maxiter
        info.status=2;
        info.niter=maxiter;
        break;
     end
     
     %Si le point optimal est atteint
     if s~=0
        if (supl<=eps1) && (supce<=eps2) && (supci_lm<=eps3)
            info.status=0;
            info.niter=nbiter;
         break;
         end
     else
         if (supl<=eps1) && (supce<=eps2) 
            info.status=0;
            info.niter=nbiter;
         break;
         end
     end
      [d, obj,~,~,lambda] = quadprog (M,g,ai,-ci,ae,-ce);
      %Calcul de delta_k
      delta_k = d;
      
      %Mise à jour de lm
      lme=lme + alpha*(lambda.eqlin - lme);
      lmi=lmi + alpha*(lambda.ineqlin - lmi);
      
      %Calcul du gradient en x_k et lm_k+1
      [e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi);
       
      if s~=0
        grad_l_current=g+ae'*lme+ai'*lmi;

      else
        grad_l_current=g+ae'*lme;

      end
      %Calcul de x_k+1
      xy=xy+alpha*delta_k;   
      
      %Calcul du gradient en x_k+1 et lm_k+1
      [e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi);
      
      if s~=0
        grad_l_next=g+ae'*lme+ai'*lmi;

      else
        grad_l_next=g+ae'*lme;

      end
      
      %Calcul gamma_k_l
      gamma_k_l = grad_l_next - grad_l_current;
      
      %Calcul de theta
      theta=1;
      if (gamma_k_l'*delta_k < 0.2 * delta_k'*M*delta_k)
          theta=0.8*(delta_k'*M*delta_k)/(delta_k'*M*delta_k-gamma_k_l'*delta_k);
          
      end
      
      %Calcul de gamma_k
      gamma_k = (1-theta)*M*delta_k+theta*gamma_k_l;
      
      %Mise à jour de M1
      if nbiter==1
        eta=norm(gamma_k)^2/(gamma_k'*delta_k);
        M=eta*eye(n,n);
      end
      
      %Formule BFGS:
      M=M-(M*delta_k*delta_k'*M/(delta_k'*M*delta_k)) + (gamma_k*gamma_k')/(gamma_k'*delta_k);
      
      %Nouvel itéré
      nbiter=nbiter+1;
      
    end
    
   end
   
   %Conditionnement de la matrice M 
   condM=cond(M);
end
return

