%%%%%%%%%%%% Variables globales %%%%%%%%%%%%
global A B L R S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHOIX DU TP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%TP= 4 ou bien TP=5 %%%%%
TP=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHOIX DU CAS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CAS TESTS TP 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Choix Cas test%%%%%
%%%%choix=1:Cas test 4.a
%%%%choix=2:Cas test 4.b
%%%%choix=3:Cas test 4.c
choix=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CAS TESTS TP 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Choix Cas test%%%%%
%%%%choice=1:Cas test 5.a
%%%%choice=2:Cas test 5.b
%%%%choice=3:Cas test 5.c
%%%%choice=4:Cas test 5.d
choice=3;



v=[0.00001 0.00001 0.00001];
options.tol=v;
options.maxit=1000; 
options.deriv=2;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TP 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if TP==4
    
if choix==1
    L = [0.7 0.5 0.3 0.2 0.5]';
    A= 1;
    B= -1;
    xy = [ 0.2 0.4 0.6 0.8 ...
           1 1.5 1.5 1.3]';
    R=[];
    S=[];
    lme=0.1*ones(length(L),1);
    [e,ce,~,g,ae,~,~,indic]=chs(4,xy,lme,[]);
    lme=-(ae*ae')\ae*g;
    %lmi=zeros(2,1);
    [e,ce,~,g,ae,~,~,indic]=chs(4,xy,lme,[]);
    [xy,lme,~,info] = sq (xy,lme,[],options);         %Optimiseur
    [~,~,~,~,~,~,~,~] = chs(1,xy,lme,0)               %traçage chaîne

else
       A= 1;
       B= 0;
       
       L = [0.2 0.2 0.2 0.3 0.3 0.5 0.2 0.2 0.3 0.1]';
       xy = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
            -0.5 -0.9 -1.2 -1.4 -1.5 -1.4 -1.2 -0.9 -0.5]';
       if choix==2
           R = -0.25 ;
           S = -0.5  ;
       
          
       else
           R = [-0.25; -0.5];
           S = [-0.5; 0];
  
       
       end;
       lme=ones(length(L),1);
       lmi=ones(length(R)*length(xy)/2,1);
       [e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi);
       Amc=[ae;ai];   
       lm=-(Amc*Amc')\Amc*g;
       lme=lm(1:size(ce,1));
       lmi=lm((size(ce,1)+1):(size(ci,1)+size(ce,1)));
       %[e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi); 
       [e,~,~,~,~,~,hl,indic]=chs(5,xy,lme,lmi);
       [xy,lme,lmi,info] = sq (xy,lme,lmi,options);       %Optimiseur
       [~,~,~,~,~,~,~,~] = chs (1,xy,lme,lmi);             %traçage chaîne
end;





else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TP 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if choice==1
    L = [0.5 0.3 0.4 1.2 0.3 0.3]';
    A = 0;
    B = 0;
    xy = [ 0.2 0.5 0.8 1.0 1.2 ...
           -0.4 -0.6 -0.4 -0.2 0.]';
    R = [-1]';
    S = [-0.1]';
    lme=ones(length(L),1);
    lmi=ones(length(R)*length(xy)/2,1);
    [e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi);
    Amc=[ae;ai];   
    lm=-(Amc*Amc')\Amc*g;
    lme=lm(1:size(ce,1));
    lmi=lm((size(ce,1)+1):(size(ci,1)+size(ce,1)));
    [e,~,~,~,~,~,hl,indic]=chs(5,xy,lme,lmi);
    [xy,lme,lmi,info] = sq (xy,lme,lmi,options);       %Optimiseur
    [~,~,~,~,~,~,~,~] = chs (1,xy,lme,lmi);
end

if choice==2
    L = [3 2.5 2.5]';
    A = 0;
    B = -4;
    xy = [ -2 0 ...
           1 -2]';
    R = [-6 -10]';
    S = [-2 100]';
    lme=ones(length(L),1);
    lmi=ones(length(R)*length(xy)/2,1);
    [e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi);
    Amc=[ae;ai];   
    lm=-(Amc*Amc')\Amc*g;
    lme=lm(1:size(ce,1));
    lmi=lm((size(ce,1)+1):(size(ci,1)+size(ce,1)));
    [e,~,~,~,~,~,hl,indic]=chs(5,xy,lme,lmi);
    [xy,lme,lmi,info] = sq (xy,lme,lmi,options);       %Optimiseur
    [~,~,~,~,~,~,~,~] = chs (1,xy,lme,lmi);
end


if choice==3
    L = [0.1 0.2 0.3 0.4 0.5 0.4 0.3 0.1]';
    xy= [0.1 0.15 0.2 0.25 0.2 0.15 0.1 ...
         0.3 0.5 0.7 0.5 0.4 0.3 0.2]';
    A = 0;
    B = 0;
    R = [-1.0 -0.2 -1]';
    S = [-7.0  0.0  7]';
    lme=ones(length(L),1);
    lmi=ones(length(R)*length(xy)/2,1);
    [e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi);
    %Amc=[ae;ai];   
    %lm=-(Amc*Amc')\Amc*g;
    %lme=lm(1:size(ce,1));
    %lmi=lm((size(ce,1)+1):(size(ci,1)+size(ce,1)));
    [e,~,~,~,~,~,hl,indic]=chs(5,xy,lme,lmi);
    [xy,lme,lmi,info] = sq (xy,lme,lmi,options);       %Optimiseur
    [~,~,~,~,~,~,~,~] = chs (1,xy,lme,lmi);
end

if choice==4
    A=1;
    B=-1;
    L = [1.5 1.7 0.8 2.2 1]' ;
    xy=[0.2 4 0.6 0.8 ...
        -1 -1.5 -1.5 -1.3]';
    R = [-0.8; -1.5] ;
    S =[-1.5; -0.5] ;
    lme=ones(length(L),1);
    lmi=ones(length(R)*length(xy)/2,1);
    [e,ce,ci,g,ae,ai,~,indic]=chs(4,xy,lme,lmi);
    %Amc=[ae;ai];   
    %lm=-(Amc*Amc')\Amc*g;
    %lme=lm(1:size(ce,1));
    %lmi=lm((size(ce,1)+1):(size(ci,1)+size(ce,1)));
    [e,~,~,~,~,~,hl,indic]=chs(5,xy,lme,lmi);
    [xy,lme,lmi,info] = sq (xy,lme,lmi,options);       %Optimiseur
    [~,~,~,~,~,~,~,~] = chs (1,xy,lme,lmi);
end

end






