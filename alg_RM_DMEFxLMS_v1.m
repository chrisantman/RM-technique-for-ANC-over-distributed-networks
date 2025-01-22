%% RM_DMEFxLMS ALGORTIHM v1: 
%  Versión donde calculamos las señales d_v de manera distribuida incremental 
%  Igual a la versión centralizada

function out =alg_RM_DMEFxLMS_v1(conf, IN, OUT)

%% IN
I         = length(conf.ref);   
J         = length(conf.alt);   
L         = conf.L;
ITE       = conf.ITE;
mu_M      = conf.mu_RM;
num_nodos = conf.num_nodos;
CSEC_m    = IN.CSEC_m;
CSEC_v    = IN.CSEC_v;
CPRI_m    = IN.CPRI_m;
csec_v    = IN.csec_v;
in        = IN.x;
O         = OUT.mopt;
%% 

[M, Km]=size(CPRI_m);
[Kv, O_size]=size(O);
P=O_size/Km;

%%  Cálculo de señal deseada
d=zeros(ITE,Km);
for k=1:Km
    d(:,k)=filter(CPRI_m(:,k),1,in);
end

%% Inicialización de variables
w=zeros(L*J,num_nodos);	
w_eval=zeros(L,num_nodos,ITE);	
buff_L=zeros(L,num_nodos);
error=zeros(ITE,num_nodos);
deseada=zeros(ITE,num_nodos); 
e=zeros(1,num_nodos);	
buff_xf=zeros(M,num_nodos);
v=zeros(L*J,num_nodos); 
buff_y=zeros(M,J);
y_hist=zeros(ITE,num_nodos);


yf=zeros(1,Km);
e_m=zeros(1,Km);
yf_v=zeros(1,Kv);
e_v=zeros(1,Kv);
d_v=zeros(1,Kv);
buff_y_nodo=zeros(M,num_nodos,num_nodos);
buff_dm=zeros(O_size,1); 

if Km == 4
    Km_j=ones(1,num_nodos);
elseif Km == 8
    Km_j=2*ones(1,num_nodos);
end

Kv_j=ones(1,num_nodos);
O_size_nodo=P*Km_j;

buff_dm_nodo=cell(num_nodos,1);
for k=1:num_nodos
    buff_dm_nodo{k}=zeros(O_size_nodo(k),1);
end

for cont=1:ITE

   y=zeros(1,J);  
   x=in(cont,I);
   D_nodo=zeros(Kv,num_nodos);
     
   for k=1:num_nodos
       
        %% ESTRATEGIA INCREMENTAL .........................................
        if k==1
            D_nodo_ant=D_nodo(num_nodos,:);
        else
            D_nodo_ant=D_nodo(k-1,:);
        end
         
        % Filtrado remote technique
        Kaux=sum(Km_j(1:k))-Km_j(1);
        d_m_nodo=e_m(1+Kaux:Km_j(k)+Kaux)-yf(1+Kaux:Km_j(k)+Kaux);          
        for km=1:Km_j(k)
              buff_dm_nodo{k}(1+P*(km-1):km*P,:)=[d_m_nodo(km);   buff_dm_nodo{k}(1+P*(km-1):km*P-1,:)];
        end
        O_nodo=O(:,1+O_size_nodo(k)*(k-1):k*O_size_nodo(k));
        D_nodo(k,:)=D_nodo_ant+(O_nodo*buff_dm_nodo{k})';
   end  
   d_v_nodo=D_nodo(end,:);         

  for nodo=1:num_nodos
      
        y_est_nodo=zeros(num_nodos,1);
        for k=1:num_nodos
            y_est_nodo(k)=w(1+L*(k-1):L*k,nodo)'*buff_L(:,nodo); 
        end
        buff_y_nodo(:,:,nodo)=[y_est_nodo'; buff_y_nodo(1:M-1,:,nodo)]; 
        yf_v_nodo=0;
        for k=1:Kv
             yf_v_nodo=yf_v_nodo+CSEC_v(:,nodo+num_nodos*(k-1))'*buff_y_nodo(:,k,nodo);
        end  
        e_v_nodo=d_v_nodo(nodo)+yf_v_nodo;
        e_v(nodo)=e_v_nodo;
        d_v(nodo)=d_v_nodo(nodo);
        yf_v(nodo)= yf_v_nodo;

        
        %% ESTRATEGIA INCREMENTAL .........................................
        w_nodo=w(1+L*(nodo-1):L*nodo,nodo);     
        if nodo==1
            w_ant=w(:,num_nodos);
        else
            w_ant=w(:,nodo-1);
        end
        
        %% FILTRADO ADAPTATIVO ............................................
        buff_L(:,nodo)=[x; buff_L(1:L-1,nodo)]; 
        y(nodo)=y(nodo)+w_nodo'*buff_L(:,nodo); 

        %% ACTUALIZACIÓN COEFICIENTES .....................................       
        w(:,nodo)=w_ant-mu_M(nodo)*v(:,nodo)*e_v_nodo;
        w_eval(:,nodo,cont)=w(1+L*(nodo-1):L*nodo,nodo);          

        %% FILTRADO ESTIMA ................................................
        buff_xf(:,nodo)=[x; buff_xf(1:M-1,nodo)]; 
        vijk=buff_xf(:,nodo)'*csec_v(:,:,nodo);
        vaux=reshape(v(:,nodo),L,Kv);
        vaux(2:L,:)=vaux(1:L-1,:);
        vaux(1,:)=vijk; 
        v(:,nodo)=vaux(:);

   end
   
    % Difusión de coeficientes
    w=w(:,end)*ones(1,num_nodos);
 


    
    %% PARTE ACÚSTICA .........................................................
    buff_y=[y; buff_y(1:M-1,:)];
    yf=zeros(1,Km);
    for k=1:Km
        for j=1:J
          yf(k)=yf(k)+CSEC_m(:,k+Km*(j-1))'*buff_y(:,j);
        end
    end
    e_m=d(cont,:)+yf;     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   


    error(cont,:)=e_v;
    y_hist(cont,:)=y;
    deseada(cont,:)=d_v;

      
end


%% OUT [error_est1,d_v_est1]
 %        [error  deseada] 
out.altavoces  = y_hist;
out.error_mv   = error;
out.deseada_mv = deseada;


