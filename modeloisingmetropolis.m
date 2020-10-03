clear all
L=200; %Lado de la red 
K=500; %Numero de iteraciones
J=1;
E(1)=0;
spin=ceil(rand(L,L)*1.5)*2-3; % red de espines
Tc=2/log(1+sqrt(2)); %Temperatura critica
T=(0.1:0.1:1.9)'*Tc; %Temperaturas alrededor de Tc
M(size(T))=0;
for t=1:size(T,1) %iteracion para las temperaturas
    p=[1 1 1 exp(-4*J/T(t)) exp(-8*J/T(t))]; %Probabilidad 
    s(1:K,1)=0;
    for k=1:K %iteraciones de metropolis
    r0=ceil(rand(L^2,2)*L); rn=mod(r0-2,L)+1; rp=mod(r0,L)+1;
    r=rand(L^2,1);
    for n=1:L^2
         Ener(n)=-log(p(((spin(rn(n,1),r0(n,2)) + spin(rp(n,1),r0(n,2)) + spin(r0(n,1),rn(n,2)) + spin(r0(n,1),rp(n,2)))*spin(r0(n,1),r0(n,2))/2 + 3)))*T(t);
        
        if (r(n) < p(((spin(rn(n,1),r0(n,2)) + spin(rp(n,1),r0(n,2)) + spin(r0(n,1),rn(n,2)) + spin(r0(n,1),rp(n,2)))*spin(r0(n,1),r0(n,2))/2 + 3)))
            spin(r0(n,1),r0(n,2))=-spin(r0(n,1),r0(n,2));
        end  %algoritmo de metropolis
    end
    s(k)=sum(sum(spin))/L^2;
    Enertot=sum(Ener)/L^2;
    end
    Ener=log(p(((spin(rn(n,1),r0(n,2)) + spin(rp(n,1),r0(n,2)) + spin(r0(n,1),rn(n,2)) + spin(r0(n,1),rp(n,2)))*spin(r0(n,1),r0(n,2))/2 + 3)));
    E(t)=Enertot; %Energia
    M(t)=sum(s(K-99:K))/100; %magnetizacion
    imagesc(spin);
    axis equal off
    drawnow;
    [t M(t)]
    [t E(t)]
end
figure;
plot(T/Tc,E,'*','DisplayName','Simulacion')
ylabel('Energia <|E|>');
xlabel('T/Tc');
grid on
figure
plot(T/Tc,abs(M),'*','DisplayName','Simulacion');
ylabel('Magnetización |M|');
xlabel('T/Tc');
hold on;
Tf=(0.5:0.001:1)*Tc';
plot(Tf/Tc,real((1-sinh(2*J*Tf.^-1).^-4).^(1/8)),'DisplayName','Teoria');
grid on