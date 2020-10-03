clear
L=200; %Lado de la red 
Q=300; %Numero de iteraciones
J=1; 
kbT=0.1; %constante de boltzman =1 *T
spin=ceil(rand(L,L)*2)*2-3; %red de espines
p=[1 1 1 exp(-4*J/kbT) exp(-8*J/kbT)]; %probabilidades
for q=1:Q
    s(1:Q,1)=0;
    r0=ceil(rand(L^2,2)*L); rn=mod(r0-2,L)+1; rp=mod(r0,L)+1;
    r=rand(L^2,1);
    for n=1:L^2
        Ener(n)=-log(p(((spin(rn(n,1),r0(n,2)) + spin(rp(n,1),r0(n,2)) + spin(r0(n,1),rn(n,2)) + spin(r0(n,1),rp(n,2)))*spin(r0(n,1),r0(n,2))/2 + 3)))*kbT;
        
        if (r(n) < p(((spin(rn(n,1),r0(n,2)) + spin(rp(n,1),r0(n,2)) + spin(r0(n,1),rn(n,2)) + spin(r0(n,1),rp(n,2)))*spin(r0(n,1),r0(n,2))/2 + 3)))
            spin(r0(n,1),r0(n,2))=-spin(r0(n,1),r0(n,2));
        end  %algoritmo metropolis
    end
    Enertot(q)=sum(Ener)/L^2;
     s(q)=sum(sum(spin))/L^2;
    %Ener(n)=log(p(((spin(rn(n,1),r0(n,2)) + spin(rp(n,1),r0(n,2)) + spin(r0(n,1),rn(n,2)) + spin(r0(n,1),rp(n,2)))*spin(r0(n,1),r0(n,2))/2 + 3)))*kbT;
    imagesc(spin);
    axis equal off
    drawnow;
end
M=sum(s(Q-99:Q))/100;
figure
plot(Enertot)
ylabel('Energia <|E|>');
xlabel('numero de iteraciones ');
grid on
%figure
%plot(M)
    