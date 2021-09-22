clear all
clc
q=dlmread('Solucionflujo.txt');
coo=dlmread('Coordenadas.txt');
n=length(q)
cont1=1
cont2=1
for i=1:n
  if rem(i,2)==0 
    u(cont1)=q(i);
    cont1=cont1+1;
  else 
    v(cont2)=q(i);
    cont2=cont2+1;
  end
end
 