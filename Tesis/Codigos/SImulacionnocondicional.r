# Script en R que a partir de la biblioteca gstat calcula los centroides de cada elemento finito y crea el conjunto de valores de simulacion no condicional y genera un conjunto de valores de carga hidraulica a partir de un variograma exponencial

library("gstat")
xy <- expand.grid(1:5, 1:5)

################ Busqueda de centroides ###########################
b=10
h=10
x=matrix(0,1,40)
y=matrix(0,1,20)
yy=matrix(0,1,400)
# Ciclo que obtiene la coordenada x de los centroides
contx=1
for(i in 1:20){
 for(j in 1:2){
  if (contx%%2==0) 
 { 
  x[contx]<-((2*b)/3)+b*(i-1)
 }
  else 
 {
  x[contx]<-(b/3)+b*(i-1)
 }
 contx=contx+1
}
}

# Ciclo que obtiene la coordenada y de los centroides
conty=1
for(m in 1:10){
 for(n in 1:2){
  if (conty%%2==0) 
 { 
  y[conty]<-((h)/3)+h*(m-1)
 }
  else 
 {
  y[conty]<-(2*h/3)+h*(m-1)
 }
 conty=conty+1
}
}

# Se crea el vector de coordenadas x por elemento
xx=c(x,x,x,x,x,x,x,x,x,x)
# Se crea el vector de coordenadas y por elemento
contyy=1
aux=1
for(i in 1:10){
 for(j in 1:40){
  if(contyy%%2==0){
   yy[contyy]<-y[aux+1]}
  else{
   yy[contyy]<-y[aux]}
  contyy=contyy+1}
 aux=aux+2}
yy=as.vector(yy)
# Se crea un data frame para los centroides
xy=data.frame(xx,yy)
names(xy) <- c("x","y")

##################Simulacion no condicional############################

# Se realizan las simulaciones no condicionales
vgm1=vgm(1,"Exp",15)
g.dummy <- gstat(formula = z~1, locations = ~x+y, dummy = TRUE, beta = 0,model = vgm1, nmax = 20)
yy <- predict(g.dummy, newdata = xy, nsim = 4)
write.csv(yy, file="Pruebayy.csv")


