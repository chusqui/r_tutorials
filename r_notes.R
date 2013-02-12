
###################### help R ########################

# funcion rm (borrar)
rm(X) 

# crear vectores (c)
c(x1,x2,x3)

# modo de objeto (mode)
mode(x)

# estadisticas basicas
mean(x) , median(x), sd(x), var(x), cor(x,y), cov(x,y)

# secuencias
seq(from=x,to=y,by=z)
rep(x,times=y)

# operadores de comparacion
 (==, !=, <, >, <=, >=)

# seleccion de indices en vectores
vector1[1:3]
vector1[c(x1,x2,x3)]

# Operadores
Indexing
Indexing > [ [[
Access variables in a name space > :: :::
Component extraction, slot extraction > $ @
Exponentiation (right to left) > ^
Unary minus and plus > - +
Sequence creation > :
Special operators > %any%
Multiplication, division > */
Addition, subtraction > +-
Comparison > == != < > <= >= 
Logical negation > !
Logical “and”, short-circuit “and” > & &&
Logical “or”, short-circuit “or” > | ||
Formula > ~
Rightward assignment > -> ->>
Assignment (right to left) > =
Assignment (right to left) > <- <<-
Help > ?

# definicion de funciones
function(param1, ...., paramN) expr

# working directory
setwd()
getwd()
chooseCRANmirror()

# salvar workspace
save.image()

# historia de comandos
history()
hstory(100)

# ayuda
help()
?nombreabuscar

# cierre de R
q()

# comentarios
signo "#"


# busqueda de recursos R
search()

# instalacion y carga de librerias
install.packages("nombre libreria")
library(nombre libreria)
installed.packages()
library()


# carga de paquetes default
data(dsname, package="pkgname")
data(package="nombre libreria")
data()

# ficheros con ancho fijo
records <- read.fwf("filename", widths=c(w1, w2, ..., wn))

# importacion de ficheros
scan()
read.table()
read.csv()
read.delim()

# clase de elemento
class(x)

# funcion stack() para crecomponer datos
Status   Age    V1     V2     V3    V4
           P 23646 45190  50333  55166 56271
          CC 26174 35535  38227  37911 41184
          CC 27723 25691  25712  26144 26398
          CC 27193 30949  29693  29754 30772
          CC 24370 50542  51966  54341 54273
          CC 28359 58591  58803  59435 61292
          CC 25136 45801  45389  47197 47126

zz <- read.csv("mr.csv", strip.white = TRUE)
zzz <- cbind(zz[gl(nrow(zz), 1, 4*nrow(zz)), 1:2], stack(zz[, 3:6]))

      Status   Age values ind
     X1         P 23646  45190  V1
     X2        CC 26174  35535  V1
     X3        CC 27723  25691  V1
     X4        CC 27193  30949  V1
     X5        CC 24370  50542  V1
     X6        CC 28359  58591  V1
     X7        CC 25136  45801  V1
     X11        P 23646  50333  V2

combinar 3 vectores con la funcion stack:
comb <- stack(list(v1=v1, v2=v2, v3=v3)) # Combine 3 vectors
combinar 3 listas:
comb <- stack(list(fresh=freshmen, soph=sophomores, jrs=juniors))

# exportacion de ficheros
write.csv(x,file="x.csv")
write.table(x,file="x.txt")

# lectura de informacion HTML

> library(XML)
> url <- "http://www.example.com/data/table.html"
> tbls <- readHTMLTable(url)
leer la tercera tabla de la pagina:
> tbl <- readHTMLTable(url, which=3)

# funcion para leer filas de ficheros
readLines()
readLines(con = stdin(), n = -1L, ok = TRUE, warn = TRUE, encoding = "unknown")

# lectura de datos de bases de datos MySQL:
library(RMySQL)
con <- dbConnect(MySQL(), user="userid", password="pswd", host="hostname", client.flag=CLIENT_MULTI_RESULTS)
sql <- "SELECT * from SurveyResults WHERE City = 'Chicago'"
rows <- dbGetQuery(con, sql)

dbConnect(MySQL(),user='alberto_quadri',password='',dbname='quadrigram',host='mysql.ekonlab.com')
dbListTables(mydb)
x <- dbGetQuery(mydb,"select * from Control_Module where example='yes'")
x <- dbGetQuery(mydb,"select * from Control_Module where example='yes' && library='resources'")

# loops

for (i in 1:n){
  do something
}

# busqueda de datos y segmentacion de data frames
childrenno <- a[a$children=="no", ]
childrenno <- a[a$children=="no" & a$yearsmarried==4, ]


# trabajo con vectores
v <- c(10, 20, 30)
names(v) <- c("Moe", "Larry", "Curly")
print(v)
  Moe Larry Curly 
   10    20    30 

funcion append:

append(vec, newvalues, after=n)
append(1:10, 99, after=5)
[1]  1  2  3  4  5 99  6  7  8  9 10

# matrices
A <- 1:6
print(A)
[1] 1 2 3 4 5 6

dim(A) <- c(2,3)
print(A)
     [,1] [,2] [,3]
[1,]    1    3    5
[2,]    2    4    6

operaciones con matrices:
transposicion
t(A)
inversa
solve(A)
multiplicacion
A %*% B
diagonalizar:
diag(n)

# arrays
D <- 1:12
dim(D) <- c(2,3,2)
print(D)
, , 1

     [,1] [,2] [,3]
[1,]    1    3    5
[2,]    2    4    6

, , 2

     [,1] [,2] [,3]
[1,]    7    9   11
[2,]    8   10   12

# funcion cbind
cbind(1:3)
     [,1]
[1,]    1
[2,]    2
[3,]    3
cbind(1:6, 1:3)
     [,1] [,2]
[1,]    1    1
[2,]    2    2
[3,]    3    3
[4,]    4    1
[5,]    5    2
[6,]    6    3

# convertir a factor
f <- factor(x$id)     donde: f es variable, x es dat frame y id es columna


# crear listas
lst <- list(x, y, z)
lst <- list(0.5, 0.841, 0.977)
lst <- list(mid=0.5, right=0.841, far.right=0.977)

# seleccionar elementos por posicion
seleccionar el elemento "n"
lst[[n]]
devuelve lista de elementos seleccionados por su posicion
lst[c(n1, n2, ..., nk)]

years <- list(1960, 1964, 1976, 1994)

years[[1]]
[1] 1960

years[c(1,2)]
[[1]]
[1] 1960
[[2]]
[1] 1964

# eliminar nulls:
lst[sapply(lst, is.null)] <- NULL

eliminar nulls condicionados:
lst[lst < 0] <- NULL
lst[is.na(lst)] <- NULL

#  nombres de filas y columnas
rownames(mat) <- c("rowname1", "rowname2", ..., "rownamem")
colnames(mat) <- c("colname1", "colname2", ..., "colnamen")

cambio de nombres en columnas:
colnames(dfrm) <- newnames

# seleccion de filas/columnas en una matrix
primera fila
vec <- mat[1,]
tercera columna
vec <- mat[,3]

seleccion de columnas por posicion:
la columna "n"
dfrm[[n]]
la columnas "na, n2..."
dfrm[c(n1, n2, ..., nk)]

seleccion con la funcion subset:
subset(dfrm, select=colname)
subset(dfrm, select=c(colname1, ..., colnameN))
subset(dfrm, select=c(predictor,response), subset=(response > 0))
subset(dfrm, select = -badboy)

# conversion a data frame
dfrm <- data.frame(v1, v2, v3, f1, f2)

conversion de datasets de una fila:
dfrm <- do.call(rbind, obs)

# introducir una fila en una columna:
newRow <- data.frame(city="West Dundee", county="Kane", state="IL", pop=5428)
suburbs <- rbind(suburbs, newRow)

# conversion a data frame
as.data.frame(x)

# eliminar na:
clean <- na.omit(dfrm)

# coeficiente de correlacion:
cor(x,y)
cor(subset(patient.data, select = c(-patient.id,-dosage)))

# combinacion de data frames:
all.cols <- cbind(dfrm1,dfrm2)
all.rows <- rbind(dfrm1,dfrm2)

# fusion de data frames con columnas en comun:
m <- merge(df1, df2, by="name")

# conversion de valores atomicos:
as.character(x)
as.complex(x)
as.numeric(x) or as.double(x)
as.integer(x)
as.logical(x)

as.data.frame(x)
as.list(x)
as.matrix(x)
as.vector(x)

# separacion de subconjuntos:
groups <- split(x, f)
groups <- unstack(data.frame(x,f))
split(x$example,x$issues_example)

# aplicar una funcion a cada elemento de una lista:
lst <- lapply(lst, fun)
vec <- sapply(lst, fun)
lapply(scores, length)
sapply(scores, length)

aplicar funcion a elementros de una fila:
results <- apply(mat, 1, fun)
apply(long, 1, mean)

aplicar funcion a elementos de una columna:
results <- apply(mat, 2, fun)
lst <- lapply(dfrm, fun)
vec <- sapply(dfrm, fun)

aplicar una funcion a grupos de datos:
tapply(x, f, fun)

# strings y dates
funcion nchar()
> nchar("Moe") [1] 3
> nchar("Curly") [1] 5

concatenacion:
paste("Everybody", "loves", "stats.") 
[1] "Everybody loves stats."
stooges <- c("Moe", "Larry", "Curly")

extraccion de substrings:
substr(string,start,end)
strsplit(string, delimiter)

reemplazo de substrings
sub(old, new, string)
gsub(old, new, string)

generacion de combinaciones de pares de valores:
m <- outer(strings1, strings2, paste, sep="")

locations <- c("NY", "LA", "CHI", "HOU") 
treatments <- c("T1", "T2", "T3")
outer(locations, treatments, paste, sep="-")
[,1] [,2] [,3]
[1,] "NY-T1" "NY-T2" "NY-T3" 
[2,] "LA-T1" "LA-T2" "LA-T3" 
[3,] "CHI-T1" "CHI-T2" "CHI-T3" 
[4,] "HOU-T1" "HOU-T2" "HOU-T3"

# conversiones:
string to date
as.Date("2010-12-31")
date to string
format(Sys.Date())
as.character(Sys.Date())
year, month and day to date
ISOdate(year, month, day)
as.Date(ISOdate(year, month, day))

# probabilidad
distribuciones 
dnorm Normal density
pnorm Normal distribution function 
qnorm Normal quantile function 
rnorm Normal random variates

Discrete distribution
Binomial                          binom
Geometric                         geom
Hypergeometric                    hyper
Negative binomial (NegBinomial)   nbinom
Poisson                           pois

Continuous distribution
Beta                              beta
Cauchy                            cauchy
Chi-squared (Chisquare)           chisq
Exponential                       exp
F                                 f
Gamma                             gamma
Log-normal (Lognormal)            lnorm
Logistic                          logis
Normal                            norm
Student´s t (Tdist)               t
Uniform                           unif
Weibull                           weibull
Wilcoxon                          wilcox


contar numero de combinaciones
choose(n, k)

generar combinaciones
combn(items, k)

generar numeros aleatorios
runif(1, min=-3, max=3)
rnorm(1, mean=100, sd=15
set.seed(165)
runif(10)

generar muestras aleatorias:
sample(vec, n)
sample(world.series$year, 10)

generar secuencias aleatorias:
sample(set, n, replace=TRUE)

sample(c("H","T"), 10, replace=TRUE)
[1] "H" "H" "H" "T" "T" "H" "T" "H" "H" "T"

# cuestiones generales

resumen de datos:
summary(x)

calculo de frecuencias relativas:
mean(x > 0)
mean(lab == “NJ”)
mean1 <- mean(af$education=="18")

calculo de cuantiles:
quantile(vec, f)
quantile(af$education)

test de la media de una muestra (t test)
t.test(x, mu=m)

intervalos de confianza:
t.test(x)
t.test(af$yearsmarried)

test de muestra:
prop.test(x, n, p)
prop.test(11, 20, 0.5, alternative="greater")

test de normalidad (shapiro)
shapiro.test(x)

# regresion lineal
lm(y ~ x)
lm(formula = y ~ x)
lm(y ~ x, data=dfrm)
lm(r~o,data=af)
regresion multiple
lm(y ~ u + v + w)
lm(y ~ u + v + w, data=dfrm)
modelos de regresion
m <- lm(y ~ u + v + w)
anova(m)
coefficients(m)
coef(m)
confint(m)
deviance(m)
effects(m)
fitted(m)
residuals(m)
resid(m)
vcov(m)

# generalidades
head(x)
tail(x)
rowSums(x)
colSums(x)
breaks <- c(-3,-2,-1,0,1,2,3)
match(80, vec)
seq_along(v) %% n == 0
pmin(1:5, 5:1)
pmax(1:5, 5:1)
str(x)
names(x)
unique(x)

# sorting
dfrm <- dfrm[order(dfrm$key),]
order(dfrm$month)
dfrm <- dfrm[order(dfrm$key1,dfrm$key2),]

# clusters:
d <- dist(x)
hc <- hclust(d)
clust <- cutree(hc, k=n)

# reshape
library(reshape)
sales<-melt(citysales)
sales$color[sales[,2]=="ProductA"] <- "red"
sales$color[sales[,2]=="ProductB"] <- "blue"
sales$color[sales[,2]=="ProductC"] <- "violet"
dotchart(sales[,3],labels=sales$City,groups=sales[,2],
   col=sales$color,pch=19,
   main="Sales Figures",
   xlab="Sales Revenue (1,000's of USD)")


############ Graficos ###############

plot(x)
plot(x, main="Forecast Results", xlab="Month", ylab="Production", + col=c("red", "black", "green"))
plot(x, y, xlim=c(lower,upper), ylim=c(lower,upper), ...)

plot = funcion generica de ploteado
boxplot = funcion de box plot
un boxplot por cada nivel de factor
boxplot(x ~ f)

hist = funcion de histograma
qqnorm = funcion de cuantil - cuantil plot
curve = grafico de funcion
points = añadir puntos
lines = añadir lineas
abline = añadir linea recta
segments = añadir linea de segmentos
polygon = añadir poligono
text = añadir texto

plot (x,y)
plot(cars,
+ main="cars: Speed vs. Stopping Distance (1920)", + xlab="Speed (MPH)",
+ ylab="Stopping Distance (ft)")
grid()
points(cars)

leyendas:
legend(x, y, labels, pch=c(pointtype1, pointtype2, ...))
legend(x, y, labels, lty=c(linetype1, linetype2, ...))
legend(x, y, labels, lwd=c(width1, width2, ...))
legend(x, y, labels, col=c(color1, color2, ...))

un scatter por cada nivel de factor
coplot(y ~ x | f)

barchart:
barplot(c(height1, height2, ..., heightn))
barplot(heights, col=colors)
barplot(c(3,5,4), col=c("red","white","blue"))
barplot(sales$ProductA,names.arg= sales$City,horiz=TRUE,col="black")

plot de lineas desde puntos x,y
plot(x, y, type="l")
plot(dfrm, type="l")
plot(x, y, type="l", lty="dashed")
plot(x, y.democr, type="l", col="blue")

histogramas:
hist(rnorm(1000))
hist(islands)

plot density:
plot(density(rnorm(1000)))

boxplot:
boxplot(metals,xlab="Metals",ylab="Atmospheric Concentration in ng per cubic metre",main="Atmospheric Metal Concentrations in London")
boxplot(copper$Cu~copper$Source,xlab="Measurement Site",ylab="Atmospheric Concentration of Copper in ng per cubic metre",main="Atmospheric Copper Concentrations in London")
boxplot(air,las=1)
boxplot(air,boxwex=0.2,las=1)
boxplot(metals[,-1],outline=FALSE)
boxplot(metals[,-1],horizontal=TRUE)



heatmaps:
heatmap(as.matrix(mtcars),Rowv=NA,Colv=NA,col = heat.colors(256),scale="column",margins=c(2,8),main = "Car characteristics by Model")
subset1 <- subset(af,select=c(age,yearmarried,education,occupation,rating))
heatmap(as.matrix(subset1),Rowv=NA,Colv=NA,col = heat.colors(256),scale="column",margins=c(2,8),main = "Car characteristics by Model")
contour(x=10*1:nrow(volcano), y=10*1:ncol(volcano), z=volcano)
filled.contour(x = 10*1:nrow(volcano),y = 10*1:ncol(volcano),z = volcano, color.palette = terrain.colors)



matriz de correlacion:
rownames(genes)<-colnames(genes)
image(x=1:ncol(genes), y=1:nrow(genes),z=t(as.matrix(genes)),axes=FALSE,xlab="",ylab="" , main="Gene Correlation Matrix")
axis(1,at=1:ncol(genes),labels=colnames(genes),col="white",las=2,cex.axis=0.8)
axis(2,at=1:nrow(genes),labels=rownames(genes),col="white",las=1,cex.axis=0.8)

pairs plots:
pairs(iris[,1:4])
pairs(subset1[,1:3])

matrix plots:
par(mfrow=c(2,3)) (filas y columnas)
plot(rnorm(100),col="blue",main="Plot No.1")
plot(rnorm(100),col="blue",main="Plot No.2")
plot(rnorm(100),col="green",main="Plot No.3")
plot(rnorm(100),col="black",main="Plot No.4")
plot(rnorm(100),col="green",main="Plot No.5")
plot(rnorm(100),col="orange",main="Plot No.6")

añadir leyendas:
plot(rain$Tokyo,type="l",col="red",
   ylim=c(0,300),
   main="Monthly Rainfall in major cities",
   xlab="Month of Year",
   ylab="Rainfall (mm)",
   lwd=2)
   lines(rain$NewYork,type="l",col="blue",lwd=2)
   lines(rain$London,type="l",col="green",lwd=2)
   lines(rain$Berlin,type="l",col="orange",lwd=2)
legend("topright",
   legend=c("Tokyo","NewYork","London","Berlin"),
   col=c("red","blue","green","orange"),
   lty=1,lwd=2)

mapas (library(maps))
map()
map('world', fill = TRUE,col=heat.colors(10))

library(sp)
load(url("http://gadm.org/data/rda/GBR_adm1.RData"))
spplot(gadm,"Shape_Area")

ajustes de parametros:

colores:
plot(rnorm(1000),col="red")
plot(sales$units~as.Date(sales$date,"%d/%m/%y"),
type="l", #Specify type of plot as l for line
col="blue")
heat.colors(5)
barplot(as.matrix(sales[,2:4]), beside=T,
   legend=sales$City,
   col=c("red","blue","green","orange","pink"),
   border="white")
barplot(as.matrix(sales[,2:4]), beside=T,
   legend=sales$City,
   col=heat.colors(length(sales$City)),
   border="white")
colores background:
par(bg="gray")
   plot(rnorm(100))
plot(rnorm(1000),type="n")
   x<-par("usr")
   rect(x[1],x[3],x[2],x[4],col="lightgray ")
   points(rnorm(1000))
otros elementos graficos:
plot(rnorm(100),
   main="Plot Title",
   col.axis="blue",
   col.lab="red",
   col.main="darkblue")

color brewer:
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(7,"YlOrRd")

plots con diferentes tamaños y estilos:
plot(rnorm(100),pch=19,cex=2)
par(oma=c(1,1,1,1))
plot(rnorm(100),bty="l")
box(which="figure")
plot(rnorm(100),xaxp=c(0,100,10))
plot(10^c(1:5),log="y",type="b")

scatter plots:
library(lattice)
library(ggplot2)
xyplot(mpg~disp,
   data=mtcars,
   groups=cyl,
   auto.key=list(corner=c(1,1)))
qplot(disp,mpg,data=mtcars,col= as.factor(cyl))
qplot(disp,mpg,data=mtcars,shape=as.factor(cyl))
qplot(disp,mpg,data=mtcars,size=as.factor(cyl))

mostrar textos de puntos dentro de un plot:
plot(mpg~disp, data=mtcars)
   text(258,22,"Hornet")

plot(health$Expenditure,health$Life_Expectancy,type="n")
text(health$Expenditure,health$Life_Expectancy,health$Country)

jitter:
plot(jitter(x), jitter(y))

linear model:
plot(mtcars$mpg~mtcars$disp)
lmfit<-lm(mtcars$mpg~mtcars$disp)
abline(lmfit)

x <- -(1:100)/10
y <- 100 + 10 * exp(x / 2) + rnorm(x)/10
nlmod <- nls(y ~  Const + A * exp(B * x), trace=TRUE)
plot(x,y)
lines(x, predict(nlmod), col="red")

scatter 3D:
library(scatterplot3d)

cuantil y cuartil plots:
qqnorm(mtcars$mpg)
qqline(mtcars$mpg)
lmfit<-lm(mtcars$mpg~mtcars$disp)
par(mfrow=c(2,2))
plot(lmfit)

rug:
x<-rnorm(1000)
plot(density(x))
rug(x)

metals<-read.csv("metals.csv")
plot(Ba~Cu,data=metals,xlim=c(0,100))
rug(metals$Cu)
rug(metals$Ba,side=2,col="red",ticksize=0.02)

smooth density:
n <- 10000
x  <- matrix(rnorm(n), ncol=2)
y  <- matrix(rnorm(n, mean=3, sd=1.5), ncol=2)
smoothScatter(x,y)


añadir marcadores en line charts:
rain <- read.csv("cityrain.csv")
plot(rain$Tokyo,type="b",lwd=2,
xaxt="n",ylim=c(0,300),col="black",
xlab="Month",ylab="Rainfall (mm)",
main="Monthly Rainfall in Tokyo")
axis(1,at=1:length(rain$Month),labels=rain$Month)
abline(v=9)
abline(h=150,col="red",lty=2)

sparklines:
rain <- read.csv("cityrain.csv")
par(mfrow=c(4,1),mar=c(5,7,4,2),omi=c(0.2,2,0.2,2))
for(i in 2:5)
   {
       plot(rain[,i],ann=FALSE,axes=FALSE,type="l",
       col="gray",lwd=2)
       mtext(side=2,at=mean(rain[,i]),names(rain[i]),
       las=2,col="black")
       mtext(side=4,at=mean(rain[,i]),mean(rain[i]),
       las=2,col="black")
       points(which.min(rain[,i]),min(rain[,i]),pch=19,col="blue")
       points(which.max(rain[,i]),max(rain[,i]),pch=19,col="red")
   }

line charts:
plot(sales$units~as.Date(sales$date,"%d/%m/%y"),type="l")
library(zoo)
plot(air$nox~as.Date(air$date,"%d/%m/%Y %H:%M"),type="l")


bar charts:
barplot(as.matrix(citysales[,2:4])
barplot(as.matrix(citysales[,2:4]), beside=TRUE,horiz=TRUE)

par(mar=c(5,4,4,8),xpd=T)
barplot(as.matrix(citysalesperc[,2:4]), horiz=TRUE,
   col=brewer.pal(5,"Set1"),border="white",
   xlab="Percentage of Sales",
   main="Perecentage Sales Figures")
legend("right",legend=citysalesperc$City,bty="n",
inset=c(-0.3,0),fill=brewer.pal(5,"Set1"))

barplot(as.matrix(citysales[,2:4]), beside=T,
   legend.text=citysales$City,args.legend=list(bty="n",horiz=T),
   ylim=c(0,100),ylab="Sales Revenue (1,000's of USD)",
   main="Sales Figures")

mapas
library(maps)
library(WDI)
library(RColorBrewer)
map("county", "new york")
map("state", region = c("california", "oregon", "nevada"))
map('italy', fill = TRUE, col = brewer.pal(7,"Set1"))

library(rgdal)
library(RgoogleMaps)
air<-read.csv("londonair.csv")
london<-GetMap(center=c(51.51,-0.116),zoom =10, destfile = "London.png",maptype = "mobile")
PlotOnStaticMap(london,lat = air$lat, lon = air$lon,cex=2,pch=19,col=as.character(air$color))
london<-GetMap(center=c(51.51,-0.116),zoom =13,destfile = "London_satellite.png",maptype = "satellite")
PlotOnStaticMap(london,lat = air$lat, lon = air$lon,cex=2,pch=19,col=as.character(air$color)

crear y leer KML:
library(rgdal)
cities <- readOGR(system.file("vectors",package = "rgdal")[1],"cities")
writeOGR(cities, "cities.kml", "cities", driver="KML")
df <- readOGR("cities.kml", "cities")

trabajar con ESRI:
library(maptools)
sfdata <- readShapeSpatial(system.file("shapes/sids.shp",package="maptools")[1], proj4string=CRS("+proj=longlat"))
plot(sfdata, col="orange", border="white", axes=TRUE)

shapefiles:
library(shapefiles)
sf<-system.file("shapes/sids.shp", package="maptools")[1]
sf<-substr(sf,1,nchar(sf)-4)
sfdata <- read.shapefile(sf)
write.shapefile(sfdata, "newsf")

exportacion de graficos:
png("cars.png",res=200,height=600,width=600)
   plot(cars$dist~cars$speed,
   main="Relationship between car distance and speed",
   xlab="Speed (miles per hour)",ylab="Distance travelled (miles)",
   xlim=c(0,30),ylim=c(0,140),
   xaxs="i",yaxs="i",col="red",pch=19)
dev.off()

png("cars.png",res=200,height=600,width=600)
   par(mar=c(4,4,3,1),omi=c(0.1,0.1,0.1,0.1),mgp=c(3,0.5,0),
   las=1,mex=0.5,cex.main=0.6,cex.lab=0.5,cex.axis=0.5)
   plot(cars$dist~cars$speed,
   main="Relationship between car distance and speed",
   xlab="Speed (miles per hour)",ylab="Distance travelled (miles)",
   xlim=c(0,30),ylim=c(0,140),
   xaxs="i",yaxs="i",
   col="red",pch=19,cex=0.5)
dev.off()

pdf("cars.pdf")
   plot(cars$dist~cars$speed,
   main="Relationship between car distance and speed",
   xlab="Speed (miles per hour)",ylab="Distance travelled (miles)",
   xlim=c(0,30),ylim=c(0,140),
   xaxs="i",yaxs="i",
   col="red",pch=19,cex=0.5)
dev.off()

svg("3067_10_03.svg")
   #plot command here
   dev.off()
postscript("3067_10_03.ps")
   #plot command here
   dev.off()

pdf("multiple.pdf")
   for(i in 1:3)
     plot(cars,pch=19,col=i)
dev.off()

 pdf("multiple.pdf",colormodel="cmyk")
   for(i in 1:3)
     plot(cars,pch=19,col=i)
dev.off()


############## ggplot2 #######################

# geoms

geomabline:
p <- qplot(wt, mpg, data = mtcars)
p + geom_abline()
p + geom_abline(intercept = 20)

geombar:
qplot(factor(cyl), data=mtcars, geom="bar")

geom_bin2d:
d <- ggplot(diamonds, aes(x = x, y = y)) + xlim(4,10) + ylim(4,10)
d + geom_bin2d()

geom_boxplot:
qplot(factor(cyl), mpg, data = mtcars, geom = "boxplot")

geom_dotplot:
ggplot(mtcars, aes(x = mpg)) + geom_dotplot()

geom_histogram:
qplot(rating, data=movies, geom="histogram")

geom_hline:
p <- ggplot(mtcars, aes(x = wt, y=mpg)) + geom_point()
p + geom_hline(aes(yintercept=mpg))

geom_jitter:
qplot(displ, hwy, data = mpg, geom = "jitter")

geom_line:
qplot(year, number, data=mry, group=rating, geom="line")

geom_rug:
p <- ggplot(mtcars, aes(x=wt, y=mpg))
p + geom_point() + geom_rug()

# statistics

stat_bin:
m <- ggplot(movies, aes(x=rating))
m + stat_bin()

stat_unique:
ggplot(mtcars, aes(vs, am)) + geom_point(alpha = 0.1, stat="unique")

# scales

guide_legend:
p1 + scale_fill_continuous(guide = "legend")

scale_alpha:
p <- qplot(mpg, cyl, data = mtcars, alpha = cyl))

scale_colour_brewer:
d <- qplot(carat, price, data=dsamp, colour=clarity))
d + scale_colour_brewer()
d + scale_colour_brewer("clarity")

scale_colour_gradient:
miss <- sample(c(NA, 1:5), nrow(mtcars), rep = TRUE)
qplot(mpg, wt, data = mtcars, colour = miss)

scale_colour_grey:
p <- qplot(mpg, wt, data=mtcars, colour=factor(cyl))
p + scale_colour_grey()

scale_colour_hue:
miss <- factor(sample(c(NA, 1:5), nrow(mtcars), rep = TRUE))
qplot(mpg, wt, data = mtcars, colour = miss)

scale_linetype:
qplot(date, value, data=ecm, geom="line", linetype=variable)

xlim / ylim:
qplot(mpg, wt, data=mtcars) + xlim(15, 20)
qplot(mpg, wt, data=mtcars) + ylim(15, 20)

# facet

facet_grid:
qplot(mpg, wt, data=mtcars, facets = . ~ vs + am)
qplot(mpg, wt, data=mtcars) + facet_grid(cyl ~ vs)
qplot(cty, model, data=mpg) + facet_grid(manufacturer ~ ., scales = "free", space = "free")

facet_wrap:
p <- qplot(price, data = diamonds, geom = "histogram", binwidth = 1000
p + facet_wrap(~ color)

position_jitter:
qplot(am, vs, data = mtcars, position = "jitter")

theme:
p + theme(panel.background = element_rect(colour = "pink"))
p + theme_bw()
p + xlab("Vehicle Weight") + ylab("Miles per Gallon")
m + theme(axis.text = element_text(colour = "blue"))


############### Notas del 25/9/2012 ##############


built in functions in R
http://www.statmethods.net/management/functions.html

order
deuda <- deuda2[order(-deudanum), ]
el negativo representa orden descendente


# notas del 16/10/2012
Charla sobre adquiring data with R
http://www.slideshare.net/jeffreybreen/tapping-the-data-deluge-with-r
https://github.com/jeffreybreen/talk-201210-data-deluge


# lectura de csv en web
ej: http://ichart.finance.yahoo.com/table.csv?s=S&d=9&e=16&f=2012&g=d&a=10&b=8&c=1984&ignore=.csv
url = "http://ichart.finance.yahoo.com/table.csv?s=S&d=9&e=16&f=2012&g=d&a=10&b=8&c=1984&ignore=.csv"
data = read.csv(url)



######## Notas de 17/10/2012

Linear Regression with R
http://msenux.redwoods.edu/math/R/regression.php

age=18:29
height=c(76.1,77,78.1,78.2,78.8,79.7,79.9,81.1,81.2,81.8,82.8,83.5)
plot(age,height)
res=lm(height~age)
height = 0.635 age + 64.928
res
Call:
lm(formula = height ~ age)

Coefficients:
(Intercept)          age  
     64.928        0.635  

abline(res)

# Extraer data frames en funcion de elementos duplicados de una de sus columnas

a <- duplicated(impure$email_address)
b <- impure[!a, ]


# función de suma acumulada cumsun()

# escribir un csv:
diamonds <- write.csv(file="diamonds.csv")



