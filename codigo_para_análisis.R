# Cargamos las librerías necesarias
library(SummarizedExperiment)
library(readr)

# Descargamos el *dataset* del repositorio de *github*.
url <- "https://raw.githubusercontent.com/nutrimetabolomics/metaboData/main/Datasets/2023-UGrX-4MetaboAnalystTutorial/ST000002_AN000002.txt"
download.file(url, destfile = "dataset.txt")

# Leemos el archivo completo
lines <- read_lines("dataset.txt")

# Separamos la información del estudio
study_info <- lines[1:70]

# Y los valores del estudio
data_values <- read_delim("dataset.txt", skip = 71, n_max = 143, delim = "\t", 
                          show_col_types = FALSE)

#Obtenemos el nombre de las columnas
colnames(data_values,)

# Separamos las últimas líneas que tienen información de las variables (metabolitos)
info_metab <- suppressWarnings(read_delim("dataset.txt", skip = 218, n_max = 142,
                                          delim = "\t", show_col_types = FALSE))

# Transformamos los datos en una matriz
assay_data_values <- data_values
assay_data_values[] <- suppressWarnings(lapply(assay_data_values, function(x) as.numeric(as.character(x))))
assay_data <- as.matrix(assay_data_values[, -1])

#Asignamos como nombre de las filas el nombre de cada metabolito
rownames(assay_data) <- data_values[[1]]

# Definimos los nombres de las columnas
colnames(assay_data) <- paste0(ifelse(1:12 <= 6, "A", "B"),colnames(assay_data))

# Suprimimos la fila correspondiente al factor
assay_data <- assay_data[-1, ]

#Creamos el objeto col_data
col_data <- DataFrame(sampleID = colnames(assay_data))

# Creamos el objeto se de clase SummarizedExperiment
se <- SummarizedExperiment(assays = list(counts = assay_data), colData = col_data)

# Y le añadimos los metadatos
metadata(se) <- list(study_info = study_info)

#Le asignamos la información de los metabolitos
rowData(se) <- DataFrame(info_metab)

# Examinamos el objeto se
se

# Nombre de las columnas de se (muestras)
colnames(se)

# Nombre de las filas de se (metabolitos)
rownames(se)

# Resumen básico de los datos
apply(assay_data, 2, summary)

# Histograma para cada muestra con la distribución de las concentraciones de los metabolitos
opt <- par(mfrow = c(3, 4))
for (i in 1:ncol(assay_data)) {
  hist(assay_data[, i], main = colnames(assay_data)[i], xlim = c(0, 1e4), breaks = 1e4)
}
par(opt)

# Boxplot
groupColors <- c(rep("red", 6), rep("blue", 6))
boxplot(assay_data, col=groupColors, main="Niveles de metabolitos en 12 muestras, dos grupos",
        xlab="Muestras",
        ylab="Concentración de metabolitos", las=2, cex.axis=0.7, cex.main=0.7)

# Corrección logarítmica de los datos
log_assay_data <- log1p(assay_data)

# Boxplot con los datos transformados
groupColors <- c(rep("red", 6), rep("blue", 6))
boxplot(log_assay_data, col=groupColors, main="Niveles de metabolitos en 12 muestras, dos grupos",
        xlab="Muestras",
        ylab="Concentración de metabolitos", las=2, cex.axis=0.7, cex.main=0.7)

# Cálculo de los componentes principales
pc<-prcomp(t(log_assay_data), scale=FALSE)
loads<- round(pc$sdev^2/sum(pc$sdev^2)*100,1)

# Representación gráfica de las dos primeras componentes
xlab<-c(paste("PC1",loads[1],"%"))
ylab<-c(paste("PC2",loads[2],"%"))
plot(pc$x[,1:2],xlab=xlab,ylab=ylab, col=groupColors, 
     main ="Componentes principales")
names2plot<-paste0(rep(c("A", "B"), each = 6), 1:6)
text(pc$x[,1],pc$x[,2],names2plot, pos=3, cex=.6)

# Clúster jerárquico
colnames(log_assay_data) <- names2plot
clust.euclid.average <- hclust(dist(t(log_assay_data)),method="average")
plot(clust.euclid.average, hang=-1)

# Definición del test t
ttest=function(x){
  tt = t.test(x[1:6], x[7:12])
  return(c(tt$statistic, tt$p.value, tt$estimate[1] - tt$estimate[2]))
}

# Aplicación del t-test a nuestra matriz
ans <- apply(log_assay_data, 1, ttest)
ts <- ans[1,]
pvals <- ans[2,]
fc <- ans[3,]

# Número de metabolitos con un p-valor significativo
for (i in c(0.05, 0.01))
  print(paste("El número de metabolitos con un p-valor inferior a", i, " es de", length(which(pvals < i))))

# Nombre de los metabolitos con un p-valor significativo
metabolitos_significativos <- rownames(assay_data)[which(pvals < 0.05)]
metabolitos_significativos

# Nombre de los metabolitos con un p-valor < 0.01
metabolitos_muy_significativos <- rownames(assay_data)[which(pvals < 0.01)]
metabolitos_muy_significativos

# Creación del objeto contenedor en formato .Rda
save(se, file = "objeto_se.Rda")