ANALISIS CON PRS
================
Lea Vilanova Mañá
January 18, 2023

``` r
rm(list = ls())
```

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/v9/9nmtx3b97392bfcshr7srp1c0000gn/T//RtmpuIhRu4/downloaded_packages

Este estudio de control de calidad y cálculo de PRS se ha realizado a
partir del código presente en el tutorial presente en este enlace:

<https://choishingwan.github.io/PRS-Tutorial/>

Dicho tutorial está asociado al articulo “P.F. Tutorial: a guide to
performing polygenic risk score analyses.” (Choi, S.W., Mak, T.S. &
O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score
analyses. Nat Protoc (2020).
<https://doi.org/10.1038/s41596-020-0353-1>)

Como herramientas se han utilizado dos softwares:

-RStudio (RStudio Team (2020). RStudio: Integrated Development for R.
RStudio, PBC, Boston, MA URL <http://www.rstudio.com/>.)

-Plink 1.9 (Shaun Purcell, Christopher Chang, URL:
www.cog-genomics.org/plink/1.9/) (Chang CC, Chow CC, Tellier LCAM,
Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising
to the challenge of larger and richer datasets. GigaScience, 4.)

El código ha sido utilizado y adaptado según conveniencia, siendo los
autores los mencionados anteriormente.

# 1 ORGANIZACION DE LOS DIRECTORIOS

Para llevar a cabo este estudio, en primer lugar se han creado algunos
directorios para, por un lado, tener los ficheros organizados y, por
otro lado, tener acceso a los softwares y sus dependencias, ya que estos
solo pueden funcionar si se incluyen en el directorio de trabajo.

Para ello se ha trabajado des del terminal con las instrucciones
siguientes:

- **mkdir tfm**

- **cd tfm**

- **mkdir prs**

Una vez creados, definimos el directorio de trabajo para tener acceso a
los ficheros. En primer lugar miramos cuál es el directorio de trabajo
por defecto de este documento.

``` r
getwd()
```

    ## [1] "/Users/leavilanovamana/tfm/prs"

En nuestro caso hemos predefinido nuestro directorio de trabajo en el
chunk de configuración con la instrucción
**knitr::opts_knit\$set(root.dir = “/Users/leavilanovamana/tfm/prs”)**,
de este modo no debemos preocuparnos en todo el análisis.

# 2 DESCOMPRESION Y ACCESO A LOS ARCHIVOS

Los archivos fueron descargados directamente des de la web,
<https://choishingwan.github.io/PRS-Tutorial/base/>, y almacenados en el
directorio “prs” creado anteriormente. Descomprimimos los archivos y
miramos que encontramos en las carpetas.

Utilizamos las siguientes instrucciones des del terminal:

- **unzip EUR.zip**

- **gunzip Height.gwas.txt.gz**

# 3 CONTROL DE CALIDAD

Antes de calcular los PRS, debemos hacer el QC de los datos, tanto base
como target.

## 3.1 CONTROL DE CALIDAD DE LOS DATOS BASE

Los datos base necesarios para calcular PRS son los obtenidos en un
estudio GWAS (resúmen estadístico, o *summary statistics* en inglés).

### 3.1.1 LECTURA DE LOS DATOS

Mostramos las primeras líneas del archivo **Height.gwas.txt** donde se
encuentran los datos que necesitamos.

``` bash
head -5 Height.gwas.txt
```

    ## CHR  BP  SNP A1  A2  N   SE  P   OR  INFO    MAF
    ## 1    756604  rs3131962   A   G   388028  0.00301666  0.483171    0.997886915712657   0.890557941364774   0.369389592764921
    ## 1    768448  rs12562034  A   G   388028  0.00329472  0.834808    1.00068731609353    0.895893511351165   0.336845754096289
    ## 1    779322  rs4040617   G   A   388028  0.00303344  0.42897 0.997603556067569   0.897508290615237   0.377368010940814
    ## 1    801536  rs79373928  G   T   388028  0.00841324  0.808999    1.00203569922793    0.908962856432993   0.483212245374095

Observamos como se organiza el documento, destacamos el número de
cromosoma (**CHR**), el nombre del SNP de interés(**SNP**), los dos
alelos (**A1** y **A2**) y el p-valor (**P**).

Además de estos estadísticos, para nuestro análisis resaltamos el *Minor
Allele Frequency* (**MAF**), el valor de imputación (**INFO**) y el
estimador del tamaño del efecto, que es **OR** en nuestro caso.

### 3.1.2 FILTRADO DE LOS DATOS

Filtramos los datos para obtener el número de SNPs que nos interesan, en
nuestro caso con un valor de imputación superior a 0.8 y un minor allele
frequency superior a 0.01.

``` r
d <- read.table(gzfile("Height.gwas.txt.gz"),sep="\t", header = T)
snp_filt <- d[d$INFO > 0.8 & d$MAF > 0.01,]
snp_filt_<-nrow(snp_filt)
```

    ## [1] "Nos quedamos con 529495 SNPs"

A continuación eliminamos los SNPs duplicados, si los hubiese.

``` r
d2 <- read.table(gzfile("Height.gz"), sep = "\t", header = T)
dup <- d2[duplicated(d2[c('SNP')]), ]
no_dup <- d2[!duplicated(d2[c('SNP')]), ]
```

    ## [1] "Hay 2 SNPs duplicados"

Para finalizar eliminamos aquellas variantes que pudiesen ser ambigüas.

``` r
nodup_data <- read.table(gzfile("Height.nodup.gz"), sep = "\t", header = T) 
no_amb <- nodup_data[!nodup_data$A1 =="A" & nodup_data$A2 =="T" |! nodup_data$A1 =="T" & nodup_data$A2=="A" |! nodup_data$A1 =="C" & nodup_data$A2=="G" |! nodup_data$A1=="G" & nodup_data$A2=="C", ]
```

    ## [1] "Quedan 499617 SNPs no ambigüos"

## 3.2 CONTROL DE CALIDAD DE LOS DATOS TARGET

Los datos target necesarios para calcular PRS son, por lo general, los
obtenidos en el laboratorio, de tipo genotipo-fenotipo.

### 3.2.1 LECTURA DE LOS DATOS

El archivo descomprimido **EUR.zip** consta de varios archivos:

- .bim
- .fam
- .cov
- .bed
- .height

### 3.2.2 QC ESTANDAR DE UN ESTUDIO GWAS

Realizamos un control de calidad en estos datos con los mismos
parámetros que un estudio GWAS en cuánto a filtrado de SNPs y otros se
refiere (discutido en la memoria en el apartado de control de calidad de
los datos).

Para llevarlo a cabo utilizamos la herramienta de software **plink**.

``` bash
./plink --bfile EUR --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.01 --write-snplist --make-just-fam --out EUR.QC
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EUR.QC.log.
    ## Options in effect:
    ##   --bfile EUR
    ##   --geno 0.01
    ##   --hwe 1e-6
    ##   --maf 0.01
    ##   --make-just-fam
    ##   --mind 0.01
    ##   --out EUR.QC
    ##   --write-snplist
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 551892 variants loaded from .bim file.
    ## 503 people (240 males, 263 females) loaded from .fam.
    ## 14 people removed due to missing genotype data (--mind).
    ## IDs written to EUR.QC.irem .
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 489 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 1806 het. haploid genotypes present (see EUR.QC.hh ); many commands
    ## treat these as missing.
    ## Total genotyping rate in remaining samples is 0.999816.
    ## 5353 variants removed due to missing genotype data (--geno).
    ## Warning: --hwe observation counts vary by more than 10%, due to the X
    ## chromosome.  You may want to use a less stringent --hwe p-value threshold for X
    ## chromosome variants.
    ## --hwe: 944 variants removed due to Hardy-Weinberg exact test.
    ## 5061 variants removed due to minor allele threshold(s)
    ## (--maf/--max-maf/--mac/--max-mac).
    ## 540534 variants and 489 people pass filters and QC.
    ## Note: No phenotypes present.
    ## List of variant IDs written to EUR.QC.snplist .
    ## --make-just-fam to EUR.QC.fam ... done.

Se especifican los parámetros y el umbral que queramos mediante los
siguientes *flags*:

- el nombre del archivo con **bfile**, en nuestro caso todos los que
  empiecen con EUR;
- la frecuencia del alelo menos común con **maf**;
- el p-valor resultante del test (Fisher o chi cuadrado) de la ley de
  Hardy-Weinberg con **hwe**;
- los SNPs poco representados, ausentes en gran parte de los sujetos con
  **geno**;
- los sujetos con gran parte del genotipo ausente con **mind**.

Las *flags* restantes las utilizamos para evitar crear los ficheros
**.bed** asociados a los archivos resultantes del análisis y para
especificar el nombre del archivo de salida.

| **genotipo ausente** | **SNPs poco representados** | **SNPs fuera de HWE** | **SNPs bajo MAF** |          **RESULTADO**           |
|:--------------------:|:---------------------------:|:---------------------:|:-----------------:|:--------------------------------:|
|          14          |            5353             |          944          |       5061        | Quedan 540534 SNPs y 489 sujetos |

A continuación filtramos los SNPs correlacionados significativamente.

``` bash
./plink --bfile EUR --keep EUR.QC.fam --extract EUR.QC.snplist --indep-pairwise 200 50 0.25 --out EUR.QC
```

Se ha especificado con las siguientes *flags*:

- los sujetos que queremos incluir en el análisis con **keep**;
- los SNPs que queremos utilizar para el análisis con **extract**;
- los parámetros de filtrado con **indep-pairwise**;

De este análisis resultan dos ficheros; Encontraremos los SNPs que nos
interesan, con un coeficiente de determinación $r^2$ de desquilibrio de
ligamiento (LD) menor a 0.25 en **EUR.QC.prune.in** Al finalizar el
análisis se eliminan 272077 SNPs de los 540534 totales. Quedan 268457
SNPs.

Calculamos el índice de heterocigosidad para filtrar aquellos sujetos
con valores demasiado elevados o leves. Nos quedamos con los SNPs
presentes en **EUR.QC.prune.in** y añadimos la *flag* –het que nos
permite hacer este cálculo.

``` bash
./plink --bfile EUR --extract EUR.QC.prune.in --keep EUR.QC.fam --het --out EUR.QC
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EUR.QC.log.
    ## Options in effect:
    ##   --bfile EUR
    ##   --extract EUR.QC.prune.in
    ##   --het
    ##   --keep EUR.QC.fam
    ##   --out EUR.QC
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 551892 variants loaded from .bim file.
    ## 503 people (240 males, 263 females) loaded from .fam.
    ## --extract: 268457 variants remaining.
    ## --keep: 489 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 489 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 558 het. haploid genotypes present (see EUR.QC.hh ); many commands
    ## treat these as missing.
    ## Total genotyping rate in remaining samples is 0.999958.
    ## 268457 variants and 489 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --het: 263228 variants scanned, report written to EUR.QC.het .

De este análisis resulta el archivo **EUR.QC.het** con los coeficientes
de estimación para la heterocigosidad.

Eliminamos aquellos que superen cierto distanciamiento de la media.

``` r
d_het <- read.table("EUR.QC.het", header = T)
media <- mean(d_het$F)
st_dev <- sd(d_het$F)
res <- subset(d_het, F <= media+3*st_dev & F >= media-3*st_dev)
write.table(res[,c(1,2)], "EUR.valid.sample", quote=F, row.names = F)
```

    ## [1] "Nos quedamos con 487 sujetos"

### 3.2.3 SNPs DESAJUSTADOS

Daremos la vuelta a las hebras de los SNPs con alelos no
complementarios. Para ello utilizaremos el archivo **.bim**, el GWAS con
los estadísticos y el QC realizado **Height.QC.gz** y la lista de SNPs
que han pasado el QC **EUR.QC.snplist**.

``` r
bim <- read.table("EUR.bim")
head(bim)
```

    ##   V1         V2       V3     V4 V5 V6
    ## 1  1  rs3131962 0.490722 756604  A  G
    ## 2  1 rs12562034 0.495714 768448  0  0
    ## 3  1  rs4040617 0.500708 779322  G  A
    ## 4  1 rs79373928 0.587220 801536  G  T
    ## 5  1 rs11240779 0.620827 808631  G  A
    ## 6  1 rs57181708 0.620827 809876  G  A

Vemos que este archivo no contiene cabecera. Le damos nombre a las
columnas para poder cruzar los datos con los demás archivos.

``` r
colnames(bim)<-c("CHR","SNP","CM","BP","B.A1","B.A2")
head(bim)
```

    ##   CHR        SNP       CM     BP B.A1 B.A2
    ## 1   1  rs3131962 0.490722 756604    A    G
    ## 2   1 rs12562034 0.495714 768448    0    0
    ## 3   1  rs4040617 0.500708 779322    G    A
    ## 4   1 rs79373928 0.587220 801536    G    T
    ## 5   1 rs11240779 0.620827 808631    G    A
    ## 6   1 rs57181708 0.620827 809876    G    A

``` r
qc <- read.table("EUR.QC.snplist", header = F, stringsAsFactors = F)
head(qc)
```

    ##           V1
    ## 1  rs3131962
    ## 2  rs4040617
    ## 3 rs79373928
    ## 4 rs11240779
    ## 5 rs57181708
    ## 6  rs4422948

``` r
height <- read.table(gzfile("Height.QC.gz"),sep="\t", header = T, stringsAsFactors = F)
head(height)
```

    ##   CHR     BP        SNP A1 A2      N         SE        P        OR      INFO
    ## 1   1 756604  rs3131962  A  G 388028 0.00301666 0.483171 0.9978869 0.8905579
    ## 2   1 768448 rs12562034  A  G 388028 0.00329472 0.834808 1.0006873 0.8958935
    ## 3   1 779322  rs4040617  G  A 388028 0.00303344 0.428970 0.9976036 0.8975083
    ## 4   1 801536 rs79373928  G  T 388028 0.00841324 0.808999 1.0020357 0.9089629
    ## 5   1 808631 rs11240779  G  A 388028 0.00242821 0.590265 1.0013083 0.8932125
    ## 6   1 809876 rs57181708  G  A 388028 0.00336785 0.714750 1.0012317 0.9235576
    ##         MAF
    ## 1 0.3693896
    ## 2 0.3368458
    ## 3 0.3773680
    ## 4 0.4832122
    ## 5 0.4504096
    ## 6 0.4997439

Pasamos todos los alelos en mayúsculas para facilitar el análisis
posterior.

``` r
library(stringr)
str_upp<- str_detect(bim$B.A1,"[[:upper:]]")
table(str_upp)
```

    ## str_upp
    ##  FALSE   TRUE 
    ##     74 551818

``` r
bim$B.A1[2:2]
```

    ## [1] "0"

Comprobamos que algunos alelos no están en mayúsculas, sin embargo
algunos simplemente no tienen alelo por lo que se detectan y devuelven
valor FALSO.

``` r
height$A1 <- toupper(height$A1)
height$A2 <- toupper(height$A2)
bim$B.A1 <- toupper(bim$B.A1)
bim$B.A2 <- toupper(bim$B.A2)
```

Buscamos los SNPs que no tengan alelos complementarios entre los datos
base y target, o sea entre el GWAS de los datos base y el archivo bim de
los datos target. Fusionamos ambos datos por las columnas que nos
interesan y filtramos para quedarnos con los SNPs que han pasado todo el
QC (lista de SNPs resultante del análisis en los datos target).

``` r
info <- merge(bim, height, by = c("CHR", "BP", "SNP"))
info <- info[info$SNP %in% qc$V1,]
```

Buscamos los SNPs que tengas el mismo alelo en ambos datos.

``` r
info.match <- subset(info, A1==B.A1 & A2==B.A2)
match<-nrow(info.match)
```

    ## [1] "485884 SNPs tienen el mismo alelo en ambos datos"

Buscamos los alelos que son complementarios entre ambos datos y los
intercambiamos cuando sea el caso, guardamos en una nueva columna.

``` r
compl <- function(x) {
    switch (x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
}

info$C.A1 <- sapply(info$B.A1, compl)
info$C.A2 <- sapply(info$B.A2, compl)
info.compl <- subset(info, A1 == C.A1 & A2 == C.A2)
compl_snp<-nrow(info.compl)
```

    ## [1] "Hay 0 SNPs complementarios"

Buscamos los alelos que necesiten recodificarse, como por ejemplo A/C y
C/A.

``` r
info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
rec <- nrow(info.recode)
```

    ## [1] "Hay 3921 SNPs que deben recodificarse"

``` r
recode.snps <- bim$SNP %in% info.recode$SNP
table(recode.snps)
```

    ## recode.snps
    ##  FALSE   TRUE 
    ## 547971   3921

Intercambiamos el orden de los alelos en los datos target.

``` r
tmp_rec <- bim[recode.snps,]$B.A1
bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
bim[recode.snps,]$B.A2 <- tmp_rec
```

``` r
info.comp.recode <- subset(info, A1 == C.A2 & A2 == C.A1)

cmp.snps <- bim$SNP %in% info.comp.recode$SNP
table(cmp.snps)
```

    ## cmp.snps
    ##  FALSE 
    ## 551892

Sigue sin haber SNPs complementarios que recodificar.

Finalmente eliminamos los SNPs que tengan alelos diferentes en ambos
datos, aún habiendo realizado los pasos anteriores.

``` r
mis <- bim$SNP[!(bim$SNP %in% info.match$SNP | bim$SNP %in% info.compl$SNP | bim$SNP %in% info.recode$SNP | bim$SNP %in% info.comp.recode$SNP)]
mis_ <- length(mis)
```

    ## [1] "Hay 62087 SNPs incompatibles"

### 3.2.4 SNPs DUPLICADOS

Miramos que no haya duplicados.

``` r
dup2 <- bim[duplicated(bim[c('SNP')]), ]
nrow(dup2)
```

    ## [1] 0

### 3.2.5 CROMOSOMAS SEXUALES

Verificamos el sexo de los sujetos, fijando un umbral del estimador de
homocigosidad en el cromosoma X.

Primero filtramos según los SNPs de interés contenidos en el archivo
**EUR.QC.prune.in** y los sujetos validados en **EUR.valid.sample**.
Añadimos la *flag* –check-sex que creara un archivo con los índices de
fijación (o estadístico F) para cada individuo.

``` bash
./plink --bfile EUR --extract EUR.QC.prune.in --keep EUR.valid.sample --check-sex --out EUR.QC
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EUR.QC.log.
    ## Options in effect:
    ##   --bfile EUR
    ##   --check-sex
    ##   --extract EUR.QC.prune.in
    ##   --keep EUR.valid.sample
    ##   --out EUR.QC
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 551892 variants loaded from .bim file.
    ## 503 people (240 males, 263 females) loaded from .fam.
    ## --extract: 268457 variants remaining.
    ## --keep: 487 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 487 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 558 het. haploid genotypes present (see EUR.QC.hh ); many commands
    ## treat these as missing.
    ## Total genotyping rate in remaining samples is 0.999958.
    ## 268457 variants and 487 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --check-sex: 5229 Xchr and 0 Ychr variant(s) scanned, 4 problems detected.
    ## Report written to EUR.QC.sexcheck .

En el informe de salida de plink vemos que hay 4 sujetos problemáticos
que tenemos que extraer de los datos. Para ello los quitaremos de la
lista de sujetos válidos cruzando las tablas de datos. En el fichero
creado por plink aparece una columna **STATUS** donde se menciona OK o
PROBLEM según el índice de fijación calculado.

``` r
validos <- read.table("EUR.valid.sample", header = T)
sex_dat <- read.table("EUr.QC.sexcheck", header = T)
validos <- subset(sex_dat, STATUS=="OK" & FID %in% validos$FID)
```

    ## [1] "Quedan 483 sujetos"

### 3.2.6 PARENTESCO

Verificamos que no hay sujetos entre ambos datos que tengan cierto grado
de parentesco. Filtramos como siempre según los SNPs y sujetos de
interés que hayan pasado los QC y añadimos la *flag* rel-cutoff poniendo
un límite de 0.125.

``` bash
./plink --bfile EUR --extract EUR.QC.prune.in --keep EUR.QC.valid --rel-cutoff 0.125 --out EUR.QC
```

No se elimina ningún sujeto. Encontramos la lista de los sujetos
validados en el archivo EUR.QC.rel.id.

### 3.2.7 ARCHIVO TARGET FINAL

Creamos el archivo final necesario par el cálculo de PRS a partir de los
sujetos y SNPs validados en las diferentes etapas del QC que acabamos de
realizar. Nos quedamos con los sujetos presentes en el archivo
EUR.QC.rel.id y los SNPs de la lista EUR.QC.snplist. Extraemos aquellos
SNPs incompatibles presentes en el archivo EUR.mismatch. Utilizamos el
archivo EUR.a1 para asegurarnos que los alelos codificantes sean los
validados en dicho archivo con la *flag* a1-allele.

``` bash
./plink --bfile EUR --make-bed --keep EUR.QC.rel.id  --extract EUR.QC.snplist --exclude EUR.mismatch --a1-allele EUR.a1 --out EUR.QC
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EUR.QC.log.
    ## Options in effect:
    ##   --a1-allele EUR.a1
    ##   --bfile EUR
    ##   --exclude EUR.mismatch
    ##   --extract EUR.QC.snplist
    ##   --keep EUR.QC.rel.id
    ##   --make-bed
    ##   --out EUR.QC
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 551892 variants loaded from .bim file.
    ## 503 people (240 males, 263 females) loaded from .fam.
    ## --extract: 540534 variants remaining.
    ## --exclude: 489805 variants remaining.
    ## --keep: 483 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate in remaining samples is exactly 1.
    ## --a1-allele: 489805 assignments made.
    ## 489805 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to EUR.QC.bed + EUR.QC.bim + EUR.QC.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

# 4 CALCULO DE PRS

Trabajamos con el software de análisis genético PLINK-1.9, probando
distintos parametrajes y validando los datos con

## 4.1 PLINK

Necesitamos los siguientes archivos, obtenidos en gran parte del QC
realizado en el apartado anterior:

En nuestro caso empezamos por cambiar el valor de OR (*odds ratio*,
estimador del tamaño de efecto para casos binarios) a logaritmo para
facilitar su cálculo y lo añadimos a una nueva columna BETA.

``` r
#Aunque podríamos usar la variable no_amb donde se encuentran los datos con el QC finalizado, preferimos cargar el archivo por si hubiese errores durante el análisis posterior.
dat_or <- read.table(gzfile("Height.QC.gz"), header = T)
dat_or$BETA <- log(dat_or$OR)
write.table(dat_or,"Height.QC.transformed", quote = F, row.names = F, sep = "\t")
```

### 4.1.1 AGLOMERACION/CLUMPING

Filtramos los SNPs de manera a crear aglomeraciones (*clumps*) para
quedarnos con aquellos mas propensos a ser los causantes del fenotipo
estudiado.

``` bash
./plink --bfile EUR.QC --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump Height.QC.transformed --clump-snp-field SNP --clump-field P --out EUR
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EUR.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --clump Height.QC.transformed
    ##   --clump-field P
    ##   --clump-kb 250
    ##   --clump-p1 1
    ##   --clump-r2 0.1
    ##   --clump-snp-field SNP
    ##   --out EUR
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 489805 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 'rs1076829' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3129818' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3118359' is missing from the main dataset, and is a top variant.
    ## 9809 more top variant IDs missing; see log file.
    ## --clump: 193758 clumps formed from 489805 top variants.
    ## Results written to EUR.clumped .

Utilizamos los siguientes parámetros

- 1 como índex del p-valor con –clump-p1;
- 0.1 como umbral del coeficiente de determinación con –clump-r2;
- 250 kilobases de radio de un SNP indexado con –clumb-kb; Se
  especifican las columnas donde se encuentran los valores de interés
  mediante –clump-snp-field y –clump-field.

Según el output se han formado 193758 aglomeraciones para 489805
variantes.

Nos aparecen algunos mensajes de error.

Se extraen los SNPs indexados para su uso posterior.

``` r
clump <- read.table("EUR.clumped", header = T, stringsAsFactors = F)
valid_clump <- clump[,c("SNP")]

write.table(valid_clump,"EUR.valid.snp", quote = F, row.names = F, col.names = F)
```

### 4.1.2 CALCULO DE PRS

Para el cálculo necesitamos los datos base y los p-valores de los SNPs
de interés. Cargamos el fichero y miramos las columnas que nos interesan
para aislarlas.

``` r
height_qc_tr <- read.table("Height.QC.transformed", header = T, sep = "\t")
head(height_qc_tr)
```

    ##   CHR     BP        SNP A1 A2      N         SE        P        OR      INFO
    ## 1   1 756604  rs3131962  A  G 388028 0.00301666 0.483171 0.9978869 0.8905579
    ## 2   1 768448 rs12562034  A  G 388028 0.00329472 0.834808 1.0006873 0.8958935
    ## 3   1 779322  rs4040617  G  A 388028 0.00303344 0.428970 0.9976036 0.8975083
    ## 4   1 801536 rs79373928  G  T 388028 0.00841324 0.808999 1.0020357 0.9089629
    ## 5   1 808631 rs11240779  G  A 388028 0.00242821 0.590265 1.0013083 0.8932125
    ## 6   1 809876 rs57181708  G  A 388028 0.00336785 0.714750 1.0012317 0.9235576
    ##         MAF        BETA
    ## 1 0.3693896 -0.00211532
    ## 2 0.3368458  0.00068708
    ## 3 0.3773680 -0.00239932
    ## 4 0.4832122  0.00203363
    ## 5 0.4504096  0.00130747
    ## 6 0.4997439  0.00123090

``` r
snp.val <- height_qc_tr[, c("SNP","P")]
```

``` bash
echo "0.001 0 0.001" > range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list
```

Con estos dos archivos creados calculamos los PRS.

``` bash
./plink --bfile EUR.QC --score Height.QC.transformed 3 4 12 header --q-score-range range_list SNP.value --extract EUR.valid.snp --out EUR
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EUR.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --extract EUR.valid.snp
    ##   --out EUR
    ##   --q-score-range range_list SNP.value
    ##   --score Height.QC.transformed 3 4 12 header
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## --extract: 193758 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 193758 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 305859 lines skipped in --score file (305859 due to variant ID
    ## mismatch, 0 due to allele code mismatch); see EUR.nopred for details.
    ## --score: 193758 valid predictors loaded.
    ## Warning: 305860 lines skipped in --q-score-range data file.
    ## --score: 7 ranges processed.
    ## Results written to EUR.*.profile.

Utilizamos los datos del archivo Height.QC.transformed para extraer los
datos que nos interesan: los identificadores de los SNPs, el alelo
efectivo (A1) y la ultima columna con el estimador (Beta). Con
q-score-range definimos los diferentes umbrales del p-valor e indicamos
el archivo donde se encuentran estos (SNP.value). Se generan 5 listas,
con los umbrales como nombre de archivo, en los que se guardan los
diferentes cálculos.

Vemos que aparece un error en la lectura de líneas. Probamos tomando 10
sujetos random.

``` r
set.seed(111)
height_prueba <- height_qc_tr[sample(nrow(height_qc_tr), size=10), ]
write.table(height_prueba, "prueba", quote = F, row.names = F, sep = "\t")
```

``` bash
./plink --bfile EUR.QC --score prueba 3 4 12 header --q-score-range range_list SNP.value --extract EUR.valid.snp --out EUR.prueba
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EUR.prueba.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --extract EUR.valid.snp
    ##   --out EUR.prueba
    ##   --q-score-range range_list SNP.value
    ##   --score prueba 3 4 12 header
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## --extract: 193758 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 193758 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 4 lines skipped in --score file (4 due to variant ID mismatch, 0 due
    ## to allele code mismatch); see EUR.prueba.nopred for details.
    ## --score: 6 valid predictors loaded.
    ## Warning: 499612 lines skipped in --q-score-range data file.
    ## --score: 7 ranges processed.
    ## Results written to EUR.prueba.*.profile.

Siguen habiendo 4 líneas saltadas. Comprobamos que todos los ID de los
SNPs estén presentes en ambos datos. Vemos que los 10 SNPs ID coinciden.

``` r
height_prueba %>%select(which((height_prueba$SNP %in% height_qc_tr$SNP)))
```

    ##        CHR        BP         SNP A1 A2      N         SE           P        OR
    ## 375291  13 109681102   rs9521132  G  A 388028 0.00282889 8.27832e-01 0.9993850
    ## 230355   7 117252977  rs73215927  G  A 388028 0.00453377 6.75064e-01 0.9981012
    ## 420633  16  72188376 rs117523217  C  T 388028 0.00652443 2.59761e-02 1.0146332
    ## 165957   5 125473558   rs3842965  T  C 388028 0.00249917 1.05315e-04 1.0097390
    ## 137212   4 144892048  rs17019056  A  C 388028 0.00243842 1.04786e-08 0.9861420
    ## 497434  22  44737788   rs9626539  T  G 388028 0.00335470 9.06833e-01 0.9996075
    ## 3117     1  14737314 rs116851506  A  G 388028 0.00707766 9.47370e-02 1.0118964
    ## 291402  10  21283327  rs77314108  T  G 388028 0.00755266 3.78840e-01 0.9933754
    ## 259951   8 127076935    rs983956  A  C 388028 0.00651095 1.73785e-01 0.9911833
    ## 499022  22  49438931  rs77779902  A  C 388028 0.00606291 4.78583e-01 0.9957131
    ##             INFO
    ## 375291 0.9201578
    ## 230355 0.9092357
    ## 420633 0.8845898
    ## 165957 0.9191940
    ## 137212 0.8862319
    ## 497434 0.9073121
    ## 3117   0.8960337
    ## 291402 0.9058828
    ## 259951 0.9004413
    ## 499022 0.8992535

### 4.1.3 ESTRATIFICACION POBLACIONAL

Hacemos un análisis de componentes principales para añadirlo como
covariable. Filtramos según los parámetros que se decidan y calculamos
los 6 componentes principales con la *flag* –pca.

``` bash
./plink --bfile EUR.QC --indep-pairwise 200 50 0.25 --out EUR2
```

``` bash
./plink --bfile EUR.QC --extract EUR2.prune.in --pca 6 --out EUR
```

Encontramos los resultados en el archivo EUR.eigenvec.

### 4.1.4 PRS MEJOR AJUSTADO

Buscamos el umbral del p-valor que proporciona el PRS mejor ajustado
para el fenotipo estudiado. Para ello hacemos un análisis de regresión.

Necesitamos el archivo con los fenotipos **EUR.height**,los PCs
calculados **EUR.eigenvec** y las covariables.

``` r
feno_h <- read.table("EUR.height", header = T)
pcs <- read.table("EUR.eigenvec", header = F)
head(pcs)
```

    ##        V1      V2           V3        V4           V5          V6         V7
    ## 1 HG00096 HG00096  0.000643305 0.0664323 -1.47374e-02 -0.03599560 -0.0163573
    ## 2 HG00097 HG00097  0.001413780 0.0736016  8.81567e-03 -0.02058210 -0.0116844
    ## 3 HG00099 HG00099  0.002646810 0.0717702 -2.09757e-02 -0.00608687 -0.0141427
    ## 4 HG00101 HG00101  0.001697520 0.0854445 -1.56991e-02 -0.00289015 -0.0335167
    ## 5 HG00102 HG00102  0.004411350 0.0696362  1.75584e-06 -0.02643000 -0.0477571
    ## 6 HG00103 HG00103 -0.004312160 0.0571794 -8.19476e-03 -0.01349190  0.0159968
    ##            V8
    ## 1 -0.02093710
    ## 2  0.02240870
    ## 3 -0.00713368
    ## 4 -0.01411770
    ## 5 -0.03144630
    ## 6 -0.01206940

``` r
colnames(pcs)<- c("FID","IID","PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
```

``` r
cov_sex <- read.table("EUR.cov", header = T)
```

Fusionamos el fenotipo junto con las dos covariables, el sexo y los pca.

``` r
feno <- merge(feno_h, cov_sex, by = c("FID","IID"))
feno <- merge(feno, pcs, by = c("FID","IID"))
colnames(feno)
```

    ##  [1] "FID"    "IID"    "Height" "Sex"    "PC1"    "PC2"    "PC3"    "PC4"   
    ##  [9] "PC5"    "PC6"

Hacemos una regresión linear para calcular el modelo con PRS, con la
altura como variable dependiente. Quitamos los identificadores para que
no sean tomados como covariables.

``` r
model <- lm(Height~., feno[,!colnames(feno) %in% c("FID","IID")])
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = Height ~ ., data = feno[, !colnames(feno) %in% c("FID", 
    ##     "IID")])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.21867 -0.58069 -0.02043  0.53206  2.63620 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept) 168.83050    0.12425 1358.768  < 2e-16 ***
    ## Sex           0.86527    0.07769   11.137  < 2e-16 ***
    ## PC1          -0.67849    0.84485   -0.803  0.42234    
    ## PC2          -2.55995    0.85715   -2.987  0.00297 ** 
    ## PC3          -0.36689    0.84544   -0.434  0.66452    
    ## PC4          -0.11097    1.57358   -0.071  0.94381    
    ## PC5           0.14327    0.92594    0.155  0.87710    
    ## PC6          -0.01685    0.88817   -0.019  0.98488    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8378 on 464 degrees of freedom
    ## Multiple R-squared:  0.2253, Adjusted R-squared:  0.2137 
    ## F-statistic: 19.28 on 7 and 464 DF,  p-value: < 2.2e-16

``` r
r2.prs<-summary(model)$r.squared
```

Cargamos los ficheros con los resultados de PRS según el umbral de
p-valor impuesto.

``` r
p.val.th <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
prs.result <- NULL
for(i in p.val.th){
  prs <- read.table(paste0("EUR.",i,".profile"), header = T)
  feno_prs <- merge(feno, prs[,c("FID","IID","SCORE")], by = c("FID", "IID"))
  model_prs <- lm(Height~., feno_prs[,!colnames(feno_prs) %in% c("FID","IID")])
  r2_model_prs <- summary(model_prs)$r.squared
  prs_r2 <- r2_model_prs-r2.prs
  prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs_r2))
}
```

``` r
prs.result[which.max(prs.result$R2),]
```

    ##   Threshold        R2
    ## 5       0.3 0.1612372

El PRS mejor ajustado se encuentra en el umbral 0.3 del p-valor y
explica un 0.1612372 la variación fenotípica.

## 4.2 CONFIGURACION DE PARAMETROS

Probamos distintas configuraciones de los parámetroa para ver si existen diferencias y
variantes raras que se excluyen cuando fijamos los umbrales estándar.

### 4.2.1 RADIO DE AGLOMERACION/CLUMPING

Probamos ampliando el radio de aglomeración pero manteniendo el $r^2_c$
a 0.1. Partiendo de una ventada de kb de entre {50, 100, 250}, creamos 3
modelos. Dividimos el radio de aglomeración usado anteriormente por el
valor del coeficiente de determinación para determinar el nuevo radio.

#### 4.2.1.1 500kb

``` bash
./plink --bfile EUR.QC --clump-p1 1 --clump-r2 0.1 --clump-kb 500 --clump Height.QC.transformed --clump-snp-field SNP --clump-field P --out EURcl500
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl500.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --clump Height.QC.transformed
    ##   --clump-field P
    ##   --clump-kb 500
    ##   --clump-p1 1
    ##   --clump-r2 0.1
    ##   --clump-snp-field SNP
    ##   --out EURcl500
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 489805 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 'rs1076829' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3129818' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3118359' is missing from the main dataset, and is a top variant.
    ## 9809 more top variant IDs missing; see log file.
    ## --clump: 181545 clumps formed from 489805 top variants.
    ## Results written to EURcl500.clumped .

Vemos que hay 181545 aglomeraciones.

``` r
clump_500 <- read.table("EURcl500.clumped", header = T, stringsAsFactors = F)
valid_clump_500 <- clump_500[,c("SNP")]

write.table(valid_clump_500,"EURcl500.valid.snp", quote = F, row.names = F, col.names = F)
```

``` bash
./plink --bfile EUR.QC --score Height.QC.transformed 3 4 12 header --q-score-range range_list SNP.value --extract EURcl500.valid.snp --out EURcl500
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl500.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --extract EURcl500.valid.snp
    ##   --out EURcl500
    ##   --q-score-range range_list SNP.value
    ##   --score Height.QC.transformed 3 4 12 header
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## --extract: 181545 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 181545 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 318072 lines skipped in --score file (318072 due to variant ID
    ## mismatch, 0 due to allele code mismatch); see EURcl500.nopred for details.
    ## --score: 181545 valid predictors loaded.
    ## Warning: 318073 lines skipped in --q-score-range data file.
    ## --score: 7 ranges processed.
    ## Results written to EURcl500.*.profile.

``` r
p.val.th <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
prs.result.cl500 <- NULL
for(i in p.val.th){
  prs.cl500 <- read.table(paste0("EURcl500.",i,".profile"), header = T)
  feno_prs.cl500 <- merge(feno, prs.cl500[,c("FID","IID","SCORE")], by = c("FID", "IID"))
  model_prs.cl500 <- lm(Height~., feno_prs.cl500[,!colnames(feno_prs.cl500) %in% c("FID","IID")])
  r2_model_prs.cl500 <- summary(model_prs.cl500)$r.squared
  prs_r2.cl500 <- r2_model_prs.cl500-r2.prs
  prs.result.cl500 <- rbind(prs.result.cl500, data.frame(Threshold=i, R2=prs_r2.cl500))
}
```

``` r
prs.result.cl500[which.max(prs.result.cl500$R2),]
```

    ##   Threshold        R2
    ## 5       0.3 0.1826487

El PRS mejor ajustado en este caso se encuentra en el umbral 03. del
p-valor y explica un 0.1826487 la variación fenotípica.

\####1000kb

``` bash
./plink --bfile EUR.QC --clump-p1 1 --clump-r2 0.1 --clump-kb 1000 --clump Height.QC.transformed --clump-snp-field SNP --clump-field P --out EURcl1000
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl1000.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --clump Height.QC.transformed
    ##   --clump-field P
    ##   --clump-kb 1000
    ##   --clump-p1 1
    ##   --clump-r2 0.1
    ##   --clump-snp-field SNP
    ##   --out EURcl1000
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 489805 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 'rs1076829' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3129818' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3118359' is missing from the main dataset, and is a top variant.
    ## 9809 more top variant IDs missing; see log file.
    ## --clump: 177471 clumps formed from 489805 top variants.
    ## Results written to EURcl1000.clumped .

Se han creado 177471 aglomeraciones.

``` r
clump_1000 <- read.table("EURcl1000.clumped", header = T, stringsAsFactors = F)
valid_clump_1000 <- clump_1000[,c("SNP")]

write.table(valid_clump_1000,"EURcl1000.valid.snp", quote = F, row.names = F, col.names = F)
```

``` bash
./plink --bfile EUR.QC --score Height.QC.transformed 3 4 12 header --q-score-range range_list SNP.value --extract EURcl1000.valid.snp --out EURcl1000
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl1000.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --extract EURcl1000.valid.snp
    ##   --out EURcl1000
    ##   --q-score-range range_list SNP.value
    ##   --score Height.QC.transformed 3 4 12 header
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## --extract: 177471 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 177471 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 322146 lines skipped in --score file (322146 due to variant ID
    ## mismatch, 0 due to allele code mismatch); see EURcl1000.nopred for details.
    ## --score: 177471 valid predictors loaded.
    ## Warning: 322147 lines skipped in --q-score-range data file.
    ## --score: 7 ranges processed.
    ## Results written to EURcl1000.*.profile.

``` r
p.val.th <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
prs.result.cl1000 <- NULL
for(i in p.val.th){
  prs.cl1000 <- read.table(paste0("EURcl1000.",i,".profile"), header = T)
  feno_prs.cl1000 <- merge(feno, prs.cl1000[,c("FID","IID","SCORE")], by = c("FID", "IID"))
  model_prs.cl1000 <- lm(Height~., feno_prs.cl1000[,!colnames(feno_prs.cl1000) %in% c("FID","IID")])
  r2_model_prs.cl1000 <- summary(model_prs.cl1000)$r.squared
  prs_r2.cl1000 <- r2_model_prs.cl1000-r2.prs
  prs.result.cl1000 <- rbind(prs.result.cl1000, data.frame(Threshold=i, R2=prs_r2.cl1000))
}
```

``` r
prs.result.cl1000[which.max(prs.result.cl1000$R2),]
```

    ##   Threshold        R2
    ## 5       0.3 0.1971698

El PRS mejor ajustado en este caso se encuentra en el umbral 0.3 del
p-valor y explica un 0.1971695 la variación fenotípica.

#### 4.2.1.2 2500kb

``` bash
./plink --bfile EUR.QC --clump-p1 1 --clump-r2 0.1 --clump-kb 2500 --clump Height.QC.transformed --clump-snp-field SNP --clump-field P --out EURcl2500
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl2500.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --clump Height.QC.transformed
    ##   --clump-field P
    ##   --clump-kb 2500
    ##   --clump-p1 1
    ##   --clump-r2 0.1
    ##   --clump-snp-field SNP
    ##   --out EURcl2500
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 489805 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 'rs1076829' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3129818' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3118359' is missing from the main dataset, and is a top variant.
    ## 9809 more top variant IDs missing; see log file.
    ## --clump: 176577 clumps formed from 489805 top variants.
    ## Results written to EURcl2500.clumped .

Vemos que se han creado 176577 aglomeraciones. Continuamos el análisis
entero.

``` r
clump_2500 <- read.table("EURcl2500.clumped", header = T, stringsAsFactors = F)
valid_clump_2500 <- clump_2500[,c("SNP")]

write.table(valid_clump_2500,"EURcl2500.valid.snp", quote = F, row.names = F, col.names = F)
```

``` bash
./plink --bfile EUR.QC --score Height.QC.transformed 3 4 12 header --q-score-range range_list SNP.value --extract EURcl2500.valid.snp --out EURcl2500
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl2500.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --extract EURcl2500.valid.snp
    ##   --out EURcl2500
    ##   --q-score-range range_list SNP.value
    ##   --score Height.QC.transformed 3 4 12 header
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## --extract: 176577 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 176577 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 323040 lines skipped in --score file (323040 due to variant ID
    ## mismatch, 0 due to allele code mismatch); see EURcl2500.nopred for details.
    ## --score: 176577 valid predictors loaded.
    ## Warning: 323041 lines skipped in --q-score-range data file.
    ## --score: 7 ranges processed.
    ## Results written to EURcl2500.*.profile.

``` r
p.val.th <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
prs.result.cl2500 <- NULL
for(i in p.val.th){
  prs.cl2500 <- read.table(paste0("EURcl2500.",i,".profile"), header = T)
  feno_prs.cl2500 <- merge(feno, prs.cl2500[,c("FID","IID","SCORE")], by = c("FID", "IID"))
  model_prs.cl2500 <- lm(Height~., feno_prs.cl2500[,!colnames(feno_prs.cl2500) %in% c("FID","IID")])
  r2_model_prs.cl2500 <- summary(model_prs.cl2500)$r.squared
  prs_r2.cl2500 <- r2_model_prs.cl2500-r2.prs
  prs.result.cl2500 <- rbind(prs.result.cl2500, data.frame(Threshold=i, R2=prs_r2.cl2500))
}
```

``` r
prs.result.cl2500[which.max(prs.result.cl2500$R2),]
```

    ##   Threshold        R2
    ## 5       0.3 0.1980358

El PRS mejor ajustado en este caso se encuentra en el umbral 0.3 del
p-valor y explica un 0.1980359 la variación fenotípica.

### 4.2.2 UMBRAL DE AGLOMERACION/CLUMPING + RADIO

Probamos en primer lugar de cambiar el umbral del $r^2_c$ a 0.2, el cuál
es menos restrictivo.

``` bash
./plink --bfile EUR.QC --clump-p1 1 --clump-r2 0.2 --clump-kb 250 --clump Height.QC.transformed --clump-snp-field SNP --clump-field P --out EURcl2
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl2.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --clump Height.QC.transformed
    ##   --clump-field P
    ##   --clump-kb 250
    ##   --clump-p1 1
    ##   --clump-r2 0.2
    ##   --clump-snp-field SNP
    ##   --out EURcl2
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 489805 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 'rs1076829' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3129818' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3118359' is missing from the main dataset, and is a top variant.
    ## 9809 more top variant IDs missing; see log file.
    ## --clump: 256119 clumps formed from 489805 top variants.
    ## Results written to EURcl2.clumped .

Vemos que se han creado 256119 aglomeraciones. Continuamos el análisis
entero.

``` r
clump_2 <- read.table("EURcl2.clumped", header = T, stringsAsFactors = F)
valid_clump_2 <- clump_2[,c("SNP")]

write.table(valid_clump_2,"EURcl2.valid.snp", quote = F, row.names = F, col.names = F)
```

``` bash
./plink --bfile EUR.QC --score Height.QC.transformed 3 4 12 header --q-score-range range_list SNP.value --extract EURcl2.valid.snp --out EURcl2
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl2.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --extract EURcl2.valid.snp
    ##   --out EURcl2
    ##   --q-score-range range_list SNP.value
    ##   --score Height.QC.transformed 3 4 12 header
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## --extract: 256119 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 256119 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 243498 lines skipped in --score file (243498 due to variant ID
    ## mismatch, 0 due to allele code mismatch); see EURcl2.nopred for details.
    ## --score: 256119 valid predictors loaded.
    ## Warning: 243499 lines skipped in --q-score-range data file.
    ## --score: 7 ranges processed.
    ## Results written to EURcl2.*.profile.

``` r
p.val.th <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
prs.result.cl2 <- NULL
for(i in p.val.th){
  prs.cl2 <- read.table(paste0("EURcl2.",i,".profile"), header = T)
  feno_prs.cl2 <- merge(feno, prs.cl2[,c("FID","IID","SCORE")], by = c("FID", "IID"))
  model_prs.cl2 <- lm(Height~., feno_prs.cl2[,!colnames(feno_prs.cl2) %in% c("FID","IID")])
  r2_model_prs.cl2 <- summary(model_prs.cl2)$r.squared
  prs_r2.cl2 <- r2_model_prs.cl2-r2.prs
  prs.result.cl2 <- rbind(prs.result.cl2, data.frame(Threshold=i, R2=prs_r2.cl2))
}
```

``` r
prs.result.cl2[which.max(prs.result.cl2$R2),]
```

    ##   Threshold        R2
    ## 6       0.4 0.1643983

El PRS mejor ajustado en este caso se encuentra en el umbral 0.4 del
p-valor y explica un 0.1643983 la variación fenotípica.

#### 4.2.2.1 500

``` bash
./plink --bfile EUR.QC --clump-p1 1 --clump-r2 0.2 --clump-kb 500 --clump Height.QC.transformed --clump-snp-field SNP --clump-field P --out EURcl500.2
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl500.2.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --clump Height.QC.transformed
    ##   --clump-field P
    ##   --clump-kb 500
    ##   --clump-p1 1
    ##   --clump-r2 0.2
    ##   --clump-snp-field SNP
    ##   --out EURcl500.2
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 489805 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 'rs1076829' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3129818' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3118359' is missing from the main dataset, and is a top variant.
    ## 9809 more top variant IDs missing; see log file.
    ## --clump: 247175 clumps formed from 489805 top variants.
    ## Results written to EURcl500.2.clumped .

Vemos que se han creado 247175 aglomeraciones. Continuamos el análisis
entero.

``` r
clump_500.2 <- read.table("EURcl500.2.clumped", header = T, stringsAsFactors = F)
valid_clump_500.2 <- clump_500.2[,c("SNP")]

write.table(valid_clump_500.2,"EURcl500.2.valid.snp", quote = F, row.names = F, col.names = F)
```

``` bash
./plink --bfile EUR.QC --score Height.QC.transformed 3 4 12 header --q-score-range range_list SNP.value --extract EURcl500.2.valid.snp --out EURcl500.2
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl500.2.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --extract EURcl500.2.valid.snp
    ##   --out EURcl500.2
    ##   --q-score-range range_list SNP.value
    ##   --score Height.QC.transformed 3 4 12 header
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## --extract: 247175 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 247175 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 252442 lines skipped in --score file (252442 due to variant ID
    ## mismatch, 0 due to allele code mismatch); see EURcl500.2.nopred for details.
    ## --score: 247175 valid predictors loaded.
    ## Warning: 252443 lines skipped in --q-score-range data file.
    ## --score: 7 ranges processed.
    ## Results written to EURcl500.2.*.profile.

``` r
p.val.th <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
prs.result.cl500.2 <- NULL
for(i in p.val.th){
  prs.cl500.2 <- read.table(paste0("EURcl500.2.",i,".profile"), header = T)
  feno_prs.cl500.2 <- merge(feno, prs.cl500.2[,c("FID","IID","SCORE")], by = c("FID", "IID"))
  model_prs.cl500.2 <- lm(Height~., feno_prs.cl500.2[,!colnames(feno_prs.cl500.2) %in% c("FID","IID")])
  r2_model_prs.cl500.2 <- summary(model_prs.cl500.2)$r.squared
  prs_r2.cl500.2 <- r2_model_prs.cl500.2-r2.prs
  prs.result.cl500.2 <- rbind(prs.result.cl500.2, data.frame(Threshold=i, R2=prs_r2.cl500.2))
}
```

``` r
prs.result.cl500.2[which.max(prs.result.cl500.2$R2),]
```

    ##   Threshold        R2
    ## 6       0.4 0.1792011

El PRS mejor ajustado en este caso se encuentra en el umbral 0.4 del
p-valor y explica un 0.179201 la variación fenotípica

#### 4.2.2.2 1250

``` bash
./plink --bfile EUR.QC --clump-p1 1 --clump-r2 0.2 --clump-kb 1250 --clump Height.QC.transformed --clump-snp-field SNP --clump-field P --out EURcl1250
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl1250.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --clump Height.QC.transformed
    ##   --clump-field P
    ##   --clump-kb 1250
    ##   --clump-p1 1
    ##   --clump-r2 0.2
    ##   --clump-snp-field SNP
    ##   --out EURcl1250
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 489805 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 'rs1076829' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3129818' is missing from the main dataset, and is a top variant.
    ## Warning: 'rs3118359' is missing from the main dataset, and is a top variant.
    ## 9809 more top variant IDs missing; see log file.
    ## --clump: 244358 clumps formed from 489805 top variants.
    ## Results written to EURcl1250.clumped .

Vemos que se han creado 244358 aglomeraciones. Continuamos el análisis
entero.

``` r
clump_1250 <- read.table("EURcl1250.clumped", header = T, stringsAsFactors = F)
valid_clump_1250 <- clump_1250[,c("SNP")]

write.table(valid_clump_1250,"EURcl1250.valid.snp", quote = F, row.names = F, col.names = F)
```

``` bash
./plink --bfile EUR.QC --score Height.QC.transformed 3 4 12 header --q-score-range range_list SNP.value --extract EURcl1250.valid.snp --out EURcl1250
```

    ## PLINK v1.90b6.26 64-bit (2 Apr 2022)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to EURcl1250.log.
    ## Options in effect:
    ##   --bfile EUR.QC
    ##   --extract EURcl1250.valid.snp
    ##   --out EURcl1250
    ##   --q-score-range range_list SNP.value
    ##   --score Height.QC.transformed 3 4 12 header
    ## 
    ## 16384 MB RAM detected; reserving 8192 MB for main workspace.
    ## 489805 variants loaded from .bim file.
    ## 483 people (232 males, 251 females) loaded from .fam.
    ## --extract: 244358 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 483 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is exactly 1.
    ## 244358 variants and 483 people pass filters and QC.
    ## Note: No phenotypes present.
    ## Warning: 255259 lines skipped in --score file (255259 due to variant ID
    ## mismatch, 0 due to allele code mismatch); see EURcl1250.nopred for details.
    ## --score: 244358 valid predictors loaded.
    ## Warning: 255260 lines skipped in --q-score-range data file.
    ## --score: 7 ranges processed.
    ## Results written to EURcl1250.*.profile.

``` r
p.val.th <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
prs.result.cl1250<- NULL
for(i in p.val.th){
  prs.cl1250 <- read.table(paste0("EURcl1250.",i,".profile"), header = T)
  feno_prs.cl1250 <- merge(feno, prs.cl1250[,c("FID","IID","SCORE")], by = c("FID", "IID"))
  model_prs.cl1250 <- lm(Height~., feno_prs.cl1250[,!colnames(feno_prs.cl1250) %in% c("FID","IID")])
  r2_model_prs.cl1250 <- summary(model_prs.cl1250)$r.squared
  prs_r2.cl1250 <- r2_model_prs.cl1250-r2.prs
  prs.result.cl1250 <- rbind(prs.result.cl1250, data.frame(Threshold=i, R2=prs_r2.cl1250))
}
```

``` r
prs.result.cl1250[which.max(prs.result.cl1250$R2),]
```

    ##   Threshold       R2
    ## 6       0.4 0.187588

El PRS mejor ajustado en este caso se encuentra en el umbral 0.4 del
p-valor y explica un 0.187588 la variación fenotípica

### 4.2.3 GRAFICOS

En primer lugar comparamos los gráficos obtenidos para $R_c^2$ = 0.1 con
distintos $w_c$ {500, 1000, 2500} además del modelo inicial con $w_c$
{250}.

![](AnalisisPRS_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

En segundo lugar comparamos los gráficos obtenidos para $R_c^2$ = 0.2
con distintos $w_c$ {250, 500, 1250}.

![](AnalisisPRS_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

En tercer lugar realizamos los boxplots, comparando el número de SNPs
utilizados por cada sujeto en el cálculo de los PRS según el $w_c$
utilizado.

![](AnalisisPRS_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

![](AnalisisPRS_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

Por último comparamos los gráficos de los PRS mejor ajustados entre
modelos de distinto umbral $R_c^2$ pero con mismo $w_c$ {250, 500}.

![](AnalisisPRS_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

![](AnalisisPRS_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->
