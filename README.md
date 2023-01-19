# AnalisisPrs
Análisis de PRS y comparación de modelos con distintos parámetros.

Las enfermedades complejas humanas son hoy día uno de los temas principales de investigación, por lo que se ha profundizado en ciertas líneas de estudio como es el del posible efecto que surge de la combinación de múltiples genes y factores ambientales. Por ello, los estudios de asociación de genoma completo (en inglés, GWAS) han adquirido notoriedad en la última década, apareciendo en 2010 el primer estudio publicado. Dichos estudios pretenden identificar variantes genéticas, o polimorfismos de nucleótido único (en inglés, SNP), asociadas a un fenotipo o enfermedad en particular, para más tarde poder identificar otras variantes cercanas que contribuyen al desarrollo del rasgo en cuestión, tal como explica el National Human Genome Research Institute (NIH).

Los SNP son variaciones que afectan a una sola base en una región precisa del genoma. Estas suelen ocurrir con una frecuencia de uno por cada 1000 pares de bases y deben darse en el 1% de la población para ser considerados como tal, y no simples mutaciones puntuales.

Dichas variantes genéticas pueden ser combinadas en Polygenic Risk Scores (PRS) para dar una puntuación que registra la vulnerabilidad de un individuo de enfermar. Estas puntuaciones, que podrían verse como una recopilación de las variantes presentes en un individuo, permiten calcular el riesgo de padecer una enfermedad o desarrollar un fenotipo(fuente NIH). 

Hoy en día esta línea de estudio ha demostrado ser valiosa en el campo de la investigación, quedando todavía por establecerse su uso clínico.

Mediante este trabajo no se pretende innovar, si no crear varios modelos de PRS cambiando los distintos parámetros y compararlos a posteriori. Los parámetros no serán establecidos de forma arbitraria y se justificará su elección, partiendo de la lectura científica y del estado del arte en la materia. Para realizar el análisis se utilizará la herramienta de software de código abierto para análisis de genoma PLINK-1.9[Shaun Purcell, Christopher Chang, URL: www.cog-genomics.org/plink/1.9/) (Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4.)]

Esta herramienta lleva implementado el método estándar C+T (Clumping and Thresholding), lo que implica que los cálculos de PRSs se basan en un subconjunto de SNPs (Clumps) que tienen un p-valor superior al límite específico (Thresholding) en el GWAS.

Este estudio de control de calidad y cálculo de PRS se ha realizado a partir del código presente en el tutorial de este enlace:

https://choishingwan.github.io/PRS-Tutorial/

Dicho tutorial está asociado al articulo "P.F. Tutorial: a guide to performing polygenic risk score analyses." (Choi, S.W., Mak, T.S. & O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score analyses. Nat Protoc (2020). https://doi.org/10.1038/s41596-020-0353-1)

Como herramientas se han utilizado dos softwares:

-RStudio (RStudio Team (2020). RStudio: Integrated Development for R. RStudio, PBC, Boston, MA URL http://www.rstudio.com/.)

-Plink 1.9 (Shaun Purcell, Christopher Chang, URL: www.cog-genomics.org/plink/1.9/) (Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4.)

El código ha sido utilizado y adaptado según conveniencia, siendo los autores los mencionados anteriormente.

