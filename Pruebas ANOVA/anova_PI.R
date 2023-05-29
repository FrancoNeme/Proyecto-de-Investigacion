###################################
#                                 #
# EVALUACION DE EFECTOS SP Y DAP  #
# SOBRE CRECIMIENTO DE LEVADURAS  #
# PROYECTO DE INVESTIGACION       #
#                                 #
###################################

## CARGA DE LIBRERIAS ----
# = = = = = = = = = = = = =

library(ggplot2)
library(effects)
library(ggeffects)
library(gvlma)
library(ggplotify)
library("ggpubr")
library(nlme)
library(multcomp)
library(emmeans)
library(car)
library(tidyverse)
library("nortest")



## CARGA DE DATOS ----
# = = = = = = = = = = = 

# ___ Carga per se 

archivo_datos = read.csv('datos_gral.csv', dec = ",")

archivo_datos <- archivo_datos[order(archivo_datos$t),]

# ___ Casteo de variables

# Convertir la columna "Temperatura" a tipo factor
archivo_datos$sp <- as.factor(archivo_datos$sp)

# Convertir la columna "N0" a tipo factor
archivo_datos$DAP <- as.factor(archivo_datos$DAP)

# Convertir la columna "Tiempo" a tipo factor
archivo_datos$t <- as.factor(as.character(archivo_datos$t))

# Convertir la columna "C_B" a tipo numérico

archivo_datos$C_B <- as.numeric(archivo_datos$C_B)
#archivo_datos$C_B <- log10(as.numeric(archivo_datos$C_B))

# Convertir la columna "Corrida_T" a tipo factor
archivo_datos$Corrida_sp <- as.factor(archivo_datos$Corrida_sp)

# Convertir la columna "Corrida_N0" a tipo factor
archivo_datos$Corrida_DAP <- as.factor(archivo_datos$Corrida_DAP)



## POTENCIALES MODELOS DE ESTRUCTURA DE MATRIZ VARIANZA-COVARIANZA ----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# ___ Modelo 0 -----

mlm.modelo.000_C_B_REML <- gls(C_B ~ 1 + sp + DAP + t + sp:DAP + sp:t + DAP:t + sp:DAP:t
                               ,method = "REML"
                               ,na.action = na.omit
                               ,data = archivo_datos)

# Tabla de resumen del modelo

resumen_modelo0 <- summary(mlm.modelo.000_C_B_REML)

print('')
print('---------------------------------------- MODELO 0 ----------------------------------------')
print('------------------------------------------------------------------------------------------')

# Imprimir medidas de ajuste

cat(paste("N:", length(resid(mlm.modelo.000_C_B_REML)), "\n"))
cat(paste("AIC:", AIC(mlm.modelo.000_C_B_REML), "\n"))
cat(paste("BIC:", BIC(mlm.modelo.000_C_B_REML), "\n"))
cat(paste("logLik:", logLik(mlm.modelo.000_C_B_REML), "\n"))
cat(paste("sigma:", summary(mlm.modelo.000_C_B_REML)$sigma, "\n"))



# ___ Modelo 1 -----

mlm.modelo.001_C_B_REML <- gls(C_B ~ 1 + sp + DAP + t + sp:DAP + sp:t + DAP:t + sp:DAP:t
                               ,correlation = corCompSymm(form = ~1 | Corrida_sp / Corrida_DAP)
                               ,method = "REML"
                               ,na.action = na.omit
                               ,data = archivo_datos)

# Tabla de resumen del modelo

resumen_modelo1 <- summary(mlm.modelo.001_C_B_REML)

print('')
print('---------------------------------------- MODELO 1 ----------------------------------------')
print('------------------------------------------------------------------------------------------')

# Imprimir medidas de ajuste

cat(paste("N:", length(resid(mlm.modelo.001_C_B_REML)), "\n"))
cat(paste("AIC:", AIC(mlm.modelo.001_C_B_REML), "\n"))
cat(paste("BIC:", BIC(mlm.modelo.001_C_B_REML), "\n"))
cat(paste("logLik:", logLik(mlm.modelo.001_C_B_REML), "\n"))
cat(paste("sigma:", summary(mlm.modelo.001_C_B_REML)$sigma, "\n"))

# Matriz de varianzas-covarianzas

VarCovar_M1 = getVarCov(mlm.modelo.001_C_B_REML)
print(VarCovar_M1)



# ___ Modelo 2 ----

mlm.modelo.002_C_B_REML <- gls(C_B ~ 1 + sp + DAP + t + sp:DAP + sp:t + DAP:t + sp:DAP:t
                               ,correlation = corAR1(form = ~as.integer(as.numeric(t)) | Corrida_sp / Corrida_DAP)
                               ,method = "REML"
                               ,na.action = na.omit
                               ,data = archivo_datos)


# Tabla de resumen del modelo

resumen_modelo2 <- summary(mlm.modelo.002_C_B_REML)

print('')
print('---------------------------------------- MODELO 2 ----------------------------------------')
print('------------------------------------------------------------------------------------------')

# Imprimir medidas de ajuste

cat(paste("N:", length(resid(mlm.modelo.002_C_B_REML)), "\n"))
cat(paste("AIC:", AIC(mlm.modelo.002_C_B_REML), "\n"))
cat(paste("BIC:", BIC(mlm.modelo.002_C_B_REML), "\n"))
cat(paste("logLik:", logLik(mlm.modelo.002_C_B_REML), "\n"))
cat(paste("sigma:", summary(mlm.modelo.002_C_B_REML)$sigma, "\n"))

# Matriz de varianzas-covarianzas

VarCovar_M2 = getVarCov(mlm.modelo.002_C_B_REML)
print(VarCovar_M2)



# ___ Modelo 3 ----

mlm.modelo.003_C_B_REML <- gls(C_B ~ 1 + sp + DAP + t + sp:DAP + sp:t+ DAP:t+ sp:DAP:t
                               ,weights = varComb(varIdent(form = ~1 | t))
                               ,correlation = corAR1(form = ~as.integer(as.numeric(t)) | Corrida_sp / Corrida_DAP)
                               ,method = "REML"
                               ,na.action = na.omit
                               ,data = archivo_datos
                               ,control = glsControl(msMaxIter = 3000))


# Tabla de resumen del modelo

resumen_modelo3 <- summary(mlm.modelo.003_C_B_REML)

print('')
print('---------------------------------------- MODELO 3 ----------------------------------------')
print('------------------------------------------------------------------------------------------')

# Imprimir medidas de ajuste

cat(paste("N:", length(resid(mlm.modelo.003_C_B_REML)), "\n"))
cat(paste("AIC:", AIC(mlm.modelo.003_C_B_REML), "\n"))
cat(paste("BIC:", BIC(mlm.modelo.003_C_B_REML), "\n"))
cat(paste("logLik:", logLik(mlm.modelo.003_C_B_REML), "\n"))
cat(paste("sigma:", summary(mlm.modelo.003_C_B_REML)$sigma, "\n"))

# Matriz de varianzas-covarianzas

VarCovar_M3 <- getVarCov(mlm.modelo.003_C_B_REML)
print(VarCovar_M3)



# ___ Modelo 4 ----

mlm.modelo.004_C_B_REML <- gls(C_B ~ 1 + sp + DAP + t + sp:DAP + sp:t + DAP:t + sp:DAP:t
                               ,weights = varComb(varIdent(form = ~1 | t))
                               ,correlation = corSymm(form = ~as.integer(as.numeric(t)) | Corrida_sp / Corrida_DAP)
                               ,method = "REML"
                               ,na.action = na.omit
                               ,data = archivo_datos
                               ,control = glsControl(msMaxIter = 3000))

# Tabla de resumen del modelo

resumen_modelo4 <- summary(mlm.modelo.004_C_B_REML)

print('---------------------------------------- MODELO 4 ----------------------------------------')
print('------------------------------------------------------------------------------------------')

# Imprimir medidas de ajuste

cat(paste("N:", length(resid(mlm.modelo.004_C_B_REML)), "\n"))
cat(paste("AIC:", AIC(mlm.modelo.004_C_B_REML), "\n"))
cat(paste("BIC:", BIC(mlm.modelo.004_C_B_REML), "\n"))
cat(paste("logLik:", logLik(mlm.modelo.004_C_B_REML), "\n"))
cat(paste("sigma:", summary(mlm.modelo.004_C_B_REML)$sigma, "\n"))
# cat(paste("R cuadrado ajustado:", resumen_modelo4$r.squared, "\n"))

# Matriz de varianzas-covarianzas

VarCovar_M4 = getVarCov(mlm.modelo.004_C_B_REML)
print(VarCovar_M4)



## COMPARACION Y SELECCION DE MODELOS SEGUN AIC Y BIC ----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

comp_1_2 = anova(mlm.modelo.001_C_B_REML, mlm.modelo.002_C_B_REML)
print(comp_1_2)

comp_2_3 = anova(mlm.modelo.002_C_B_REML, mlm.modelo.003_C_B_REML)
print(comp_2_3)

comp_2_4 = anova(mlm.modelo.002_C_B_REML, mlm.modelo.004_C_B_REML)
print(comp_2_4)

comp_3_4 = anova(mlm.modelo.003_C_B_REML, mlm.modelo.004_C_B_REML)
print(comp_3_4)


# ___ Modelo Elegido ----

modelo_elegido = mlm.modelo.003_C_B_REML



## VERIFICACION DE SUPUESTOS ESTADISTICOS DEL ANOVA ----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# ___ Normalidad ----


# ______ Visualizacion ----

# QQPlot 

plot_norm <- qqnorm(modelo_elegido, abline = c(0,1))
ggsave("QQPlot.png", as.ggplot(plot_norm), width = 6, height = 6, dpi = 300, bg = "white")

# Histograma de residuos

png("hist_resid.png", width = 6, height = 6, units = 'in', res = 300)
hist(modelo_elegido$residuals, freq = F, border = "gray50", main = "Histograma de residuales")

lines(density(modelo_elegido$residuals), lwd = 2)

curve(dnorm(x, mean(modelo_elegido$residuals), sd(modelo_elegido$residuals)), lwd = 2, col = "blue", add = T)

legend("topleft", c("curva observada", "curva (normal) teórica"),
       lty = 1, lwd = 2, col = c("black", "blue"), bty = "n",
       cex = 0.8)  

dev.off()

# ______Test de Shapiro-Wilks ----

shapiro = shapiro.test(modelo_elegido$residuals)

# ______Test de Lilliefords ----

lillie = lillie.test(x = modelo_elegido$residuals )



# ___ Homocedasticidad ----

# Grafica de residuos vs predichos

plot_residuals <- plot(modelo_elegido, abline = c(0))
ggsave("residuos_vs_predichos.png", as.ggplot(plot_residuals), width = 6, height = 6, dpi = 300, bg = "white")



## ANOVA DEL MODELO SELECCIONADO ----
# = = = = = = = = = = = = = = = = = = 

anova_modelo_elegido <- anova(modelo_elegido, test = "F")



# PRUEBAS A POSTERIORI ----
# = = = = = = = = = = = = =

# ___ Tukey DAP:sp

# Calcular medias ajustadas y prueba de tukey, segun DAP y sp 

medias <- emmeans(modelo_elegido, ~ DAP + sp)
tukey <- cld(medias, alpha = 0.05, Letters = letters)

# Convertir medias y letras en un data.frame

medias_df <- as.data.frame(medias)
tukey_df <- as.data.frame(tukey)

# Graficar barras y agregar letras de tukey

ggplot(medias_df, aes(x = factor(DAP), y = emmean, fill = factor(sp))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(data = tukey_df, aes(x = factor(DAP), y = emmean, label = .group),
            position = position_dodge(width = 0.9), vjust = 1.5, size = 4) +
  scale_fill_discrete(name = "sp", labels = c("c", "e", "u")) +
  labs(x = "DAP", y = "C_B") +
  theme_bw()

ggsave("tukey_general_labo_DAPsp.png", plot = last_plot(), dpi = 300, width = 6, height = 4, units = "in")



# ___ Tukey sp:t

# Calcular medias ajustadas y prueba de tukey, segun t y sp

medias <- emmeans(modelo_elegido, ~ t + sp)
tukey <- cld(medias, alpha = 0.05, Letters = letters)

# Convertir medias y letras en un data.frame

medias_df <- as.data.frame(medias)
tukey_df <- as.data.frame(tukey)

# Graficar barras y agregar letras de tukey

ggplot(medias_df, aes(x = factor(t), y = emmean, fill = factor(sp))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(data = tukey_df, aes(x = factor(t), y = emmean, label = .group),
            position = position_dodge(width = 0.9), vjust = 1.5, size = 4) +
  scale_fill_discrete(name = "sp", labels = c("c", "e", "u")) +
  labs(x = "t", y = "C_B") +
  theme_bw()

ggsave("tukey_general_labo_tsp.png", plot = last_plot(), dpi = 300, width = 6, height = 4, units = "in")



# ___ Gráfica de interacción particionando sp 

png("C_B_particion_sp.png", width = 6, height = 5, units = 'in', res = 300)#, dpi = 300)

interaction.plot(x.factor = archivo_datos$t,
                 trace.factor = archivo_datos$sp,
                 response = archivo_datos$C_B,
                 type = "l",
                 xlab = "t",
                 ylab = "C_B",
                 col = c("blue", "red", "green"),
                 lty = 1,
                 lwd = 2,
                 trace.label = "sp")

dev.off()



# ___ Gráfica de interacción particionando DAP

png("C_B_particion_DAP.png", width = 6, height = 5, units = 'in', res = 300)#, dpi = 300)

interaction.plot(x.factor = archivo_datos$t,
                 trace.factor = archivo_datos$DAP,
                 response = archivo_datos$C_B,
                 type = "l",
                 xlab = "t",
                 ylab = "C_B",
                 col = c("blue", "red", "green"),
                 lty = 1,
                 lwd = 2,
                 trace.label = "DAP")

dev.off()



# ___ Grafico3 de interacción 

interaccion_gral <- plot(effect("sp*DAP*t",modelo_elegido,confidence.level=0.95))

ggp <- as.ggplot(interaccion_gral)

ggsave("interaccion_gral.png", ggp, width = 8, height = 6, dpi = 300, bg = "white")



## SALIDAS .TXT ----
# = = = = = = = = =

# ___ Salida modelo_lineal.txt ----


output <- c(
  "MODELO ELEGIDO:",
  "---------------",
  "",
  as.character(modelo_elegido$call),
  "",
  "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =",
  "",
  "NORMALIDAD:",
  "-----------",
  "",
  "Test de Lilliefors:",
  "-------------------",
  "",
  paste("Estadístico D =", lillie$statistic[["D"]]),
  paste("p-valor =", lillie$p.value),
  "",
  "Test de Shapiro-Wilks:",
  "----------------------",
  "",
  paste("Estadístico W =", shapiro$statistic[["W"]]),
  paste("p-valor =", shapiro$p.value)

)

# Especificar el nombre del archivo de salida
nombre_archivo <- "modelo_lineal.txt"

# Escribir el vector 'output' en el archivo de salida
writeLines(output, nombre_archivo)



# ___ Salidas anova y tukey .csv ----

write.csv(anova_modelo_elegido, 'anova.csv')
write.csv(tukey,'tukey.csv', row.names = FALSE)



# ___ Tabla comparacion de modelos ----

# Especificar el nombre del archivo de salida

nombre_archivo <- "comp_modelos.csv"

# Abrir el archivo en modo de escritura
archivo <- file(nombre_archivo, "w")

writeLines("", con = archivo)
writeLines("COMPARACION 1_2", con = archivo)
#writeLines("", con = archivo)

# Escribir la tabla 1 en el archivo CSV
write.csv(comp_1_2, file = archivo)
writeLines("", con = archivo)

writeLines("COMPARACION 2_3", con = archivo)
#writeLines("", con = archivo)

# Escribir la tabla 2 en el archivo CSV
write.csv(comp_2_3, file = archivo)
writeLines("", con = archivo)

# Escribir la tabla 3 en el archivo CSV
writeLines("COMPARACION 2_4", con = archivo)

write.csv(comp_2_4, file = archivo)
writeLines("", con = archivo)

# Escribir la tabla 4 en el archivo CSV
writeLines("COMPARACION 3_4", con = archivo)

write.csv(comp_3_4, file = archivo)

close(archivo)