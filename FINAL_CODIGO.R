# ============================================================
# TFM - ANÁLISIS TRANSCRIPTÓMICO DE RESPUESTA A NEOADYUVANCIA
# CÁNCER DE MAMA 
#
# OBJETIVO GENERAL
# Evaluar la asociación entre los perfiles transcriptómicos y la
# respuesta al tratamiento neoadyuvante en pacientes con cáncer
# de mama, con el fin de identificar biomarcadores candidatos y
# comprender los mecanismos moleculares implicados.
#
# OBJETIVOS ESPECÍFICOS
# 1. Analizar la actividad de rutas moleculares según respuesta.
# 2. Construir modelos de machine learning a partir de rutas.
# 3. Identificar firmas moleculares y redes transcripcionales.
# 4. Identificar biomarcadores transcriptómicos candidatos.
#
# NOTAS METODOLÓGICAS IMPORTANTES
# - Se utiliza TODA la cohorte con RNA-seq y respuesta conocida.
# - Se realizan DESeq2 y limma-voom como análisis complementarios.
# # - El filtrado de expresión combina filterByExpr() y un umbral
#   de expresión mínima (CPM).
# - El análisis diferencial y GSVA se realizan sobre todos los genes
#   que superan el filtro de expresión.
# - Los ~8.000 genes más variables se utilizan únicamente para
#   análisis exploratorios como el PCA.
# - El ML utiliza las puntuaciones de rutas, no miles de genes.
# - La validación es cruzada y repetida por el reducido tamaño muestral.
# - Los biomarcadores se consideran CANDIDATOS exploratorios y no
#   biomarcadores clínicamente validados.
# ============================================================


# ============================================================
# 0. INSTALACIÓN Y CARGA DE PAQUETES
# ============================================================
cran_packages <- c(
  "readxl",
  "tidyverse",
  "janitor",
  "caret",
  "pROC",
  "ranger",
  "glmnet",
  "e1071",
  "pheatmap",
  "corrplot",
  "igraph",
  "ggrepel",
  "msigdbr",
  "kernlab"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

# Instalar CRAN
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

# Paquetes Bioconductor
bioc_packages <- c(
  "edgeR",
  "limma",
  "DESeq2",
  "GSVA"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

# ============================================================
# 1. CARGAR PAQUETES
# ============================================================

library(ggplot2)
library(dplyr)

library(readxl)
library(tidyverse)
library(janitor)

library(edgeR)
library(limma)
library(DESeq2)
library(GSVA)

library(msigdbr)

library(caret)
library(pROC)
library(ranger)
library(glmnet)
library(e1071)

library(pheatmap)
library(corrplot)
library(igraph)
library(ggrepel)

library(kernlab)
# ============================================================
# 2. SEMILLA Y OPCIONES
# ============================================================

set.seed(12345)
options(stringsAsFactors = FALSE)
# ============================================================
# 3. RUTAS DE ARCHIVOS
# ============================================================

counts_file <- file.path("datos", "RNAseq-rawcounts.tsv.gz")

clinical_file <- file.path("datos", "Supplementary Tables 1-6.xlsx")

outdir <- "resultados_TFM"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

figdir <- file.path(outdir, "figuras")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

tabdir <- file.path(outdir, "tablas")
dir.create(tabdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 4. PARÁMETROS PRINCIPALES DEL ANÁLISIS
# ============================================================

# Expresión mínima:
# al menos 1 CPM en >= 20% de las muestras.
# filterByExpr() para que el filtrado tenga en cuenta el diseño/grupos.

min_cpm <- 1

min_prop_samples <- 0.20

# Objetivo aproximado de genes variables utilizados para PCA.
target_genes <- 8000

fdr_threshold <- 0.05
log2fc_threshold <- 0.5

set.seed(12345)
# ============================================================
# 5. IMPORTAR RNA-seq
# ============================================================

cat("\n============================================\n")
cat("IMPORTANDO MATRIZ RNA-seq\n")
cat("============================================\n")

raw_counts <- read.delim(
  gzfile(counts_file),
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

raw_counts <- as.matrix(raw_counts)

mode(raw_counts) <- "numeric"

cat(
  "Genes:",
  nrow(raw_counts),
  "\n"
)

cat(
  "Muestras:",
  ncol(raw_counts),
  "\n"
)

sample_ids <- colnames(raw_counts)

# Comprobar si existen genes duplicados
duplicated_genes <- duplicated(rownames(raw_counts))

cat(
  "Genes duplicados:",
  sum(duplicated_genes),
  "\n"
)

# Si hay IDs duplicados, sumar counts de los mismos genes
if (any(duplicated_genes)) {
  
  raw_counts <- rowsum(
    raw_counts,
    group = rownames(raw_counts)
  )
  
}

cat(
  "Genes después de resolver duplicados:",
  nrow(raw_counts),
  "\n"
)


# ============================================================
# 6. IMPORTAR DATOS CLÍNICOS
# ============================================================

cat("\n============================================\n")
cat("IMPORTANDO DATOS CLÍNICOS\n")
cat("============================================\n")

print(
  excel_sheets(clinical_file)
)

clinical <- read_excel(
  clinical_file,
  sheet = "Supplementary Table 1"
) %>%
  clean_names()


# Comprobar variables
cat("\nVariables clínicas disponibles:\n")
print(names(clinical))


# ============================================================
# 7. PREPARAR VARIABLES CLÍNICAS
# ============================================================

clinical <- clinical %>%
  mutate(
    
    donor_id = as.character(donor_id),
    
    age = suppressWarnings(
      as.numeric(age)
    ),
    
    rcb_score = suppressWarnings(
      as.numeric(rcb_score)
    ),
    
    response = case_when(
      p_cr_rd == "pCR" ~ "pCR",
      p_cr_rd == "RD" ~ "RD",
      TRUE ~ NA_character_
    ),
    
    response = factor(
      response,
      levels = c("RD", "pCR")
    ),
    
    t_stage = factor(t_stage),
    
    ln_status_at_diagnosis =
      factor(ln_status_at_diagnosis),
    
    er_status =
      factor(er_status),
    
    her2_status =
      factor(her2_status),
    
    grade_pre_nat =
      factor(grade_pre_nat),
    
    pam50 =
      factor(pam50),
    
    # Clasificación RCB para visualizar la gravedad de la enfermedad residual.
    # Umbrales estándar del RCB score: I <= 1.36; II > 1.36 y <= 3.28; III > 3.28.
    rcb_class = case_when(
      response == "pCR" ~ "pCR",
      !is.na(rcb_score) & rcb_score > 0 & rcb_score <= 1.36 ~ "RCB-I",
      !is.na(rcb_score) & rcb_score > 1.36 & rcb_score <= 3.28 ~ "RCB-II",
      !is.na(rcb_score) & rcb_score > 3.28 ~ "RCB-III",
      TRUE ~ NA_character_
    ),
    rcb_class = factor(
      rcb_class,
      levels = c("pCR", "RCB-I", "RCB-II", "RCB-III"),
      ordered = TRUE
    )
  )


# ============================================================
# 8. RESUMEN CLÍNICO INICIAL
# ============================================================

cat("\n============================================\n")
cat("RESUMEN CLÍNICO\n")
cat("============================================\n")

cat(
  "Pacientes con información clínica:",
  nrow(clinical),
  "\n"
)

cat(
  "pCR:",
  sum(clinical$response == "pCR", na.rm = TRUE),
  "\n"
)

cat(
  "RD:",
  sum(clinical$response == "RD", na.rm = TRUE),
  "\n"
)


# ============================================================
# 9. INTEGRACIÓN RNA-seq + CLÍNICA
# ============================================================

# Se utilizan todas las muestras disponibles que:
# 1. tengan RNA-seq
# 2. tengan donor_id coincidente
# 3. tengan respuesta pCR/RD conocida.

metadata <- clinical %>%
  filter(
    donor_id %in% sample_ids,
    !is.na(response)
  ) %>%
  arrange(
    match(donor_id, sample_ids)
  )


cat("\n============================================\n")
cat("COHORTE FINAL\n")
cat("============================================\n")

cat(
  "Muestras con RNA-seq + respuesta:",
  nrow(metadata),
  "\n"
)

cat(
  "pCR:",
  sum(metadata$response == "pCR"),
  "\n"
)

cat(
  "RD:",
  sum(metadata$response == "RD"),
  "\n"
)


# Comprobar IDs
common_ids <- intersect(
  sample_ids,
  metadata$donor_id
)

if (length(common_ids) == 0) {
  stop(
    "No hay IDs coincidentes entre RNA-seq y clínica."
  )
}

if (length(common_ids) != nrow(metadata)) {
  warning(
    "Hay pacientes clínicos con IDs duplicados o no coincidentes."
  )
}


# ============================================================
# 10. EXTRAER MATRIZ FINAL
# ============================================================

counts <- raw_counts[
  ,
  metadata$donor_id,
  drop = FALSE
]

# Asegurar correspondencia
stopifnot(
  identical(
    colnames(counts),
    metadata$donor_id
  )
)

cat(
  "Dimensiones matriz final:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "muestras\n"
)


# Guardar metadata
write.csv(
  metadata,
  file.path(
    tabdir,
    "01_metadata_cohorte_completa.csv"
  ),
  row.names = FALSE
)

# Comprobaciones de coherencia de la matriz definitiva

stopifnot(!anyDuplicated(colnames(counts)))
stopifnot(!anyDuplicated(metadata$donor_id))

stopifnot(
  all(colnames(counts) %in% metadata$donor_id)
)

cat("Dimensiones de la matriz definitiva:\n")
print(dim(counts))

cat("Número de pacientes por respuesta:\n")
print(table(metadata$response, useNA = "ifany"))

# ============================================================
# 11. CONTROL DE CALIDAD
# ============================================================

lib_sizes <- colSums(counts)

detected_genes <- colSums(
  counts > 0
)

qc_df <- metadata %>%
  mutate(
    library_size =
      lib_sizes[donor_id],
    
    detected_genes =
      detected_genes[donor_id]
  )


# Guardar
write.csv(
  qc_df,
  file.path(
    tabdir,
    "02_QC_summary.csv"
  ),
  row.names = FALSE
)


# Library size
p_lib <- ggplot(
  qc_df,
  aes(
    x = response,
    y = library_size
  )
) +
  
  geom_boxplot(
    aes(fill = response)
  ) +
  
  geom_jitter(
    width = 0.15
  ) +
  
  theme_bw() +
  
  labs(
    title = "Tamaño de librería según respuesta",
    x = "Respuesta",
    y = "Raw counts"
  )


ggsave(
  file.path(
    figdir,
    "01_QC_library_size.png"
  ),
  p_lib,
  width = 7,
  height = 5,
  dpi = 300
)


# Genes detectados
p_detected <- ggplot(
  qc_df,
  aes(
    x = response,
    y = detected_genes
  )
) +
  
  geom_boxplot(
    aes(fill = response)
  ) +
  
  geom_jitter(
    width = 0.15
  ) +
  
  theme_bw() +
  
  labs(
    title = "Genes detectados según respuesta",
    x = "Respuesta",
    y = "Genes con counts > 0"
  )


ggsave(
  file.path(
    figdir,
    "02_QC_detected_genes.png"
  ),
  p_detected,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# 11. CONTROL DE CALIDAD (QC) DE LAS MUESTRAS RNA-seq
# ============================================================

# Comprobar que los recuentos son válidos
stopifnot(
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == round(counts))
)

# ------------------------------------------------------------
# 11.1. Tamaño de las bibliotecas
# ------------------------------------------------------------

library_size <- colSums(counts)

# ------------------------------------------------------------
# 11.2. Número de genes detectados por muestra
# ------------------------------------------------------------
# Se considera detectado un gen con al menos 1 recuento.
# Este criterio es descriptivo para QC, no el filtro definitivo.

detected_genes <- colSums(counts >= 1)

# ------------------------------------------------------------
# 11.3. Crear tabla de QC
# ------------------------------------------------------------

qc_table <- data.frame(
  donor_id = colnames(counts),
  library_size = library_size,
  detected_genes = detected_genes,
  stringsAsFactors = FALSE
)

# Añadir respuesta clínica
qc_table <- qc_table %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(donor_id, response),
    by = "donor_id"
  )

# Ordenar por tamaño de biblioteca
qc_table <- qc_table %>%
  dplyr::arrange(library_size)

# Mostrar resumen
cat("\nRESUMEN DEL TAMAÑO DE LAS BIBLIOTECAS\n")
print(summary(qc_table$library_size))

cat("\nRESUMEN DE GENES DETECTADOS POR MUESTRA\n")
print(summary(qc_table$detected_genes))

cat("\nMuestras con menor tamaño de biblioteca:\n")
print(head(qc_table, 10))

cat("\nMuestras con mayor tamaño de biblioteca:\n")
print(tail(qc_table, 10))

# Guardar tabla completa
write.csv(
  qc_table,
  file.path(tabdir, "QC_tamano_bibliotecas_genes_detectados.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11.4. Gráfico del tamaño de las bibliotecas
# ------------------------------------------------------------

p_library <- ggplot(
  qc_table,
  aes(
    x = reorder(donor_id, library_size),
    y = library_size,
    fill = response
  )
) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title = "Tamaño de las bibliotecas RNA-seq",
    x = "Paciente",
    y = "Número total de recuentos",
    fill = "Respuesta"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

print(p_library)

ggsave(
  file.path(figdir, "QC_tamano_bibliotecas.png"),
  plot = p_library,
  width = 9,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 11.5. Gráfico del número de genes detectados
# ------------------------------------------------------------

p_detected <- ggplot(
  qc_table,
  aes(
    x = reorder(donor_id, detected_genes),
    y = detected_genes,
    fill = response
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Genes detectados por muestra",
    x = "Paciente",
    y = "Número de genes con >= 1 recuento",
    fill = "Respuesta"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

print(p_detected)

ggsave(
  file.path(figdir, "QC_genes_detectados.png"),
  plot = p_detected,
  width = 9,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 11.6. Detectar muestras potencialmente atípicas
# ------------------------------------------------------------
# Se utiliza el criterio exploratorio de 1.5 veces el rango intercuartílico (IQR). 

iqr_outlier <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr_value <- IQR(x, na.rm = TRUE)
  
  (x < q1 - 1.5 * iqr_value) |
    (x > q3 + 1.5 * iqr_value)
}

qc_table <- qc_table %>%
  dplyr::mutate(
    outlier_library = iqr_outlier(library_size),
    outlier_detected_genes = iqr_outlier(detected_genes),
    potential_outlier = outlier_library |
      outlier_detected_genes
  )

cat("\nMuestras potencialmente atípicas según el criterio IQR:\n")

print(
  qc_table %>%
    dplyr::filter(potential_outlier) %>%
    dplyr::select(
      donor_id,
      response,
      library_size,
      detected_genes,
      outlier_library,
      outlier_detected_genes
    )
)

cat(
  "\nNúmero de muestras potencialmente atípicas:",
  sum(qc_table$potential_outlier),
  "\n"
)

write.csv(
  qc_table,
  file.path(tabdir, "QC_tamano_bibliotecas_genes_detectados.csv"),
  row.names = FALSE
)


# Crear el gráfico del tamaño de bibliotecas
png(
  filename = file.path(figdir, "QC_boxplot_bibliotecas.png"),
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  qc_table$library_size / 1e6,
  main = "Tamaño de las bibliotecas RNA-seq",
  ylab = "Millones de recuentos",
  col = "lightblue"
)

dev.off()

# Crear el gráfico de genes detectados
png(
  filename = file.path(figdir, "QC_boxplot_genes_detectados.png"),
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  qc_table$detected_genes,
  main = "Genes detectados por muestra",
  ylab = "Número de genes",
  col = "lightgreen"
)

dev.off()

# Comprobar que se han guardado
file.exists(
  file.path(figdir, "QC_boxplot_bibliotecas.png")
)

file.exists(
  file.path(figdir, "QC_boxplot_genes_detectados.png")
)

# ============================================================
# REVISIÓN DETALLADA DE MUESTRAS POTENCIALMENTE ATÍPICAS
# ============================================================

# Límites IQR del tamaño de las bibliotecas
q1_library <- quantile(qc_table$library_size, 0.25, na.rm = TRUE)
q3_library <- quantile(qc_table$library_size, 0.75, na.rm = TRUE)
iqr_library <- IQR(qc_table$library_size, na.rm = TRUE)

lower_library <- q1_library - 1.5 * iqr_library
upper_library <- q3_library + 1.5 * iqr_library

# Límites IQR de genes detectados
q1_genes <- quantile(qc_table$detected_genes, 0.25, na.rm = TRUE)
q3_genes <- quantile(qc_table$detected_genes, 0.75, na.rm = TRUE)
iqr_genes <- IQR(qc_table$detected_genes, na.rm = TRUE)

lower_genes <- q1_genes - 1.5 * iqr_genes
upper_genes <- q3_genes + 1.5 * iqr_genes

cat("Límites de tamaño de biblioteca:\n")
print(c(
  inferior = lower_library,
  superior = upper_library
))

cat("\nLímites de genes detectados:\n")
print(c(
  inferior = lower_genes,
  superior = upper_genes
))

# Tabla detallada de las muestras señaladas
qc_outliers_detail <- qc_table %>%
  dplyr::filter(potential_outlier) %>%
  dplyr::mutate(
    motivo = dplyr::case_when(
      library_size < lower_library ~ "Biblioteca baja",
      library_size > upper_library ~ "Biblioteca alta",
      detected_genes < lower_genes ~ "Pocos genes detectados",
      detected_genes > upper_genes ~ "Muchos genes detectados",
      TRUE ~ "Revisar"
    )
  ) %>%
  dplyr::select(
    donor_id,
    response,
    library_size,
    detected_genes,
    motivo
  )

print(qc_outliers_detail)

cat("Resumen del tamaño de las bibliotecas:\n")
print(summary(qc_table$library_size))

cat("\nResumen de genes detectados:\n")
print(summary(qc_table$detected_genes))

cat("\nLímites IQR para genes detectados:\n")
print(c(
  inferior = quantile(qc_table$detected_genes, 0.25) -
    1.5 * IQR(qc_table$detected_genes),
  superior = quantile(qc_table$detected_genes, 0.75) +
    1.5 * IQR(qc_table$detected_genes)
))

# Boxplots de ambas métricas
par(mfrow = c(1, 2))

boxplot(
  qc_table$library_size / 1e6,
  main = "Tamaño de las bibliotecas",
  ylab = "Millones de recuentos",
  col = "lightblue"
)

boxplot(
  qc_table$detected_genes,
  main = "Genes detectados",
  ylab = "Número de genes",
  col = "lightgreen"
)

par(mfrow = c(1, 1))


# ============================================================
# 12. FILTRADO DE EXPRESIÓN
# ============================================================

cat("\n============================================\n")
cat("FILTRADO DE GENES\n")
cat("============================================\n")

# Comprobar que las columnas de expresión corresponden
# al orden de los pacientes del metadata
stopifnot(
  identical(colnames(counts), metadata$donor_id)
)

# Crear objeto DGEList con los recuentos originales
dge_raw <- edgeR::DGEList(
  counts = counts,
  group = metadata$response
)

# ------------------------------------------------------------
# 12.1. Filtrado estándar con filterByExpr()
# ------------------------------------------------------------

keep_filterByExpr <- edgeR::filterByExpr(
  dge_raw,
  group = metadata$response
)

# ------------------------------------------------------------
# 12.2. Filtrado adicional por CPM
# Criterio: >= 1 CPM en al menos el 20% de las muestras
# ------------------------------------------------------------

cpm_raw <- edgeR::cpm(
  dge_raw,
  log = FALSE
)

min_samples <- max(
  2,
  ceiling(min_prop_samples * ncol(counts))
)

keep_cpm <- rowSums(
  cpm_raw >= min_cpm
) >= min_samples

# ------------------------------------------------------------
# 12.3. Combinar ambos criterios
# Se conservan genes que superan AMBOS filtros
# ------------------------------------------------------------

keep_expression <- keep_filterByExpr & keep_cpm

cat("Genes iniciales:", nrow(counts), "\n")
cat("Genes tras filterByExpr:", sum(keep_filterByExpr), "\n")
cat("Genes que cumplen el criterio CPM:", sum(keep_cpm), "\n")
cat("Genes tras el filtrado combinado:", sum(keep_expression), "\n")

# Crear matriz de recuentos filtrada
counts_expr <- counts[
  keep_expression,
  ,
  drop = FALSE
]

cat(
  "Dimensiones de la matriz filtrada:",
  nrow(counts_expr), "genes x",
  ncol(counts_expr), "muestras\n"
)

# Comprobar que no se han perdido muestras
stopifnot(
  identical(colnames(counts_expr), colnames(counts))
)

# ============================================================
# 13. NORMALIZACIÓN TMM
# ============================================================

cat("\n============================================\n")
cat("NORMALIZACIÓN TMM\n")
cat("============================================\n")

# Crear objeto con los genes filtrados
dge_expr <- edgeR::DGEList(
  counts = counts_expr,
  group = metadata$response
)

# Calcular factores de normalización TMM
dge_expr <- edgeR::calcNormFactors(
  dge_expr,
  method = "TMM"
)

# Obtener expresión logCPM normalizada
logcpm_expr <- edgeR::cpm(
  dge_expr,
  log = TRUE,
  prior.count = 1
)

# Comprobaciones finales
stopifnot(
  all(is.finite(logcpm_expr)),
  identical(colnames(logcpm_expr), metadata$donor_id)
)

cat(
  "Dimensiones de logCPM:",
  nrow(logcpm_expr), "genes x",
  ncol(logcpm_expr), "muestras\n"
)

# ============================================================
# 14. FILTRADO POR VARIABILIDAD
# ============================================================

# Desviación estándar por gen sobre logCPM normalizado
gene_sd <- apply(
  logcpm_expr,
  1,
  sd,
  na.rm = TRUE
)

# Crear tabla de variabilidad
variability_df <- data.frame(
  gene = rownames(logcpm_expr),
  sd_logcpm = gene_sd
) %>%
  arrange(
    desc(sd_logcpm)
  )

# Número de genes después del filtro de expresión
n_expressed <- nrow(variability_df)

# Seleccionar aproximadamente los 8.000 genes más variables para el PCA
if (n_expressed > target_genes) {
  
  selected_genes <- variability_df %>%
    slice_head(
      n = target_genes
    ) %>%
    pull(gene)
  
} else {
  
  selected_genes <- variability_df %>%
    pull(gene)
  
}

cat(
  "Genes tras filtro de expresión:",
  n_expressed,
  "\n"
)

cat(
  "Genes seleccionados por variabilidad:",
  length(selected_genes),
  "\n"
)


# Guardar información de variabilidad
write.csv(
  variability_df,
  file.path(
    tabdir,
    "03_gene_filter_variability.csv"
  ),
  row.names = FALSE
)

# Resumen del filtrado
write.csv(
  data.frame(
    genes_iniciales = nrow(counts),
    genes_filterByExpr =
      sum(keep_filterByExpr),
    genes_CPM =
      sum(keep_cpm),
    genes_expression =
      sum(keep_expression),
    genes_variables_PCA =
      length(selected_genes)
  ),
  file.path(
    tabdir,
    "04_gene_filter_summary.csv"
  ),
  row.names = FALSE
)

# Matriz específica para análisis exploratorio
logcpm_variable <- logcpm_expr[
  selected_genes,
  ,
  drop = FALSE
]
# ============================================================
# 15. MATRICES PARA LOS DIFERENTES ANÁLISIS
# ============================================================

# ------------------------------------------------------------
# Expresión diferencial
# ------------------------------------------------------------
# Se utilizan TODOS los genes que han superado el filtro de expresión.

counts_de <- counts_expr

# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------
# Se utilizan únicamente los ~8.000 genes más variables.

counts_variable <- counts_expr[
  selected_genes,
  ,
  drop = FALSE
]

cat(
  "\nMatriz para expresión diferencial:",
  nrow(counts_de),
  "genes x",
  ncol(counts_de),
  "muestras\n"
)

cat(
  "Matriz para PCA:",
  nrow(counts_variable),
  "genes x",
  ncol(counts_variable),
  "muestras\n"
)

# ============================================================
# 16. PCA TRANSCRIPTÓMICO
# ============================================================

# Genes seleccionados únicamente para exploración de la estructura global de los datos.

pca <- prcomp(
  t(logcpm_variable),
  scale. = FALSE
)

percent_var <-
  100 *
  pca$sdev^2 /
  sum(pca$sdev^2)


pca_df <- as.data.frame(
  pca$x[, 1:2]
) %>%
  
  mutate(
    donor_id =
      rownames(pca$x),
    
    response =
      metadata$response[
        match(
          rownames(pca$x),
          metadata$donor_id
        )
      ]
  )


p_pca <- ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    color = response,
    label = donor_id
  )
) +
  
  geom_point(
    size = 3
  ) +
  
  geom_text_repel(
    size = 2.8,
    max.overlaps = 30
  ) +
  
  theme_bw() +
  
  labs(
    title = "PCA transcriptómico - cohorte completa",
    x = paste0(
      "PC1 (",
      round(percent_var[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(percent_var[2], 1),
      "%)"
    )
  )


ggsave(
  file.path(
    figdir,
    "03F_PCA_transcriptomico.png"
  ),
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 17. DISEÑO CLÍNICO PARA EXPRESIÓN DIFERENCIAL
# ============================================================

# El análisis principal compara pCR frente a RD en toda la cohorte, se utilizan las variables: 
#response + age + ER + HER2 + estado ganglionar.
# Si alguna covariable tiene demasiados datos ausentes, se elimina automáticamente del modelo ajustado.

metadata_de <- metadata %>%
  
  mutate(
    age = as.numeric(age)
  )


candidate_covariates <- c(
  "age",
  "er_status",
  "her2_status",
  "ln_status_at_diagnosis"
)


# Determinar covariables utilizables
usable_covariates <- c()

for (v in candidate_covariates) {
  
  if (v %in% names(metadata_de)) {
    
    complete_fraction <-
      mean(
        !is.na(
          metadata_de[[v]]
        )
      )
    
    if (complete_fraction >= 0.80) {
      usable_covariates <-
        c(
          usable_covariates,
          v
        )
    }
  }
}


cat(
  "\nCovariables utilizadas:",
  paste(
    usable_covariates,
    collapse = ", "
  ),
  "\n"
)


# Para el modelo ajustado utilizamos únicamente pacientes completos para las variables necesarias.
model_vars <- c(
  "response",
  usable_covariates
)

metadata_model <- metadata_de %>%
  select(
    donor_id,
    all_of(model_vars)
  ) %>%
  drop_na()


# Matriz de expresión para análisis diferencial
# Se utilizan TODOS los genes que han superado el filtro de expresión.

counts_model <- counts_de[
  ,
  metadata_model$donor_id,
  drop = FALSE
]

cat("Pacientes incluidos en el modelo diferencial:\n")
print(nrow(metadata_de))

cat("Distribución de respuesta en el modelo:\n")
print(table(metadata_de$response))

# ============================================================
# 18. LIMMA-VOOM: MODELO AJUSTADO
# ============================================================

formula_limma <- as.formula(
  paste(
    "~",
    paste(
      model_vars,
      collapse = " + "
    )
  )
)

design_limma <- model.matrix(
  formula_limma,
  data = metadata_model
)

cat(
  "\nDiseño limma:\n"
)

print(
  colnames(design_limma)
)


dge_model <- DGEList(
  counts = counts_model,
  group = metadata_model$response
)

dge_model <- calcNormFactors(
  dge_model,
  method = "TMM"
)


v <- voom(
  dge_model,
  design_limma,
  plot = TRUE
)


fit_limma <- lmFit(
  v,
  design_limma
)

fit_limma <- eBayes(
  fit_limma
)


response_coef <- grep(
  "^response",
  colnames(design_limma),
  value = TRUE
)[1]


if (is.na(response_coef)) {
  stop(
    "No se ha encontrado el coeficiente de respuesta en limma."
  )
}


limma_res <- topTable(
  fit_limma,
  coef = response_coef,
  number = Inf,
  sort.by = "P"
) %>%
  
  rownames_to_column(
    "ensembl_id"
  ) %>%
  
  mutate(
    method = "limma_voom"
  )


write.csv(
  limma_res,
  file.path(
    tabdir,
    "05_DEG_limma_voom.csv"
  ),
  row.names = FALSE
)
# ============================================================
# 19. DESEQ2
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "EJECUTANDO DESeq2\n"
)

cat(
  "============================================\n"
)


counts_model_integer <-
  round(
    counts_model
  )


formula_deseq <- formula_limma


dds <- DESeqDataSetFromMatrix(
  countData =
    counts_model_integer,
  colData =
    as.data.frame(
      metadata_model %>%
        select(
          -donor_id
        )
    ),
  design =
    formula_deseq
)


# Filtro adicional de DESeq2
keep_deseq <- rowSums(
  counts(dds) >= 10
) >= max(
  2,
  ceiling(
    0.10 *
      ncol(dds)
  )
)

dds <- dds[
  keep_deseq,
]


dds <- DESeq(
  dds
)


# Resultado pCR vs RD
res_deseq <- results(
  dds,
  contrast = c(
    "response",
    "pCR",
    "RD"
  ),
  alpha = 0.05
)


res_deseq_df <- as.data.frame(
  res_deseq
) %>%
  
  rownames_to_column(
    "ensembl_id"
  ) %>%
  
  arrange(
    padj
  )


write.csv(
  res_deseq_df,
  file.path(
    tabdir,
    "06_DEG_DESeq2.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 20. COMPARACIÓN DE DESeq2 Y LIMMA-VOOM
# ============================================================

# Genes significativos por limma-voom
limma_sig_genes <- limma_res %>%
  dplyr::filter(
    !is.na(adj.P.Val),
    adj.P.Val < 0.05
  ) %>%
  dplyr::pull(ensembl_id)

# Genes significativos por DESeq2
deseq_sig_genes <- res_deseq_df %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  dplyr::pull(ensembl_id)

# Intersección de genes significativos
overlap_deg <- intersect(
  limma_sig_genes,
  deseq_sig_genes
)

# Resumen de la comparación
deg_comparison <- data.frame(
  limma_FDR05 = length(limma_sig_genes),
  DESeq2_FDR05 = length(deseq_sig_genes),
  overlap = length(overlap_deg)
)

print(deg_comparison)

# Guardar resumen
write.csv(
  deg_comparison,
  file.path(tabdir, "07_DEG_methods_comparison.csv"),
  row.names = FALSE
)

# Guardar genes compartidos
write.csv(
  data.frame(ensembl_id = overlap_deg),
  file.path(tabdir, "08_DEG_overlap_limma_DESeq2.csv"),
  row.names = FALSE
)

# ============================================================
# 21. ANÁLISIS DIFERENCIAL COMPLEMENTARIO: pCR vs RCB-III
# ============================================================
# RCB-III representa enfermedad residual extensa. 
# Se compara con pCR para explorar diferencias entre los extremos de respuesta.
# El RCB score > 3.28 define RCB-III; las muestras pCR se conservan
# según la variable clínica original, aunque el score esté ausente.

# ------------------------------------------------------------
# 21.1. Seleccionar pacientes
# ------------------------------------------------------------

metadata_extreme <- metadata_de %>%
  dplyr::filter(
    response == "pCR" |
      (!is.na(rcb_score) & rcb_score > 3.28)
  ) %>%
  dplyr::mutate(
    response_extreme = dplyr::if_else(
      response == "pCR",
      "pCR",
      "RCB-III"
    ),
    response_extreme = factor(
      response_extreme,
      levels = c("RCB-III", "pCR")
    )
  )

# ------------------------------------------------------------
# 21.2. Seleccionar covariables disponibles
# ------------------------------------------------------------

covariates_extreme <- usable_covariates[
  vapply(
    usable_covariates,
    function(v) {
      mean(!is.na(metadata_extreme[[v]])) >= 0.80
    },
    logical(1)
  )
]

vars_extreme <- c(
  "response_extreme",
  covariates_extreme
)

metadata_extreme_model <- metadata_extreme %>%
  dplyr::select(
    donor_id,
    dplyr::all_of(vars_extreme)
  ) %>%
  tidyr::drop_na()

# Resumen de la cohorte
cat("\nComparación pCR vs RCB-III\n")
cat("Muestras seleccionadas:", nrow(metadata_extreme), "\n")

cat(
  "pCR:",
  sum(metadata_extreme$response_extreme == "pCR"),
  "\n"
)

cat(
  "RCB-III:",
  sum(metadata_extreme$response_extreme == "RCB-III"),
  "\n"
)

cat(
  "Muestras completas para el modelo:",
  nrow(metadata_extreme_model),
  "\n"
)

cat(
  "Covariables utilizadas:",
  paste(covariates_extreme, collapse = ", "),
  "\n"
)

print(table(metadata_extreme_model$response_extreme))

# Guardar resumen de la cohorte
write.csv(
  data.frame(
    n_total = nrow(metadata_extreme),
    n_pCR = sum(metadata_extreme$response_extreme == "pCR"),
    n_RCB_III = sum(metadata_extreme$response_extreme == "RCB-III"),
    n_model = nrow(metadata_extreme_model),
    covariables = paste(covariates_extreme, collapse = ";")
  ),
  file.path(tabdir, "07b_cohorte_pCR_vs_RCBIII.csv"),
  row.names = FALSE
)

if (
  sum(metadata_extreme_model$response_extreme == "pCR") < 5 ||
  sum(metadata_extreme_model$response_extreme == "RCB-III") < 5
) {
  warning(
    "Hay menos de 5 muestras en uno de los grupos. ",
    "Interpreta los resultados con especial cautela."
  )
}

# ------------------------------------------------------------
# 21.3. Preparar matriz y diseño
# ------------------------------------------------------------

counts_extreme <- counts_de[
  ,
  metadata_extreme_model$donor_id,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(counts_extreme),
    metadata_extreme_model$donor_id
  )
)

formula_extreme <- stats::as.formula(
  paste("~", paste(vars_extreme, collapse = " + "))
)

# ------------------------------------------------------------
# 21.4. Análisis diferencial con limma-voom
# ------------------------------------------------------------

design_extreme <- stats::model.matrix(
  formula_extreme,
  data = metadata_extreme_model
)

dge_extreme <- edgeR::DGEList(
  counts = counts_extreme
)

dge_extreme <- edgeR::calcNormFactors(
  dge_extreme,
  method = "TMM"
)

voom_extreme <- limma::voom(
  dge_extreme,
  design_extreme,
  plot = FALSE
)

fit_extreme <- limma::eBayes(
  limma::lmFit(
    voom_extreme,
    design_extreme
  )
)

coef_extreme <- grep(
  "^response_extremepCR$",
  colnames(design_extreme),
  value = TRUE
)

if (length(coef_extreme) != 1) {
  stop(
    "No se pudo identificar el coeficiente pCR ",
    "en el diseño pCR vs RCB-III."
  )
}

limma_extreme <- limma::topTable(
  fit_extreme,
  coef = coef_extreme,
  number = Inf,
  sort.by = "P"
) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  dplyr::mutate(
    method = "limma_voom_pCR_vs_RCBIII"
  )

write.csv(
  limma_extreme,
  file.path(tabdir, "07c_DEG_limma_pCR_vs_RCBIII.csv"),
  row.names = FALSE
)

# Genes priorizados por limma-voom
limma_extreme_prioritized <- limma_extreme %>%
  dplyr::filter(
    !is.na(adj.P.Val),
    adj.P.Val < 0.05,
    abs(logFC) >= 0.5
  )

write.csv(
  limma_extreme_prioritized,
  file.path(
    tabdir,
    "07h_DEG_limma_priorizados_pCR_vs_RCBIII_FDR05_LFC05.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 21.5. Análisis diferencial con DESeq2
# ------------------------------------------------------------

dds_extreme <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(counts_extreme),
  colData = as.data.frame(
    dplyr::select(
      metadata_extreme_model,
      -donor_id
    )
  ),
  design = formula_extreme
)

keep_deseq_extreme <- rowSums(
  DESeq2::counts(dds_extreme) >= 10
) >= max(
  2,
  ceiling(0.10 * ncol(dds_extreme))
)

dds_extreme <- dds_extreme[
  keep_deseq_extreme,
]

dds_extreme <- DESeq2::DESeq(
  dds_extreme,
  quiet = TRUE
)

res_extreme <- DESeq2::results(
  dds_extreme,
  contrast = c(
    "response_extreme",
    "pCR",
    "RCB-III"
  ),
  alpha = 0.05
)

res_extreme_df <- as.data.frame(res_extreme) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  dplyr::arrange(padj)

write.csv(
  res_extreme_df,
  file.path(tabdir, "07d_DEG_DESeq2_pCR_vs_RCBIII.csv"),
  row.names = FALSE
)


# ============================================================
# 22. VOLCANO PLOT DESeq2: pCR vs RCB-III
# ============================================================

fdr_threshold <- 0.05
log2fc_threshold <- 0.5

volcano_extreme <- res_extreme_df %>%
  dplyr::mutate(
    neg_log10_fdr = -log10(
      pmax(padj, .Machine$double.xmin)
    ),
    significance = dplyr::case_when(
      !is.na(padj) &
        padj < fdr_threshold &
        log2FoldChange >= log2fc_threshold ~ "Mayor en pCR",
      
      !is.na(padj) &
        padj < fdr_threshold &
        log2FoldChange <= -log2fc_threshold ~ "Mayor en RCB-III",
      
      TRUE ~ "No priorizado"
    )
  )

p_volcano_extreme <- ggplot2::ggplot(
  volcano_extreme,
  ggplot2::aes(
    x = log2FoldChange,
    y = neg_log10_fdr,
    color = significance
  )
) +
  ggplot2::geom_point(
    alpha = 0.65,
    size = 1.4
  ) +
  ggplot2::geom_vline(
    xintercept = c(-log2fc_threshold, log2fc_threshold),
    linetype = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(fdr_threshold),
    linetype = "dashed"
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "DESeq2: pCR frente a RCB-III",
    subtitle = "FDR < 0,05 y |log2FC| >= 0,5",
    x = "log2 fold change (pCR frente a RCB-III)",
    y = "-log10(FDR)",
    color = "Clasificación"
  )

print(p_volcano_extreme)

ggplot2::ggsave(
  file.path(figdir, "04b_Volcano_pCR_vs_RCBIII.png"),
  plot = p_volcano_extreme,
  width = 8,
  height = 6,
  dpi = 300
)

# Guardar genes priorizados
volcano_extreme_prioritized <- volcano_extreme %>%
  dplyr::filter(
    !is.na(padj),
    padj < fdr_threshold,
    abs(log2FoldChange) >= log2fc_threshold
  )

write.csv(
  volcano_extreme_prioritized,
  file.path(
    tabdir,
    "07e_DEG_priorizados_pCR_vs_RCBIII_FDR05_LFC05.csv"
  ),
  row.names = FALSE
)

cat(
  "Genes priorizados en pCR vs RCB-III:",
  nrow(volcano_extreme_prioritized),
  "\n"
)


# ============================================================
# 23. VOLCANO PLOT DESeq2: pCR vs RD
# ============================================================

# Umbrales de significación y tamaño del efecto
fdr_threshold <- 0.05
log2fc_threshold <- 0.5

# ------------------------------------------------------------
# 23.1. Crear el volcano plot
# ------------------------------------------------------------

volcano_deseq <- res_deseq_df %>%
  dplyr::mutate(
    neg_log10_fdr = -log10(
      pmax(padj, .Machine$double.xmin)
    ),
    significance = dplyr::case_when(
      !is.na(padj) &
        padj < fdr_threshold &
        log2FoldChange >= log2fc_threshold ~ "Mayor en pCR",
      
      !is.na(padj) &
        padj < fdr_threshold &
        log2FoldChange <= -log2fc_threshold ~ "Mayor en RD",
      
      TRUE ~ "No priorizado"
    )
  )

p_volcano <- ggplot2::ggplot(
  volcano_deseq,
  ggplot2::aes(
    x = log2FoldChange,
    y = neg_log10_fdr,
    color = significance
  )
) +
  ggplot2::geom_point(
    alpha = 0.65,
    size = 1.4
  ) +
  ggplot2::geom_vline(
    xintercept = c(-log2fc_threshold, log2fc_threshold),
    linetype = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(fdr_threshold),
    linetype = "dashed"
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "DESeq2: pCR frente a RD",
    subtitle = "Criterio de priorización: FDR < 0,05 y |log2FC| >= 0,5",
    x = "log2 fold change (pCR frente a RD)",
    y = "-log10(FDR)",
    color = "Clasificación"
  )

# Mostrar el gráfico
print(p_volcano)

# Guardar el gráfico
ggplot2::ggsave(
  filename = file.path(figdir, "04_Volcano_DESeq2.png"),
  plot = p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 23.2. Guardar genes priorizados por DESeq2
# ------------------------------------------------------------

volcano_prioritized <- volcano_deseq %>%
  dplyr::filter(
    !is.na(padj),
    padj < fdr_threshold,
    abs(log2FoldChange) >= log2fc_threshold
  )

write.csv(
  volcano_prioritized,
  file.path(
    tabdir,
    "07f_DEG_priorizados_pCR_vs_RD_FDR05_LFC05.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 23.3. Guardar genes priorizados por limma-voom
# ------------------------------------------------------------

limma_prioritized <- limma_res %>%
  dplyr::filter(
    !is.na(adj.P.Val),
    adj.P.Val < fdr_threshold,
    abs(logFC) >= log2fc_threshold
  )

write.csv(
  limma_prioritized,
  file.path(
    tabdir,
    "07g_DEG_limma_priorizados_pCR_vs_RD_FDR05_LFC05.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 23.4. Mostrar resumen de resultados
# ------------------------------------------------------------

cat(
  "\nGenes priorizados por DESeq2 (pCR vs RD):",
  nrow(volcano_prioritized),
  "\n"
)

cat(
  "Genes priorizados por limma-voom (pCR vs RD):",
  nrow(limma_prioritized),
  "\n"
)

cat("\nClasificación de los genes priorizados por DESeq2:\n")
print(table(volcano_prioritized$significance))

# ============================================================
# 24. TOP GENES DIFERENCIALMENTE EXPRESADOS CONSISTENTES
# ============================================================

top_limma <- limma_res %>%
  dplyr::filter(
    adj.P.Val < 0.05
  ) %>%
  dplyr::select(
    ensembl_id,
    logFC,
    P.Value,
    adj.P.Val
  )

top_deseq <- res_deseq_df %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  dplyr::select(
    ensembl_id,
    log2FoldChange,
    pvalue,
    padj
  )


consistent_genes <- dplyr::inner_join(
  top_limma,
  top_deseq,
  by = "ensembl_id"
)


write.csv(
  consistent_genes,
  file.path(
    tabdir,
    "09_consistent_DEGs_limma_DESeq2.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 25. PREPARAR EXPRESIÓN PARA GSVA
# ============================================================

# Para GSVA trabajamos con genes Ensembl.
# msigdbr permite obtener Hallmark con símbolos.
# Para hacer la correspondencia de forma robusta se obtiene también la equivalencia Ensembl -> símbolo desde msigdbr.

hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

# Comprobar columnas disponibles
cat(
  "\nColumnas msigdbr:\n"
)

print(
  names(hallmark)
)

# ============================================================
# 26. CONVERSIÓN ENSEMBL -> SYMBOL USANDO MSIGDB
# ============================================================
gene_annotation <- hallmark %>%
  dplyr::filter(
    !is.na(ensembl_gene),
    !is.na(gene_symbol),
    ensembl_gene != "",
    gene_symbol != ""
  ) %>%
  dplyr::distinct(
    ensembl_gene,
    .keep_all = TRUE
  )

clean_ensembl <- function(x) {
  sub("\\..*$", "", x)
}

expr_ensembl <- logcpm_expr

rownames(expr_ensembl) <-
  clean_ensembl(
    rownames(expr_ensembl)
  )

gene_annotation$ensembl_gene <-
  clean_ensembl(
    gene_annotation$ensembl_gene
  )

annotation_present <- gene_annotation %>%
  dplyr::filter(
    ensembl_gene %in% rownames(expr_ensembl)
  )

# ============================================================
# 27. PREPARAR CONJUNTOS DE GENES HALLMARK 
# ============================================================

# Asociar cada gen Ensembl con su símbolo
gene_symbols <- annotation_present$gene_symbol[
  match(
    rownames(expr_ensembl),
    annotation_present$ensembl_gene
  )
]

# Identificar los genes que tienen símbolo
valid_genes <- !is.na(gene_symbols) & gene_symbols != ""

# Crear la matriz conservando únicamente los genes anotados
expr_symbol <- expr_ensembl[valid_genes, , drop = FALSE]

# Asignar los símbolos como nombres de fila
rownames(expr_symbol) <- gene_symbols[valid_genes]

# Comprobar las dimensiones
cat("Genes en la matriz Ensembl:", nrow(expr_ensembl), "\n")
cat("Genes con símbolo asignado:", nrow(expr_symbol), "\n")
cat("Muestras:", ncol(expr_symbol), "\n")

# Comprobar los nombres
head(rownames(expr_symbol))

gene_sets <- split(
  hallmark$gene_symbol,
  hallmark$gs_name
)


gene_sets <- lapply(
  gene_sets,
  function(x) {
    
    intersect(
      unique(x),
      rownames(expr_symbol)
    )
    
  }
)


# Eliminar sets demasiado pequeños
gene_sets <- gene_sets[
  sapply(
    gene_sets,
    length
  ) >= 10
]


cat(
  "Número de Hallmarks:",
  length(gene_sets),
  "\n"
)


# ============================================================
# 28. CALCULAR GSVA
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "CALCULANDO GSVA\n"
)

cat(
  "============================================\n"
)


gsva_param <- gsvaParam(
  exprData =
    expr_symbol,
  
  geneSets =
    gene_sets,
  
  kcdf =
    "Gaussian"
)


pathway_scores <- gsva(
  gsva_param
)


pathway_scores <-
  as.matrix(
    pathway_scores
  )


cat(
  "Pathways calculados:",
  nrow(pathway_scores),
  "\n"
)

# ============================================================
# 29. GUARDAR TABLA DE PUNTUACIONES GSVA
# ============================================================

scores_df <- t(
  pathway_scores
) %>%
  
  as.data.frame() %>%
  
  tibble::rownames_to_column(
    "donor_id"
  ) %>%
  
  dplyr::left_join(
    metadata %>%
      dplyr::select(
        donor_id,
        response,
        age,
        t_stage,
        ln_status_at_diagnosis,
        er_status,
        her2_status,
        pam50
      ),
    by = "donor_id"
  )


write.csv(
  scores_df,
  file.path(
    tabdir,
    "10_GSVA_Hallmark_scores.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 30. COMPARACIÓN DE RUTAS pCR vs RD
# ============================================================

pathway_names <- rownames(
  pathway_scores
)


pathway_tests <- map_dfr(
  pathway_names,
  function(pathway) {
    
    tmp <- data.frame(
      
      score =
        pathway_scores[
          pathway,
          metadata$donor_id
        ],
      
      response =
        metadata$response[
          match(
            metadata$donor_id,
            colnames(pathway_scores)
          )
        ]
      
    )
    
    wt <- wilcox.test(
      score ~ response,
      data = tmp
    )
    
    median_pcr <- median(
      tmp$score[
        tmp$response == "pCR"
      ],
      na.rm = TRUE
    )
    
    median_rd <- median(
      tmp$score[
        tmp$response == "RD"
      ],
      na.rm = TRUE
    )
    
    tibble(
      pathway = pathway,
      median_pCR = median_pcr,
      median_RD = median_rd,
      difference_median =
        median_pcr -
        median_rd,
      p_value =
        wt$p.value
    )
  }
) %>%
  
  mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  
  arrange(
    FDR,
    p_value
  )


write.csv(
  pathway_tests,
  file.path(
    tabdir,
    "11_GSVA_pathway_tests.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 31. TOP PATHWAYS
# ============================================================

# Seleccion de rutas mas relevantes 
top_pathways <- pathway_tests %>%
  
  slice_head(
    n = min(
      12,
      nrow(pathway_tests)
    )
  ) %>%
  
  pull(
    pathway
  )


scores_long <- scores_df %>%
  
  pivot_longer(
    
    cols =
      all_of(
        setdiff(
          colnames(scores_df),
          c(
            "donor_id",
            "response",
            "age",
            "t_stage",
            "ln_status_at_diagnosis",
            "er_status",
            "her2_status",
            "pam50"
          )
        )
      ),
    
    names_to =
      "pathway",
    
    values_to =
      "score"
  )


# ============================================================
# 32. BOXPLOTS TOP PATHWAYS
# ============================================================

p_pathways <- scores_long %>%
  
  filter(
    pathway %in%
      top_pathways
  ) %>%
  
  ggplot(
    aes(
      x = response,
      y = score,
      fill = response
    )
  ) +
  
  geom_boxplot(
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15
  ) +
  
  facet_wrap(
    ~ pathway,
    scales = "free_y"
  ) +
  
  theme_bw() +
  
  labs(
    title =
      "Actividad de rutas Hallmark según respuesta",
    x =
      "Respuesta",
    y =
      "GSVA score"
  )


ggsave(
  file.path(
    figdir,
    "05_GSVA_top_pathways.png"
  ),
  p_pathways,
  width = 12,
  height = 8,
  dpi = 300
)

# ============================================================
# 33. HEATMAP TOP PATHWAYS
# ============================================================

heatmap_matrix <-
  scores_df %>%
  dplyr::select(
    donor_id,
    dplyr::all_of(top_pathways)
  ) %>%
  tibble::column_to_rownames("donor_id") %>%
  as.matrix()


# Pacientes en columnas
heatmap_matrix <- t(heatmap_matrix)


# Z-score por pathway
heatmap_scaled <- t(
  scale(
    t(heatmap_matrix)
  )
)


# Orden de gravedad: pCR, RCB-I, RCB-II, RCB-III
heatmap_order <- metadata %>%
  dplyr::mutate(
    rcb_order = match(
      as.character(rcb_class),
      c("pCR", "RCB-I", "RCB-II", "RCB-III")
    ),
    rcb_order = dplyr::if_else(
      is.na(rcb_order),
      99L,
      as.integer(rcb_order)
    )
  ) %>%
  dplyr::arrange(rcb_order, rcb_score, donor_id) %>%
  dplyr::pull(donor_id)

heatmap_order <- heatmap_order[heatmap_order %in% colnames(heatmap_matrix)]
heatmap_matrix <- heatmap_matrix[, heatmap_order, drop = FALSE]

annotation_col <- metadata %>%
  dplyr::select(donor_id, rcb_class, rcb_score) %>%
  dplyr::filter(donor_id %in% colnames(heatmap_matrix)) %>%
  dplyr::mutate(
    donor_id = factor(
      donor_id,
      levels = colnames(heatmap_matrix)
    )
  ) %>%
  dplyr::arrange(donor_id) %>%
  dplyr::mutate(
    donor_id = as.character(donor_id)
  ) %>%
  tibble::column_to_rownames("donor_id")
annotation_col <- annotation_col[colnames(heatmap_matrix), , drop = FALSE]


# ============================================================
# GUARDAR HEATMAP
# ============================================================

# Crear la carpeta si no existe
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# Guardar el heatmap
heatmap_file <- file.path(
  figdir,
  "06_GSVA_heatmap_v2.png"
)

pheatmap::pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  main = "Actividad de rutas Hallmark (muestras ordenadas por RCB)",
  filename = heatmap_file,
  width = 10,
  height = 8
)

# EN GRANDE
pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  main = "Actividad de rutas Hallmark (muestras ordenadas por RCB)",
  filename = file.path(figdir, "06_GSVA_heatmap_grande.png"),
  width = 16,
  height = 10
)

# ============================================================
# 34. PCA DE LAS RUTAS
# ============================================================

pca_pathways <- prcomp(
  t(pathway_scores),
  scale. = FALSE
)


pathway_percent_var <-
  100 *
  pca_pathways$sdev^2 /
  sum(
    pca_pathways$sdev^2
  )


pca_pathway_df <-
  as.data.frame(
    pca_pathways$x[, 1:2]
  ) %>%
  mutate(
    donor_id =
      rownames(
        pca_pathways$x
      ),
    
    response =
      metadata$response[
        match(
          rownames(
            pca_pathways$x
          ),
          metadata$donor_id
        )
      ]
  )


p_pca_pathways <- ggplot(
  pca_pathway_df,
  aes(
    PC1,
    PC2,
    color = response,
    label = donor_id
  )
) +
  geom_point(
    size = 3
  ) +
  geom_text_repel(
    size = 2.8,
    max.overlaps = 30
  ) +
  theme_bw() +
  labs(
    title =
      "PCA de actividad de rutas",
    x = paste0(
      "PC1 (",
      round(
        pathway_percent_var[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        pathway_percent_var[2],
        1
      ),
      "%)"
    )
  )


ggsave(
  file.path(
    figdir,
    "07_PCA_pathways.png"
  ),
  p_pca_pathways,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 35. CORRELACIÓN ENTRE RUTAS
# ============================================================

pathway_cor <- cor(
  t(pathway_scores),
  method = "spearman",
  use = "pairwise.complete.obs"
)


write.csv(
  pathway_cor,
  file.path(
    tabdir,
    "12_pathway_correlation_matrix.csv"
  )
)


# ============================================================
# HEATMAP DE CORRELACIÓN
# ============================================================

cor_file <- file.path(
  figdir,
  "08_pathway_correlation_heatmap.png"
)

pheatmap(
  pathway_cor,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Correlación entre rutas Hallmark",
  filename = cor_file,
  width = 12,
  height = 12
)

# Para visualizarlo en el plot
pheatmap(
  pathway_cor,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Correlación entre rutas Hallmark"
)

# ============================================================
# 36. ANÁLISIS COMPLEMENTARIO DE RUTAS HALLMARK: pCR vs RCB-III
# ============================================================
# Objetivo:
# Evaluar si las diferencias en actividad relativa de las rutas
# Hallmark observadas en pCR vs RD se mantienen al comparar
# pCR frente a enfermedad residual extrema (RCB-III).
#
# Método:
# - Prueba de Wilcoxon para cada ruta.
# - Corrección por comparaciones múltiples: Benjamini-Hochberg.
# - Diferencia de puntuación: media pCR - media RCB-III.
#
# Nota: análisis exploratorio no ajustado por covariables.

# ------------------------------------------------------------
# 36.1. Preparar las puntuaciones GSVA
# ------------------------------------------------------------

scores_mat <- as.matrix(pathway_scores)

if (nrow(scores_mat) != 50 && ncol(scores_mat) == 50) {
  scores_mat <- t(scores_mat)
}

if (is.null(rownames(scores_mat)) ||
    is.null(colnames(scores_mat))) {
  stop("pathway_scores debe tener nombres de rutas y de muestras.")
}

cat("\nDimensiones de las puntuaciones GSVA:\n")
print(dim(scores_mat))

# Convertir a formato: una fila por paciente
scores_df_extreme <- as.data.frame(t(scores_mat)) %>%
  rownames_to_column("donor_id")

# Comprobar que las muestras están disponibles
missing_samples <- setdiff(
  metadata_extreme_model$donor_id,
  scores_df_extreme$donor_id
)

if (length(missing_samples) > 0) {
  stop(
    "Faltan muestras en pathway_scores: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Añadir grupo de respuesta y conservar solo pCR y RCB-III
scores_extreme <- metadata_extreme_model %>%
  dplyr::select(donor_id, response_extreme) %>%
  dplyr::inner_join(scores_df_extreme, by = "donor_id") %>%
  dplyr::mutate(
    response_extreme = factor(
      response_extreme,
      levels = c("RCB-III", "pCR")
    )
  )

cat("\nNúmero de pacientes incluidos:\n")
print(table(scores_extreme$response_extreme))

# ------------------------------------------------------------
# 36.2. Comparar las rutas mediante Wilcoxon
# ------------------------------------------------------------

hallmark_names <- setdiff(
  colnames(scores_extreme),
  c("donor_id", "response_extreme")
)


results_gsva_extreme <- lapply(hallmark_names, function(pathway) {
  
  data_pathway <- scores_extreme %>%
    dplyr::select(response_extreme, dplyr::all_of(pathway)) %>%
    dplyr::filter(!is.na(.data[[pathway]]))
  
  scores_rcb <- data_pathway %>%
    dplyr::filter(response_extreme == "RCB-III") %>%
    dplyr::pull(dplyr::all_of(pathway))
  
  scores_pcr <- data_pathway %>%
    dplyr::filter(response_extreme == "pCR") %>%
    dplyr::pull(dplyr::all_of(pathway))
  
  test <- stats::wilcox.test(
    scores_pcr,
    scores_rcb,
    exact = FALSE
  )
  
  tibble::tibble(
    pathway = pathway,
    n_pCR = length(scores_pcr),
    n_RCB_III = length(scores_rcb),
    mean_pCR = mean(scores_pcr, na.rm = TRUE),
    mean_RCB_III = mean(scores_rcb, na.rm = TRUE),
    median_pCR = median(scores_pcr, na.rm = TRUE),
    median_RCB_III = median(scores_rcb, na.rm = TRUE),
    difference_mean = mean(scores_pcr, na.rm = TRUE) -
      mean(scores_rcb, na.rm = TRUE),
    p_value = test$p.value
  )
  
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    FDR = stats::p.adjust(p_value, method = "BH"),
    direction = dplyr::case_when(
      difference_mean > 0 ~ "Mayor puntuación en pCR",
      difference_mean < 0 ~ "Mayor puntuación en RCB-III",
      TRUE ~ "Sin diferencia en las medias"
    ),
    significant_FDR05 = !is.na(FDR) & FDR < 0.05
  ) %>%
  dplyr::arrange(FDR)


# Mostrar resultados en consola
cat("\nResultados GSVA: pCR frente a RCB-III\n")
print(results_gsva_extreme, n = Inf)

cat(
  "\nRutas con FDR < 0,05:",
  sum(results_gsva_extreme$FDR < 0.05, na.rm = TRUE),
  "\n"
)

# Guardar tabla completa
write.csv(
  results_gsva_extreme,
  file.path(tabdir, "10_GSVA_pCR_vs_RCBIII.csv"),
  row.names = FALSE
)

# Guardar únicamente rutas significativas
write.csv(
  results_gsva_extreme %>%
    filter(!is.na(FDR), FDR < 0.05),
  file.path(tabdir, "10b_GSVA_significativas_pCR_vs_RCBIII.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 36.3. Gráfico de diferencias entre grupos
# ------------------------------------------------------------

plot_data_extreme <- results_gsva_extreme %>%
  mutate(
    pathway_label = gsub("^HALLMARK_", "", pathway),
    pathway_label = gsub("_", " ", pathway_label),
    pathway_label = factor(
      pathway_label,
      levels = pathway_label[order(difference_mean)]
    ),
    significance = if_else(
      !is.na(FDR) & FDR < 0.05,
      "FDR < 0,05",
      "FDR >= 0,05"
    )
  )

p_gsva_extreme <- ggplot(
  plot_data_extreme,
  aes(x = difference_mean, y = pathway_label, color = significance)
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = c(
      "FDR < 0,05" = "#C0392B",
      "FDR >= 0,05" = "grey60"
    )
  ) +
  theme_bw() +
  labs(
    title = "Diferencias en rutas Hallmark: pCR frente a RCB-III",
    subtitle = "Diferencia de medias de puntuación GSVA (pCR - RCB-III)",
    x = "Diferencia media de puntuación GSVA",
    y = "Ruta Hallmark",
    color = "Significación"
  )

ggsave(
  file.path(figdir, "10_GSVA_pCR_vs_RCBIII.png"),
  p_gsva_extreme,
  width = 10,
  height = 12,
  dpi = 300
)

print(p_gsva_extreme)

# ------------------------------------------------------------
# 36.4. Comparar el analisis complementario con pCR frente a RD
# ------------------------------------------------------------
if (
  exists("metadata") &&
  all(c("donor_id", "response") %in% names(metadata))
) {
  
  # Seleccionar las pacientes de la comparación principal pCR vs RD
  metadata_main_gsva <- metadata %>%
    dplyr::filter(response %in% c("pCR", "RD")) %>%
    dplyr::select(donor_id, response)
  
  # Unir metadatos y puntuaciones GSVA
  scores_main <- metadata_main_gsva %>%
    dplyr::inner_join(scores_df_extreme, by = "donor_id")
  
  # Comparar las puntuaciones de cada ruta entre pCR y RD
  main_pathway_results <- lapply(hallmark_names, function(pathway) {
    
    dat <- scores_main %>%
      dplyr::select(response, dplyr::all_of(pathway)) %>%
      dplyr::filter(!is.na(.data[[pathway]]))
    
    x <- dat %>%
      dplyr::filter(response == "pCR") %>%
      dplyr::pull(dplyr::all_of(pathway))
    
    y <- dat %>%
      dplyr::filter(response == "RD") %>%
      dplyr::pull(dplyr::all_of(pathway))
    
    if (length(x) == 0 || length(y) == 0) {
      return(tibble::tibble(
        pathway = pathway,
        mean_difference_pCR_RD = NA_real_,
        p_value_main = NA_real_
      ))
    }
    
    test <- stats::wilcox.test(x, y, exact = FALSE)
    
    tibble::tibble(
      pathway = pathway,
      mean_difference_pCR_RD = mean(x, na.rm = TRUE) -
        mean(y, na.rm = TRUE),
      p_value_main = test$p.value
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      FDR_main = stats::p.adjust(p_value_main, method = "BH")
    )
  
  # Unir los resultados de ambos análisis
  comparison_gsva <- results_gsva_extreme %>%
    dplyr::select(
      pathway,
      difference_mean,
      p_value,
      FDR
    ) %>%
    dplyr::rename(
      difference_mean_pCR_RCBIII = difference_mean,
      p_value_extreme = p_value,
      FDR_extreme = FDR
    ) %>%
    dplyr::left_join(main_pathway_results, by = "pathway") %>%
    dplyr::mutate(
      significant_main = !is.na(FDR_main) & FDR_main < 0.05,
      significant_extreme = !is.na(FDR_extreme) & FDR_extreme < 0.05,
      
      direction_concordant = dplyr::case_when(
        is.na(mean_difference_pCR_RD) |
          is.na(difference_mean_pCR_RCBIII) ~ NA,
        sign(mean_difference_pCR_RD) ==
          sign(difference_mean_pCR_RCBIII) ~ TRUE,
        TRUE ~ FALSE
      ),
      
      interpretation = dplyr::case_when(
        significant_main & significant_extreme &
          direction_concordant ~
          "Significativa en ambos análisis y dirección concordante",
        
        significant_main & significant_extreme &
          !direction_concordant ~
          "Significativa en ambos, dirección discordante",
        
        significant_main & !significant_extreme &
          direction_concordant ~
          "Significativa solo en pCR vs RD; dirección concordante",
        
        !significant_main & significant_extreme ~
          "Significativa solo en pCR vs RCB-III",
        
        TRUE ~ "Sin significación en ambos análisis"
      )
    )
  
  # Guardar la comparación
  write.csv(
    comparison_gsva,
    file.path(tabdir, "10c_comparacion_GSVA_RD_vs_RCBIII.csv"),
    row.names = FALSE
  )
  
  # Mostrar los resultados ordenados por FDR del análisis extremo
  cat("\nComparación entre los dos análisis de rutas:\n")
  
  print(
    comparison_gsva %>%
      dplyr::select(
        pathway,
        FDR_main,
        FDR_extreme,
        direction_concordant,
        interpretation
      ) %>%
      dplyr::arrange(FDR_extreme),
    n = Inf
  )
  
} else {
  message(
    "No se generó la comparación principal: ",
    "metadata no contiene las columnas donor_id y response. ",
    "El análisis pCR vs RCB-III sí se ha realizado."
  )
}

cat("\nAnálisis GSVA pCR vs RCB-III finalizado.\n")

# ============================================================
# 37. MACHINE LEARNING
# ============================================================

# OBJETIVO:
# Evaluar la capacidad de diferentes conjuntos de variables
# para discriminar entre pacientes con respuesta patológica
# completa (pCR) y enfermedad residual (RD).
#
# Se comparan tres tipos de modelos:
#
# 1. Modelo clínico
# 2. Modelo basado en rutas Hallmark
# 3. Modelo combinado: variables clínicas + rutas Hallmark
#
# Se comparan tres algoritmos:
#
# - Elastic Net
# - Random Forest
# - Support Vector Machine (SVM)
#
# CARACTERÍSTICAS DE LA COHORTE:
#
# - 156 pacientes
# - 41 pacientes pCR
# - 115 pacientes RD
#
# Existe un desbalance moderado entre las clases,con predominio de pacientes RD.
#
# Para evitar que el desbalance distorsione la interpretación:
#
# - Se define pCR como clase positiva.
# - Se utiliza ROC-AUC como métrica principal.
# - Se incluyen PR-AUC y balanced accuracy.
# - Se calculan sensibilidad y especificidad.
#
# No se aplica SMOTE ni downsampling en el análisis principal.
#
# VALIDACIÓN:
#
# - Validación cruzada estratificada.
# - 5 particiones.
# - 10 repeticiones.
# - Mismas particiones para los tres tipos de modelos.

# ============================================================
# 37. PREPARAR DATOS PARA MACHINE LEARNING
# ============================================================


# ------------------------------------------------------------
# 37.1. MODELO CLÍNICO
# ------------------------------------------------------------

clinical_ml <- scores_df %>%
  
  select(
    donor_id,
    response,
    age,
    t_stage,
    ln_status_at_diagnosis,
    er_status,
    her2_status,
    pam50
  )


# ------------------------------------------------------------
# 37.2. MODELO DE RUTAS HALLMARK
# ------------------------------------------------------------

# Se utilizan todas las rutas Hallmark.

pathway_ml <- scores_df %>%
  
  select(
    donor_id,
    response,
    all_of(pathway_names)
  )


# ------------------------------------------------------------
# 37.3. MODELO COMBINADO
# ------------------------------------------------------------

combined_ml <- scores_df %>%
  
  select(
    donor_id,
    response,
    age,
    t_stage,
    ln_status_at_diagnosis,
    er_status,
    her2_status,
    pam50,
    all_of(pathway_names)
  )


# ------------------------------------------------------------
# 37.4. COMPROBAR EL NÚMERO DE PACIENTES
# ------------------------------------------------------------

cat(
  "\nPacientes en modelo clínico:",
  nrow(clinical_ml),
  "\n"
)

cat(
  "Pacientes en modelo de rutas:",
  nrow(pathway_ml),
  "\n"
)

cat(
  "Pacientes en modelo combinado:",
  nrow(combined_ml),
  "\n"
)


# Comprobar que las tres matrices contienen los mismos pacientes en el mismo orden.

stopifnot(
  identical(
    clinical_ml$donor_id,
    pathway_ml$donor_id
  )
)

stopifnot(
  identical(
    clinical_ml$donor_id,
    combined_ml$donor_id
  )
)


# Comprobar que la respuesta coincide entre modelos.

stopifnot(
  identical(
    as.character(clinical_ml$response),
    as.character(pathway_ml$response)
  )
)

stopifnot(
  identical(
    as.character(clinical_ml$response),
    as.character(combined_ml$response)
  )
)


# Mostrar distribución de clases.

print(
  table(clinical_ml$response)
)

print(
  prop.table(
    table(clinical_ml$response)
  )
)

# ============================================================
# 38. FUNCIÓN DE PREPROCESAMIENTO
# ============================================================

prepare_ml_data <- function(df) {
  df <- df %>%
    select(
      -donor_id
    )
  
  # ----------------------------------------------------------
  # DEFINIR LA VARIABLE DE RESPUESTA
  # ----------------------------------------------------------
  
  # pCR = clase positiva
  # RD  = clase negativa
  
  df$response <- factor(
    df$response,
    levels = c(
      "pCR",
      "RD"
    )
  )
  
  # ----------------------------------------------------------
  # CONVERTIR VARIABLES CLÍNICAS A FACTORES
  # ----------------------------------------------------------
  
  if ("t_stage" %in% names(df)) {
    
    df$t_stage <- factor(
      df$t_stage
    )
    
  }
  
  
  if ("ln_status_at_diagnosis" %in% names(df)) {
    
    df$ln_status_at_diagnosis <- factor(
      df$ln_status_at_diagnosis
    )
    
  }
  
  
  if ("er_status" %in% names(df)) {
    
    df$er_status <- factor(
      df$er_status
    )
    
  }
  
  
  if ("her2_status" %in% names(df)) {
    
    df$her2_status <- factor(
      df$her2_status
    )
    
  }
  
  
  if ("pam50" %in% names(df)) {
    
    df$pam50 <- factor(
      df$pam50
    )
    
  }
  
  
  # ----------------------------------------------------------
  # ELIMINAR FILAS SIN RESPUESTA
  # ----------------------------------------------------------
  
  df <- df %>%
    
    filter(
      !is.na(response)
    )
  
  return(df)
  
}

# Preparar los tres conjuntos de datos.

ml_clinical <- prepare_ml_data(
  clinical_ml
)

ml_pathways <- prepare_ml_data(
  pathway_ml
)

ml_combined <- prepare_ml_data(
  combined_ml
)


# Comprobar el número de pacientes tras el preprocesamiento.

stopifnot(
  nrow(ml_clinical) == nrow(ml_pathways),
  nrow(ml_clinical) == nrow(ml_combined)
)

# Comprobar los niveles de respuesta.

levels(
  ml_clinical$response
)
# ============================================================
# 39. CONTROL DE VALIDACIÓN CRUZADA
# ============================================================

# Validación cruzada estratificada repetida.
#
# - 5 folds.
# - 10 repeticiones.

stopifnot(
  identical(as.character(ml_clinical$response), as.character(ml_pathways$response)),
  identical(as.character(ml_clinical$response), as.character(ml_combined$response)),
  nrow(ml_clinical) == nrow(ml_pathways),
  nrow(ml_clinical) == nrow(ml_combined)
)

set.seed(12345)

# Crear las particiones a partir de la variable respuesta.

cv_index <- caret::createMultiFolds(
  ml_clinical$response,
  k = 5,
  times = 10
)

# createMultiFolds devuelve los índices de entrenamiento.
# indexOut contiene los pacientes reservados para validación en cada fold.
cv_index_out <- lapply(
  cv_index,
  function(train_rows) setdiff(seq_len(nrow(ml_clinical)), train_rows)
)

# ============================================================
# CONTROL DE VALIDACIÓN CRUZADA
# ============================================================

ctrl <- caret::trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  index = cv_index,
  indexOut = cv_index_out,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  savePredictions = "final",
  verboseIter = FALSE,
  allowParallel = TRUE
)

# ============================================================
# FUNCIÓN PARA ENTRENAR LOS MODELOS
# ============================================================

train_models <- function(data, nombre_modelo) {
  
  cat("\n============================================\n")
  cat("MODELO:", nombre_modelo, "\n")
  cat("============================================\n")

  # ----------------------------------------------------------
  # Preparar variables predictoras
  # ----------------------------------------------------------
  
  # Separar variable respuesta
  x <- data %>%
    dplyr::select(-response)
  
  y <- data$response
  
  # ----------------------------------------------------------
  # Convertir variables categóricas a variables dummy
  # ----------------------------------------------------------
  
  dummy_model <- caret::dummyVars(
    ~ .,
    data = x,
    fullRank = TRUE
  )
  
  x_dummy <- predict(
    dummy_model,
    newdata = x
  ) %>%
    as.data.frame()
  
  # ----------------------------------------------------------
  # Eliminar variables con varianza cero
  # ----------------------------------------------------------
  
  nzv <- caret::nearZeroVar(
    x_dummy
  )
  
  if (length(nzv) > 0) {
    x_dummy <- x_dummy[, -nzv, drop = FALSE]
  }
  
  # Volver a añadir la respuesta
  
  data_model <- x_dummy %>%
    dplyr::mutate(
      response = y
    )
  
  data_model$response <- factor(
    data_model$response,
    levels = c(
      "pCR",
      "RD"
    )
  )

    
  # ----------------------------------------------------------
  # 39.1. ELASTIC NET
  # ----------------------------------------------------------
  
  set.seed(12345)
  
  model_glmnet <- caret::train(
    
    response ~ .,
    
    data = data_model,
    
    method = "glmnet",
    
    metric = "ROC",
    
    trControl = ctrl,
    
    preProcess = c(
      "medianImpute",
      "center",
      "scale"
    ),
    
    tuneLength = 10
    
  )
  
  
  # ----------------------------------------------------------
  # 39.2. RANDOM FOREST
  # ----------------------------------------------------------
  
  model_rf <- caret::train(
    
    response ~ .,
    
    data = data_model,
    
    method = "ranger",
    
    metric = "ROC",
    
    trControl = ctrl,
    
    preProcess = c(
      "medianImpute"
    ),
    
    tuneLength = 5,
    
    importance = "permutation",
    
    num.trees = 1000
    
  )
  
  # ----------------------------------------------------------
  # 39.3. SUPPORT VECTOR MACHINE
  # ----------------------------------------------------------
  
  model_svm <- caret::train(
    
    response ~ .,
    
    data = data_model,
    
    method = "svmRadial",
    
    metric = "ROC",
    
    trControl = ctrl,
    
    preProcess = c(
      "medianImpute",
      "center",
      "scale"
    ),
    
    tuneLength = 5
    
  )
  
  # Devolver los modelos entrenados
  
  return(
    list(
      glmnet = model_glmnet,
      rf = model_rf,
      svm = model_svm
    )
  )
  
}

# ============================================================
# 40. MODELO CLÍNICO
# ============================================================

models_clinical <- train_models(
  ml_clinical,
  "Clínica"
)

# ============================================================
# COMPROBAR MODELOS
# ============================================================

names(models_clinical)

models_clinical$svm

models_clinical$glmnet

models_clinical$rf

# ============================================================
# MÉTRICAS DEL MODELO CLÍNICO
# ============================================================

# -------------------------
# GLMNET
# -------------------------

res_glmnet <- models_clinical$glmnet$results

best_glmnet <- res_glmnet[
  res_glmnet$alpha == models_clinical$glmnet$bestTune$alpha &
    res_glmnet$lambda == models_clinical$glmnet$bestTune$lambda,
]

best_glmnet$BalancedAccuracy <- (
  best_glmnet$Sens + best_glmnet$Spec
) / 2


# -------------------------
# RANDOM FOREST
# -------------------------

res_rf <- models_clinical$rf$results

best_rf <- res_rf[
  res_rf$mtry == models_clinical$rf$bestTune$mtry &
    res_rf$splitrule == models_clinical$rf$bestTune$splitrule &
    res_rf$min.node.size == models_clinical$rf$bestTune$min.node.size,
]

best_rf$BalancedAccuracy <- (
  best_rf$Sens + best_rf$Spec
) / 2


# -------------------------
# SVM
# -------------------------

res_svm <- models_clinical$svm$results

best_svm <- res_svm[
  res_svm$C == models_clinical$svm$bestTune$C,
]

best_svm$BalancedAccuracy <- (
  best_svm$Sens + best_svm$Spec
) / 2

# ============================================================
# MOSTRAR RESULTADOS
# ============================================================

best_glmnet[, c(
  "alpha",
  "lambda",
  "ROC",
  "ROCSD",
  "Sens",
  "SensSD",
  "Spec",
  "SpecSD",
  "BalancedAccuracy"
)]

best_rf[, c(
  "mtry",
  "splitrule",
  "ROC",
  "ROCSD",
  "Sens",
  "SensSD",
  "Spec",
  "SpecSD",
  "BalancedAccuracy"
)]

best_svm[, c(
  "C",
  "ROC",
  "ROCSD",
  "Sens",
  "SensSD",
  "Spec",
  "SpecSD",
  "BalancedAccuracy"
)]

names(models_clinical$glmnet$pred)
head(models_clinical$glmnet$pred)

# ============================================================
# 41. MODELO DE RUTAS HALLMARK
# ============================================================

models_pathways <- train_models(
  ml_pathways,
  "Rutas Hallmark"
)

# ============================================================
# 42. MODELO COMBINADO
# ============================================================

models_combined <- train_models(
  ml_combined,
  "Clínica + Rutas Hallmark"
)
names(models_combined)
models_combined$glmnet$bestTune

# ============================================================
# 43. EVALUACIÓN DEL RENDIMIENTO
# ============================================================

# ------------------------------------------------------------
# 43.1. FUNCIÓN PARA CALCULAR PR-AUC
# ------------------------------------------------------------

calculate_pr_auc <- function(
    observed,
    probability,
    positive_class = "pCR"
) {
  
  observed <- as.character(
    observed
  )
  
  valid <- !is.na(observed) &
    !is.na(probability) &
    is.finite(probability)
  
  observed <- observed[valid]
  
  probability <- probability[valid]
  
  
  # Etiquetas binarias.
  
  positive <- as.integer(
    observed == positive_class
  )
  
  if (
    length(unique(positive)) < 2 ||
    sum(positive) == 0
  ) {
    
    return(NA_real_)
    
  }
  
  
  # Ordenar las probabilidades de mayor a menor.
  
  ord <- order(
    probability,
    decreasing = TRUE
  )
  
  positive <- positive[ord]
  
  probability <- probability[ord]
  
  
  # Verdaderos positivos y falsos positivos acumulados.
  
  tp <- cumsum(
    positive
  )
  
  fp <- cumsum(
    1 - positive
  )
  
  
  recall <- tp / sum(
    positive
  )
  
  precision <- tp / (
    tp + fp
  )
  
  # Conservar el último punto de cada umbral de probabilidad, para agrupar observaciones con probabilidades empatadas.
  
  keep <- !duplicated(
    probability,
    fromLast = TRUE
  )
  
  
  recall <- c(
    0,
    recall[keep]
  )
  
  precision <- c(
    1,
    precision[keep]
  )
  
  
  # Integración trapezoidal.
  
  pr_auc <- sum(
    diff(recall) *
      (
        head(precision, -1) +
          tail(precision, -1)
      ) / 2
  )
  
  
  return(
    as.numeric(pr_auc)
  )
  
}


# ------------------------------------------------------------
# 43.2. EVALUAR LAS PREDICCIONES DE VALIDACIÓN CRUZADA
# ------------------------------------------------------------

evaluate_model_cv <- function(
    model,
    nombre_modelo
) {
  
  predictions <- model$pred
  
  
  if (
    is.null(predictions) ||
    nrow(predictions) == 0
  ) {
    
    stop(
      paste(
        "No hay predicciones guardadas para:",
        nombre_modelo
      )
    )
    
  }
  
  
  # Utilizar la probabilidad de pCR.
  
  if (!"pCR" %in% names(predictions)) {
    
    stop(
      paste(
        "No se encuentra la columna de probabilidades pCR:",
        nombre_modelo
      )
    )
    
  }
  
  
  # Guardar las predicciones del mejor conjunto de hiperparámetros seleccionado por caret.
  
  resamples <- unique(
    predictions$Resample
  )
  
  
  metrics_list <- lapply(
    
    resamples,
    
    function(resample_name) {
      
      d <- predictions[
        predictions$Resample == resample_name,
        ,
        drop = FALSE
      ]
      
      
      # Eliminar filas incompletas.
      
      d <- d[
        !is.na(d$obs) &
          !is.na(d$pCR) &
          !is.na(d$pred),
        ,
        drop = FALSE
      ]
      
      d$obs <- factor(
        d$obs,
        levels = c(
          "pCR",
          "RD"
        )
      )
      
      
      if (
        length(unique(d$obs)) < 2
      ) {
        
        return(
          NULL
        )
        
      }
      
      
      # ------------------------------------------------------
      # ROC-AUC
      # ------------------------------------------------------
      
      roc_object <- pROC::roc(
        
        response = d$obs,
        
        predictor = d$pCR,
        
        levels = c(
          "RD",
          "pCR"
        ),
        
        direction = "<",
        
        quiet = TRUE
        
      )
      
      
      roc_auc <- as.numeric(
        pROC::auc(
          roc_object
        )
      )
      
      
      # ------------------------------------------------------
      # PR-AUC
      # ------------------------------------------------------
      
      pr_auc <- calculate_pr_auc(
        
        observed = d$obs,
        
        probability = d$pCR,
        
        positive_class = "pCR"
        
      )
      
      
      # ------------------------------------------------------
      # SENSIBILIDAD
      # ------------------------------------------------------
      
      true_positive <- sum(
        d$obs == "pCR" &
          d$pred == "pCR"
      )
      
      false_negative <- sum(
        d$obs == "pCR" &
          d$pred != "pCR"
      )
      
      
      sensitivity <- if (
        true_positive + false_negative > 0
      ) {
        
        true_positive /
          (
            true_positive +
              false_negative
          )
        
      } else {
        
        NA_real_
        
      }
      
      
      # ------------------------------------------------------
      # ESPECIFICIDAD
      # ------------------------------------------------------
      
      true_negative <- sum(
        d$obs == "RD" &
          d$pred == "RD"
      )
      
      false_positive <- sum(
        d$obs == "RD" &
          d$pred != "RD"
      )
      
      
      specificity <- if (
        true_negative + false_positive > 0
      ) {
        
        true_negative /
          (
            true_negative +
              false_positive
          )
        
      } else {
        
        NA_real_
        
      }
      
      
      # ------------------------------------------------------
      # BALANCED ACCURACY
      # ------------------------------------------------------
      
      balanced_accuracy <- mean(
        c(
          sensitivity,
          specificity
        ),
        na.rm = TRUE
      )
      
      
      # Devolver métricas de la partición.
      
      data.frame(
        
        modelo = nombre_modelo,
        
        resample = resample_name,
        
        ROC_AUC = roc_auc,
        
        PR_AUC = pr_auc,
        
        sensibilidad = sensitivity,
        
        especificidad = specificity,
        
        balanced_accuracy = balanced_accuracy
        
      )
      
    }
    
  )
  
  
  metrics_list <- Filter(
    Negate(is.null),
    metrics_list
  )
  
  
  if (
    length(metrics_list) == 0
  ) {
    
    stop(
      paste(
        "No se pudieron calcular métricas para:",
        nombre_modelo
      )
    )
    
  }
  
  
  dplyr::bind_rows(
    metrics_list
  )
  
}

# ============================================================
# 44. COMPARACIÓN DE MODELOS
# ============================================================


# ------------------------------------------------------------
# 44.1. MODELOS CLÍNICOS
# ------------------------------------------------------------

cv_clinical_glmnet <- evaluate_model_cv(
  models_clinical$glmnet,
  "Clínica - Elastic Net"
)

cv_clinical_rf <- evaluate_model_cv(
  models_clinical$rf,
  "Clínica - Random Forest"
)

cv_clinical_svm <- evaluate_model_cv(
  models_clinical$svm,
  "Clínica - SVM"
)


# ------------------------------------------------------------
# 44.2. MODELOS DE RUTAS
# ------------------------------------------------------------

cv_pathways_glmnet <- evaluate_model_cv(
  models_pathways$glmnet,
  "Rutas - Elastic Net"
)

cv_pathways_rf <- evaluate_model_cv(
  models_pathways$rf,
  "Rutas - Random Forest"
)

cv_pathways_svm <- evaluate_model_cv(
  models_pathways$svm,
  "Rutas - SVM"
)


# ------------------------------------------------------------
# 44.3. MODELOS COMBINADOS
# ------------------------------------------------------------

cv_combined_glmnet <- evaluate_model_cv(
  models_combined$glmnet,
  "Clínica + Rutas - Elastic Net"
)

cv_combined_rf <- evaluate_model_cv(
  models_combined$rf,
  "Clínica + Rutas - Random Forest"
)

cv_combined_svm <- evaluate_model_cv(
  models_combined$svm,
  "Clínica + Rutas - SVM"
)


# ------------------------------------------------------------
# 44.4. UNIR TODOS LOS RESULTADOS
# ------------------------------------------------------------

cv_metrics <- dplyr::bind_rows(
  
  cv_clinical_glmnet,
  cv_clinical_rf,
  cv_clinical_svm,
  
  cv_pathways_glmnet,
  cv_pathways_rf,
  cv_pathways_svm,
  
  cv_combined_glmnet,
  cv_combined_rf,
  cv_combined_svm
  
)


# Guardar las métricas individuales por partición.

write.csv(
  
  cv_metrics,
  
  file.path(
    tabdir,
    "13_machine_learning_cv_metrics.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 44.5. RESUMIR MEDIA Y DESVIACIÓN ESTÁNDAR
# ------------------------------------------------------------

results_ml <- cv_metrics %>%
  
  group_by(
    modelo
  ) %>%
  
  summarise(
    
    n_particiones = n(),
    
    ROC_AUC_mean = mean(
      ROC_AUC,
      na.rm = TRUE
    ),
    
    ROC_AUC_sd = sd(
      ROC_AUC,
      na.rm = TRUE
    ),
    
    PR_AUC_mean = mean(
      PR_AUC,
      na.rm = TRUE
    ),
    
    PR_AUC_sd = sd(
      PR_AUC,
      na.rm = TRUE
    ),
    
    sensibilidad_mean = mean(
      sensibilidad,
      na.rm = TRUE
    ),
    
    sensibilidad_sd = sd(
      sensibilidad,
      na.rm = TRUE
    ),
    
    especificidad_mean = mean(
      especificidad,
      na.rm = TRUE
    ),
    
    especificidad_sd = sd(
      especificidad,
      na.rm = TRUE
    ),
    
    balanced_accuracy_mean = mean(
      balanced_accuracy,
      na.rm = TRUE
    ),
    
    balanced_accuracy_sd = sd(
      balanced_accuracy,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(
    desc(
      ROC_AUC_mean
    )
  )


# Mostrar resultados.

print(
  results_ml
)

print(
  results_ml,
  width = Inf
)

# Guardar resumen.

write.csv(
  
  results_ml,
  
  file.path(
    tabdir,
    "14_machine_learning_model_comparison.csv"
  ),
  
  row.names = FALSE
  
)

# ============================================================
# 44.6. TABLA COMPARATIVA FINAL DEL RENDIMIENTO
# ============================================================

tabla_rendimiento <- results_ml %>%
  dplyr::select(
    modelo,
    n_particiones,
    ROC_AUC_mean,
    ROC_AUC_sd,
    PR_AUC_mean,
    PR_AUC_sd,
    sensibilidad_mean,
    especificidad_mean,
    balanced_accuracy_mean
  ) %>%
  dplyr::mutate(
    dplyr::across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )

print(tabla_rendimiento, n = Inf)

write.csv(
  tabla_rendimiento,
  file.path(
    tabdir,
    "20_tabla_comparativa_rendimiento_ML.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 44.7. CURVAS ROC COMPARATIVAS
# ============================================================

# Asociar cada resultado con el objeto de modelo correspondiente
lista_modelos <- list(
  "Clínica - Elastic Net" = models_clinical$glmnet,
  "Clínica - Random Forest" = models_clinical$rf,
  "Clínica - SVM" = models_clinical$svm,
  "Rutas - Elastic Net" = models_pathways$glmnet,
  "Rutas - Random Forest" = models_pathways$rf,
  "Rutas - SVM" = models_pathways$svm,
  "Clínica + Rutas - Elastic Net" = models_combined$glmnet,
  "Clínica + Rutas - Random Forest" = models_combined$rf,
  "Clínica + Rutas - SVM" = models_combined$svm
)

# Seleccionar el mejor algoritmo de cada conjunto
mejores_modelos <- results_ml %>%
  dplyr::mutate(
    conjunto = dplyr::case_when(
      grepl("^Clínica -", modelo) ~ "Clínico",
      grepl("^Rutas -", modelo) ~ "Hallmark",
      grepl("^Clínica \\+ Rutas -", modelo) ~ "Combinado"
    )
  ) %>%
  dplyr::group_by(conjunto) %>%
  dplyr::slice_max(
    order_by = ROC_AUC_mean,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

# Obtener puntos de las curvas ROC a partir de las predicciones de validación cruzada de cada modelo
roc_data <- list()
roc_auc_labels <- character()

for (i in seq_len(nrow(mejores_modelos))) {
  
  nombre <- mejores_modelos$modelo[i]
  conjunto <- mejores_modelos$conjunto[i]
  modelo <- lista_modelos[[nombre]]
  pred <- modelo$pred
  
  if (is.null(pred) || !("pCR" %in% names(pred))) {
    warning("No hay predicciones pCR para: ", nombre)
    next
  }
  
  pred <- pred[
    !is.na(pred$obs) &
      !is.na(pred$pCR),
    ,
    drop = FALSE
  ]
  
  roc_obj <- pROC::roc(
    response = pred$obs,
    predictor = pred$pCR,
    levels = c("RD", "pCR"),
    direction = "<",
    quiet = TRUE
  )
  
  puntos <- as.data.frame(
    pROC::coords(
      roc_obj,
      x = "all",
      ret = c("specificity", "sensitivity"),
      transpose = FALSE
    )
  )
  
  puntos$FPR <- 1 - puntos$specificity
  puntos$TPR <- puntos$sensitivity
  puntos$conjunto <- conjunto
  
  roc_data[[conjunto]] <- puntos[
    order(puntos$FPR, puntos$TPR),
  ]
  
  roc_auc_labels[conjunto] <- sprintf(
    "%s (AUC = %.3f)",
    conjunto,
    as.numeric(pROC::auc(roc_obj))
  )
}

roc_plot_data <- dplyr::bind_rows(roc_data)

roc_plot_data$conjunto <- factor(
  roc_plot_data$conjunto,
  levels = names(roc_auc_labels),
  labels = unname(roc_auc_labels)
)

p_roc_comparativa <- ggplot(
  roc_plot_data,
  aes(
    x = FPR,
    y = TPR,
    color = conjunto
  )
) +
  geom_line(linewidth = 1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Comparación de curvas ROC",
    subtitle = "Mejor algoritmo por conjunto de predictores",
    x = "Tasa de falsos positivos",
    y = "Tasa de verdaderos positivos",
    color = "Modelo"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "bottom"
  )

print(p_roc_comparativa)

ggsave(
  filename = file.path(
    figdir,
    "20_curvas_ROC_comparativas.png"
  ),
  plot = p_roc_comparativa,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# 44.8. CURVAS PRECISION-RECALL COMPARATIVAS
# ============================================================

pr_data <- list()

for (i in seq_len(nrow(mejores_modelos))) {
  
  nombre <- mejores_modelos$modelo[i]
  conjunto <- mejores_modelos$conjunto[i]
  modelo <- lista_modelos[[nombre]]
  pred <- modelo$pred
  
  if (is.null(pred) || !("pCR" %in% names(pred))) {
    next
  }
  
  pred <- pred[
    !is.na(pred$obs) &
      !is.na(pred$pCR),
    ,
    drop = FALSE
  ]
  
  # La clase positiva es pCR
  positivo <- as.integer(
    pred$obs == "pCR"
  )
  
  # Ordenar probabilidades de mayor a menor
  orden <- order(
    pred$pCR,
    decreasing = TRUE
  )
  
  positivo <- positivo[orden]
  prob <- pred$pCR[orden]
  
  # Calcular TP y FP acumulados
  tp <- cumsum(positivo)
  fp <- cumsum(1 - positivo)
  
  recall <- tp / sum(positivo)
  precision <- tp / (tp + fp)
  
  # Conservar un punto por umbral de probabilidad
  conservar <- !duplicated(
    prob,
    fromLast = TRUE
  )
  
  curva <- data.frame(
    recall = c(0, recall[conservar]),
    precision = c(1, precision[conservar]),
    conjunto = conjunto
  )
  
  pr_data[[conjunto]] <- curva
}

pr_plot_data <- dplyr::bind_rows(pr_data)

pr_plot_data$conjunto <- factor(
  pr_plot_data$conjunto,
  levels = unique(pr_plot_data$conjunto)
)

# Prevalencia de pCR como referencia
prevalencia_pcr <- mean(
  as.character(ml_clinical$response) == "pCR"
)

p_pr_comparativa <- ggplot(
  pr_plot_data,
  aes(
    x = recall,
    y = precision,
    color = conjunto
  )
) +
  geom_step(linewidth = 1) +
  geom_hline(
    yintercept = prevalencia_pcr,
    linetype = "dashed",
    color = "grey50"
  ) +
  scale_x_continuous(
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    title = "Comparación de curvas Precision–Recall",
    subtitle = "Clase positiva: respuesta patológica completa (pCR)",
    x = "Recall (sensibilidad)",
    y = "Precisión",
    color = "Modelo",
    caption = paste0(
      "Línea discontinua: prevalencia de pCR en la cohorte (",
      round(prevalencia_pcr, 3),
      ")."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "bottom"
  )

print(p_pr_comparativa)

ggsave(
  filename = file.path(
    figdir,
    "21_curvas_PR_comparativas.png"
  ),
  plot = p_pr_comparativa,
  width = 9,
  height = 7,
  dpi = 300
)
# ============================================================
# 45. IMPORTANCIA DE VARIABLES
# ============================================================

# ------------------------------------------------------------
# 45.1. IMPORTANCIA DE RUTAS EN RANDOM FOREST
# ------------------------------------------------------------

rf_imp <- caret::varImp(
  models_pathways$rf,
  scale = FALSE
)

rf_importance <- rf_imp$importance %>%
  
  tibble::rownames_to_column(
    "feature"
  ) %>%
  
  arrange(
    desc(
      Overall
    )
  )

# Guardar tabla.

write.csv(
  
  rf_importance,
  
  file.path(
    tabdir,
    "15_RF_pathway_importance.csv"
  ),
  
  row.names = FALSE
  
)

# ------------------------------------------------------------
# REPRESENTACIÓN GRÁFICA
# ------------------------------------------------------------

png(
  filename = file.path(
    figdir,
    "09_RF_pathway_importance.png"
  ),
  width = 1800,
  height = 1500,
  res = 250
)

plot(
  rf_imp,
  top = min(20, nrow(rf_imp$importance))
)

dev.off()
# ============================================================
# 45.2. IMPORTANCIA DE RUTAS EN ELASTIC NET
# ============================================================

glmnet_importance <- caret::varImp(
  
  models_pathways$glmnet,
  
  scale = FALSE
  
)$importance %>%
  
  tibble::rownames_to_column(
    "feature"
  ) %>%
  
  arrange(
    desc(
      Overall
    )
  )

# Guardar tabla.

write.csv(
  
  glmnet_importance,
  
  file.path(
    tabdir,
    "16_ElasticNet_pathway_importance.csv"
  ),
  
  row.names = FALSE
  
)

# ============================================================
# 45.3. IMPORTANCIA DE VARIABLES DEL MODELO COMBINADO
# ============================================================

rf_combined_imp <- caret::varImp(
  models_combined$rf,
  scale = FALSE
)

importance_combined <- rf_combined_imp$importance %>%
  tibble::rownames_to_column("variable") %>%
  dplyr::arrange(
    dplyr::desc(Overall)
  )

# Guardar tabla completa
write.csv(
  importance_combined,
  file.path(
    tabdir,
    "17_RF_combined_variable_importance.csv"
  ),
  row.names = FALSE
)

# Seleccionar las 20 variables de mayor importancia
top20_importance <- importance_combined %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::mutate(
    variable = factor(
      variable,
      levels = rev(variable)
    )
  )

# Gráfico horizontal
p_importance_combined <- ggplot(
  top20_importance,
  aes(
    x = Overall,
    y = variable
  )
) +
  geom_col() +
  labs(
    title = "Importancia de las variables del modelo combinado",
    subtitle = "Random Forest: variables clínicas y rutas Hallmark",
    x = "Importancia relativa",
    y = "Variable"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

print(p_importance_combined)

ggsave(
  filename = file.path(
    figdir,
    "22_RF_combined_variable_importance.png"
  ),
  plot = p_importance_combined,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# 46. BIOMARCADORES GÉNICOS CANDIDATOS
# ============================================================

# Criterios de selección:
# 1. El gen debe tener p-valor y FDR disponibles en limma.
# 2. El gen debe tener p-valor y FDR disponibles en DESeq2.
# 3. Debe presentar FDR < 0.10 en al menos uno de los métodos.
# 4. Se priorizan como máximo 200 genes.

# Extraer resultados estadísticos de limma
limma_priority <- limma_res %>%
  dplyr::select(
    ensembl_id,
    limma_p = P.Value,
    limma_fdr = adj.P.Val
  )

# Extraer resultados estadísticos de DESeq2
deseq_priority <- res_deseq_df %>%
  dplyr::select(
    ensembl_id,
    deseq_p = pvalue,
    deseq_fdr = padj
  )

# Conservar únicamente genes presentes en ambos métodos
gene_priority <- dplyr::inner_join(
  limma_priority,
  deseq_priority,
  by = "ensembl_id"
) %>%
  # Exigir que los cuatro estadísticos estén disponibles
  dplyr::filter(
    !is.na(limma_p),
    !is.na(limma_fdr),
    !is.na(deseq_p),
    !is.na(deseq_fdr)
  ) %>%
  # FDR < 0.10 en al menos uno de los dos métodos
  dplyr::filter(
    limma_fdr < 0.10 |
      deseq_fdr < 0.10
  ) %>%
  # Mantener únicamente genes presentes en la matriz de expresión
  dplyr::filter(
    ensembl_id %in% rownames(logcpm_expr)
  ) %>%
  # Prioridad basada en ambos p-valores
  dplyr::mutate(
    combined_rank = rowMeans(
      cbind(
        -log10(pmax(limma_p, 1e-300)),
        -log10(pmax(deseq_p, 1e-300))
      )
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(combined_rank)
  ) %>%
  dplyr::slice_head(n = 200)

# Lista definitiva de genes candidatos
top_candidate_genes <- gene_priority$ensembl_id

# Mantener este objeto para el resumen y el guardado final
candidate_gene_pool <- gene_priority$ensembl_id

# Comprobar el resultado
cat(
  "Genes candidatos seleccionados:",
  nrow(gene_priority),
  "\n"
)

cat(
  "Genes con FDR < 0.10 en limma:",
  sum(gene_priority$limma_fdr < 0.10),
  "\n"
)

cat(
  "Genes con FDR < 0.10 en DESeq2:",
  sum(gene_priority$deseq_fdr < 0.10),
  "\n"
)

head(gene_priority)

# ============================================================
# ANOTACIÓN HUMANA: DESCARGA DE ENSEMBL
# ============================================================

gtf_url <- paste0(
  "https://ftp.ensembl.org/pub/release-116/gtf/",
  "homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz"
)

gtf_file <- file.path(
  outdir,
  "Homo_sapiens.GRCh38.116.gtf.gz"
)

# Comprobar que existe
file.exists(gtf_file)

download.file(
  url = gtf_url,
  destfile = gtf_file,
  mode = "wb",
  method = "libcurl"
)

file.exists(gtf_file)

# ============================================================
# ANOTACIÓN ENSEMBL -> SÍMBOLO DE GEN DESDE GTF LOCAL
# ============================================================

# Leer el archivo GTF comprimido
gtf <- read.delim(
  gzfile(gtf_file),
  header = FALSE,
  comment.char = "#",
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE
)

# Seleccionar únicamente los registros correspondientes a genes
gtf_genes <- gtf[gtf$V3 == "gene", ]

# Extraer identificadores Ensembl
ensembl_gene <- sub(
  '.*gene_id "([^"]+)".*',
  '\\1',
  gtf_genes$V9
)

# Eliminar la versión del identificador Ensembl
ensembl_gene <- sub(
  "\\.[0-9]+$",
  "",
  ensembl_gene
)

# Extraer símbolos de genes
gene_symbol <- ifelse(
  grepl('gene_name "', gtf_genes$V9),
  sub(
    '.*gene_name "([^"]+)".*',
    '\\1',
    gtf_genes$V9
  ),
  NA_character_
)

# Crear la tabla de anotación
gene_annotation_complete <- data.frame(
  ensembl_gene = ensembl_gene,
  gene_symbol = gene_symbol,
  stringsAsFactors = FALSE
)

# Eliminar duplicados
gene_annotation_complete <- gene_annotation_complete %>%
  dplyr::distinct(
    ensembl_gene,
    .keep_all = TRUE
  )

# Guardar la anotación completa
write.csv(
  gene_annotation_complete,
  file.path(
    outdir,
    "gene_annotation_complete.csv"
  ),
  row.names = FALSE
)

# Comprobar los resultados
print(head(gene_annotation_complete))

cat(
  "Identificadores Ensembl anotados:",
  nrow(gene_annotation_complete),
  "\n"
)

cat(
  "Genes con símbolo disponible:",
  sum(!is.na(gene_annotation_complete$gene_symbol)),
  "\n"
)

# ============================================================
# ANOTAR LOS GENES CANDIDATOS PRIORIZADOS
# ============================================================

# Normalizar los identificadores para asegurar la coincidencia
gene_priority <- gene_priority %>%
  dplyr::mutate(
    ensembl_id = sub("\\.[0-9]+$", "", ensembl_id)
  )

gene_annotation_complete <- gene_annotation_complete %>%
  dplyr::mutate(
    ensembl_gene = sub("\\.[0-9]+$", "", ensembl_gene)
  )

# Añadir el símbolo a los genes candidatos
gene_priority_annotated <- gene_priority %>%
  dplyr::left_join(
    gene_annotation_complete,
    by = c("ensembl_id" = "ensembl_gene")
  )

# Crear la tabla final de candidatos
candidate_genes_table <- gene_priority_annotated %>%
  dplyr::filter(
    ensembl_id %in% top_candidate_genes
  ) %>%
  dplyr::select(
    ensembl_id,
    gene_symbol,
    dplyr::everything()
  )

# Guardar la tabla final
write.csv(
  candidate_genes_table,
  file.path(
    outdir,
    "candidate_genes_annotated.csv"
  ),
  row.names = FALSE
)

# Comprobar la anotación de los candidatos
cat(
  "Genes candidatos:",
  nrow(candidate_genes_table),
  "\n"
)

cat(
  "Candidatos con símbolo:",
  sum(!is.na(candidate_genes_table$gene_symbol)),
  "\n"
)

cat(
  "Candidatos sin símbolo:",
  sum(is.na(candidate_genes_table$gene_symbol)),
  "\n"
)

print(head(candidate_genes_table))

# ============================================================
# 47. AUC INDIVIDUAL DE GENES
# ============================================================

gene_auc <- map_dfr(
  
  top_candidate_genes,
  
  function(g) {
    
    tmp <- tibble(
      
      response =
        metadata$response[
          match(
            colnames(logcpm_expr),
            metadata$donor_id
          )
        ],
      
      expression =
        as.numeric(
          logcpm_expr[
            g,
          ]
        )
    )
    
    
    roc_obj <- roc(
      
      response =
        tmp$response,
      
      predictor =
        tmp$expression,
      
      levels =
        c(
          "RD",
          "pCR"
        ),
      
      quiet =
        TRUE
    )
    
    
    tibble(
      
      ensembl_id =
        g,
      
      auc =
        as.numeric(
          auc(
            roc_obj
          )
        ),
      
      auc_distance =
        abs(
          as.numeric(
            auc(
              roc_obj
            )
          ) -
            0.5
        )
    )
    
  }
) %>%
  
  left_join(
    gene_priority,
    by =
      "ensembl_id"
  ) %>%
  
  arrange(
    desc(
      auc_distance
    )
  )


write.csv(
  gene_auc,
  file.path(
    tabdir,
    "17_candidate_gene_biomarkers.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 47.1. GRÁFICO AUC DE LOS 20 GENES CANDIDATOS PRIORITARIOS
# ============================================================

# ------------------------------------------------------------
# Seleccionar los 20 genes con mayor distancia respecto a 0.5
# ------------------------------------------------------------

top20_auc <- gene_auc %>%
  dplyr::filter(!is.na(auc)) %>%
  dplyr::arrange(dplyr::desc(auc_distance)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::left_join(
    candidate_genes_table %>%
      dplyr::select(ensembl_id, gene_symbol),
    by = "ensembl_id"
  ) %>%
  dplyr::mutate(
    gene_label = dplyr::if_else(
      !is.na(gene_symbol) & gene_symbol != "",
      gene_symbol,
      ensembl_id
    )
  )

# Comprobar resultados
nrow(top20_auc)
head(top20_auc)

# ------------------------------------------------------------
# Ordenar los genes para la representación
# ------------------------------------------------------------

top20_auc <- top20_auc %>%
  dplyr::arrange(auc) %>%
  dplyr::mutate(
    gene_label = factor(
      gene_label,
      levels = gene_label
    )
  )

# ------------------------------------------------------------
# Generar el gráfico
# ------------------------------------------------------------

p_top20_auc <- ggplot(
  top20_auc,
  aes(
    x = auc,
    y = gene_label
  )
) +
  geom_point(
    size = 3.5,
    color = "black"
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  labs(
    title = "AUC de los 20 genes candidatos prioritarios",
    x = "Área bajo la curva ROC (AUC)",
    y = "Gen"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(
      size = 13
    ),
    axis.text.y = element_text(
      size = 10
    ),
    panel.grid.minor = element_line(
      color = "grey90"
    )
  )

# Mostrar figura
print(p_top20_auc)

# ------------------------------------------------------------
# Guardar figura en alta resolución
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figdir,
    "10_top20_genes_auc.png"
  ),
  plot = p_top20_auc,
  width = 10,
  height = 7.5,
  dpi = 300
)

# Guardar en PDF 
ggsave(
  filename = file.path(
    figdir,
    "10_top20_genes_auc.pdf"
  ),
  plot = p_top20_auc,
  width = 10,
  height = 7.5
)

# ============================================================
# 48. ANÁLISIS DE REDES Y CORRELACIÓN ENTRE GENES
# ============================================================

# ------------------------------------------------------------
# 48.1. SELECCIÓN DE LOS 200 GENES CANDIDATOS
# ------------------------------------------------------------

top_network_genes <- unique(
  top_candidate_genes
)

top_network_genes <- intersect(
  top_network_genes,
  rownames(logcpm_expr)
)

cat(
  "Genes candidatos para el análisis de red:",
  length(top_network_genes),
  "\n"
)


# ------------------------------------------------------------
# 48.2. MATRIZ DE EXPRESIÓN DE LOS 200 GENES
# ------------------------------------------------------------

gene_expression_matrix <- t(
  logcpm_expr[
    top_network_genes,
    metadata$donor_id,
    drop = FALSE
  ]
)


cat(
  "Dimensiones de la matriz:",
  nrow(gene_expression_matrix),
  "muestras x",
  ncol(gene_expression_matrix),
  "genes\n"
)


# ============================================================
# 48.3. CORRELACIÓN ENTRE LOS 200 GENES
# ============================================================

gene_cor_200 <- cor(
  gene_expression_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)


# ------------------------------------------------------------
# Guardar matriz de correlación
# ------------------------------------------------------------

saveRDS(
  gene_cor_200,
  file = file.path(
    outdir,
    "gene_correlation_200.rds"
  )
)

write.csv(
  gene_cor_200,
  file = file.path(
    outdir,
    "gene_correlation_200.csv"
  ),
  row.names = TRUE
)
# ------------------------------------------------------------
# HEATMAP DE CORRELACIÓN - 200 GENES
# ------------------------------------------------------------

p_cor_200 <- pheatmap::pheatmap(
  gene_cor_200,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Correlación entre 200 genes candidatos",
  silent = TRUE
)

png(
  filename = file.path(
    figdir,
    "09_gene_correlation_200.png"
  ),
  width = 2600,
  height = 2400,
  res = 250
)

grid::grid.newpage()
grid::grid.draw(p_cor_200$gtable)

dev.off()

# ============================================================
# 48.4. SELECCIÓN DE LOS 20 GENES PRINCIPALES
# ============================================================

top20_genes <- gene_priority %>%
  slice_head(n = 20) %>%
  pull(ensembl_id)

top20_genes <- intersect(
  top20_genes,
  rownames(logcpm_expr)
)

cat(
  "Genes prioritarios para el análisis de 20 genes:",
  length(top20_genes),
  "\n"
)


# ------------------------------------------------------------
# MATRIZ DE EXPRESIÓN DE LOS 20 GENES
# ------------------------------------------------------------

gene_expression_matrix_20 <- t(
  logcpm_expr[
    top20_genes,
    metadata$donor_id,
    drop = FALSE
  ]
)

# ============================================================
# 48.5. CORRELACIÓN ENTRE LOS 20 GENES
# ============================================================

gene_cor_20 <- cor(
  gene_expression_matrix_20,
  method = "spearman",
  use = "pairwise.complete.obs"
)


# ------------------------------------------------------------
# Guardar matriz
# ------------------------------------------------------------

saveRDS(
  gene_cor_20,
  file = file.path(
    outdir,
    "gene_correlation_20.rds"
  )
)

write.csv(
  gene_cor_20,
  file = file.path(
    outdir,
    "gene_correlation_20.csv"
  ),
  row.names = TRUE
)

# ============================================================
# CONVERTIR LOS ID ENSEMBL A SÍMBOLOS GÉNICOS
# ============================================================

# Normalizar los identificadores Ensembl
gene_annotation_plot <- gene_annotation %>%
  dplyr::mutate(
    ensembl_gene = clean_ensembl(ensembl_gene)
  ) %>%
  dplyr::filter(
    !is.na(ensembl_gene),
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  dplyr::distinct(
    ensembl_gene,
    .keep_all = TRUE
  )

# Crear correspondencia entre Ensembl y símbolo
gene_symbols <- gene_annotation_plot$gene_symbol[
  match(
    clean_ensembl(colnames(gene_cor_20)),
    gene_annotation_plot$ensembl_gene
  )
]

# Mantener el identificador Ensembl si no hay símbolo disponible
missing_symbols <- is.na(gene_symbols) | gene_symbols == ""

gene_symbols[missing_symbols] <-
  colnames(gene_cor_20)[missing_symbols]

# Evitar nombres duplicados en la matriz
gene_symbols <- make.unique(gene_symbols)

# Sustituir los nombres de filas y columnas
rownames(gene_cor_20) <- gene_symbols
colnames(gene_cor_20) <- gene_symbols

# Comprobar la conversión
print(
  data.frame(
    Ensembl = top20_genes,
    Simbolo = gene_symbols
  )
)

# ------------------------------------------------------------
# HEATMAP DE CORRELACIÓN - 20 GENES
# ------------------------------------------------------------

p_cor_20 <- pheatmap::pheatmap(
  gene_cor_20,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Correlación entre 20 genes candidatos prioritarios",
  silent = TRUE
)

png(
  filename = file.path(
    figdir,
    "09_gene_correlation_20.png"
  ),
  width = 2200,
  height = 2000,
  res = 250
)

grid::grid.newpage()
grid::grid.draw(p_cor_20$gtable)

dev.off()

# ============================================================
# 49. RED TRANSCRIPTÓMICA - 200 GENES
# ============================================================

# ------------------------------------------------------------
# Crear edges a partir de la matriz de correlación
# ------------------------------------------------------------

edges_200 <- as.data.frame(
  as.table(gene_cor_200),
  stringsAsFactors = FALSE
)

colnames(edges_200) <- c(
  "from",
  "to",
  "rho"
)

edges_200 <- edges_200 %>%
  dplyr::filter(
    from != to,
    abs(rho) >= 0.50
  ) %>%
  dplyr::mutate(
    pair = purrr::map2_chr(
      from,
      to,
      ~ paste(
        sort(c(.x, .y)),
        collapse = "_"
      )
    )
  ) %>%
  dplyr::distinct(
    pair,
    .keep_all = TRUE
  ) %>%
  dplyr::select(
    -pair
  )


# ------------------------------------------------------------
# Guardar edges
# ------------------------------------------------------------

write.csv(
  edges_200,
  file.path(
    tabdir,
    "18_network_edges_200.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 49.1. RED TRANSCRIPTÓMICA - 20 GENES
# ============================================================

# ------------------------------------------------------------
# Crear edges
# ------------------------------------------------------------

edges_20 <- as.data.frame(
  as.table(gene_cor_20),
  stringsAsFactors = FALSE
)

colnames(edges_20) <- c(
  "from",
  "to",
  "rho"
)

edges_20 <- edges_20 %>%
  dplyr::filter(
    from != to,
    abs(rho) >= 0.50
  ) %>%
  dplyr::mutate(
    pair = purrr::map2_chr(
      from,
      to,
      ~ paste(
        sort(c(.x, .y)),
        collapse = "_"
      )
    )
  ) %>%
  dplyr::distinct(
    pair,
    .keep_all = TRUE
  ) %>%
  dplyr::select(-pair)

# Comprobar el resultado
nrow(edges_20)
head(edges_20)

# ------------------------------------------------------------
# Guardar edges
# ------------------------------------------------------------

write.csv(
  edges_20,
  file.path(
    tabdir,
    "18_network_edges_20.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Generar red de 20 genes
# ------------------------------------------------------------

if (nrow(edges_20) > 0) {
  
  g_20 <- igraph::graph_from_data_frame(
    edges_20,
    directed = FALSE
  )
  
  set.seed(123)
  
  layout_20 <- igraph::layout_with_fr(
    g_20
  )
  
  png(
    filename = file.path(
      figdir,
      "10_transcriptomic_network_20.png"
    ),
    width = 2200,
    height = 1800,
    res = 250
  )
  
  plot(
    g_20,
    layout = layout_20,
    
    vertex.size = 12,
    
    vertex.color = "gold",
    
    vertex.frame.color = "black",
    
    vertex.label = V(g_20)$name,
    
    vertex.label.cex = 0.65,
    
    vertex.label.color = "black",
    
    vertex.label.dist = 1.3,
    
    edge.color = "grey60",
    
    edge.width = 1 + abs(E(g_20)$rho) * 3,
    
    main = "Red transcriptomica exploratoria de genes candidatos prioritarios"
  )
  
  dev.off()
}
# ============================================================
# 50. RESUMEN FINAL
# ============================================================

n_total <- nrow(metadata)

n_pcr <- sum(
  metadata$response == "pCR",
  na.rm = TRUE
)

n_rd <- sum(
  metadata$response == "RD",
  na.rm = TRUE
)


# ------------------------------------------------------------
# Mejor modelo según ROC-AUC media
# ------------------------------------------------------------

results_ml_ordered <-
  results_ml %>%
  arrange(
    desc(
      ROC_AUC_mean
    )
  )


best_model <-
  results_ml_ordered$modelo[1]

best_roc_auc <-
  results_ml_ordered$ROC_AUC_mean[1]

best_pr_auc <-
  results_ml_ordered$PR_AUC_mean[1]


# ------------------------------------------------------------
# Resumen de resultados
# ------------------------------------------------------------

summary_results <- data.frame(
  
  metric = c(
    
    "Muestras RNA-seq iniciales",
    "Muestras con respuesta",
    "pCR",
    "RD",
    
    "Genes iniciales",
    "Genes tras expresión",
    "Genes variables para PCA",
    
    "DEGs limma FDR < 0.05",
    "DEGs DESeq2 FDR < 0.05",
    "DEGs comunes",
    
    "Hallmark pathways",
    "Pathways FDR < 0.05",
    
    "Mejor modelo ML",
    "Mejor ROC-AUC media",
    "Mejor PR-AUC media"
    
  ),
  
  value = c(
    
    ncol(raw_counts),
    n_total,
    n_pcr,
    n_rd,
    
    nrow(counts),
    nrow(counts_expr),
    nrow(counts_variable),
    
    sum(
      limma_res$adj.P.Val < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      res_deseq_df$padj < 0.05,
      na.rm = TRUE
    ),
    
    length(overlap_deg),
    
    length(pathway_names),
    
    sum(
      pathway_tests$FDR < 0.05,
      na.rm = TRUE
    ),
    
    best_model,
    best_roc_auc,
    best_pr_auc
    
  )
)


write.csv(
  summary_results,
  file.path(
    tabdir,
    "19_TFM_analysis_summary.csv"
  ),
  row.names = FALSE
)


print(
  summary_results
)

# ============================================================
# Comparación de genes diferencialmente expresados entre limma y DESeq2
# ============================================================


overlap_table <- limma_res %>%
  dplyr::filter(
    ensembl_id %in% overlap_deg
  ) %>%
  dplyr::select(
    ensembl_id,
    logFC,
    P.Value,
    adj.P.Val
  ) %>%
  dplyr::left_join(
    res_deseq_df %>%
      dplyr::filter(
        ensembl_id %in% overlap_deg
      ) %>%
      dplyr::select(
        ensembl_id,
        log2FoldChange,
        pvalue,
        padj
      ),
    by = "ensembl_id"
  ) %>%
  dplyr::arrange(adj.P.Val)

# Comprobar los resultados
nrow(overlap_table)
head(overlap_table)

overlap_table

write.csv(
  overlap_table,
  file.path(outdir, "DEGs_comunes_limma_DESeq2.csv"),
  row.names = FALSE
)

# ============================================================
# 51. GUARDAR OBJETOS COMPLETOS
# ============================================================
saveRDS(
  
  list(
    
    # ========================================================
    # DATOS Y MATRICES
    # ========================================================
    
    metadata =
      metadata,
    
    # Matriz RNA-seq original
    counts =
      counts,
    
    # Matriz tras el filtro de expresión
    # (~17.334 genes)
    counts_expression_filtered =
      counts_expr,
    
    # Matriz de los genes más variables
    # (~8.000 genes), utilizada para PCA
    counts_variable =
      counts_variable,
    
    # DGEList de la matriz filtrada por expresión
    dge_expression_filtered =
      dge_expr,
    
    # logCPM normalizado de todos los genes
    logcpm_expression_filtered =
      logcpm_expr,
    
    # logCPM de los genes más variables
    logcpm_variable =
      logcpm_variable,
    
    
    # ========================================================
    # ANÁLISIS EXPLORATORIO
    # ========================================================
    
    pca =
      pca,
    
    
    # ========================================================
    # EXPRESIÓN DIFERENCIAL: pCR vs RD
    # ========================================================
    
    limma =
      limma_res,
    
    DESeq2 =
      res_deseq_df,
    
    consistent_DEGs =
      consistent_genes,
    
    
    # ========================================================
    # NUEVO: EXPRESIÓN DIFERENCIAL pCR vs RCB-III
    # ========================================================
    
    metadata_extreme =
      if (exists("metadata_extreme"))
        metadata_extreme else NULL,
    
    metadata_extreme_model =
      if (exists("metadata_extreme_model"))
        metadata_extreme_model else NULL,
    
    counts_extreme =
      if (exists("counts_extreme"))
        counts_extreme else NULL,
    
    covariates_extreme =
      if (exists("covariates_extreme"))
        covariates_extreme else NULL,
    
    design_extreme =
      if (exists("design_extreme"))
        design_extreme else NULL,
    
    limma_extreme =
      if (exists("limma_extreme"))
        limma_extreme else NULL,
    
    DESeq2_extreme =
      if (exists("res_extreme_df"))
        res_extreme_df else NULL,
    
    
    # ========================================================
    # ANÁLISIS DE RUTAS HALLMARK: pCR vs RD
    # ========================================================
    
    pathway_scores =
      scores_df,
    
    pathway_tests =
      pathway_tests,
    
    pathway_correlation =
      pathway_cor,
    
    
    # ========================================================
    # NUEVO: ANÁLISIS GSVA pCR vs RCB-III
    # ========================================================
    
    pathway_scores_matrix =
      if (exists("pathway_scores"))
        pathway_scores else NULL,
    
    pathway_scores_extreme =
      if (exists("scores_extreme"))
        scores_extreme else NULL,
    
    pathway_tests_extreme =
      if (exists("results_gsva_extreme"))
        results_gsva_extreme else NULL,
    
    pathway_comparison_RD_RCBIII =
      if (exists("comparison_gsva"))
        comparison_gsva else NULL,
    
    
    # ========================================================
    # NUEVO: HEATMAP ORDENADO POR RCB
    # ========================================================
    
    heatmap_matrix =
      if (exists("heatmap_matrix"))
        heatmap_matrix else NULL,
    
    heatmap_scaled =
      if (exists("heatmap_scaled"))
        heatmap_scaled else NULL,
    
    heatmap_order =
      if (exists("heatmap_order"))
        heatmap_order else NULL,
    
    heatmap_annotation =
      if (exists("annotation_col"))
        annotation_col else NULL,
    
    
    # ========================================================
    # MACHINE LEARNING
    # ========================================================
    
    machine_learning =
      list(
        
        clinical =
          models_clinical,
        
        pathways =
          models_pathways,
        
        combined =
          models_combined,
        
        comparison =
          results_ml,
        
        summary =
          summary_results
      ),
    
    
    # ========================================================
    # BIOMARCADORES GÉNICOS
    # ========================================================
    
    candidate_gene_pool =
      candidate_gene_pool,
    
    gene_priority =
      gene_priority,
    
    top_candidate_genes =
      top_candidate_genes,
    
    gene_biomarkers =
      gene_auc,
    
    
    # ========================================================
    # REDES DE COEXPRESIÓN
    # ========================================================
    
    network_edges_200 =
      edges_200,
    
    network_edges_20 =
      if (exists("edges_20"))
        edges_20 else NULL
    
  ),
  
  file.path(
    outdir,
    "TFM_full_analysis_COMPLETE.rds"
  )
  
)


# ============================================================
# 52. GUARDAR ENTORNO COMPLETO DE R
# ============================================================

save.image(
  file.path(
    outdir,
    "TFM_entorno_completo.RData"
  )
)


# ============================================================
# 53. COMPROBAR QUE LOS ARCHIVOS SE HAN GUARDADO
# ============================================================

file.exists(
  file.path(outdir, "TFM_full_analysis_COMPLETE.rds")
)

file.exists(
  file.path(outdir, "TFM_entorno_completo.RData")
)


# Versiones del entorno computacional
R.version.string

paquetes <- c(
  "edgeR", "limma", "DESeq2", "GSVA", "msigdbr",
  "caret", "pROC", "glmnet", "ranger", "e1071",
  "pheatmap", "corrplot", "igraph"
)

versiones <- data.frame(
  paquete = paquetes,
  version = vapply(
    paquetes,
    function(x) {
      if (requireNamespace(x, quietly = TRUE)) {
        as.character(packageVersion(x))
      } else {
        NA_character_
      }
    },
    character(1)
  )
)

write.csv(
  versiones,
  file.path(tabdir, "versiones_paquetes.csv"),
  row.names = FALSE
)

print(versiones)
