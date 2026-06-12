# ============================================================================
# ANÁLISE ESTATÍSTICA AUTOMATIZADA EM SAÚDE
# Descritiva + Tabela Cruzada + Regressão Logística
# Autor: Cruz S Sebastião (https://orcid.org/0000-0003-1232-0119)
# Versão: 1.0
# ============================================================================


# ----------------------------------------------------------------------------
# 1. INSTALAÇÃO E CARREGAMENTO DE PACOTES
# ----------------------------------------------------------------------------
pacotes <- c("tidyverse", "openxlsx", "gt", "haven")
novos   <- pacotes[!pacotes %in% installed.packages()[, "Package"]]
if (length(novos)) install.packages(novos)
invisible(lapply(pacotes, library, character.only = TRUE))


# ----------------------------------------------------------------------------
# 2. IMPORTAÇÃO DOS DADOS (DETEÇÃO AUTOMÁTICA DA EXTENSÃO)
# ----------------------------------------------------------------------------
getwd()
if (interactive()) {
  
  diretorio <- readline("Indique o caminho completo do diretório onde está o ficheiro: ")
  if (!dir.exists(diretorio)) stop("Diretório não existe.")
  setwd(diretorio)
  
  message("Ficheiros disponíveis:")
  print(list.files())
  
  ficheiro <- readline("Indique o nome completo do ficheiro (ex: dados.xlsx ou dados.sav): ")
  if (!file.exists(ficheiro)) stop("Ficheiro não encontrado.")
  
  extensao <- tools::file_ext(ficheiro)
  
  if (tolower(extensao) == "xlsx") {
    
    dados <- read.xlsx(ficheiro)
    
  } else if (tolower(extensao) == "sav") {
    
    dados <- haven::read_sav(ficheiro)
    
  } else {
    stop("Extensão não suportada. Utilize .xlsx ou .sav")
  }
  
} else stop("Execução interativa necessária.")

message("Ficheiro importado com sucesso.")
message("Observações: ", nrow(dados),
        " | Variáveis: ", ncol(dados))

# Leitura automática da estrutura
str(dados)


# ----------------------------------------------------------------------------
# 3. PREPARAÇÃO DOS DADOS
# ----------------------------------------------------------------------------

# 3.1 Converter character em factor (exceto identificadores iniciados por id)
colunas_converter <- names(dados)[
  sapply(dados, is.character) &
    !grepl("^id", names(dados), ignore.case = TRUE)
]

dados[colunas_converter] <- lapply(dados[colunas_converter], factor)

message("Variáveis character convertidas para factor (exceto identificadores).")


# 3.2 Criação interativa de novas categorias
repeat {
  
  criar <- menu(c("Sim", "Não"),
                title = "Pretende criar categorias em alguma variável?")
  
  if (criar == 2) break
  
  tipo_var <- menu(c("Variável Quantitativa",
                     "Variável Qualitativa"),
                   title = "Tipo de variável base?")
  
  ordered_choice <- menu(c("Sim", "Não"),
                         title = "A nova variável deve ser ordinal (ordered = TRUE)?")
  
  is_ordered <- ordered_choice == 1
  
  
# 3.3 Identificação e Tratamento de Valores Omissos
  
  # Percentagem de omissos por variável (colunas)
  missing_col <- round(colMeans(is.na(dados)) * 100, 2)
  
  tabela_missing_col <- data.frame(
    Variavel = names(missing_col),
    Percentagem_Omissos = missing_col,
    row.names = NULL
  )
  
  print(tabela_missing_col)
  
  # Percentagem de omissos por observação (linhas)
  missing_row <- round(rowMeans(is.na(dados)) * 100, 2)
  
  print(summary(missing_row))
  
  # Remover variáveis com mais de 5% de omissos?
  remover_col <- menu(c("Sim", "Não"),
                      title = "Deseja remover variáveis com mais de 5% de omissos?")
  
  if (remover_col == 1) {
    
    vars_remover <- names(missing_col[missing_col > 5])
    
    if (length(vars_remover) > 0) {
      dados <- dados[, !names(dados) %in% vars_remover]
      message("Variáveis removidas: ", paste(vars_remover, collapse = ", "))
    } else {
      message("Nenhuma variável ultrapassa 5% de omissos.")
    }
  }
  
  # Remover observações com mais de 5% de omissos?
  remover_row <- menu(c("Sim", "Não"),
                      title = "Deseja remover observações com mais de 5% de omissos?")
  
  if (remover_row == 1) {
    
    linhas_remover <- which(missing_row > 5)
    
    if (length(linhas_remover) > 0) {
      dados <- dados[-linhas_remover, ]
      message(length(linhas_remover), " observações removidas.")
    } else {
      message("Nenhuma observação ultrapassa 5% de omissos.")
    }
  }
  
  print(dim(dados))
  
  # -------------------------------------------------------------------------
  # A PARTIR DE VARIÁVEL QUANTITATIVA
  # -------------------------------------------------------------------------
  if (tipo_var == 1) {
    
    vars_num <- names(dados)[sapply(dados, is.numeric)]
    if (length(vars_num) == 0) {
      message("Não existem variáveis quantitativas disponíveis.")
      next
    }
    
    var_base <- select.list(vars_num,
                            title = "Selecione a variável quantitativa:")
    
    cortes_input <- readline("Indique os pontos de corte separados por vírgula (ex: 20,40): ")
    cortes <- as.numeric(strsplit(cortes_input, ",")[[1]])
    cortes <- c(-Inf, cortes, Inf)
    
    labels_input <- readline("Indique os rótulos separados por vírgula (ex: <20,20-40,>40): ")
    labels <- strsplit(labels_input, ",")[[1]]
    
    if (length(labels) != (length(cortes) - 1))
      stop("Número de rótulos incompatível com número de intervalos.")
    
    nome_nova <- readline("Indique o nome da nova variável: ")
    
    if (nome_nova %in% names(dados)) {
      message("Já existe uma variável com esse nome. Operação cancelada.")
      next
    }
    
    nova_var <- cut(
      dados[[var_base]],
      breaks = cortes,
      labels = labels,
      right = FALSE
    )
    
    if (is_ordered)
      nova_var <- factor(nova_var, ordered = TRUE)
    
    dados[[nome_nova]] <- nova_var
    
    message("Variável criada com sucesso: ", nome_nova)
    print(table(dados[[nome_nova]]))
    
  }
  
  # -------------------------------------------------------------------------
  # A PARTIR DE VARIÁVEL QUALITATIVA
  # -------------------------------------------------------------------------
  if (tipo_var == 2) {
    
    vars_cat <- names(dados)[sapply(dados, is.factor)]
    if (length(vars_cat) == 0) {
      message("Não existem variáveis qualitativas disponíveis.")
      next
    }
    
    var_base <- select.list(vars_cat,
                            title = "Selecione a variável qualitativa:")
    
    niveis_orig <- levels(dados[[var_base]])
    print(niveis_orig)
    
    nome_nova <- readline("Indique o nome da nova variável: ")
    
    if (nome_nova %in% names(dados)) {
      message("Já existe uma variável com esse nome. Operação cancelada.")
      next
    }
    
    nova_var <- dados[[var_base]]
    
    for (nivel in niveis_orig) {
      novo_nome <- readline(paste("Novo grupo para categoria", nivel, ": "))
      levels(nova_var)[levels(nova_var) == nivel] <- novo_nome
    }
    
    if (is_ordered)
      nova_var <- factor(nova_var, ordered = TRUE)
    
    dados[[nome_nova]] <- nova_var
    
    message("Variável recategorizada criada com sucesso: ", nome_nova)
    print(table(dados[[nome_nova]]))
  }
  
}

message("Finalizada a etapa de criação de categorias.")


################### 44444444444444444 NEWWWWEEEEEEEEEE
# ----------------------------------------------------------------------------
# 4. DEFINIÇÃO DA VARIÁVEL DEPENDENTE (BINÁRIA)
# ----------------------------------------------------------------------------

variaveis_factor <- names(dados)[sapply(dados, is.factor)]

if (length(variaveis_factor) == 0)
  stop("Não existem variáveis qualitativas disponíveis.")

VD <- select.list(
  variaveis_factor,
  multiple = FALSE,
  title = "Selecione a variável dependente (binária):"
)

if (length(VD) == 0)
  stop("Nenhuma variável dependente selecionada.")

if (nlevels(dados[[VD]]) != 2)
  stop("A regressão logística requer variável dependente binária.")

cat("\nVariável dependente:", VD, "\n")


# ----------------------------------------------------------------------------
# 5. DEFINIÇÃO DAS VARIÁVEIS INDEPENDENTES
# ----------------------------------------------------------------------------

variaveis_disponiveis <- setdiff(names(dados), VD)

VI <- select.list(
  variaveis_disponiveis,
  multiple = TRUE,
  title = "Selecione as variáveis independentes:"
)

if (length(VI) == 0)
  stop("Nenhuma variável independente selecionada.")

VI_factors <- VI[sapply(dados[VI], is.factor)]
vars_quant <- VI[sapply(dados[VI], is.numeric)]

cat("\nVariáveis independentes selecionadas:\n")
print(VI)


# ----------------------------------------------------------------------------
# 6.1 TABELA DE FREQUÊNCIAS (APENAS CATEGÓRICAS)
# ----------------------------------------------------------------------------

if (length(VI_factors) > 0) {
  
  Frequencias <- do.call(rbind, lapply(VI_factors, function(var) {
    
    tab <- table(dados[[var]])
    
    resultado <- data.frame(
      Categoria = names(tab),
      N_perc = paste0(tab, " (", round(prop.table(tab)*100,1), ")"),
      stringsAsFactors = FALSE
    )
    
    linha_var <- data.frame(
      Categoria = var,
      N_perc = "",
      stringsAsFactors = FALSE
    )
    
    rbind(linha_var, resultado)
  }))
  
  colnames(Frequencias)[2] <- "N (%)"
  
  cat("\nTabela de Frequências:\n")
  print(Frequencias)
}


# ----------------------------------------------------------------------------
# 6.2 TABELA CRUZADA
# ----------------------------------------------------------------------------

library(gt)
library(dplyr)

if (length(VI_factors) == 0) {
  stop("Nenhuma variável independente categórica (factor) foi selecionada.")
}

tipo_percent <- menu(
  c("Percentagem por linha", "Percentagem por coluna"),
  title = "Como pretende calcular as percentagens?"
)

# Níveis fixos da VD para garantir mesmas colunas em todos os blocos
niveis_vd <- levels(dados[[VD]])

# ----------------------------------------------------------------------------
# CONSTRUÇÃO DOS BLOCOS CATEGÓRICOS
# ----------------------------------------------------------------------------
blocos_cat <- lapply(VI_factors, function(var) {
  
  dados_temp <- dados[, c(var, VD)]
  names(dados_temp) <- c("VI", "VD")
  dados_temp <- na.omit(dados_temp)
  
  # Fixar níveis para manter estrutura consistente
  dados_temp$VI <- factor(dados_temp$VI, levels = levels(dados[[var]]))
  dados_temp$VD <- factor(dados_temp$VD, levels = niveis_vd)
  
  tab <- table(dados_temp$VI, dados_temp$VD, dnn = NULL)
  
  # Remover categorias sem observações
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  
  if (nrow(tab) == 0) return(NULL)
  
  # ------------------------------------------------------------
  # p-value
  # ------------------------------------------------------------
  p_val_fmt <- ""
  
  if (nrow(tab) >= 2 && ncol(tab) == 2 && all(colSums(tab) > 0)) {
    qui_obj <- suppressWarnings(chisq.test(tab, correct = FALSE))
    
    teste <- if (any(qui_obj$expected < 5)) {
      fisher.test(tab)
    } else {
      qui_obj
    }
    
    p_val <- teste$p.value
    p_val_fmt <- ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))
  }
  
  # ------------------------------------------------------------
  # N (%) total por categoria da VI
  # ------------------------------------------------------------
  freq_total <- rowSums(tab)
  N_total_var <- paste0(
    freq_total, " (",
    round(freq_total / sum(freq_total) * 100, 1), ")"
  )
  
  # ------------------------------------------------------------
  # Percentagem por linha ou coluna
  # ------------------------------------------------------------
  if (tipo_percent == 1) {
    percent_matrix <- prop.table(tab, 1) * 100
  } else {
    percent_matrix <- prop.table(tab, 2) * 100
  }
  
  # ------------------------------------------------------------
  # Colunas da VD no formato N (%)
  # ------------------------------------------------------------
  vd_formatado <- as.data.frame(
    lapply(seq_len(ncol(tab)), function(j) {
      paste0(tab[, j], " (", round(percent_matrix[, j], 1), ")")
    }),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  colnames(vd_formatado) <- colnames(tab)
  
  # ------------------------------------------------------------
  # Tabela da variável
  # ------------------------------------------------------------
  df_tab <- cbind(
    data.frame(
      Categoria = rownames(tab),
      `N (%)` = N_total_var,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    vd_formatado
  )
  
  df_tab$p <- ""
  df_tab$p[1] <- p_val_fmt
  
  # ------------------------------------------------------------
  # Linha com nome da variável
  # ------------------------------------------------------------
  linha_var <- as.data.frame(
    setNames(as.list(rep("", ncol(df_tab))), names(df_tab)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  linha_var$Categoria <- var
  
  bind_rows(linha_var, df_tab)
})

# Remover blocos nulos
blocos_cat <- blocos_cat[!sapply(blocos_cat, is.null)]

if (length(blocos_cat) == 0) {
  stop("Não foi possível construir a tabela cruzada.")
}

# ----------------------------------------------------------------------------
# JUNTAR TODOS OS BLOCOS
# ----------------------------------------------------------------------------
tabela_cruzada <- bind_rows(blocos_cat)
rownames(tabela_cruzada) <- NULL

# ----------------------------------------------------------------------------
# FORMATAR COM gt
# ----------------------------------------------------------------------------
colunas_vd <- setdiff(colnames(tabela_cruzada), c("Categoria", "N (%)", "p"))

linhas_p_sig <- tabela_cruzada$p == "<0.001" |
  suppressWarnings(as.numeric(tabela_cruzada$p)) < 0.05

tabela_cruzada_gt <- tabela_cruzada %>%
  gt() %>%
  tab_header(
    title = md("**Tabela Cruzada**")
  ) %>%
  tab_spanner(
    label = VD,
    columns = all_of(colunas_vd)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Categoria,
      rows = Categoria %in% VI_factors
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p,
      rows = linhas_p_sig
    )
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_options(
    table.font.size = 12
  )

# ----------------------------------------------------------------------------
# PRINT FINAL
# ----------------------------------------------------------------------------
print(tabela_cruzada_gt)


# ----------------------------------------------------------------------------
# 6.3 REGRESSÃO LOGÍSTICA UNIVARIADA (REFERÊNCIA = 1ª CATEGORIA)
# ----------------------------------------------------------------------------

resp_reg <- readline("Deseja realizar regressão logística univariada? (s/n): ")

if (tolower(resp_reg) == "s") {
  
  resp_var <- readline("Usar apenas variáveis com p<0.05 na tabela cruzada? (s/n): ")
  
  p_numeric <- suppressWarnings(as.numeric(tabela_cruzada$p))
  
  variaveis_sig <- tabela_cruzada$Categoria[
    tabela_cruzada$p == "<0.001" | p_numeric < 0.05
  ]
  
  variaveis_sig <- intersect(variaveis_sig, VI_factors)
  
  if (tolower(resp_var) == "s") {
    variaveis_reg <- variaveis_sig
  } else {
    variaveis_reg <- VI_factors
  }
  
  # Criar colunas
  tabela_cruzada$`OR (95% CI)` <- ""
  tabela_cruzada$p_reg <- ""
  
  if (length(variaveis_reg) > 0) {
    
    for (var in variaveis_reg) {
      
      # Garantir que a primeira categoria é referência
      dados[[var]] <- relevel(dados[[var]], ref = levels(dados[[var]])[1])
      
      modelo <- glm(
        as.formula(paste(VD, "~", var)),
        data = dados,
        family = binomial
      )
      
      resumo <- summary(modelo)
      coeficientes <- resumo$coefficients
      
      categorias <- levels(dados[[var]])
      
      # Linha da variável (título)
      linha_var <- which(tabela_cruzada$Categoria == var)[1]
      
      # Primeira categoria = referência
      linha_ref <- which(tabela_cruzada$Categoria == categorias[1])
      if (length(linha_ref) > 0) {
        tabela_cruzada$`OR (95% CI)`[linha_ref] <- "Referência"
        tabela_cruzada$p_reg[linha_ref] <- ""
      }
      
      # Outras categorias
      for (i in 2:length(categorias)) {
        
        nome_coef <- paste0(var, categorias[i])
        
        if (nome_coef %in% rownames(coeficientes)) {
          
          beta <- coef(modelo)[nome_coef]
          OR   <- exp(beta)
          IC   <- exp(confint.default(modelo)[nome_coef, ])
          p    <- coeficientes[nome_coef, 4]
          
          linha_cat <- which(tabela_cruzada$Categoria == categorias[i])
          
          if (length(linha_cat) > 0) {
            
            tabela_cruzada$`OR (95% CI)`[linha_cat] <-
              paste0(round(OR, 2), " (",
                     round(IC[1], 2), "-",
                     round(IC[2], 2), ")")
            
            tabela_cruzada$p_reg[linha_cat] <-
              ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
          }
        }
      }
    }
  }
  
  cat("\nRegressão univariada concluída.\n")
  
} else {
  
  tabela_cruzada$`OR (95% CI)` <- ""
  tabela_cruzada$p_reg <- ""
  
  cat("\nRegressão não realizada.\n")
}

# ----------------------------------------------------------------------------
# REIMPRIMIR TABELA COM REGRESSÃO (MESMA ESTRUTURA)
# ----------------------------------------------------------------------------

colunas_vd <- setdiff(colnames(tabela_cruzada),
                      c("Categoria", "N (%)", "p",
                        "OR (95% CI)", "p_reg"))

linhas_p_sig <- tabela_cruzada$p == "<0.001" |
  suppressWarnings(as.numeric(tabela_cruzada$p)) < 0.05

linhas_reg_sig <- tabela_cruzada$p_reg == "<0.001" |
  suppressWarnings(as.numeric(tabela_cruzada$p_reg)) < 0.05

tabela_regressao_gt <- tabela_cruzada %>%
  gt() %>%
  tab_header(
    title = md("**Tabela Cruzada com Regressão Logística Univariada**")
  ) %>%
  tab_spanner(
    label = VD,
    columns = all_of(colunas_vd)
  ) %>%
  tab_spanner(
    label = "Regressão Logística",
    columns = c("OR (95% CI)", "p_reg")
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Categoria,
      rows = Categoria %in% VI_factors
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p,
      rows = linhas_p_sig
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p_reg,
      rows = linhas_reg_sig
    )
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_options(
    table.font.size = 12
  )

print(tabela_regressao_gt)

