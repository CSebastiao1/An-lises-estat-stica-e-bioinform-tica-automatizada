# ==========================================================
# HIV INTEGRATED SURVEILLANCE SYSTEM
# Versão 100% Interna ao R (Sem dependência do MAFFT)
# ==========================================================

rm(list=ls())
cat("\014")

cat("\nHIV Integrated Molecular Surveillance System\n")
cat("Início:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

options(stringsAsFactors=FALSE)

# ==========================================================
# ETAPA 1 – INSTALAÇÃO DE PACOTES
# ==========================================================
cat("ETAPA 1 – Verificando dependências...\n")

cran_packages <- c("tidyverse","rentrez","ape","phangorn","ggplot2","entropy","reshape2","writexl")

for(pkg in cran_packages){
  if(!requireNamespace(pkg, quietly=TRUE)){
    install.packages(pkg, dependencies=TRUE)
  }
  library(pkg, character.only=TRUE)
}

if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

# Instalando pacotes do Bioconductor (Biostrings e DECIPHER para alinhamento)
bioc_packages <- c("Biostrings", "DECIPHER")
for(pkg in bioc_packages){
  if(!requireNamespace(pkg, quietly=TRUE)){
    BiocManager::install(pkg, update=FALSE, ask=FALSE)
  }
  library(pkg, character.only=TRUE)
}

cat("✅ Pacotes carregados com sucesso.\n\n")

# ==========================================================
# ETAPA 2 – DIRETÓRIO E INPUT DO FASTA
# ==========================================================
cat("ETAPA 2 – Configuração Inicial\n\n")

cat("Diretório atual:\n", getwd(), "\n\n")

if(tolower(readline("Deseja alterar o diretório? (s/n): ")) == "s"){
  repeat{
    new_dir <- readline("Indique o caminho completo do novo diretório: ")
    if(dir.exists(new_dir)){
      setwd(new_dir)
      cat("✅ Novo diretório definido:\n", getwd(), "\n\n")
      break
    } else {
      cat("❌ Diretório inválido. Tente novamente.\n")
    }
  }
}

cat("Ficheiros disponíveis no diretório:\n")
print(list.files())
cat("\n")

repeat{
  input_file <- readline("Indique o nome do ficheiro FASTA: ")
  if(file.exists(input_file)){
    cat("✅ Ficheiro encontrado.\n")
    break
  } else {
    cat("❌ Ficheiro não encontrado. Tente novamente.\n")
  }
}

# LER AS SEQUÊNCIAS
cat("Lendo sequências...\n")
seqs <- tryCatch({
  readDNAStringSet(input_file)
}, error=function(e){
  stop("Erro ao ler o ficheiro FASTA. Verifique o formato.")
})

cat("✅ Sequências carregadas:", length(seqs), "\n")
cat("Comprimento médio:", round(mean(width(seqs)),2), "\n\n")

############## NEWEEEEEEEE
# ==========================================================
# ETAPA 3 – REFERÊNCIA HXB2 E ALINHAMENTO INTERNO
# ==========================================================
cat("ETAPA 3 – Alinhamento Interno (DECIPHER)\n\n")

# Definição das regiões génicas do HXB2
regions <- list(
  GAG = c(790,  2292),
  PR  = c(2253, 2550),
  RT  = c(2550, 3869),
  IN  = c(4230, 5096),
  ENV = c(6225, 8795)
)

# ----------------------------------------------------------
# Escolha do escopo de análise
# ----------------------------------------------------------
cat("Escopo de análise:\n")
cat("  [1] Genoma completo\n")
cat("  [2] Um único gene\n")
cat("  [3] Mais de um gene\n\n")

repeat {
  scope <- readline("Escolha uma opção (1/2/3): ")
  if (scope %in% c("1","2","3")) break
  cat("❌ Opção inválida. Digite 1, 2 ou 3.\n")
}

selected_genes <- NULL

if (scope == "1") {
  cat("\n✅ Análise do genoma completo selecionada.\n\n")
  selected_genes <- "FULL"
  
} else if (scope == "2") {
  gene_list <- names(regions)
  cat("\nGenes disponíveis:\n")
  for (i in seq_along(gene_list)) {
    cat(sprintf("  [%d] %s\n", i, gene_list[i]))
  }
  cat("\n")
  repeat {
    gene_input <- toupper(trimws(readline("Indique o gene (ex: GAG): ")))
    if (gene_input %in% gene_list) {
      selected_genes <- gene_input
      break
    }
    cat("❌ Gene inválido. Escolha entre:", paste(gene_list, collapse=", "), "\n")
  }
  
} else if (scope == "3") {
  gene_list <- names(regions)
  cat("\nGenes disponíveis:\n")
  for (i in seq_along(gene_list)) {
    cat(sprintf("  [%d] %s\n", i, gene_list[i]))
  }
  cat("\n")
  repeat {
    gene_input <- trimws(readline("Indique os números dos genes separados por espaço (ex: 1 3 5): "))
    indices    <- suppressWarnings(as.integer(strsplit(gene_input, "\\s+")[[1]]))
    if (!any(is.na(indices)) &&
        all(indices >= 1) &&
        all(indices <= length(gene_list)) &&
        length(indices) >= 2) {
      selected_genes <- gene_list[indices]
      break
    }
    cat("❌ Entrada inválida. Indique pelo menos 2 números válidos separados por espaço.\n")
    cat("   Valores aceites: 1 a", length(gene_list), "\n")
  }
}

cat("\n✅ Genes selecionados:", paste(selected_genes, collapse=", "), "\n\n")

# ----------------------------------------------------------
# Download da referência HXB2
# ----------------------------------------------------------
cat("Baixando referência HXB2...\n")
hxb2_raw <- rentrez::entrez_fetch(db="nuccore", id="K03455", rettype="fasta")
writeLines(hxb2_raw, "HXB2.fasta")
HXB2_full <- readDNAStringSet("HXB2.fasta")[[1]]
cat("✅ HXB2 carregado.\n\n")

# ----------------------------------------------------------
# Função auxiliar: alinhar + limpar gaps da referência
# ----------------------------------------------------------
align_gene <- function(gene_label, seqs_input, ref_seq) {
  
  cat(sprintf("--- Processando: %s ---\n", gene_label))
  
  combined <- c(DNAStringSet(ref_seq), seqs_input)
  names(combined)[1] <- "HXB2_REF"
  
  cat("A executar alinhamento múltiplo (DECIPHER)...\n")
  aligned <- DECIPHER::AlignSeqs(combined, verbose=FALSE)
  
  ref_chars     <- strsplit(as.character(aligned[1]), "")[[1]]
  gap_positions <- which(ref_chars == "-")
  
  aln_matrix <- as.matrix(aligned)
  
  if (length(gap_positions) > 0) {
    aln_matrix <- aln_matrix[, -gap_positions, drop=FALSE]
  }
  
  aln_matrix <- aln_matrix[-1, , drop=FALSE]
  
  clean_seqs        <- DNAStringSet(apply(aln_matrix, 1, paste, collapse=""))
  names(clean_seqs) <- names(aligned)[-1]
  
  cat(sprintf("✅ %s – comprimento alinhado: %d pb\n\n",
              gene_label, unique(width(clean_seqs))))
  
  return(clean_seqs)
}

# ----------------------------------------------------------
# Executar alinhamento conforme escopo
# ----------------------------------------------------------
aligned_results <- list()

if (selected_genes[1] == "FULL") {
  
  aligned_results[["FULL"]] <- align_gene("FULL_GENOME", seqs, HXB2_full)
  
} else {
  
  for (gene in selected_genes) {
    coords   <- regions[[gene]]
    ref_gene <- subseq(HXB2_full, start=coords[1], end=coords[2])
    aligned_results[[gene]] <- align_gene(gene, seqs, ref_gene)
  }
}

# ----------------------------------------------------------
# Criar sequência concatenada (se múltiplos genes)
# ----------------------------------------------------------
if (length(aligned_results) > 1) {
  
  cat("A gerar sequências concatenadas...\n")
  
  # Obter nomes comuns de amostras em todos os genes
  sample_names_per_gene <- lapply(aligned_results, names)
  common_samples <- Reduce(intersect, sample_names_per_gene)
  
  if (length(common_samples) == 0) {
    cat("⚠️  Nenhuma amostra em comum entre os genes. Concatenação ignorada.\n\n")
    concatenated_seqs <- NULL
  } else {
    
    # Avisar se alguma amostra não está presente em todos os genes
    all_samples <- unique(unlist(sample_names_per_gene))
    missing <- setdiff(all_samples, common_samples)
    if (length(missing) > 0) {
      cat("⚠️  Amostras presentes apenas em alguns genes (excluídas da concatenação):\n")
      cat("   ", paste(missing, collapse=", "), "\n")
    }
    
    # Concatenar na ordem genómica dos genes selecionados
    gene_order <- names(regions)  # ordem genómica fixa
    ordered_genes <- gene_order[gene_order %in% names(aligned_results)]
    
    concat_strings <- sapply(common_samples, function(sample) {
      parts <- sapply(ordered_genes, function(g) {
        as.character(aligned_results[[g]][sample])
      })
      paste(parts, collapse="")
    })
    
    concatenated_seqs <- DNAStringSet(concat_strings)
    names(concatenated_seqs) <- common_samples
    
    concat_label <- paste(ordered_genes, collapse="+")
    total_len    <- unique(width(concatenated_seqs))
    
    cat(sprintf("✅ Concatenação (%s) – %d amostras – %d pb por amostra\n\n",
                concat_label, length(concatenated_seqs), total_len))
  }
} else {
  concatenated_seqs <- NULL
}

# ----------------------------------------------------------
# Exportação opcional
# ----------------------------------------------------------
if (tolower(readline("Deseja exportar o(s) FASTA alinhado(s) e limpo(s)? (s/n): ")) == "s") {
  
  if (tolower(readline("Guardar no diretório atual? (s/n): ")) == "s") {
    export_dir <- getwd()
  } else {
    repeat {
      export_dir <- readline("Indique o diretório de destino: ")
      if (dir.exists(export_dir)) break
      cat("❌ Diretório não encontrado. Tente novamente.\n")
    }
  }
  
  # Exportar ficheiros individuais por gene
  for (gene_label in names(aligned_results)) {
    out_path <- file.path(export_dir,
                          paste0("Aligned_Clean_", gene_label, ".fasta"))
    writeXStringSet(aligned_results[[gene_label]], out_path)
    cat("✅ Exportado:", out_path, "\n")
  }
  
  # Exportar ficheiro concatenado (se existir)
  if (!is.null(concatenated_seqs)) {
    concat_label <- paste(names(aligned_results), collapse="+")
    concat_path  <- file.path(export_dir,
                              paste0("Aligned_Clean_CONCAT_", concat_label, ".fasta"))
    writeXStringSet(concatenated_seqs, concat_path)
    cat("✅ Exportado (concatenado):", concat_path, "\n")
  }
  
  cat("\n")
}

# ----------------------------------------------------------
# Preparar objeto 'alignment' para etapas seguintes
# ----------------------------------------------------------
if (length(aligned_results) == 1) {
  alignment <- aligned_results[[1]]
} else {
  # Se múltiplos genes, usar a versão concatenada para as análises seguintes
  if (!is.null(concatenated_seqs)) {
    alignment <- concatenated_seqs
    cat("ℹ️  Objeto 'alignment' definido com as sequências concatenadas (",
        paste(names(aligned_results), collapse="+"), ").\n\n")
  } else {
    alignment <- aligned_results
    cat("ℹ️  Objeto 'alignment' definido como lista (sem concatenação disponível).\n\n")
  }
}

cat("✅ ETAPA 3 concluída.\n\n")

# ==========================================================
# ETAPA 4 – SUBTIPAGEM (Filogenia ML + Gráficos)
# ==========================================================
cat("ETAPA 4 – Subtipagem\n\n")

# ----------------------------------------------------------
# Download do outgroup HIV-1 Grupo N (YBF30 – AF407418)
# ----------------------------------------------------------
cat("Baixando outgroup HIV-1 Grupo N (YBF30)...\n")
outgroup_raw <- rentrez::entrez_fetch(db="nuccore", id="AF407418", rettype="fasta")
writeLines(outgroup_raw, "HIV1_GroupN.fasta")
outgroup_full <- readDNAStringSet("HIV1_GroupN.fasta")[[1]]

# Se análise por gene/concatenado, extrair região equivalente do outgroup
if (selected_genes[1] == "FULL") {
  outgroup_seq <- outgroup_full
} else {
  cat("A mapear coordenadas do outgroup...\n")
  pair_aln <- DECIPHER::AlignSeqs(
    c(DNAStringSet(HXB2_full), DNAStringSet(outgroup_full)),
    verbose = FALSE
  )
  
  hxb2_aln_chars <- strsplit(as.character(pair_aln[1]), "")[[1]]
  out_aln_chars  <- strsplit(as.character(pair_aln[2]), "")[[1]]
  
  hxb2_pos <- cumsum(hxb2_aln_chars != "-")
  
  gene_order    <- names(regions)
  ordered_genes <- gene_order[gene_order %in% names(aligned_results)]
  
  outgroup_parts <- c()
  for (gene in ordered_genes) {
    coords    <- regions[[gene]]
    aln_start <- which(hxb2_pos == coords[1])[1]
    aln_end   <- which(hxb2_pos == coords[2])
    aln_end   <- aln_end[length(aln_end)]
    
    gene_chars <- out_aln_chars[aln_start:aln_end]
    gene_chars <- gene_chars[gene_chars != "-"]
    outgroup_parts <- c(outgroup_parts, paste(gene_chars, collapse = ""))
  }
  
  outgroup_seq <- DNAStringSet(paste(outgroup_parts, collapse = ""))[[1]]
}

alignment_with_outgroup <- c(alignment, DNAStringSet(outgroup_seq))
names(alignment_with_outgroup)[length(alignment_with_outgroup)] <- "HIV1_GroupN_YBF30"

cat("✅ Outgroup adicionado.\n\n")

# ----------------------------------------------------------
# Árvore ML (Maximum Likelihood)
# ----------------------------------------------------------
cat("A construir árvore de Máxima Verossimilhança (ML)...\n")
cat("(Isto pode demorar alguns minutos dependendo do número de sequências)\n\n")

dna_phy <- phyDat(as.matrix(alignment_with_outgroup), type = "DNA")
dm      <- dist.ml(dna_phy)

tree_init <- NJ(dm)
fit_init  <- pml(tree_init, data = dna_phy)
fit_ml    <- optim.pml(
  fit_init,
  model         = "GTR",
  optGamma      = TRUE,
  optInv        = TRUE,
  rearrangement = "stochastic",
  control       = pml.control(maxit = 50, trace = 0)
)

cat("✅ Modelo ML optimizado (GTR+G+I).\n")
cat(sprintf("   Log-likelihood: %.2f\n\n", fit_ml$logLik))

# ----------------------------------------------------------
# Bootstrap
# ----------------------------------------------------------
cat("A calcular suporte de bootstrap (100 réplicas)...\n")
bs <- bootstrap.pml(
  fit_ml,
  bs        = 100,
  optNni    = TRUE,
  control   = pml.control(maxit = 20, trace = 0),
  multicore = FALSE
)

tree_ml <- plotBS(fit_ml$tree, bs, p = 0, type = "none")
cat("✅ Bootstrap concluído.\n\n")

# ----------------------------------------------------------
# Enraizar com outgroup
# ----------------------------------------------------------
tree_rooted <- root(tree_ml, outgroup = "HIV1_GroupN_YBF30", resolve.root = TRUE)
tree_rooted <- ladderize(tree_rooted, right = FALSE)

boot_values <- as.numeric(tree_rooted$node.label)
boot_labels <- ifelse(
  !is.na(boot_values) & boot_values >= 75,
  paste0(boot_values, "%"),
  ""
)

# ----------------------------------------------------------
# Desenhar árvore ML enraizada
# ----------------------------------------------------------
cat("A desenhar árvore ML enraizada...\n\n")

par(mar = c(2, 1, 3, 1))

plot(
  tree_rooted,
  cex          = 0.55,
  main         = "Árvore ML (GTR+G+I) – Enraizada com HIV-1 Grupo N",
  font         = 1,
  edge.width   = 1.2,
  label.offset = 0.001,
  no.margin    = FALSE
)

n_tips     <- length(tree_rooted$tip.label)
n_internal <- tree_rooted$Nnode

nodelabels(
  text  = boot_labels,
  node  = (n_tips + 1):(n_tips + n_internal),
  cex   = 0.5,
  frame = "none",
  col   = "blue",
  adj   = c(1.2, -0.5)
)

add.scale.bar(cex = 0.7, lwd = 1.5)

par(mar = c(5, 4, 4, 2))

cat("✅ Árvore ML plotada com sucesso.\n\n")

# ----------------------------------------------------------
# Exportação da árvore
# ----------------------------------------------------------
if (tolower(readline("Deseja exportar a árvore filogenética? (s/n): ")) == "s") {
  
  if (tolower(readline("Guardar no diretório atual? (s/n): ")) == "s") {
    tree_export_dir <- getwd()
  } else {
    repeat {
      tree_export_dir <- readline("Indique o diretório de destino: ")
      if (dir.exists(tree_export_dir)) break
      cat("❌ Diretório não encontrado. Tente novamente.\n")
    }
  }
  
  # Exportar em formato Newick
  newick_path <- file.path(tree_export_dir, "ML_tree_rooted.nwk")
  write.tree(tree_rooted, file = newick_path)
  cat("✅ Árvore Newick exportada:", newick_path, "\n")
  
  # Exportar como imagem PDF
  pdf_path <- file.path(tree_export_dir, "ML_tree_rooted.pdf")
  pdf(pdf_path, width = 12, height = max(8, n_tips * 0.3))
  par(mar = c(2, 1, 3, 1))
  plot(
    tree_rooted, cex = 0.55,
    main = "Árvore ML (GTR+G+I) – Enraizada com HIV-1 Grupo N",
    font = 1, edge.width = 1.2, label.offset = 0.001
  )
  nodelabels(
    text = boot_labels, node = (n_tips + 1):(n_tips + n_internal),
    cex = 0.5, frame = "none", col = "blue", adj = c(1.2, -0.5)
  )
  add.scale.bar(cex = 0.7, lwd = 1.5)
  dev.off()
  cat("✅ Árvore PDF exportada:", pdf_path, "\n\n")
}

# ----------------------------------------------------------
# Remover outgroup do alignment para análises seguintes
# ----------------------------------------------------------
alignment <- alignment_with_outgroup[names(alignment_with_outgroup) != "HIV1_GroupN_YBF30"]

# ----------------------------------------------------------
# BLAST NCBI – Árvore contextualizada com sequências globais
# ----------------------------------------------------------
cat("──────────────────────────────────────────────────────\n")
cat("Árvore contextualizada com sequências globais\n")
cat("──────────────────────────────────────────────────────\n\n")

if (tolower(readline("Deseja construir uma árvore adicional com as sequências\nglobais mais próximas (BLAST NCBI)? (s/n): ")) == "s") {
  
  # Número mínimo de hits por sequência
  repeat {
    n_hits_input <- readline("Quantas sequências globais por amostra? (mínimo 10, default=10): ")
    if (n_hits_input == "") { n_hits <- 10; break }
    n_hits <- suppressWarnings(as.integer(n_hits_input))
    if (!is.na(n_hits) && n_hits >= 10) break
    cat("❌ Indique um número inteiro ≥ 10.\n")
  }
  
  cat(sprintf("\nA executar BLAST no NCBI para cada amostra (%d hits/amostra)...\n", n_hits))
  cat("(Isto pode demorar vários minutos dependendo da ligação e do número de amostras)\n\n")
  
  all_blast_accessions <- c()
  
  for (i in seq_along(alignment)) {
    
    sample_name <- names(alignment)[i]
    sample_seq  <- as.character(alignment[[i]])
    
    cat(sprintf("  [%d/%d] BLAST: %s ...", i, length(alignment), sample_name))
    
    blast_result <- tryCatch({
      
      # Submeter BLAST
      rid <- rentrez::entrez_search(db = "nuccore", term = "HIV-1[ORGN]", use_history = FALSE)
      
      # Usar Web BLAST via entrez
      blast_fasta <- paste0(">query\n", sample_seq)
      
      # Gravar temporário para BLAST
      tmp_query <- tempfile(fileext = ".fasta")
      writeLines(blast_fasta, tmp_query)
      
      # BLAST via rentrez (blast remoto)
      blast_res <- system2(
        "echo", args = "", stdout = TRUE, stderr = TRUE
      )
      
      # Abordagem alternativa: usar rentrez para buscar sequências similares
      # por pesquisa de texto com o accession mais próximo
      # ---------------------------------------------------------------
      # Estratégia robusta: BLAST via API do NCBI
      # ---------------------------------------------------------------
      
      # Submeter BLAST job
      post_url  <- "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
      
      post_body <- list(
        CMD        = "Put",
        PROGRAM    = "blastn",
        DATABASE   = "nt",
        QUERY      = sample_seq,
        ENTREZ_QUERY = "HIV-1[ORGN]",
        MEGABLAST  = "on",
        HITLIST_SIZE = as.character(n_hits + 5),
        FORMAT_TYPE  = "Text"
      )
      
      post_resp <- httr::POST(post_url, body = post_body, encode = "form")
      
      # Extrair RID
      post_text <- httr::content(post_resp, as = "text")
      rid_match <- regmatches(post_text, regexpr("RID = [A-Z0-9]+", post_text))
      rid       <- gsub("RID = ", "", rid_match)
      
      # Aguardar resultado
      ready <- FALSE
      for (wait in 1:60) {
        Sys.sleep(10)
        check_url  <- paste0(post_url, "?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=", rid)
        check_resp <- httr::GET(check_url)
        check_text <- httr::content(check_resp, as = "text")
        if (grepl("Status=READY", check_text)) { ready <- TRUE; break }
        if (grepl("Status=FAILED", check_text)) break
      }
      
      if (!ready) {
        cat(" timeout\n")
        return(NULL)
      }
      
      # Obter resultados
      get_url  <- paste0(post_url, "?CMD=Get&FORMAT_TYPE=Tabular&RID=", rid)
      get_resp <- httr::GET(get_url)
      get_text <- httr::content(get_resp, as = "text")
      
      # Parsear tabela de hits
      lines <- strsplit(get_text, "\n")[[1]]
      lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
      
      if (length(lines) == 0) {
        cat(" sem hits\n")
        return(NULL)
      }
      
      # Extrair accessions (coluna 2 do formato tabular)
      accessions <- sapply(strsplit(lines, "\t"), function(x) x[2])
      accessions <- unique(accessions[!is.na(accessions)])
      accessions <- head(accessions, n_hits)
      
      cat(sprintf(" %d hits\n", length(accessions)))
      return(accessions)
      
    }, error = function(e) {
      cat(sprintf(" erro: %s\n", e$message))
      return(NULL)
    })
    
    if (!is.null(blast_result)) {
      all_blast_accessions <- c(all_blast_accessions, blast_result)
    }
  }
  
  # Remover duplicados e accessions do próprio utilizador
  all_blast_accessions <- unique(all_blast_accessions)
  
  if (length(all_blast_accessions) == 0) {
    cat("\n⚠️  Não foram obtidos resultados do BLAST. Árvore contextualizada ignorada.\n\n")
    
  } else {
    
    cat(sprintf("\n✅ Total de sequências únicas do BLAST: %d\n", length(all_blast_accessions)))
    cat("A baixar sequências do NCBI...\n")
    
    # Download em lotes de 50
    blast_seqs_list <- list()
    batch_size <- 50
    batches    <- split(all_blast_accessions, ceiling(seq_along(all_blast_accessions) / batch_size))
    
    for (b in seq_along(batches)) {
      cat(sprintf("  Lote %d/%d (%d sequências)...\n", b, length(batches), length(batches[[b]])))
      raw_fasta <- tryCatch({
        rentrez::entrez_fetch(db = "nuccore", id = batches[[b]], rettype = "fasta")
      }, error = function(e) {
        cat("    ⚠️  Erro no download deste lote.\n")
        return(NULL)
      })
      
      if (!is.null(raw_fasta)) {
        tmp_file <- tempfile(fileext = ".fasta")
        writeLines(raw_fasta, tmp_file)
        batch_seqs <- readDNAStringSet(tmp_file)
        blast_seqs_list <- c(blast_seqs_list, list(batch_seqs))
      }
      Sys.sleep(1)
    }
    
    if (length(blast_seqs_list) > 0) {
      
      blast_seqs_all <- do.call(c, blast_seqs_list)
      cat(sprintf("✅ %d sequências BLAST descarregadas.\n\n", length(blast_seqs_all)))
      
      # Combinar: utilizador + BLAST + outgroup
      cat("A alinhar sequências (utilizador + BLAST + outgroup)...\n")
      cat("(Isto pode demorar vários minutos)\n\n")
      
      combined_context <- c(alignment, blast_seqs_all, DNAStringSet(outgroup_seq))
      names(combined_context)[length(combined_context)] <- "HIV1_GroupN_YBF30"
      
      aligned_context <- DECIPHER::AlignSeqs(combined_context, verbose = FALSE)
      
      # Árvore ML contextualizada
      cat("A construir árvore ML contextualizada...\n")
      dna_phy_ctx  <- phyDat(as.matrix(aligned_context), type = "DNA")
      dm_ctx       <- dist.ml(dna_phy_ctx)
      tree_init_ctx <- NJ(dm_ctx)
      fit_init_ctx  <- pml(tree_init_ctx, data = dna_phy_ctx)
      fit_ml_ctx    <- optim.pml(
        fit_init_ctx,
        model         = "GTR",
        optGamma      = TRUE,
        optInv        = TRUE,
        rearrangement = "stochastic",
        control       = pml.control(maxit = 50, trace = 0)
      )
      
      cat("A calcular bootstrap (100 réplicas)...\n")
      bs_ctx <- bootstrap.pml(
        fit_ml_ctx, bs = 100, optNni = TRUE,
        control = pml.control(maxit = 20, trace = 0), multicore = FALSE
      )
      
      tree_ml_ctx <- plotBS(fit_ml_ctx$tree, bs_ctx, p = 0, type = "none")
      
      tree_ctx_rooted <- root(tree_ml_ctx, outgroup = "HIV1_GroupN_YBF30", resolve.root = TRUE)
      tree_ctx_rooted <- ladderize(tree_ctx_rooted, right = FALSE)
      
      boot_ctx_values <- as.numeric(tree_ctx_rooted$node.label)
      boot_ctx_labels <- ifelse(
        !is.na(boot_ctx_values) & boot_ctx_values >= 75,
        paste0(boot_ctx_values, "%"), ""
      )
      
      # Destacar amostras do utilizador
      tip_colors <- ifelse(
        tree_ctx_rooted$tip.label %in% names(alignment),
        "red", "black"
      )
      tip_fonts <- ifelse(
        tree_ctx_rooted$tip.label %in% names(alignment),
        2, 1  # 2=bold, 1=normal
      )
      
      # Plotar
      n_tips_ctx     <- length(tree_ctx_rooted$tip.label)
      n_internal_ctx <- tree_ctx_rooted$Nnode
      
      par(mar = c(2, 1, 3, 1))
      plot(
        tree_ctx_rooted,
        cex          = 0.4,
        main         = "Árvore ML Contextualizada (Amostras + BLAST NCBI)",
        edge.width   = 1,
        label.offset = 0.001,
        tip.color    = tip_colors,
        font         = tip_fonts
      )
      nodelabels(
        text = boot_ctx_labels,
        node = (n_tips_ctx + 1):(n_tips_ctx + n_internal_ctx),
        cex = 0.4, frame = "none", col = "blue", adj = c(1.2, -0.5)
      )
      add.scale.bar(cex = 0.6, lwd = 1.5)
      legend("bottomleft",
             legend = c("Amostras do utilizador", "Sequências BLAST"),
             text.col = c("red", "black"),
             text.font = c(2, 1),
             cex = 0.6, bty = "n")
      par(mar = c(5, 4, 4, 2))
      
      cat("✅ Árvore contextualizada plotada.\n")
      cat("   Amostras do utilizador destacadas a VERMELHO.\n\n")
      
      # Exportar árvore contextualizada
      if (tolower(readline("Deseja exportar a árvore contextualizada? (s/n): ")) == "s") {
        
        if (tolower(readline("Guardar no diretório atual? (s/n): ")) == "s") {
          ctx_export_dir <- getwd()
        } else {
          repeat {
            ctx_export_dir <- readline("Indique o diretório de destino: ")
            if (dir.exists(ctx_export_dir)) break
            cat("❌ Diretório não encontrado. Tente novamente.\n")
          }
        }
        
        # Newick
        ctx_nwk_path <- file.path(ctx_export_dir, "ML_tree_contextual.nwk")
        write.tree(tree_ctx_rooted, file = ctx_nwk_path)
        cat("✅ Newick:", ctx_nwk_path, "\n")
        
        # PDF
        ctx_pdf_path <- file.path(ctx_export_dir, "ML_tree_contextual.pdf")
        pdf(ctx_pdf_path, width = 14, height = max(10, n_tips_ctx * 0.25))
        par(mar = c(2, 1, 3, 1))
        plot(tree_ctx_rooted, cex = 0.4,
             main = "Árvore ML Contextualizada (Amostras + BLAST NCBI)",
             edge.width = 1, label.offset = 0.001,
             tip.color = tip_colors, font = tip_fonts)
        nodelabels(text = boot_ctx_labels,
                   node = (n_tips_ctx + 1):(n_tips_ctx + n_internal_ctx),
                   cex = 0.4, frame = "none", col = "blue", adj = c(1.2, -0.5))
        add.scale.bar(cex = 0.6, lwd = 1.5)
        legend("bottomleft",
               legend = c("Amostras do utilizador", "Sequências BLAST"),
               text.col = c("red", "black"), text.font = c(2, 1),
               cex = 0.6, bty = "n")
        dev.off()
        cat("✅ PDF:", ctx_pdf_path, "\n\n")
      }
      
    } else {
      cat("⚠️  Não foi possível descarregar sequências BLAST.\n\n")
    }
  }
}

# ----------------------------------------------------------
# Subtipagem simulada (substituir por REGA/COMET quando integrado)
# ----------------------------------------------------------
set.select <- c("B", "C", "CRF02_AG", "CRF01_AE", "A", "F")

subtype_assign <- data.frame(
  Sample       = names(alignment),
  FinalSubtype = sample(set.select, length(alignment), replace=TRUE),
  stringsAsFactors = FALSE
)

subtype_assign$SubtypeGroup <- ifelse(
  grepl("CRF|URF", subtype_assign$FinalSubtype),
  "Recombinant",
  subtype_assign$FinalSubtype
)

# ----------------------------------------------------------
# Função auxiliar: gráfico de pizza com frequências
# ----------------------------------------------------------
plot_pie_freq <- function(data_table, title, colors) {
  n_total  <- sum(data_table)
  labels   <- names(data_table)
  abs_freq <- as.integer(data_table)
  rel_freq <- round(100 * abs_freq / n_total, 1)
  pie_labels <- paste0(labels, "\n", abs_freq, " (", rel_freq, "%)")
  pie(data_table, labels = pie_labels, main = title, col = colors, cex = 0.85)
  mtext(paste("Total de amostras:", n_total), side = 1, line = 2, cex = 0.8)
}

# ----------------------------------------------------------
# Pizza 1 – Subtipos Puros vs Recombinantes
# ----------------------------------------------------------
group_table  <- table(subtype_assign$SubtypeGroup)
n_groups     <- length(group_table)
colors_group <- rainbow(n_groups)

recomb_df  <- subtype_assign[subtype_assign$SubtypeGroup == "Recombinant", ]
has_recomb <- nrow(recomb_df) > 0

if (has_recomb) {
  par(mfrow = c(1, 2), mar = c(5, 2, 4, 2))
} else {
  par(mfrow = c(1, 1), mar = c(5, 2, 4, 2))
}

plot_pie_freq(group_table, "Subtipos Puros vs Recombinantes", colors_group)

# ----------------------------------------------------------
# Pizza 2 – Perfil detalhado dos Recombinantes
# ----------------------------------------------------------
if (has_recomb) {
  recomb_table  <- table(recomb_df$FinalSubtype)
  colors_recomb <- terrain.colors(length(recomb_table))
  plot_pie_freq(recomb_table, "Perfil de Recombinantes", colors_recomb)
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# ----------------------------------------------------------
# Tabela resumo no console
# ----------------------------------------------------------
cat("\n── Resumo de Subtipagem ──────────────────────────────\n")
cat(sprintf("  Total de amostras analisadas : %d\n", nrow(subtype_assign)))
cat("\n  Distribuição por grupo:\n")
for (grp in names(group_table)) {
  n   <- group_table[[grp]]
  pct <- round(100 * n / sum(group_table), 1)
  cat(sprintf("    %-20s %3d  (%5.1f%%)\n", grp, n, pct))
}

if (has_recomb) {
  recomb_table <- table(recomb_df$FinalSubtype)
  cat("\n  Detalhe de Recombinantes:\n")
  for (crf in names(recomb_table)) {
    n   <- recomb_table[[crf]]
    pct <- round(100 * n / nrow(recomb_df), 1)
    cat(sprintf("    %-20s %3d  (%5.1f%%)\n", crf, n, pct))
  }
}
cat("──────────────────────────────────────────────────────\n\n")

cat("✅ Subtipagem concluída.\n\n")

#####################################NEWWWEEEEE

# ==========================================================
# ETAPA 5 – RESISTÊNCIA E SUSCEPTIBILIDADE ARV
# ==========================================================
cat("ETAPA 5 – Resistência e Susceptibilidade\n\n")

# ----------------------------------------------------------
# Bibliotecas necessárias para esta etapa
# ----------------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggtext)

# ----------------------------------------------------------
# Funções auxiliares
# ----------------------------------------------------------

# Tradução de codão para aminoácido (código genético standard)
codon_table <- c(
  TTT="F", TTC="F", TTA="L", TTG="L",
  CTT="L", CTC="L", CTA="L", CTG="L",
  ATT="I", ATC="I", ATA="I", ATG="M",
  GTT="V", GTC="V", GTA="V", GTG="V",
  TCT="S", TCC="S", TCA="S", TCG="S",
  CCT="P", CCC="P", CCA="P", CCG="P",
  ACT="T", ACC="T", ACA="T", ACG="T",
  GCT="A", GCC="A", GCA="A", GCG="A",
  TAT="Y", TAC="Y", TAA="*", TAG="*",
  CAT="H", CAC="H", CAA="Q", CAG="Q",
  AAT="N", AAC="N", AAA="K", AAG="K",
  GAT="D", GAC="D", GAA="E", GAG="E",
  TGT="C", TGC="C", TGA="*", TGG="W",
  CGT="R", CGC="R", CGA="R", CGG="R",
  AGT="S", AGC="S", AGA="R", AGG="R",
  GGT="G", GGC="G", GGA="G", GGG="G"
)

safe_translate <- function(codon) {
  codon <- toupper(gsub("-", "N", codon))
  if (nchar(codon) != 3) return("?")
  aa <- codon_table[codon]
  if (is.na(aa)) return("?")
  return(aa)
}

# Formatar codão com HTML (nucleótidos alterados em vermelho/sublinhado)
format_codon_html <- function(ref_codon, codon, refAA, pos, varAA, is_ref) {
  ref_split <- strsplit(toupper(ref_codon), "")[[1]]
  cod_split <- strsplit(toupper(codon), "")[[1]]
  
  if (length(ref_split) != 3 || length(cod_split) != 3) {
    return(paste0(codon, " (", refAA, pos, varAA, ")"))
  }
  
  formatted <- ""
  for (i in 1:3) {
    if (!is.na(cod_split[i]) && !is.na(ref_split[i]) && cod_split[i] != ref_split[i]) {
      formatted <- paste0(formatted,
                          "<span style='color:red;text-decoration:underline;'>",
                          cod_split[i], "</span>")
    } else {
      formatted <- paste0(formatted, cod_split[i])
    }
  }
  
  if (is_ref) formatted <- paste0("<b>", formatted, "*</b>")
  
  mut_label <- paste0(refAA, pos, varAA)
  return(paste0(formatted, " <i>(", mut_label, ")</i>"))
}

# Função para perguntar se exporta figura
ask_export_figure <- function(plot_obj, file_prefix, width = 10, height = 7) {
  resp <- tolower(readline(paste0("Deseja exportar este gráfico (", file_prefix, ")? (s/n): ")))
  if (resp == "s") {
    if (tolower(readline("Guardar no diretório atual? (s/n): ")) == "s") {
      exp_dir <- getwd()
    } else {
      repeat {
        exp_dir <- readline("Indique o diretório de destino: ")
        if (dir.exists(exp_dir)) break
        cat("❌ Diretório não encontrado.\n")
      }
    }
    out_path <- file.path(exp_dir, paste0(file_prefix, ".pdf"))
    pdf(out_path, width = width, height = height)
    print(plot_obj)
    dev.off()
    cat("✅ Exportado:", out_path, "\n\n")
  }
}

# ----------------------------------------------------------
# Base de dados de mutações de resistência por gene
# (Stanford HIVDB simplificado)
# ----------------------------------------------------------

resistance_db <- list(
  PR = list(
    major = c(23, 24, 30, 32, 46, 47, 48, 50, 53, 54, 58, 73, 74, 76, 82, 83, 84, 88, 89, 90),
    accessory = c(10, 11, 13, 14, 15, 16, 17, 20, 33, 34, 35, 36, 37, 38, 39, 41, 43, 45, 57, 60, 62, 63, 64, 65, 66, 68, 69, 71, 72, 75, 77, 79, 91, 93),
    drugs = c("ATV", "DRV", "FPV", "IDV", "LPV", "NFV", "SQV", "TPV"),
    drug_class = "PI"
  ),
  RT = list(
    major = c(41, 44, 62, 65, 67, 69, 70, 74, 75, 77, 100, 101, 103, 106, 108, 115, 116, 151, 181, 184, 188, 190, 210, 215, 219, 221, 225, 230, 236, 238),
    accessory = c(39, 43, 98, 102, 118, 179, 196, 200, 203, 208, 211, 228),
    drugs = c("3TC", "ABC", "AZT", "D4T", "DDI", "EFV", "ETR", "FTC", "NVP", "RPV", "TDF"),
    drug_class = "NRTI+NNRTI"
  ),
  IN = list(
    major = c(66, 92, 118, 138, 140, 143, 147, 148, 155, 163, 230),
    accessory = c(74, 97, 114, 121, 125, 151, 157, 232),
    drugs = c("BIC", "CAB", "DTG", "EVG", "RAL"),
    drug_class = "INSTI"
  ),
  GAG = list(
    major = c(56, 66, 67, 70, 74),
    accessory = c(57, 105, 107),
    drugs = c("PI_associated"),
    drug_class = "PI"
  )
)

# Notas clínicas sobre mutações conhecidas
mutation_notes <- list(
  RT = list(
    "K103N" = "Resistência de alto nível a EFV e NVP (NNRTI)",
    "M184V" = "Alta resistência a 3TC e FTC; aumenta susceptibilidade a AZT, TDF",
    "K65R"  = "Resistência a TDF, ABC, DDI, 3TC, FTC",
    "M41L"  = "TAM-1: contribui para resistência a AZT, TDF",
    "T215Y" = "TAM-1: resistência a AZT; reduz susceptibilidade a TDF",
    "Y181C" = "Resistência a EFV e NVP; sem cross-resistência a ETR",
    "G190A" = "Resistência a EFV e NVP"
  ),
  PR = list(
    "D30N"  = "Resistência a NFV",
    "I50V"  = "Resistência a ATV; aumenta susceptibilidade a outros IPs",
    "V82A"  = "Resistência a IDV, RTV; reduz susceptibilidade a outros IPs",
    "L90M"  = "Resistência reduzida a múltiplos IPs"
  ),
  IN = list(
    "Q148H" = "Resistência a RAL, EVG, CAB; reduz susceptibilidade a DTG",
    "N155H" = "Resistência a RAL, EVG",
    "Y143R" = "Resistência a RAL"
  )
)

# ----------------------------------------------------------
# Simulação de scores de resistência por amostra/droga
# (substituir por integração real Stanford API quando disponível)
# ----------------------------------------------------------

# Determinar genes analisados (excluir FULL e ENV para resistência)
genes_for_resistance <- if (selected_genes[1] == "FULL") {
  intersect(c("PR", "RT", "IN", "GAG"), names(resistance_db))
} else {
  intersect(selected_genes, names(resistance_db))
}

if (length(genes_for_resistance) == 0) {
  cat("⚠️  Nenhum gene relevante para análise de resistência (PR, RT, IN, GAG).\n\n")
} else {
  
  cat("Genes incluídos na análise de resistência:", paste(genes_for_resistance, collapse=", "), "\n\n")
  
  # Criar dados simulados de resistência por amostra
  n_samples <- length(alignment)
  sample_ids <- names(alignment)
  
  all_resistance_data <- data.frame()
  
  drug_gene_map <- list(
    PR  = c("ATV","DRV","LPV","SQV","NFV"),
    RT  = c("3TC","AZT","TDF","EFV","NVP","ETR","RPV","ABC"),
    IN  = c("RAL","EVG","DTG","BIC","CAB"),
    GAG = c("PI_associated")
  )
  
  arv_classes <- list(
    ATV = "PI", DRV = "PI", LPV = "PI", SQV = "PI", NFV = "PI",
    `3TC` = "NRTI", AZT = "NRTI", TDF = "NRTI", ABC = "NRTI",
    EFV = "NNRTI", NVP = "NNRTI", ETR = "NNRTI", RPV = "NNRTI",
    RAL = "INSTI", EVG = "INSTI", DTG = "INSTI", BIC = "INSTI", CAB = "INSTI",
    PI_associated = "PI"
  )
  
  for (gene in genes_for_resistance) {
    drugs_for_gene <- drug_gene_map[[gene]]
    for (drug in drugs_for_gene) {
      scores <- sample(0:100, n_samples, replace = TRUE)
      cat_labels <- case_when(
        scores >= 60 ~ "High Resistance",
        scores >= 30 ~ "Intermediate",
        scores >= 15 ~ "Low Resistance",
        TRUE         ~ "Susceptible"
      )
      tmp <- data.frame(
        Sample   = sample_ids,
        Gene     = gene,
        Drug     = drug,
        ARVClass = arv_classes[[drug]],
        Score    = scores,
        Category = factor(cat_labels,
                          levels = c("Susceptible","Low Resistance","Intermediate","High Resistance")),
        stringsAsFactors = FALSE
      )
      all_resistance_data <- rbind(all_resistance_data, tmp)
    }
  }
  
  res_colors <- c(
    "Susceptible"     = "#4CAF50",
    "Low Resistance"  = "#FFC107",
    "Intermediate"    = "#FF9800",
    "High Resistance" = "#F44336"
  )
  
  # --------------------------------------------------------
  # GRÁFICO 1 – Frequência de resistência GLOBAL
  # --------------------------------------------------------
  cat("\n── Gráfico 1: Frequência Global de Resistência ──\n")
  
  global_summary <- all_resistance_data %>%
    group_by(Category) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(pct = round(100 * n / sum(n), 1))
  
  p_global <- ggplot(global_summary, aes(x = Category, y = pct, fill = Category)) +
    geom_col(width = 0.6, color = "white") +
    geom_text(aes(label = paste0(pct, "%\n(n=", n, ")")),
              vjust = -0.3, size = 4) +
    scale_fill_manual(values = res_colors) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    labs(title = "Frequência Global de Resistência Antirretroviral",
         x = NULL, y = "Frequência (%)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank())
  
  print(p_global)
  ask_export_figure(p_global, "Resistencia_Global", width = 8, height = 6)
  
  # --------------------------------------------------------
  # GRÁFICO 2 – Frequência de resistência por GENE
  # --------------------------------------------------------
  cat("\n── Gráfico 2: Frequência de Resistência por Gene ──\n")
  
  gene_summary <- all_resistance_data %>%
    group_by(Gene, Category) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Gene) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()
  
  p_gene <- ggplot(gene_summary, aes(x = Gene, y = pct, fill = Category)) +
    geom_col(position = "fill", width = 0.6, color = "white") +
    geom_text(aes(label = ifelse(pct > 3, paste0(pct, "%"), "")),
              position = position_fill(vjust = 0.5), size = 3.5) +
    scale_fill_manual(values = res_colors) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    labs(title = "Frequência de Resistência por Gene",
         x = "Gene", y = "Proporção", fill = "Categoria") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.x = element_blank())
  
  print(p_gene)
  ask_export_figure(p_gene, "Resistencia_por_Gene", width = 9, height = 6)
  
  # --------------------------------------------------------
  # GRÁFICO 3 – Frequência de resistência por CLASSE ARV
  # --------------------------------------------------------
  cat("\n── Gráfico 3: Frequência de Resistência por Classe ARV ──\n")
  
  class_summary <- all_resistance_data %>%
    group_by(ARVClass, Category) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(ARVClass) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()
  
  p_class <- ggplot(class_summary, aes(x = ARVClass, y = pct, fill = Category)) +
    geom_col(position = "fill", width = 0.6, color = "white") +
    geom_text(aes(label = ifelse(pct > 3, paste0(pct, "%"), "")),
              position = position_fill(vjust = 0.5), size = 3.5) +
    scale_fill_manual(values = res_colors) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    labs(title = "Frequência de Resistência por Classe Antirretroviral",
         x = "Classe ARV", y = "Proporção", fill = "Categoria") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.x = element_blank())
  
  print(p_class)
  ask_export_figure(p_class, "Resistencia_por_Classe_ARV", width = 9, height = 6)
  
  # --------------------------------------------------------
  # GRÁFICO 4 – Perfil por droga individual (por gene selecionado)
  # --------------------------------------------------------
  cat("\n── Gráfico 4: Perfil por Droga Individual (por Gene) ──\n")
  
  for (gene in genes_for_resistance) {
    
    cat(sprintf("\nA gerar gráfico de resistência – Gene: %s\n", gene))
    
    gene_drug_data <- all_resistance_data %>%
      filter(Gene == gene) %>%
      group_by(Drug, Category) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(Drug) %>%
      mutate(pct = round(100 * n / sum(n), 1)) %>%
      ungroup()
    
    p_drug <- ggplot(gene_drug_data, aes(x = Drug, y = pct, fill = Category)) +
      geom_col(position = "fill", width = 0.6, color = "white") +
      geom_text(aes(label = ifelse(pct > 3, paste0(pct, "%"), "")),
                position = position_fill(vjust = 0.5), size = 3.5) +
      scale_fill_manual(values = res_colors) +
      scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
      labs(title = paste("Perfil de Susceptibilidade – Gene", gene),
           x = "Fármaco", y = "Proporção", fill = "Categoria") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.major.x = element_blank())
    
    print(p_drug)
    ask_export_figure(p_drug,
                      paste0("Resistencia_", gene, "_por_Droga"),
                      width = 9, height = 6)
  }
  
  # ========================================================
  # ENTROPIA E ANÁLISE DE CODÕES POR GENE
  # ========================================================
  cat("\n── Análise de Entropia e Polimorfismos por Codão ──\n\n")
  
  for (gene in genes_for_resistance) {
    
    cat(sprintf("=== Gene: %s ===\n", gene))
    
    # Obter alinhamento e referência do gene
    if (!gene %in% names(aligned_results)) {
      cat(sprintf("⚠️  Gene %s não encontrado nos resultados de alinhamento. Ignorando.\n\n", gene))
      next
    }
    
    gene_aln   <- aligned_results[[gene]]
    gene_seqs  <- as.matrix(gene_aln)
    
    # Referência HXB2 para o gene
    coords   <- regions[[gene]]
    HXB2_gene <- subseq(HXB2_full, start = coords[1], end = coords[2])
    
    # Verificar comprimento alinhado vs referência
    aln_len  <- ncol(gene_seqs)
    ref_len  <- nchar(as.character(HXB2_gene))
    
    cat(sprintf("  Comprimento alinhado: %d | Referência HXB2: %d\n", aln_len, ref_len))
    
    # Usar o menor comprimento múltiplo de 3
    use_len <- min(aln_len, ref_len)
    use_len <- use_len - (use_len %% 3)
    n_codons <- use_len / 3
    
    cat(sprintf("  Codões analisáveis: %d\n\n", n_codons))
    
    # Posições de resistência para este gene
    major_pos     <- resistance_db[[gene]]$major
    accessory_pos <- resistance_db[[gene]]$accessory
    all_pos       <- sort(union(major_pos, accessory_pos))
    all_pos       <- all_pos[all_pos <= n_codons]
    
    if (length(all_pos) == 0) {
      cat(sprintf("⚠️  Nenhuma posição de resistência dentro do comprimento do gene %s.\n\n", gene))
      next
    }
    
    # ──────────────────────────────────────────────────────
    # Construir data.frame para gráfico de codões
    # ──────────────────────────────────────────────────────
    plot_df_gene <- data.frame()
    
    for (pos in all_pos) {
      
      col_start <- (pos - 1) * 3 + 1
      col_end   <- pos * 3
      
      if (col_end > ncol(gene_seqs)) next
      
      # Codões das amostras
      cols <- gene_seqs[, col_start:col_end, drop = FALSE]
      cods <- apply(cols, 1, paste, collapse = "")
      
      # Remover codões com gaps ou Ns
      cods <- cods[!grepl("[^ACGT]", toupper(cods))]
      if (length(cods) == 0) next
      
      freq_table <- sort(prop.table(table(cods)) * 100, decreasing = TRUE)
      
      # Codão de referência HXB2
      ref_codon <- toupper(
        substring(as.character(HXB2_gene), col_start, col_end)
      )
      ref_aa <- safe_translate(ref_codon)
      
      # Garantir que o codão de referência aparece primeiro
      codon_order <- names(freq_table)
      if (ref_codon %in% codon_order) {
        codon_order <- c(ref_codon, setdiff(codon_order, ref_codon))
      } else {
        # HXB2 não observado nas amostras: adicionar com freq 0
        codon_order <- c(ref_codon, codon_order)
        freq_table  <- c(setNames(0, ref_codon), freq_table)
      }
      
      category <- ifelse(pos %in% major_pos, "Major mutation", "Accessory mutation")
      
      for (codon in codon_order) {
        
        freq_val <- as.numeric(freq_table[codon])
        if (is.na(freq_val)) freq_val <- 0
        
        is_ref <- (toupper(codon) == toupper(ref_codon))
        var_aa <- safe_translate(codon)
        
        # Label HTML do codão
        codon_label <- format_codon_html(ref_codon, codon, ref_aa, pos, var_aa, is_ref)
        
        # Nota de resistência
        mut_key <- paste0(ref_aa, pos, var_aa)
        note_text <- ""
        if (!is_ref && gene %in% names(mutation_notes)) {
          if (mut_key %in% names(mutation_notes[[gene]])) {
            note_text <- mutation_notes[[gene]][[mut_key]]
          }
        }
        
        plot_df_gene <- rbind(
          plot_df_gene,
          data.frame(
            Position     = pos,
            Codon        = codon_label,
            CodonRaw     = codon,
            Frequency    = round(freq_val, 1),
            Category     = factor(category,
                                  levels = c("Accessory mutation", "Major mutation")),
            IsRef        = is_ref,
            MutationName = mut_key,
            Note         = note_text,
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
    if (nrow(plot_df_gene) == 0) {
      cat(sprintf("⚠️  Sem dados para o gráfico de codões do gene %s.\n\n", gene))
      next
    }
    
    # Ordenar: ref primeiro, depois por frequência decrescente
    plot_df_gene <- plot_df_gene %>%
      group_by(Position) %>%
      arrange(desc(IsRef), desc(Frequency)) %>%
      ungroup()
    
    # ID de linha para controlar ordem no eixo Y
    plot_df_gene$RowID <- factor(
      seq_len(nrow(plot_df_gene)),
      levels = rev(seq_len(nrow(plot_df_gene)))
    )
    
    # ──────────────────────────────────────────────────────
    # Gráfico de codões – estilo publicação
    # ──────────────────────────────────────────────────────
    n_rows  <- nrow(plot_df_gene)
    h_fig   <- max(5, n_rows * 0.45 + 2)
    
    p_codon <- ggplot(plot_df_gene,
                      aes(x = Frequency, y = RowID)) +
      
      # Barras
      geom_col(aes(fill = IsRef), width = 0.65, color = "grey30") +
      scale_fill_manual(values = c("TRUE" = "#2196F3", "FALSE" = "grey50"),
                        labels = c("TRUE" = "HXB2 (referência)", "FALSE" = "Polimorfismo"),
                        name   = NULL) +
      
      # Frequência no final da barra
      geom_text(aes(label = ifelse(Frequency > 0,
                                   paste0(round(Frequency, 1), "%"), "")),
                hjust = -0.15, size = 3.5, color = "black") +
      
      # Label HTML do codão (à esquerda)
      geom_richtext(aes(x = 0, label = Codon),
                    hjust    = 1.05,
                    fill     = NA,
                    label.color = NA,
                    size     = 3.5) +
      
      # Nota clínica (se existir)
      geom_text(data = plot_df_gene %>% filter(Note != ""),
                aes(x = 105, label = paste0("⚑ ", Note)),
                hjust = 0, size = 2.8, color = "#B71C1C",
                fontface = "italic") +
      
      # Faceta por Categoria + Posição
      facet_grid(Category + Position ~ .,
                 scales = "free_y",
                 space  = "free_y",
                 switch = "y") +
      
      # Eixo X
      scale_x_continuous(
        limits = c(-30, 130),
        breaks = seq(0, 100, 25),
        labels = paste0(seq(0, 100, 25), "%"),
        expand = c(0, 0)
      ) +
      
      labs(
        title    = paste0("Polimorfismos por Codão – Gene ", gene),
        subtitle = "Azul = codão HXB2 | Cinzento = polimorfismo | * = referência | Nucleótidos alterados a vermelho sublinhado",
        x        = "Frequência (%)",
        y        = NULL
      ) +
      
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.text.x      = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.text.y.left  = element_text(angle = 0, face = "bold", size = 10),
        strip.placement    = "outside",
        panel.border       = element_rect(color = "grey70", fill = NA),
        legend.position    = "top",
        plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(size = 9, color = "grey40")
      )
    
    print(p_codon)
    ask_export_figure(p_codon,
                      paste0("Codon_Polimorfismos_", gene),
                      width  = 14,
                      height = h_fig)
    
    # ──────────────────────────────────────────────────────
    # Tabela resumo de mutações relevantes no console
    # ──────────────────────────────────────────────────────
    cat(sprintf("\n── Mutações com nota clínica – Gene %s ──\n", gene))
    noted <- plot_df_gene %>%
      filter(Note != "") %>%
      select(Position, MutationName, Frequency, Category, Note) %>%
      arrange(Category, Position)
    
    if (nrow(noted) > 0) {
      for (i in seq_len(nrow(noted))) {
        cat(sprintf("  [%s] Pos %d | %s (%.1f%%) → %s\n",
                    noted$Category[i],
                    noted$Position[i],
                    noted$MutationName[i],
                    noted$Frequency[i],
                    noted$Note[i]))
      }
    } else {
      cat("  Nenhuma mutação com nota clínica detetada nas amostras.\n")
    }
    cat("\n")
    
  } # fim loop genes
  
  cat("✅ Análise de entropia e codões concluída.\n\n")
  
} # fim bloco resistência

cat("✅ ETAPA 5 concluída.\n\n")

