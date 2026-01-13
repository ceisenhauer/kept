#remove_problem_genes <- function(genes) {
  #out <- genes |>
    #str_split(';') |>
    #map_chr(~{
      #filtered <- .x[!str_detect(.x, "[%?$]")]
      #cleaned <- str_remove_all(filtered, "[*^]")
      #paste(cleaned, collapse = ";")
    #}) |>
    #na_if("") |>
    #na_if("-")

  #return(out)
#}

remove_problem_genes <- function(genes) {
  sapply(strsplit(genes, ';'), function(items) {
    filtered <- items[!grepl("[%?$]", items)]
    cleaned <- gsub("[*^]", "", filtered)
    cleaned <- cleaned[nzchar(cleaned)]
    out <- paste(cleaned, collapse = ';')

    if (out == '' || out == '-') {
      return(NA_character_)
    }
    return(out)
  })
}

exclude_genes <- function(genes, to_exclude) {
  sapply(strsplit(genes, ';'), function(items) {
    filtered <- items[!items %in% to_exclude]
    result <- paste(filtered, collapse = ';')

    if (result == "") {
      return(NA_character_)
    }

    return(result)
  })
}

create_gene_cluster <- function(genes, cluster_genes, cluster_name) {
  sapply(strsplit(genes, ';'), function(items) {
    all_present <- all(cluster_genes %in% items)
    filtered <- items[!items %in% cluster_genes]

    if (all_present) {
      filtered <- c(filtered, cluster_name)
    }

    result <- paste(filtered, collapse = ';')

    if (result == '') {
      return(NA_character_)
    }

    return(result)
  })
}
