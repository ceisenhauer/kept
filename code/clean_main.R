# --------------------------------------------------------------------------------------------------
# cleaning of raw data
#
# author : cat eisenhauer
# date : january 2026
# --------------------------------------------------------------------------------------------------

library(tidyverse)


# import -------------------------------------------------------------------------------------------
meta <- rio::import(here::here('data', 'ref', 'meta_clean.rds'))

df_raw <- data.frame()
for (i in 1:2) {
  tmp <- rio::import(here::here('data', 'raw', paste0('48k_parte', i, '.csv')), 
                    format = 'tsv',
                    header = FALSE,
                    blank.lines.skip = TRUE,
                    comment.char = '',
                    fill = TRUE) |>
    filter(!grepl('^strain', V1))

  colnames(tmp) <- scan(here::here('data', 'raw', paste0('48k_parte', i, '.csv')),
                       what = character(), 
                       nlines = 1, 
                       sep = "\t", 
                       quiet = TRUE)

  df_raw <- bind_rows(df_raw, tmp)
}

df_raw <- left_join(df_raw, meta)



# wrangle ------------------------------------------------------------------------------------------

df <- df_raw |>
  mutate(
    total_size = as.numeric(total_size),
    n_ambiguous = case_when(
      ambiguous_bases == 'no' ~ 0,
      .default = as.numeric(str_extract(ambiguous_bases, "\\d+"))
     )
  ) |>

  # filter by exclusion criteria
  filter(
    n_ambiguous < 40000,            # 2 removed
    !is.na(year),                   # 5728 removed
    !is.na(country),                # 524 removed
  ) |>
  filter(
    year >= 2000 & year != 2023     # 518 removed
  ) |>

  # cleanup genes
  mutate(

    # remove problem genes (%, ?, $) and remove [^, *] characters and harmonize NA style
    across(ends_with('_acquired'), remove_problem_genes),
    Bla_chr = remove_problem_genes(Bla_chr)

    # clean tetR
    Tet_acquired = exclude_genes(Tet_acquired, 'tetR'),

    # create gene clusters
    Tgc_acquired = create_gene_cluster(Tgc_acquired,
                                       cluster_genes = c('tmexC1', 'tmexD1', 'toprJ1'),
                                       cluster_name = 'tmexCD1-toprJ1'),
    MLS_acquired = create_gene_cluster(MLS_acquired,
                                       cluster_genes = c('Mrx', 'mphA'),
                                       cluster_name = 'mphA-mrx')
  )


