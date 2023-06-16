# uninspired helper function to aggregate meta celltypes

aggregated_celltype = function(x) {
  if (grepl('CD8', x)) {
    return ('CD8+ T cell')
  } else if (grepl('CD4', x)) {
    return ('CD4+ T cell')
  } else if (grepl('Treg', x)) {
    return ('CD4+ T cell')
  } else if (x %in% c('Memory B', 'Naive B', 'Plasmablast')) {
    return ('B cell')
  } else if (x %in% c('NK', 'CD56 bright NK')) {
    return ('NK')
  } else if (x %in% c('CD14 Mono', 'CD16 Mono')) {
    return ('Mono')
  } else if (x %in% c('gdT', 'MAIT')) {
    return ('innate T cell')
  } else if (x %in% c('Prog_B 1', 'Prog_B 2')) {
    return ('progenitor B cell')
  } else if (x %in% c('Prog_Mk')) {
    return ('megakaryopoiesis')
  } else if (x %in% c('Prog_RBC')) {
    return ('erythropoiesis')
  }
  return (x)
}
