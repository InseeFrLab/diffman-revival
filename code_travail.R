
test_toy_example <- function(toy_example){
  
  t_ind <- create_fictive_ind_table(toy_example)
  
  res_diff_ind <- find_pbm_diff(
    t_ind = t_ind,
    threshold = 11,
    max_agregate_size = 15
  )
  
  # res_diff_ind <- res_diff_ind[!duplicated(paste0(agregat_z1,type_diff)),]
  
  link_table <- build_link_table(toy_example)
  
  res_diff_tab <- find_pbm_diff_tab(
    link_table = link_table,
    threshold = 8,
    max_agregate_size = 15
  )
  
  res_diff_tab <- do.call(rbind,lapply(res_diff_tab,return_diff_info, link_table = link_table, threshold = 11))
  
  return(list(res_diff_ind = res_diff_ind, res_diff_tab = res_diff_tab))
}


test_toy_example(toy_example_2)
