
test_toy_example <- function(toy_example){
  
  t_ind <- create_fictive_ind_table(toy_example)
  
  res_diff_ind <- find_pbm_diff(
    t_ind = t_ind,
    threshold = 11,
    max_agregate_size = 15
  )
  
  
  if (! is.null(res_diff_ind)) res_diff_ind <- res_diff_ind[!duplicated(paste0(agregat_z1,type_diff)),] else res_diff_ind <- NULL
  
  link_table <- build_link_table(toy_example)
  
  res_diff_tab <- find_pbm_diff_tab(
    link_table = link_table,
    threshold = 11,
    max_agregate_size = 15
  )
  
  res_diff_tab <- do.call(rbind,lapply(res_diff_tab,return_diff_info, link_table = link_table, threshold = 11))
  
  return(list(res_diff_ind = res_diff_ind, res_diff_tab = res_diff_tab))
}


test_toy_example(toy_example_1)
test_toy_example(toy_example_2)


toy_example_3 <- # pb dde diff dans c1 cé inter M2 sur c4 inter M5 et c3 est sous le seuil, pas détecté car pas dans le cadres de la différenciation mais c'est du secret primaire
  data.table(
    z1 = c("M1","M1","M2","M2","M3","M4","M5"),
    z2 =  c("c1","c2","c1","c2","c3","c4","c4"),
    nb_obs = c(2,10,5,3,4,13,5)
  )

usethis::use_data(toy_example_3, overwrite = TRUE)
test_toy_example(toy_example_3)


toy_example_4 <- # c1 à cheval sur M1, M2, M3
  data.table(
    z1 = c("M1","M2","M3","M2","M3"),
    z2 =  c("c1","c1","c1","c2","c2"),
    nb_obs = c(2,13,10,20,30)
  )

usethis::use_data(toy_example_4, overwrite = TRUE)
test_toy_example(toy_example_4)



