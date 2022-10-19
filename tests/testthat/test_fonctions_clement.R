#' Compare output from diffman_ind and diffman_tab function
#'
#' @param toy_example The toy_example data.table which containsd counts for each crossing z1 x z2
#'
#' @return the list containin the result of the 2 ssearch_diff_pbm functions
#'
#' @examples 
#' data("toy_example_1")
#' test_toy_example(toy_example_1)


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


data("toy_example_1")

test_that("same output toy example 1", {
  expect_true(
    unique(test_toy_example(toy_example_1)$res_diff_ind$agregat_z1) %in% 
      unique(test_toy_example(toy_example_1)$res_diff_tab$checked_area)
  )
})

