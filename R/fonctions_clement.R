
#' Cleans the initial dataset by removing rows that are not in the correct format and removes duplicated lines
#'
#' @param input_df The tabulation table (data.table). Each row
#' corresponds to the number of statistical unitsin a cross defined by a modality of the z1 nomenclature and a modality of the z2 nomenclature
#'
#' @return the cleaned tabulation table 
#'
#' @examples 
#' cleaned_dt<- clean_init_dt(input_dt)

clean_init_dt <- function(input_dt){
  #input_dt <- dt
  dt <- copy(input_dt)
  dt <- dt[,.(nb_obs = sum(nb_obs)),by = .(z1,z2)] # agregates count on duplicated lines !
  
  ## Suppress lines int the wrong format
  dt <- dt[z2 != "FR_unallocated"]
  dt <- dt[nb_obs != 0]
  
  dt
}

#' Return elements of z1 with associated cobnnected components id from a given z1 x z2 crossing data.table
#'
#' @param input_df (data.table) containing the z1xz2 crossing data with nb_obs
#'
#' @return a data table where associating each z1 to the id of its connected component 
#'
#' @examples 
#' link_table<- return_connected_components(link_table)

return_connected_components<-function(input_dt){
  
  #input_dt <- toy_example_3
  #input_dt <- toy_example_4
  #input_dt <- data.table(data_rp)
  # Construction du graph
  # https://igraph.org/r/#docs  joli !
  dt <- copy(input_dt)
  
  l <- prepare_data(dt)
  m_crois <- build_m_crois(l$intersecting_z2)
  
  adjacency_matrix <- m_crois%*%t(m_crois)
  adjacency_matrix <- as(adjacency_matrix,"lsparseMatrix")
  
  graph <- igraph::graph_from_adjacency_matrix(adjacency_matrix)
  clust <- igraph::clusters(graph)
  nodes <- data.table(from = names(clust$membership), id_comp = clust$membership)
  
  # je regarde les clusters
  colnames(nodes)<- c("z1","id_comp")
  
  nodes
}

#' Build the cross matrix 
#'
#' @param dt (data.table) containing the z1 x z2 crossing with out duplicated lines
#'
#' @return returns the crossover matrix (SparseMatrix), whose rows represent the elements of z1 and the columns the elements of z2 (possibly aggregated if they contribute to the same connection) and whose values are equal to the number of statistical units in the considered z1 x z2 crossover
#'
#' @examples 
#' m_crois<- build_mcrois(link_table)

build_m_crois <- function(input_dt){
  # dt <- ltable
  dt <- copy(input_dt)
  dt[, ":="(z1 = factor(z1), z2 = factor(z2))]
  
  m_crois <- Matrix::sparseMatrix(
    i=as.numeric(dt$z1),
    j=as.numeric(dt$z2),
    x=dt$nb_obs,
    dimnames=list(levels(dt$z1),levels(dt$z2))
  )
  m_crois
}




#' Prepare data , cleaning input_df with z1xzé crossin and extract lines containing z2 whichh are not fully contained by one element of z1  
#'
#' @param input_df (data.table) containing the z2 x z2 crossing with counts, tabulation table
#' 
#' @return data.table containing only the z2 intersecting 2 z1

prepare_data <- function(input_dt){
  # input_df <- toy_example_3
  dt <- copy(input_dt)
  dt <- clean_init_dt(dt)
  
  ## z2 intersecting only one element of z1
  dt_z2 <- dt[,.(nb_z1=length(z1)),by=.(z2)]
  dt_z2_mono_z1 <- dt[z2 %in% dt_z2[nb_z1==1]$z2]
  dt_z2_multi_z1 <- dt[z2 %in% dt_z2[nb_z1>1]$z2]
  
  out <- list(intersecting_z2 = dt_z2_multi_z1,fully_included_z2 = dt_z2_mono_z1)
  
  return(out)
}

#' Find at-risk-of-differenciation areas.
#' Build the crossover matrix defined in build_m_crois, and operates the graph reduction functions on it and then returns the z1-zones (union of elements of z1) at risk 
#'
#' @param input_df (data.table) containing the z2 x z2 crossing with counts, tabulation table
#' @param max_agregate_size Integer indicating the maximal size of agregates
#' which are tested exhaustively. If that number is too large (greater than 30), the
#' computations may not end because of the combinations number that can become very large.
#' Also the RAM can be overloaded.
#' @param save_intermediate_data_file Character indicating the suffix of the name of the saved intermediate results.
#' If is null, results are not writing on the hardware. The path root is taken from the working directory (getwd()).
#' if not null a rds file is saved containing the initial crossing matrix, the simplificated matrix after merging and the list of matrix after the splitting operation
#' @param simplify Boolean. If TRUE then the graph simplification (merging + splitting)
#' occures. Otherwise the exhaustive search is directly applied on the original graph.
#' @param verbose Boolean. If TRUE (default), the different steps of the process are 
#' notified and progress bars provide an estimation of time left.
#' @param threshold Strictly positive integer indicating the confidentiality
#' threshold. Observations are considered at risk if one can deduce information
#' on a agregate of n observations where n < threshold.
#' 
#' @return the list of  at-risk zones defined by unions of z1 elements  

find_pbm_diff_tab <- function(
    input_dt,
    max_agregate_size = 15,
    save_intermediate_data_file = NULL,
    simplify = TRUE,
    verbose = TRUE,
    threshold = 11
){
  # test sur grosses composante si ça bugg pas splitter le taf par composantes
  
  # input_df <- readRDS("data/data_rp.rds")
  # df <- setDT(input_df)
  # link_table <- build_link_table(df)
  # ltable <- unique(long_table(link_table))
  # ltable[,.(n_com = length(unique(z1))),by=.(id_comp)][rev(order(n_com)),]
  # list_z1 <- unique(ltable[id_comp == 88]$z1)
  # data_rp <- readRDS("data/data_rp.rds")
  # input_df <- setDT(data_rp)
  # input_df <- input_df[input_df$z1 %in% list_z1]
  
  # input_df <- toy_example_4
  # threshold = 11; max_agregate_size = 15;save_file = NULL; simplify = TRUE; verbose = TRUE
  dt <- copy(input_dt)
  intersecting_z2 <- prepare_data(dt)$intersecting_z2
  
  intersecting_z2[ , z2_b := paste0(sort(z1),collapse="-"),by =.(z2)][ ,.(nb_obs=sum(nb_obs)),by=.(z1,z2_b)]
  intersecting_z2[, z2:= z2_b]
  
  m_crois <- build_m_crois(intersecting_z2)
  
  if(simplify){ #one can choose to skip these steps of graph reduction if desired
    
    if(verbose) message("< --- Merging method 1 --- > ")
    m_crois_1 <- agregate(m_crois, threshold, methode = "m1", verbose = verbose)
    
    if(verbose) message("< --- Merging methods 1 and 2 --- >")
    m_crois_2 <- agregate(m_crois_1, threshold, methode = "both", verbose = verbose)
    
    if(sum(dim(m_crois_2)==0)>0) {
      message("No differentiation problems detected !")
      return(NULL)
    }
    
    if(verbose) message("< --- Splitting the graph --- >") 
    l_decomp <- decompose_m_crois(m_crois_2, max_agregate_size)
  
  }else{
    l_decomp <- comp_connexe_list(m_crois)
  }
  
  if(verbose) message("< --- Exhaustive search of differentiation problems --- >")
  l_ag <- search_diff_agregate(l_decomp, threshold, max_agregate_size) 
  
  # Contains the list of z1-zones at risk
  l_ag <- desagregate_list(l_ag) 
  
  if(!is.null(save_intermediate_data_file)) {
    # save_1 save_intermediate_data_file = "res_grosse_composante"
    saveRDS(
      list(
        m_crois = m_crois,
        m_crois_agreg_1 = m_crois_1,
        m_crois_agreg_2 = m_crois_2,
        list_splitted_m_crois = l_decomp,
        list_area_at_risk = l_ag
      ),
      paste0(save_intermediate_data_file,".RDS")
    )
  }
  
  return(l_ag)
}



#' Return information on at-risk-of-differenciation z1 zones
#' From a given list of elements of z1, outputs all the information allowing to evaluate if a differentiation problem exists on this zone (no information if no differneciation issue)
#'
#' @param list_z1  vector containing the elements of z1 constituting the area to be evaluated
#' @param link_table (data.table) containing the triplets z1-z2-z1 corresponding to a connection between 2 elements of z1 throug one element of z2 with metadata
#' @param threshold (data.table) the frequency rule threshold  below which a zone is considered at risk
#'
#' @return A data.table with the following columns
#' \enumerate{
#' \item $checked_area the area checked (z1 elements concatenated with "-")
#' \item $z1 : the z1 elements in the z1xz2 intersections at risk
#' \item $z2 : the z2 elements in the z1xz2 intersections at risk
#' \item $nb_obs : the number of statistical units in the intersection z1xZ2
#' \item $type_diff : a boolean indicating whether the zone defined by list_z1 is subject to an internal differentiation risk
#' }
#'
#' @examples 
#' ltable<- return_diff_info(c("A","B","C"),link_table,threshold)

return_diff_info <- function(list_z1,input_dt, threshold){
  # l_ag <- find_pbm_diff_tab(toy_example_3,15,threshold = 11,verbose = FALSE)
  # list_z1 <- l_ag[[1]]
  dt <- copy(input_dt)
  l <- prepare_data(dt)
  dt <- l$intersecting_z2
  
  # all of the z2 elements partially or fully included in the area defined by list_z1, (and not fully included in one element of z1, only crossing z1, z2 elements are interesting here)
  z2_target <- dt[z1 %in% list_z1,]$z2 
  
  # z2 elements crossing the list_z1 area and other elements of z1 than those in list_z1 (frontiers of the area)
  # z2 at the frontier of the area
  at_risk_crossing <- unique(dt[z2 %in% z2_target &  ! z1 %in% list_z1,"z2"])
  
  # build the internal_diff_table  (part of the intersection of at_risk_crossing inside the list_z1 area) 
  # and the external_diff_table  (part of the intersection of at_risk_crossing outside the list_z1 area) 
  diff_table <- dt[z2 %in% at_risk_crossing$z2 ]
  external_diff_table <-diff_table[ !z1 %in% list_z1]
  internal_diff_table <-diff_table[ z1 %in% list_z1]
  
  external_diff_issue <- internal_diff_issue <- FALSE
  
  # check for differenciation issue (if under the threshold)
  if (sum(internal_diff_table$nb_obs) < threshold) internal_diff_issue <- TRUE
  if (sum(external_diff_table$nb_obs) < threshold) external_diff_issue <- TRUE
  
  out <- data.table()
  
  if(internal_diff_issue){
    internal_diff_table$checked_area <- paste0(list_z1,collapse = "-")
    internal_diff_table$type_diff <- "internal"
    out <-rbind(out,internal_diff_table)
  }
  
  if(external_diff_issue){
    external_diff_table$checked_area <- paste0(list_z1,collapse = "-")
    external_diff_table$type_diff <- "external"
    out <-rbind(out,external_diff_table)
  }
  
  out
}



#' output a leaflet interactive map from a \emph{situation}, subset of lines from the link table
#'
#' @param situation_table  a subset of crossing z1 x z2 with nb_obs at the intersection
#' @param geom_z1 the sf data.frame containing geometry of z1 elements referenced in the situation table 
#' @param geom_z2 the sf data.frame containing geometry of z2 elements referenced in the situation table 
#' @param list_z1_to_color list of z1 elements whichh polygon will be colored in the output interactive map
#' @param threshold Strictly positive integer indicating the confidentiality
#' threshold. z1 x z2 intersections which number of statistical units is under the threshold are colored in red
#' @param save_name boolean, if not nul the map is saved in the diffman_results with the given name
#' 
#' @return A list with the following elements
#' 
#' @examples 
#' ltable<- return_diff_info(c("A","B","C"),link_table,threshold)

draw_situation <- function(situation_table,geom_z1,geom_z2,list_z1_to_color = NULL ,threshold = 11,save_name = NULL){
  
  # input_df <- readRDS("data/data_rp.rds")
  # df <- setDT(input_df)
  # link_table <- build_link_table(df)
  # ltable <- unique(long_table(link_table))
  # ltable[,.(n_com = length(unique(z1))),by=.(id_comp)][rev(order(n_com)),]
  # list_z1 <- unique(ltable[id_comp == 17]$z1)
  # data_rp <- readRDS("data/data_rp.rds")
  # input_df <- setDT(data_rp)
  # situation_table <- input_df[input_df$z1 %in% list_z1]

  # install.packages("btb")
  # get_centroid_carreau <- function(data, var_carreau){
  #   
  #   taille_carreau <- readr::parse_number(stringr::str_extract(data[[var_carreau]][1], "RES[0-9]*m"))
  #   centroides <- data %>% 
  #     select(carreau = all_of(var_carreau)) %>%
  #     mutate(
  #       ll_coord_x = readr::parse_number(stringr::str_extract(carreau, "E[0-9]*$")),
  #       ll_coord_y = readr::parse_number(stringr::str_extract(carreau, "N[0-9]*")),
  #       centroid_x = ll_coord_x + taille_carreau/2,
  #       centroid_y = ll_coord_y + taille_carreau/2,
  #       crs = readr::parse_number(stringr::str_extract(carreau, "CRS[0-9]*"))
  #     )
  #   return(list(df = centroides, epsg=centroides$crs[1], taille=taille_carreau))
  # }
  # 
  # carreaux_to_polygon <- function(data, var_carreau){
  #   
  #   require(btb)
  #   
  #   centroides_l <- get_centroid_carreau(data, var_carreau) 
  #   
  #   centroides_l$df %>% 
  #     select(carreau, x=centroid_x, y=centroid_y) %>% 
  #     btb::dfToGrid(sEPSG = centroides_l$epsg, iCellSize = centroides_l$taille)
  #   
  # }
  # 
  # list_z1_to_color <- NULL
  
#mc cp s3/cguillo/commune_franceentiere_2021.dbf data/commune_franceentiere_2021.dbf
#mc cp s3/cguillo/commune_franceentiere_2021.fix data/commune_franceentiere_2021.fix
#mc cp s3/cguillo/commune_franceentiere_2021.prj data/commune_franceentiere_2021.prj
#mc cp s3/cguillo/commune_franceentiere_2021.shp data/commune_franceentiere_2021.shp
#mc cp s3/cguillo/commune_franceentiere_2021.shx data/commune_franceentiere_2021.shx
#mc cp s3/cguillo/data_rp.rds data/data_rp.rds
  
  # communes <-st_read("data/commune_franceentiere_2021.shp")
  # situation_table <- prepare_data(situation_table)$intersecting_z2
  # liste_carreau <- situation_table$z2
  # 
  # polygone_carreau <-
  #   carreaux_to_polygon(data.frame(carreau = liste_carreau), var_carreau = "carreau") %>%
  #   st_transform(crs = 4326)
  # geom_z1 <- communes %>%
  #   filter(code %in% unique(c(situation_table$z1))) %>%
  #   select(code) %>%
  #   rename(z1 = code)
  # geom_z2 <- polygone_carreau %>%
  #   rename(z2 = carreau) %>%
  #   select(-x,-y)

  
  l <- prepare_data(situation_table)
  situation_table <-l$intersecting_z2
  
  z2_to_nb_obs <- situation_table[,.(nb_obs_z2=sum(nb_obs)),by = z2]
  
  # build the z1 x z2 intersection geometry
  st_agr(geom_z1) = "constant"
  st_agr(geom_z2) = "constant"
  
  inter_carreau_commune <-st_intersection(
    geom_z1 %>% select(z1),
    geom_z2 %>% select(z2)
  )
  
  inter_carreau_commune <- merge(situation_table,inter_carreau_commune,by = c("z1","z2"))
  
  geom_z2 <- merge(geom_z2,z2_to_nb_obs,by ="z2",nomatch = 0)
  
  z1_fillColor <- with(geom_z1,ifelse(z1 %in% list_z1_to_color,"orange","#3FC8FC"))
  
  highlightOptions_defaut <- highlightOptions(
    stroke = TRUE,
    weight = 6,
    color = "black",
    fillColor = "black",
    bringToFront = TRUE
  )
  
  m <- 
    leaflet() %>% 
    addProviderTiles("GeoportailFrance.orthos") %>%  
    addPolygons(
      data = geom_z1,
      color = "#3FC8FC",
      fillColor = z1_fillColor,
      weight = 2,
      fillOpacity = 0.25,
      opacity = 1,
      label = geom_z1$z1
    ) %>% 
    addPolygons(
      data = geom_z2,
      color = "red",
      label = with(geom_z2, 
                   sprintf(
                     "<b> id z2 : </b> %s  <br/> <b> Number of observations : </b>  %s", 
                     z2, round(nb_obs_z2,1)
                   ) %>% lapply(htmltools::HTML)
      ),
      weight = 2,
      fillOpacity = 0,
      opacity = 1,
      group ="z2 on two sides of one z1 area"
    ) %>% 
    addPolygons(
      data = inter_carreau_commune$geometry,
      color = ifelse(inter_carreau_commune$nb_obs < threshold,"red","#6E3DFF"),
      weight = 2,
      fillOpacity = 0.5,
      group = "intersections",
      highlightOptions = highlightOptions_defaut,
      label  =  with(inter_carreau_commune, 
                     sprintf(
                       "<b> id z1 : </b> %s  <br/> <b> id z2 : </b>  %s <br/> <b> Number of observations : </b>  %s", 
                       z1, z2, round(nb_obs,1)
                     ) %>% lapply(htmltools::HTML)
      )
    ) %>% 
    addScaleBar(position="bottomright") %>% 
    hideGroup(c("z2 on two sides of one z1 area","intersections")) %>% 
    addLayersControl(
      overlayGroups = c("z2 on two sides of one z1 area","intersections"),
      options = layersControlOptions(collapsed = FALSE)
    ) %>% 
    addScaleBar(position="bottomright")
  
  if(!is.null(save_name)) htmlwidgets::saveWidget(m, file=paste0(save_name,".html"),selfcontained = TRUE)
  
  m
}

#' create the equivalent individual data.frame from a tabulate z1 x z2 data.frame with integer (round in any case)
#'
#' @param input_df The tabulation table (data.frame or data.table). Each row
#' corresponds to the number of statistical units in a cross defined by a modality of the z1 nomenclature and a modality of the z2 nomenclature
#'
#' @return the fictive individual table from the tabulation table 
#'
#' @examples 
#' tab_table <- data.table(
#'   z1 =c("A","A","B","D"),
#'   z2 = c("a","b","b","a"),
#'   nb_obs = c(2,5,6,0)
#'   )
#'
#' create_fictive_ind_table(tab_table)

create_fictive_ind_table <- function(tab_table){
  
  l <- split(tab_table,with(tab_table,paste0(z1,"_",z2)))
  tmp <- lapply(
    1:nrow(tab_table),
    function(i){
      # i <- 1
      crossing <- tab_table[i,]
      if (crossing$nb_obs == 0) return(NULL)
      
      with(
        crossing,
        data.table(
          z1 = rep(z1,nb_obs),
          z2 = rep(z2,nb_obs),
          stringsAsFactors = FALSE
        )
      )
    }
  )
  
  out <- do.call(rbind,tmp)
  
  out$id <- seq(1,nrow(out))
  
  out[,c("id","z1","z2")]
}




### Fonction pour comparer diffman pas tab et diffman tab, à effacer par la suite
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
