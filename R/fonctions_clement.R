
#' Cleans the initial dataset by removing rows that are not in the correct format and removes duplicated lines
#'
#' @param input_df The tabulation table (data.frame or data.table). Each row
#' corresponds to the number of statistical unitsin a cross defined by a modality of the z1 nomenclature and a modality of the z2 nomenclature
#'
#' @return the cleaned tabulation table 
#'
#' @examples 
#' cleaned_df<- clean_init_df(input_df)

clean_init_df <- function(input_df){
  #input_df <- df
  df <- copy(input_df)
  df <- df[,.(nb_obs = sum(nb_obs)),by = .(z1,z2)] # agregates count on duplicated lines !
  
  ## Suppress lines int the wrong format
  df <- df[z2 != "FR_unallocated"]
  df <- df[nb_obs != 0]
  
  df
}


#' Build the link table from the initial tabulation table (data.frame or data.table) where each row
#' correspond to the number of statistical units in a cross defined by a modality of the z1 nomenclature and a modality of the z2 nomenclature. 
#' One line of the unit table represents a connection, triplets composed of 2 elements of z1 and one element of z2 connecting. 
#' The identifier of connected component in the graph of the elements of z2 is added as well as the number of observations contained 
#' in the intersection between the zone of z2 connecting and the 2 elements of z1 connected by z2
#'
#' @param input_df The tabulation table (data.frame or data.table). Each row
#' corresponds to the number of statistical unitsin a cross defined by a modality of the z1 nomenclature and a modality of the z2 nomenclature
#'
#' @return the link_table containing the triplet and metadatas (data.table)
#'
#' @examples 
#' link_table<- build_link_table(input_df)


build_link_table <- function(input_df){
  # donnees_rp <- readRDS("../data_rp.rds")
  # df <- setDT(donnees_rp)
  # input_df <- df[substr(z1,1,2) == "22"]
  
  df <- copy(input_df)
  df <- clean_init_df(df)
  
  ## z2 intersecting only one element of z1
  df_z2 <- df[,.(nb_z1=length(z1)),by=.(z2)]
  df_z2_mono_z1 <- df[z2 %in% df_z2[nb_z1==1]$z2]
  df_z2_multi_z1 <- df[z2 %in% df_z2[nb_z1>1]$z2]
  
  
  link_table <- merge( 
    df_z2_multi_z1,
    df_z2_multi_z1,
    by=c('z2'),
    allow.cartesian = TRUE,
    nomatch = 0,
    suffix = c("_from","_to")
  ) 
  
  colnames(link_table)[colnames(link_table)=="z1_from"] <- "from"
  colnames(link_table)[colnames(link_table)=="z1_to"] <- "to"
  
  link_table <- link_table[from != to]
  
  link_table$from_to_z2 <- apply(link_table[,c("from","to","z2")],1,function(x) paste0(sort(x),collapse= "-"))
  link_table <- link_table[!duplicated(from_to_z2)][,c("from","to","z2","nb_obs_from","nb_obs_to")]
  
  link_table[,n_z2 := length(z2) ,by =.(from,to)]
  
  return_connected_components(link_table)
}




#' Return elements of z1 graph nodes with the associated connected components  
#'
#' @param link_table (data.table) containing the triplets z1-z2-z1 corresponding to a connection between 2 elements of z1 throug one element of z2
#'
#' @return the node_table, each a line corresponds to an element of z1 with the identifier of the connected component to which it belongs in the graph of elements of z1
#'
#' @examples 
#' link_table<- return_connected_components(link_table)

return_connected_components<-function(link_table){
  # Construction du graph
  # https://igraph.org/r/#docs  joli !
  
  graph <- graph_from_data_frame(link_table,directed = FALSE)
  clust <- clusters(graph)
  nodes <- data.table(from = names(clust$membership), id_comp = clust$membership)
  
  # je regarde les clusters
  out <- nodes[link_table, on = "from"]
  
  out
}


#' Transform the link table containing triplets z1-z2-z1 with metadata about connection into a long table where each line correspond to  a part of the connections z1-z2 and the number of statistical unitis inside the intersection of the corresponding z1 and z2    
#'
#' @param link_table (data.table) containing the triplets z1-z2-z1 corresponding to a connection between 2 elements of z1 throug one element of z2 
#'
#' @return the long table with the metadata, identifiying the actual connected components in the z1 graph
#'
#' @examples 
#' ltable<- long_table(link_table)

long_table <- function(link_table){
  
  tab_from <- link_table[,c("from","z2","nb_obs_from","id_comp")] 
  colnames(tab_from) <- c("z1","z2","nb_obs","id_comp")
  
  tab_to <-  link_table[,c("to","z2","nb_obs_to","id_comp")] 
  colnames(tab_to) <- c("z1","z2","nb_obs","id_comp")
  
  out <- rbind(tab_from,tab_to)
  
  return(out)
  
}



#' Build the cross matrix 
#'
#' @param link_table (data.table) containing the triplets z1-z2-z1 corresponding to a connection between 2 elements of z1 throug one element of z2
#'
#' @return returns the crossover matrix (SparseMatrix), whose rows represent the elements of z1 and the columns the elements of z2 (possibly aggregated if they contribute to the same connection) and whose values are equal to the number of statistical units in the considered z1 x z2 crossover
#'
#' @examples 
#' m_crois<- build_mcrois(link_table)

build_m_crois <- function(link_table){
  
  df <- link_table[ ,
                    .(
                      nb_obs_from = sum(nb_obs_from),
                      nb_obs_to = sum(nb_obs_to),
                      id_comp = unique(id_comp),
                      z2 = paste0(z2,collapse = "-"),
                      n_z2 = unique(n_z2)
                    ),
                    by = .(from,to)
  ]
  
  ltable <-long_table(df)
  
  ltable$z1  <- as.factor(ltable$z1)
  ltable$z2  <- as.factor(ltable$z2)
  
  m_crois <- Matrix::sparseMatrix(
    i=as.numeric(ltable$z1),
    j=as.numeric(ltable$z2),
    x=ltable$nb_obs,
    dimnames=list(levels(ltable$z1),levels(ltable$z2))
  )
  
  m_crois 
}


#' Find at-risk-of-differenciation areas.
#' Build the crossover matrix defined in build_m_crossover, and operates the graph reduction functions on it and then returns the z1-zones (union of elements of z1) at risk 
#'
#' @param link_table (data.table) containing the triplets z1-z2-z1 corresponding to a connection between 2 elements of z1 throug one element of z2 with metadata
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
    link_table,
    max_agregate_size = 15,
    save_intermediate_data_file = NULL,
    simplify = TRUE,
    verbose = TRUE,
    threshold = 11
){
  
  # threshold = 7; max_agregate_size = 15;save_file = NULL; simplify = TRUE; verbose = TRUE
  
  m_crois <- build_m_crois(link_table)
  
  if(simplify){ #one can choose to skip these steps of graph reduction if desired
    
    if(verbose) message("< --- Merging method 1 --- > ")
    m_crois_1 <- diffman:::agregate(m_crois, threshold, methode = "m1", verbose = verbose)
    
    if(verbose) message("< --- Merging methods 1 and 2 --- >")
    m_crois_2 <- diffman:::agregate(m_crois_1, threshold, methode = "both", verbose = verbose)
    
    if(sum(dim(m_crois)==0)>0) {
      message("No differentiation problems detected !")
      return(NULL)
    }
    
    if(verbose) message("< --- Splitting the graph --- >") 
    l_decomp <- diffman:::decompose_m_crois(m_crois_2, max_agregate_size)
    
    
    if(!is.null(save_intermediate_data_file)) {
      saveRDS(
        list(
          m_crois = m_crois,
          m_crois_agreg_1 = m_crois_1,
          m_crois_agreg_2 = m_crois_2,
          list_splitted_m_crois = l_decomp
        ),
        paste0("diffman_results/",save_file,".RDS")
      )
    }
    
  }else{
    l_decomp <- diffman:::comp_connexe_list(m_crois)
  }

  if(verbose) message("< --- Exhaustive search of differentiation problems --- >")
  l_ag <- diffman:::search_diff_agregate(l_decomp, threshold, max_agregate_size) 
  
  # Contains the list of z1-zones at risk
  l_ag <- diffman:::desagregate_list(l_ag) 
  
  return(l_ag)
}



#' Return information on at-risk-of-differnciation z1 zones
#' From a given list of elements of z1, outputs all the information allowing to evaluate if a differentiation problem exists on this zone
#'
#' @param list_z1  vector containing the elements of z1 constituting the area to be evaluated
#' @param link_table (data.table) containing the triplets z1-z2-z1 corresponding to a connection between 2 elements of z1 throug one element of z2 with metadata
#' @param threshold (data.table) the frequency rule threshold  below which a zone is considered at risk
#'
#' @return A list with the following elements
#' \enumerate{
#' \item $checked_area the area checked (z1 elements concatenated with "-")
#' \item $internal_diff_table : the data.table containing intersection z1 x z2 at risk of internal differenciation
#' \item $external_diff_table : the data.table containing intersection z1 x z2 at risk of external differenciation
#' \item $internal_diff_issue : a boolean indicating whether the zone defined by list_z1 is subject to an internal differentiation risk
#' \item $external_diff_issue : a boolean indicating whether the zone defined by list_z1 is subject to an external differentiation risk
#' }
#'
#' @examples 
#' ltable<- return_diff_info(c("A","B","C"),link_table,threshold)

return_diff_info <- function(list_z1,link_table, threshold){
  
  dataf <- long_table(link_table) 
  
  # all of the z2 elements partially or fully included in the area defined by list_z1, (and not fully included in one element of z1, only crossing z1, z2 elements are interesting here)
  z2_target <- dataf[z1 %in% list_z1,]$z2 
  
  # z2 elements crossing the list_z1 area and other elements of z1
  at_risk_crossing <- unique(dataf[z2 %in% z2_target &  ! z1 %in% list_z1,"z2"])
  
  # build the internal_diff_table  (part of the intersection of at_risk_crossing inside the list_z1 area) 
  # and the external_diff_table  (part of the intersection of at_risk_crossing outside the list_z1 area) 
  diff_table <- dataf[z2 %in% at_risk_crossing$z2 ]
  external_diff_table <-diff_table[ !z1 %in% list_z1]
  internal_diff_table <-diff_table[ z1 %in% list_z1]
  
  external_diff_issue <- internal_diff_issue <- FALSE
  
  # check for differenciation issue (if under the threshold)
  if (sum(internal_diff_table$nb_obs) < threshold) internal_diff_issue <- TRUE
  if (sum(external_diff_table$nb_obs) < threshold) external_diff_issue <- TRUE
  
  out <- list(
    checked_area = paste0(list_z1,collapse = "-"),
    internal_diff_table = internal_diff_table,
    external_diff_table = external_diff_table,
    internal_diff_issue = internal_diff_issue, 
    external_diff_issue = external_diff_issue
  )
  
  return(out)
  
}







#' output a leaflet interactive map from a \emph{situation}, subset of lines from the link table
#'
#' @param situation_table  a subset of the link_table build with the build_link_table_function
#' @param z2_to_nb_obs data.table containing the number of statistical units inside each elements of z2
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

draw_situation <- function(situation_table,z2_to_nb_obs,geom_z1,geom_z2,list_z1_to_color = NULL ,threshold = 11,save_name = NULL){
  
  # z2_to_nb_obs <- long_table(link_table)[ ,.(nb_obs_z2 = sum(nb_obs)) ,by = .(z2)]
  
  # situation_table <- link_table[id_comp == 6] ;unique(c(situation_table$from,situation_table$to))
  # list_z1_to_color <- c("01014")
  # liste_carreau <- situation_table$z2
  # polygone_carreau <-
  #   carreaux_to_polygon(data.frame(carreau = liste_carreau), var_carreau = "carreau") %>%
  #   st_transform(crs = 4326)
  # 
  # geom_z1 <- communes %>%
  #   filter(code %in% unique(c(situation_table$from,situation_table$to))) %>%
  #   select(code) %>%
  #   rename(z1 = code)
  # 
  # geom_z2 <- polygone_carreau %>%
  #   rename(z2 = carreau) %>%
  #   select(-x,-y)
  
  
  # build the z1 x z2 intersection geometry
  st_agr(geom_z1) = "constant"
  st_agr(geom_z2) = "constant"
  
  inter_carreau_commune <-st_intersection(
    geom_z1 %>% select(z1),
    geom_z2 %>% select(z2)
  )
  
  tab_croisement<- long_table(situation_table)
  
  inter_carreau_commune <- merge(tab_croisement,inter_carreau_commune,by = c("z1","z2"))
  
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
