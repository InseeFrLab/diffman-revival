
#' Cleans the initial dataset by removing rows that are not in the correct format and removes duplicated lines
#'
#' @param input_dt The tabulation table (data.table). Each row
#' corresponds to the number of statistical unitsin a cross defined by a modality of the z1 nomenclature and a modality of the z2 nomenclature
#'
#' @return the cleaned tabulation table 
#'
#' @examples
#' input_dt<- toy_example_6
#' cleaned_dt <- diffman:::clean_init_dt(input_dt)
#' @export

clean_init_dt <- function(input_dt){
  #input_dt <- dt
  
  nb_obs <- z1 <- z2 <- NULL
  dt <- copy(input_dt)
  dt <- dt[,.(nb_obs = sum(nb_obs)),by = .(z1,z2)] # agregates count on duplicated lines !
  
  ## Suppress lines int the wrong format
  dt <- dt[z2 != "FR_unallocated"]
  dt <- dt[nb_obs != 0]
  
  dt
}

#' Return elements of z1 with associated cobnnected components id from a given z1 x z2 crossing data.table
#'
#' @param input_dt (data.table) containing the z1xz2 crossing data with nb_obs
#'
#' @return a data table where associating each z1 to the id of its connected component 
#'
#' @examples 
#' input_dt <- toy_example_6
#' link_table <- return_connected_components(input_dt)
#' @export

return_connected_components<-function(input_dt){
  
  # Construction du graph
  # https://igraph.org/r/#docs  joli !
  dt <- copy(input_dt)
  
  l <- prepare_data(dt)
  m_crois <- build_m_crois(l$intersecting_z2)
  
  adjacency_matrix <- m_crois%*%t(m_crois)
  adjacency_matrix <- methods::as(adjacency_matrix,"lsparseMatrix")
  
  graph <- igraph::graph_from_adjacency_matrix(adjacency_matrix)
  clust <- igraph::clusters(graph)
  nodes <- data.table(from = names(clust$membership), id_comp = clust$membership)
  
  # je regarde les clusters
  colnames(nodes)<- c("z1","id_comp")
  
  nodes
}

#' Build the cross matrix 
#'
#' @param input_dt (data.table) containing the z1 x z2 crossing with out duplicated lines
#'
#' @return returns the crossover matrix (SparseMatrix), whose rows represent the elements of z1 and the columns the elements of z2 (possibly aggregated if they contribute to the same connection) and whose values are equal to the number of statistical units in the considered z1 x z2 crossover
#'
#' @examples 
#' input_dt <- toy_example_6
#' m_crois<- diffman:::build_m_crois(input_dt)

build_m_crois <- function(input_dt){
  # dt <- ltable
  z1 <- z2 <- NULL
  
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
#' @param input_dt (data.table) containing the z2 x z2 crossing with counts, tabulation table
#' 
#' @return a list with the intersection on the first element and the llines corresponding to the fully included square in the second element

prepare_data <- function(input_dt){
  # input_df <- toy_example_3
  
  dt <- copy(input_dt)
  dt <- clean_init_dt(dt)
  
  out <- isolate_intersection(dt)
  
  return(out)
}

#' Isolate intersection z2xz1 from z2 lines where z2 is fully included in z1 
#' @param input_dt (data.table) containing the z2 x z2 crossing with counts, tabulation table
#' 
#' @return a list with the intersection on the first element and the llines corresponding to the fully included square in the second element
#' @export

isolate_intersection <- function(input_dt){
  
  z1 <- z2 <- nb_z1 <- NULL
  dt <- copy(input_dt)
  
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
#' @param input_dt (data.table) containing the z2 x z2 crossing with counts, tabulation table
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
#' 
#' @examples 
#' input_dt <- toy_example_6
#' find_pbm_diff_tab(input_dt)
#' 
#' @export
find_pbm_diff_tab <- function(
    input_dt,
    max_agregate_size = 15,
    save_intermediate_data_file = NULL,
    simplify = TRUE,
    verbose = TRUE,
    threshold = 11
){
  
  z2_b<- z1 <- z2 <- nb_obs <- NULL
  
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
#' @param input_dt (data.table) containing the triplets z1-z2-z1 corresponding to a connection between 2 elements of z1 throug one element of z2 with metadata
#' @param threshold (data.table) the frequency rule threshold  below which a zone is considered at risk
#'
#' @return A data.table with the diff information
#'
#' @examples 
#' list_z1 <- c("M1","M2",c("M1","M2"))
#' input_dt <- toy_example_1
#' diffman:::return_diff_info(list_z1,input_dt)
#' 
#' @export

return_diff_info <- function(list_z1,input_dt, threshold = 11){
  # l_ag <- find_pbm_diff_tab(toy_example_1,15,threshold = 11,verbose = FALSE)
  # list_z1 <- l_ag[[1]]de
  # list_z1 = area_z1_at_risk;threshold = 11;dt = input_dt
  
  z1 <- z2 <- NULL
  
  dt <- copy(input_dt)
  
  dt <- clean_init_dt(dt)
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



#' output a leaflet interactive map from a \emph{situation}, corresponding to tthe set of information needed ton interpretn differencing issue s in a zone built frolmz1 elements
#'
#' @param situation_table  a subset of crossing z1 x z2 with nb_obs at the intersection
#' @param geom_z1 the sf data.frame containing geometry of z1 elements referenced in the situation table 
#' @param geom_z2 the sf data.frame containing geometry of z2 elements referenced in the situation table 
#' @param list_z1_to_color list of z1 elements whichh polygon will be colored in the output interactive map
#' @param list_z2_to_color list of z2 elements whichh polygon will be colored in the output interactive map
#' @param threshold Strictly positive integer indicating the confidentiality
#' threshold. z1 x z2 intersections which number of statistical units is under the threshold are colored in red
#' @param save_name boolean, if not nul the map is saved in the diffman_results with the given name
#' 
#' @return A list with the following elements
#' 
#' @export

draw_situation <- function(situation_table,geom_z1,geom_z2,list_z1_to_color = NULL,list_z2_to_color = NULL, threshold = 11,save_name = NULL){
  
  nb_obs <- z2 <- z1 <- NULL
  
  dt <- copy(situation_table)
  dt <- clean_init_dt(dt)
  z2_to_nb_obs <- dt[,.(nb_obs_z2=sum(nb_obs)),by = z2]
  z1_to_nb_obs <- dt[,.(nb_obs_z1=sum(nb_obs)),by = z1]
  
  l <- isolate_intersection(dt)
  intersect_dt <-l$intersecting_z2
  inside_dt <- l$fully_included_z2
  
  geom_z2_inter <- geom_z2[geom_z2$z2 %in% intersect_dt$z2,]
  #geom_z2_inside <-  geom_z2[geom_z2$z2 %in% inside_dt$z2,]
  
  geom_z2_inter <- merge(geom_z2_inter,z2_to_nb_obs,by ="z2",nomatch = 0)
  geom_z2 <- merge(geom_z2,z2_to_nb_obs,by ="z2",nomatch = 0)
  
  # build the z1 x z2 intersection geometry
  sf::st_agr(geom_z1) = "constant"
  sf::st_agr(geom_z2_inter) = "constant"
  
  inter_z2_z1 <-sf::st_intersection(
    geom_z1[,"z1"],
    geom_z2_inter[,"z2"]
  )
  
  inter_z2_z1 <- merge(dt,inter_z2_z1,by = c("z1","z2"))
  z1_fillColor <- with(geom_z1,ifelse(z1 %in% list_z1_to_color,"orange","#3FC8FC"))
  z2_fillColor <- with(geom_z2,ifelse(z2 %in% list_z2_to_color,
                                      ifelse(nb_obs_z2 >= threshold,"#1FF064","#FFD5C7"),
                                      "#04117A"))
  z2_fillOpacity <- with(geom_z2,ifelse(z2 %in% list_z2_to_color,1,0))
  
  geom_z1 <- merge(geom_z1,z1_to_nb_obs,by ="z1",nomatch = 0)
  
  highlightOptions_defaut <- leaflet::highlightOptions(
    stroke = TRUE,
    weight = 6,
    color = "black",
    fillColor = "black",
    bringToFront = TRUE
  )
  
  m <- leaflet::leaflet()
  
  m <- leaflet::addProviderTiles(m,"GeoportailFrance.orthos") 
  m <- leaflet::addPolygons(
    m,
    data = geom_z1,
    color = "#3FC8FC",
    fillColor = z1_fillColor,
    weight = 2,
    fillOpacity = 0.25,
    opacity = 1,
    label = with(geom_z1,
                 lapply( 
                   sprintf(
                     "<b> id z1 : </b> %s  <br/> <b> Number of observations : </b>  %s", 
                     z1, round(nb_obs_z1,1)
                   ),
                   htmltools::HTML
                 )
    )
  )
  m <- leaflet::addPolygons(
    m,
    data = geom_z2_inter,
    color = "red",
    label = with(geom_z2_inter,
                 lapply( 
                   sprintf(
                     "<b> id z2 : </b> %s  <br/> <b> Number of observations : </b>  %s", 
                     z2, round(nb_obs_z2,1)
                   ),
                   htmltools::HTML
                 )
    ),
    weight = 2,
    fillOpacity = 0,
    opacity = 1,
    group ="z2 on two sides of one z1 area"
  )
  m <- 
    leaflet::addPolygons(
      m,
      data = inter_z2_z1$geometry,
      color = ifelse(inter_z2_z1$nb_obs < threshold,"red","#6E3DFF"),
      weight = 2,
      fillOpacity = 0.5,
      group = "intersections",
      highlightOptions = highlightOptions_defaut,
      label  =  with(
        inter_z2_z1, 
        lapply(
          sprintf(
            "<b> id z1 : </b> %s  <br/> <b> id z2 : </b>  %s <br/> <b> Number of observations : </b>  %s", 
            z1, z2, round(nb_obs,3)
          ),
          htmltools::HTML)
      )
    )
  
  
  
  m <-
    leaflet::addPolygons(
      m,
      data = geom_z2,
      color = "#04117A",
      weight = 2,
      fillOpacity = z2_fillOpacity,
      fillColor = z2_fillColor,
      group = "inside z2",
      label = lapply(
        sprintf(
          "<b> id z2 : </b> %s <br/> <b> Number of observations : </b>  %s",
          geom_z2$z2, round(geom_z2$nb_obs_z2,2)
        ),
        htmltools::HTML
      )
    )
  
  group_names <- c("z2 on two sides of one z1 area","intersections","inside z2")  
  
  m <- leaflet::addScaleBar(m,position="bottomright") 
  m <- leaflet::hideGroup(m,group_names) 
  m <- leaflet::addLayersControl(m,
                                 overlayGroups = group_names,
                                 options = leaflet::layersControlOptions(collapsed = FALSE)
  )
  m <- leaflet::addScaleBar(m,position="bottomright")
  
  if(!is.null(save_name)) htmlwidgets::saveWidget(m, file=paste0(save_name,".html"),selfcontained = TRUE)
  
  m
}

#' create the equivalent individual data.frame from a tabulate z1 x z2 data.frame with integer (round in any case)
#'
#' @param tab_table The tabulation table (data.frame or data.table). Each row
#' corresponds to the number of statistical units in a cross defined by a modality of the z1 nomenclature and a modality of the z2 nomenclature
#'
#' @return the fictive individual table from the tabulation table 
#'
#' @examples 
#' tab_table <- toy_example_1
#' create_fictive_ind_table(tab_table)
#' 
#' @export

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



#' return the diff info from a input data table corresponding to one connected component of z1
#' perform the search of diff issue and return the diff information fro the corresponding component 
#' @param input_dt The tabulation table (data.frame or data.table). Each row
#' corresponds to the number of statistical units in a cross defined by a modality of the z1 nomenclature and a modality of the z2 nomenclature
#' z1 has to belong to the same related component (regarding the graph construction)
#' @param max_agregate_size Integer indicating the maximal size of agregates
#' which are tested exhaustively. If that number is too large (greater than 30), the
#' computations may not end because of the combinations number that can become very large.
#' Also the RAM can be overloaded.
#' @param threshold Strictly positive integer indicating the confidentiality
#' threshold. Observations are considered at risk if one can deduce information
#' on a agregate of n observations where n < threshold.
#' @return a data.table containning the diff info for this componet (same format than the return diff info function)
#'
#' @examples 
#' input_dt <- toy_example_6
#' compo <- return_connected_components(input_dt)
#' one_component_risk_extraction(input_dt[z1 %in% compo[id_comp == 1]$z1], 11, 15)
#' 
#' @export

one_component_risk_extraction <- function(input_dt,
                                          threshold,
                                          max_agregate_size
){
  
  # list_z1_compo <- compo[id_comp == 1232]$z1
  # input_dt <- data_rp[z1 %in% list_z1_compo]
  verbose  = FALSE
  if(length(unique(input_dt$z1 > 1000))) verbose = TRUE
  
  list_area_z1_at_risk <- find_pbm_diff_tab(input_dt,
                                            verbose = verbose,
                                            threshold = threshold,
                                            max_agregate_size = max_agregate_size
  )
  
  dt <- copy(input_dt)
  l <- prepare_data(dt)
  dt <- l$intersecting_z2
  
  l_diff_info <- lapply(
    list_area_z1_at_risk,
    function(list_z1) return_diff_info(list_z1,dt,11)
  )
  
  tot_diff_info <- Reduce(rbind,l_diff_info)
  
  tot_diff_info
}

#' return the diff info from a input data table corresponding the whole graph which is a list of connected component of z1
#' perform the search of diff issue and return the diff information for all the corresponding component 
#' @param input_dt The tabulation table (data.frame or data.table). Each row
#' corresponds to the number of statistical units in a cross defined by a modality of the z1 nomenclature and a modality of the z2 nomenclature
#' z1 has to belong to the same related component (regarding the graph construction)
#' 
#' @param max_agregate_size Integer indicating the maximal size of agregates
#' which are tested exhaustively. If that number is too large (greater than 30), the
#' computations may not end because of the combinations number that can become very large.
#' Also the RAM can be overloaded.
#' @param save_dir Character indicating the folder where ntermediate results willbe saved (useful if there is too much components)
#' @param threshold Strictly positive integer indicating the confidentiality
#' threshold. Observations are considered at risk if one can deduce information
#' on a agregate of n observations where n < threshold.
#' 
#' @return a data.table containning the diff info for this componet (same format than the return diff info function)
#' 
#' @examples 
#' input_dt <- toy_example_6
#' all_component_risk_extraction(input_dt)
#' 
#' @export
all_component_risk_extraction <- function(input_dt,threshold = 11, max_agregate_size = 15, save_dir = "diff_info"){
  
  # list_z1_compo <- compo[id_comp %in% c(22),]$z1
  # input_dt <- data_rp[z1 %in% list_z1_compo]
  
  # first build the related components id
  z1_to_component <- return_connected_components(input_dt)
  
  # merge with input_dt
  input_dt_wit_comp <- z1_to_component[input_dt,on ="z1"]
  
  # 
  l_input_dt <- split(input_dt_wit_comp,input_dt_wit_comp$id_comp)
  
  # parellilise here  if  an option parametr say yes
  message(paste0(length(l_input_dt)), " components to handle")
  
  extract_info_and_save <- function(i){
    # i <- 1
    input_dt <- l_input_dt[[i]]
    id_compo <- names(l_input_dt[i])
    nz1 <- length(unique(input_dt$z1))
    s <- Sys.time()
    message("work on component number ",i," with ", nz1," elements of z1 ",appendLF = FALSE)
    tot_diff_info <- one_component_risk_extraction(
      input_dt,threshold = threshold,
      max_agregate_size = max_agregate_size
    )
    #print(tot_diff_info)
    if(!is.null(tot_diff_info) ) tot_diff_info$id_comp <- id_compo
    
    e <- Sys.time()
    message(round(e-s)," seconds")
    
    dir.create(save_dir,showWarnings = FALSE)
    saveRDS(
      tot_diff_info,
      file = paste0(save_dir,"/res_",id_compo,".RDS")
    )
    
    # tot_diff_info
  }
  
  l_risk_compo <-lapply(seq_along(l_input_dt),extract_info_and_save)
  
  return(z1_to_component)# in order to be able to link saved file name and z1 elements
}


#' Read diff_info info from files into a given directory and rbind them
#' @param save_dir the directory where the diff_info files (from all_component_risk_extraction) are included
#' 
#' @return a data.table containning all the diff info for all related components
#' 
#' @export

read_diff_info <- function(save_dir){
  # save_dir <- "diff_info"
  l_res<- paste0(save_dir,"/",list.files(save_dir))
  dt <- data.table()
  for(file in l_res){
    # print(file)
    # avoid loading all the files at the same time
    read_dt <- readRDS(file)
    dt <- data.table::rbindlist(list(dt,read_dt),fill = TRUE)
  }
  dt 
}

