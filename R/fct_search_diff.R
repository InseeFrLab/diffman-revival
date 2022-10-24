#Ensemble des fonctions pour rechercher les problemes de differenciation


#' Tester toutes les combinaisons possibles (jusqua une certaine taille)
#'
#' Permet de tester, pour une liste de matrice de croisement, les agregats
#' de zones de z1 conduisant a un probleme de differenciation.
#'
#' @param list_m_crois Liste de matrices de croisement.
#' @param threshold threshold de confidentialite.
#' @param max_agregate_size Entier indiquant la taille maximale des agregats
#' a tester
#'
#' @return On retourne une liste d'agregats.
search_diff_agregate=function(list_m_crois, threshold, max_agregate_size = 20){
  # list_m_crois <- l_decomp
  a <- 1
  list_agregat <- list()
  
  #Barres de progressions :
  pb1 <- progress::progress_bar$new(
    format = "  composantes connexes [:bar] :percent, reste :eta",
    total = length(list_m_crois), clear = FALSE, width= 60,
    complete = "*", incomplete = ".")
  
  pb1$tick(0)
  
  
  for(i in 1:length(list_m_crois)){
    # i <- 1
    m <- list_m_crois[[i]]
    base <- as.matrix(m)  
    if(ncol(base) == 1) colnames(base) = "1carreau"
    
    #### Preparation des donnees ####
    commune0 <- data.frame(z1 = rownames(base), men = rowSums(base), stringsAsFactors = FALSE)
    carreau0 <- data.frame(z2 = colnames(base), men = colSums(base), stringsAsFactors = FALSE)
    rownames(base) <- NULL
    colnames(base) <- NULL
    mat_conti_commune_carreau <- (base > 0) * 1L
    mat_conti_commune <- (mat_conti_commune_carreau %*% t(mat_conti_commune_carreau) > 0) * 1L
    diag(mat_conti_commune) <- 0L
    
    taille_test <- min(nrow(m), max_agregate_size)
    
    #### On test tous les agregats un a un ####
    for(taille_agregat in 1:taille_test){
      
      pb2 <- progress::progress_bar$new(
        format = "  agregats [:bar] :current ",
        total = taille_test, clear = TRUE, width= 50)
      
      #### appel de la detection de differenciation ####
      dfResultat <- differencierRcpp(iTailleCible = taille_agregat, 
                                     iSeuil = threshold,
                                     vNbObsTerritoire = commune0$men,
                                     vNbObsCarreaux = carreau0$men,
                                     mContiguiteT = mat_conti_commune,
                                     mContiguiteTC = mat_conti_commune_carreau)
      
      if(length(dfResultat)>0){
        for(k in 1:length(dfResultat)){
          num <- dfResultat[[k]][1:taille_agregat]
          list_agregat[[a]] <- commune0$z1[num]
          a <- a+1
        }
      }
      
      pb2$tick()
      
    }
    
    pb1$tick()
  }
  return(list_agregat)
}


#' Desagreger la liste d'agregat
#'
#' Permet de creer une liste de vecteur, chaque vecteur donnant
#' la composition en termes de zones de z1 de l'agregat.
#'
#' @param list_agregat Liste de chaînes de caractere. Chaque element de la liste
#' est une chaîne de caracteres donnant les nomes des zones de z1 composant l'agregat.
#' Chaque nom est separe dans la chaîne de caractere par le symbole ".".
#'
#' @return On retourne une liste de vecteur de chaînes de caracteres.
desagregate_list=function(list_agregat){
  lapply(list_agregat,function(x){ unlist(strsplit(x, ".", fixed = TRUE)) })
}

