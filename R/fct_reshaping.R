#Ce script rassemble les fonctions permettant de mettre en forme les données:
#



#' Simplifier la matrice de croisement
#'
#' Simplifier la matrice de croisement revient à supprimer deux types de colones
#' (et donc supprimer des zones du zonage z2).
#' 1 - on supprime les colonnes vides, c'est-à-dire ne contenant que des 0.
#' 2 - on supprime les colonnes ne contenant qu'un seul élément différent de 0.
#'
#' @param m_crois Une matrice de croisement.
#'
#' @return En matrice on a une matrice de croisement simplifiée, c'est-à-dire
#' avec moins de colonnes.
simplify_mcrois <- function(m_crois){
  # col_to_delete_1 <- which(colSums(m_crois>0) == 1)
  # col_to_delete_2 <- which(colSums(m_crois>0) == 0)
  # col_to_delete <- c(col_to_delete_1, col_to_delete_2)
  
  #nouvelle façon:
  col_to_delete <- which(colSums(m_crois>0) <= 1)
  
  if(length(col_to_delete)>0)
    m_crois <- m_crois[, -col_to_delete, drop = FALSE]
  
  return(m_crois)
}

#' Créer la matrice de liens (ou matrice de contiguïté)
#'
#' Cette fonction permet de créer une matrice carré de booléens de taille
#' égale au nombre de zones du zonage z1.
#'
#' Un élément (i,j) de cette matrice
#' vaut TRUE (ou 1) si les zones i et j du zonage z1 sont contigues, c'est-à-dire
#' s'il existe au moins une zone de z2 recouvrant à la fois i et j. Les élements
#' de la diagonales portent la valeur FALSE.
#'
#' @param m_crois Matrice de croisement.
#'
#' @return En sortie on a une matrice carré de booléens.
matrix_liens <- function(m_crois){
  
  m_liens <- m_crois %*% t(m_crois) > 0
  diag(m_liens) <- FALSE
  
  return(m_liens)
}



#' Constuire la matrice d'adjacence du graphe
#'
#' Fonction permettant à partir de la matrice de croisement
#' de déterminer la matrice d'adjacence du graphe. Cette matrice
#' de graphe est pondérée et non symmétrique (ce qui correspond
#' à un graphe orienté).
#'
#' L'option \code{multi} permet de choisir si on prend en compte
#' ou non les zones de z2 recouvrant 3 zones de z1 ou plus. En effet
#' si on les prend en compte, alors certaines observations sont comptées
#' plusieurs fois dans le graphe, ce qui peut conduire à de mauvaises
#' interprétations.
#'
#' @param m_crois Matrice de croisement.
#' @param multi Booléen indiquant s'il faut considérer les zones de z2
#' recouvrant trois zones ou plus de z1.
#'
#' @return En sortie on obtient une matrice carré d'adjacence.
matrix_graphe <- function(m_crois, multi = TRUE){
  if(multi == FALSE) { #on enlève les carreaux sur 3 communes ou plus
    car_sel <- colSums(m_crois > 0) >= 3
    m_crois <- m_crois[, car_sel, drop = FALSE]
  }
  m_graph <- m_crois %*% (t(m_crois) > 0)
  diag(m_graph) <- 0
  return(m_graph)
}


