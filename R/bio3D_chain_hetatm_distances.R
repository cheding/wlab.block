#' A Function to get distances from chain to hetatm
#'
#' @param input_file path to PDB file (required)
#' @param chain_query query chain id (default:A)
#' @param hetatm_target target hetatm id (default:GNP)
#'
#' @return data.table with hetatm to chain (side-chain) heavy atom distances
#' @export
#'
#' @import data.table
bio3D_chain_hetatm_distances <- function(
    input_file,
    chain_query = "A",
    hetatm_target = "GNP"
){

  #load PDB structure
  pdb <- bio3d::read.pdb(input_file, rm.alt = TRUE)

  ### Atom selections
  ###########################

  #Protein atoms
  sele_protein <- bio3d::atom.select(pdb, "protein", verbose=FALSE)
  #Hydrogen atoms
  sele_H <-bio3d::atom.select(pdb, "h", verbose=FALSE)
  #Water atoms
  sele_water <- bio3d::atom.select(pdb, "water", verbose=FALSE)
  #Side chain atoms
  sele_sc <- bio3d::atom.select(pdb, "sidechain", verbose=FALSE)
  #C-alpha atoms
  sele_ca <- bio3d::atom.select(pdb, "calpha", verbose=FALSE)
  #Glycine c-alpha atoms
  sele_glyca <- bio3d::atom.select(pdb, resid = "GLY", string = "calpha", verbose=FALSE)
  #HETATM
  sele_hetatm <- bio3d::atom.select(pdb, resid = hetatm_target, verbose=FALSE)
  ### Combine atom selections
  ###########################

  #Heavy atoms
  sele_prot_HA <- bio3d::combine.select(sele_protein, sele_hetatm, operator = "OR", verbose=FALSE)
  sele_HA <- bio3d::combine.select(sele_prot_HA, sele_H, sele_water, operator = "-", verbose=FALSE)
  #Side chain heavy atoms + c-alpha for glycine
  sele_prot_sc <- bio3d::combine.select(sele_protein, sele_sc, operator = "AND", verbose=FALSE)
  sele_prot_sc_glyca <- bio3d::combine.select(sele_prot_sc, sele_glyca, sele_hetatm, operator = "OR", verbose=FALSE)
  sele_scHA <- bio3d::combine.select(sele_prot_sc_glyca, sele_H, sele_water, operator = "-", verbose=FALSE)

  #List
  sele_list <- list(
    "HA" = sele_HA,
    "scHA" = sele_scHA)
  ### Calculate minimum target chain distances
  ###########################
  result_dt <- data.table()
  for(metric in names(sele_list)){
    #Distance matrix
    pdb_sub <- bio3d::trim.pdb(pdb, sele_list[[metric]])
    dist_mat <- bio3d::dm.xyz(pdb_sub$xyz, grpby=apply(pdb_sub$atom[,c("resno", "chain")], 1, paste, collapse = "_"), scut=0, mask.lower = FALSE)
    resno_sub <- unique(pdb_sub$atom[,c("resno", "chain","resid")])
    #Ligand distance matrix
    ligand_dist <- dist_mat[resno_sub[,"chain"]==chain_query&resno_sub[,"resid"]!=hetatm_target,resno_sub[,"resid"]==hetatm_target]
    #Absolute residue number
    ligand_dist_dt <- data.table(Pos = resno_sub[resno_sub[,"chain"]==chain_query&resno_sub[,"resid"]!=hetatm_target,"resno"])
    #Minimum ligand distance
    ligand_dist_dt[, min_dist := ligand_dist]
    names(ligand_dist_dt)[2] <- paste0(metric, "min_ligand")
    if(nrow(result_dt)==0){
      result_dt <- ligand_dist_dt
    }else{
      result_dt <- merge(result_dt, ligand_dist_dt, by = "Pos", all = T)
    }
  }

  #Return
  return(result_dt)

}
