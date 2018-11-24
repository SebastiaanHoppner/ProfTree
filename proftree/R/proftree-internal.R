.initializeNode <- function(mtree) {
  proftree_gid <- 1
  
  initNode <- function(id, mtree) {
    if (mtree$varType[mtree$splitV[id]] > 0) {
      split = partysplit(as.integer(mtree$splitV[id]), breaks = mtree$splitP[id], right = FALSE)
    } else {
      index <- as.integer(mtree$csplit[((id-1)*(mtree$maxCat)+1 ):((id-1)*(mtree$maxCat) + 1 + abs(mtree$varType[mtree$splitV[id]]) - 1)])
      index[index == 2] <- NA
      index[index == 3] <- as.integer(2)
      breaks = NULL
      split = partysplit(as.integer(mtree$splitV[id]), index = index)
    }
    node <- partynode(id = length(proftree_gid), split = split,
                      if (id*2 < 2^(mtree$maxdepth)) {
                        kids = list(
                          if (mtree$splitV[id*2] >= 0) {
                            initNode(id = id*2, mtree)
                          } else {
                            assign("proftree_gid", c(proftree_gid, id*2), inherits = TRUE)
                            partynode(as.integer(length(proftree_gid)))
                          },
                          if (mtree$splitV[id*2+1] >= 0) {
                            initNode(id = id*2+1, mtree)
                          } else {
                            assign("proftree_gid", c(proftree_gid, id*2+1), inherits = TRUE)
                            partynode(as.integer(length(proftree_gid)))
                          }
                        ) #kids
                      } else {
                        assign("proftree_gid", c(proftree_gid, id*2),   inherits = TRUE)
                        assign("proftree_gid", c(proftree_gid, id*2+1), inherits = TRUE)
                        kids = list(
                          partynode(as.integer(length(proftree_gid)-1)),
                          partynode(as.integer(length(proftree_gid) ))
                        )
                      } #else
    ) #node
  } #initNode()
  
  return(list(initNode(id = 1L, mtree = mtree), proftree_gid))
}