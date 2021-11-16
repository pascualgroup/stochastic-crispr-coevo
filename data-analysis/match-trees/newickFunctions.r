library(dplyr)


# nodes: data frame of node rows: id, parent_id, creation_time
# the root is the node with is.na(parent_id)
nodes_dataframe_to_newick <- function(nodes) {
  root <- nodes %>% filter(is.na(parent_id))
  stopifnot(nrow(root) == 1)
  
  nodes_dataframe_to_newick_recursive(nodes, NULL, root)
}

nodes_dataframe_to_newick_recursive <- function(nodes, parent, root) {
  children <- nodes %>% filter(parent_id == root$id)
  
  do.call(
    paste,
    as.list(c(
      if(nrow(children) == 0) "" else "(",
      do.call(paste,
              c(
                if(nrow(children) >= 1) {
                  as.list(sapply(1:nrow(children), function(i) nodes_dataframe_to_newick_recursive(nodes, root, children[i,])))
                } else {
                  list()
                },
                sep = ","
              )
      ),
      if(nrow(children) == 0) "" else ")",
      sprintf("%d", root$id),
      if(!is.null(parent)) {
        sprintf(":%f", root$creation_time - (if(is.null(parent)) 0.0 else parent$creation_time))
      } else {
        ""
      },
      sep = ""
    )
    ))
}



# #################

args <- commandArgs(trailingOnly = TRUE)
file = args[1]
name = strsplit(file, '.txt')
out = paste(name,".nwk",sep="")

Datos <- read.delim(file=file, header=FALSE, sep="\t")
names(Datos) <- c("Recordtime","id","parent_id","creation_time")
MyData = subset(Datos, select = -c(Recordtime))
newick_str <- nodes_dataframe_to_newick(MyData)

write(newick_str, file = out)