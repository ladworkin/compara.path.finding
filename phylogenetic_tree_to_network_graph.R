library(ggraph)
library(igraph)

# 1. convert xml file to dendrogram object ####

dend.obj <- DECIPHER::ReadDendrogram('./data/GALNT3_gene_tree.nh', keepRoot = T)

# 2. convert dendrogram object object to network graph object with ggraph/igraph ####
layout <- create_layout(dend.obj[[1]], "dendrogram", repel=T)
edges <- ggraph:::getEdges.layout_dendrogram(layout)
names(layout)[1:2] <- paste0("layout.", names(layout)[1:2])
layout <- cbind(node.name = seq_len(nrow(layout)), layout)
edges <- merge(edges, layout[,c("node.name","layout.y")], by.x = "to", by.y = "node.name")
names(edges)[6] <- c("weight")

igraph.obj <- igraph::graph_from_data_frame(edges, vertices = layout)

# 3. curating list of vertices with label corresponding to taxon ####
vertex.labels <- igraph::vertex.attributes(igraph.obj)$label
taxon.specific.vertices <-  list(
    human.labels   = grep("Hsap", vertex.labels),
    gorilla.labels = grep("Ggor", vertex.labels),
    chicken.labels = grep("Ggal", vertex.labels)
)

# 4. apply pathfinding algorithms to network graph object to find shortest distance ####
bfs.searches <- function(from.vertices, to.vertices, igraph.obj){
  
  setNames(
    lapply(seq_along(1:length(from.vertices)), function(bfs.root){
      
      # performing breadth first searches ####
      paths <- 
        shortest_paths(
          graph = igraph.obj, 
          from = from.vertices[bfs.root], 
          to = to.vertices, 
          mode = "all",
          weights = NA,
          output = c("vpath"))$vpath
      
      # calculating distances ####
      distances <- 
        lapply(seq_along(1:length(paths)), function(path){
          igraph.subset <- induced_subgraph(igraph.obj, vids=as.vector(paths[[path]]))
          distances(
            graph = igraph.subset,
            v = which(as.numeric(V(igraph.subset)$name) %in% from.vertices[bfs.root]),
            to = which(as.numeric(V(igraph.subset)$name) %in% to.vertices[path]),
            weights = NULL, algorithm = "dijkstra")
        })
      
      paths.distances <- setNames(
        mapply(function(p, d){return(list(path=p, distance=d))}, paths, distances, SIMPLIFY =F),
        get.vertex.attribute(igraph.obj, 'label', to.vertices))
      
    }),
    get.vertex.attribute(igraph.obj, 'label', from.vertices))
}

chicken.to.gorilla <- bfs.searches(taxon.specific.vertices$chicken.labels, taxon.specific.vertices$gorilla.labels, igraph.obj)
chicken.to.human   <- bfs.searches(taxon.specific.vertices$chicken.labels, taxon.specific.vertices$human.labels, igraph.obj)
gorilla.to.human   <- bfs.searches(taxon.specific.vertices$gorilla.labels, taxon.specific.vertices$human.labels, igraph.obj)

