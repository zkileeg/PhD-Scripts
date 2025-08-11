#Normalizing Data Frame Title - 
name_edit = function(df, base = FALSE)
{
  #mutate if take the df data.frame and applies the function list if the variable is.character and does that to each element, the ToUpper changes all characters to upper case?
  df = mutate_if(df, is.character, list(toupper))
  
  if (base == TRUE)
  {
    #if base is exactly equal to TRUE then use these names for the columns in the df
    names(df) = c("Bait","Prey","Score")
  }
  
  else
  {
    #if base is anything but TRUE (including the default FALSE), then use these names for the columns of the data frame
    names(df) = c("Gene", "Value")
  }
  
  return(df)
}

#Mapping Gene Nicknames to ATG Names - the df is the base item from the call which is the interaction data (bait-prey) and the map is a file containing two columns, one with ATGs and the other with nicknames
mapper = function(df, map)
{
  for (i in 1:nrow(map))
  {
    #check if each bait has a match in the map file and if so replace
    df$Bait[df$Bait == map$GeneName[i]] = map$GeneID[i]
    #check if each prey has a match in the map file and if so replace
    df$Prey[df$Prey == map$GeneName[i]] = map$GeneID[i]
  }
  return(df)
}

#Now build tissue specific sub-graphs for each tissue, using only nodes with expression >1
#Arguments: graph is the network that we want to limit, df_base is the base input of interaction data, df_restrict is the df with each gene and a Value column that contains the mean expression value, cutoff is the numerical value that you want VALUE to meet or exceed for inclusion in the sub-graph, and simpl should be set to T if you want to remove self loops and multiple edges 
rest_subgraph = function(graph, df_base, df_restrict, cutoff = median(df_restrict$Value), simpl = F)
{
  #find all bait nodes expressed above the threshold
  l1 = intersect(df_restrict$Gene[df_restrict$Value>=cutoff], df_base$Bait)
  #find all bait nodes expressed above the threshold
  l2 = intersect(df_restrict$Gene[df_restrict$Value>=cutoff], df_base$Prey)
  #this will build a list of all nodes that are expression and remove any duplicates
  li = unique(append(l1,l2))
  #include any included nodes/edges in the subgraph
  sg = induced_subgraph(graph,li)
  if (simpl == T)
  {
    #this removes loops and multiple edges
    return(simplify(sg))
  }
  return(sg)
}

#Generic plot used for tissue specific networks, arguments: subgraph this is just the graph object to be plotted, communities is the predicted community membership from the walktrap in general, nodesize can be standard or defined in a list such as the one produced by the pagerank calculation
visPlot = function(subgraph, communities = rep(1, length(V(subgraph))), nodesize = 0)
{
  #V is a call from igraph that creates a vertex sequence of all vertices in a given graph object, so this is a df with the nodes and their community assignment
  nodes <- data.frame(id = V(subgraph)$name, group = communities)
  #set the default font size for the node names
  nodes$font.size<-20
  #This scales the node size
  nodes$value = (nodesize)*3000 + 10
  #Pull the edge list from the graph object
  edges <- data.frame(get.edgelist(subgraph))
  #Name the columns
  colnames(edges)<-c("from","to")
  plt = visNetwork(nodes, edges, height = "600px")%>%
    visIgraphLayout(layout = "layout_nicely") %>%
    visOptions(selectedBy = "group", 
               highlightNearest = TRUE, 
               nodesIdSelection = TRUE) %>%
    visInteraction(keyboard = TRUE,
                   dragNodes = T, 
                   dragView = T, 
                   zoomView = T)
  return(plt)
}

#Custom function to retreive Genes with confidence level = 2 of immune-relatedness and in the interaction network
#Arguments: known is a df with the list of ATGs and a column called Code with 0 (not known immune-related), 1 (some evidence), 2 (known immune function), cutoff is what values you want to highlight in the network
k_list = function(graph, known, cutoff=2)
{
  #initialize a variable
  d = c()
  #look at each name in the graph file
  for(i in 1:length(V(graph)$name))
  {
    #give it a value of 2 in the variable d if it has a Code of at least the cutoff value given
    if(V(graph)$name[i] %in% known$ATG[known$Code>=cutoff])
      d = c(d,2)
    #The other nodes get a value of 1, so group 2 contains all the known immune related genes for community plotting
    else
      d = c(d,1)
  }
  return(d)
}

# Custom function that produces a list of edges present in g1 and not g2 marked by an 'a' and a list of edges present in g2 and not g1 makred by a 'b'
#this creates the objects needed to make a network plot showing the differences between two graphs
#Arguments: g1 and g2 are the two graphs, or subgraphs, to be compared
diff_plot = function(g1,g2)
{
  #use the difference function in igraph that identifies edges that are present in g1 BUT NOT g2
  dg1 = difference(g1,g2)
  #the E function in R creates an edge sequence, in this case assigning all the edges present in g1 BUT NOT g2 to the value 1
  E(dg1)$type1 = 1
  #use the difference function in igraph that identifies edges that are present in g2 BUT NOT g1
  dg2 = difference(g2,g1)
  E(dg2)$type2 = 2
  #Now combine the two lists with the built in values of 1 and 2
  dg = union(dg1,dg2)
  temp = c()
  #Turn those edges with NA for type 1 and 2 into a and b
  E(dg)$type1[is.na(E(dg)$type1)] = 'a'
  E(dg)$type2[is.na(E(dg)$type2)] = 'b'
  return(dg)
}

#Custom function that produces a vertex list for the differences
#Arguments: graph is the object returned from the diff_plot function, which is edge lists showing which are present in g1 and g2 only, list is a list of vertices in g2 that are NOT found in g1
high_nodes = function(graph, list)
{
  d = c()
  #for every vertex name in the difference graph
  for(i in 1:length(V(graph)$name))
  {
    #If that vertex is UNIQUE to the second subgraph then assign it a 2, else assign it a 1
    if(V(graph)$name[i] %in% list)
      d = c(d,2)
    else
      d = c(d,1)
  }
  return(d)
}

#Custom plotting for difference graphs
#Arguments: subgraph is the difference graph, communities default is length of edge list (assign via diff_plot), nodesize can be PR, v-comm defaults to length of vertex list (assign via high_nodes)
visPlot2 = function(subgraph, communities = rep(1, length(E(subgraph))), nodesize = 0, v_comm = rep(1, length(V(subgraph))))
{
  #the node list has the id column with the vertex names of the subgraph and group membership assigned by high_nodes, so that nodes present in one sg are highlighted
  nodes <- data.frame(id = V(subgraph)$name, group = v_comm)
  nodes$font.size<-20
  #Assign node size (PR value)
  nodes$value = (nodesize)*3000 + 10
  #the edge list taken from the subgraph, but with communities assigned via communities (ie. the diff_plot values)
  edges <- data.frame(from = get.edgelist(subgraph)[,1], to = get.edgelist(subgraph)[,2], value = communities)
  #colnames(edges)<-c("from","to")
  #assign edge color based on the community membership as assigned via diff_plot
  edges$color <- c('black','red')[edges$value]
  #plot all that data
  plt = visNetwork(nodes, edges, height = "600px")%>%
    visIgraphLayout(layout = "layout_nicely") %>%
    visOptions(highlightNearest = TRUE, 
               nodesIdSelection = TRUE) %>%
    visInteraction(keyboard = TRUE,
                   dragNodes = T, 
                   dragView = T, 
                   zoomView = T)
  return(plt)
}

combine_graphs = function(df1, df2)
{
  df3<-merge(df1, df2, all = TRUE)
  df3<-distinct(df3)
  g3<-graph_from_data_frame(df3, directed = FALSE, vertices = NULL) 
}

#Custom function that produces an edge list for the differences
#Arguments: graph is the object returned from the diff_plot function, which is edge lists showing which are present in g1 and g2 only, list is a list of vertices in g2 that are NOT found in g1
high_edges2 = function(graph, list, list2)
{
  d = c()
  #for every vertex name in the difference graph
  for(i in 1:nrow(get.edgelist(graph)))
  {
    edgestemp<-get.edgelist(graph)
    edge1<-edgestemp[,1]
    edge2<-edgestemp[,2]
    #If either vertex is UNIQUE to the second subgraph then assign it a 2, else assign it a 1
    if(edge1[i] %in% list |  edge2[i] %in% list)
      d = c(d,"#F4EF18")
    else if(edge1[i] %in% list2 |  edge2[i] %in% list2)
      d = c(d,"#E64C2B")
    else
      d = c(d,"#184FF4")
  }
  return(d)
}

#Custom function that produces a vertex list for the differences
#Arguments: graph is the object returned from the diff_plot function, which is edge lists showing which are present in g1 and g2 only, list is a list of vertices in g2 that are NOT found in g1
high_vertices = function(graph, list, list2)
{
  d = c()
  #for every vertex name in the difference graph
  for(i in 1:length(V(graph)$name))
  {
    #If either vertex is UNIQUE to the second subgraph then assign it a 2, else assign it a 1
    if(V(graph)$name[i] %in% list)
      d = c(d,"Root")
    else if(V(graph)$name[i] %in% list2)
      d = c(d,"Leaf")
    else
      d = c(d,"Shared")
  }
  return(d)
}

#Plot to show difference between two graphs. Plots the combinedf network structure using all nodes, highlights the nodes/edges that are unique to the two subgraphs
visPlot4 = function(subgraph, vcommunities = rep(1, length(V(subgraph))), edgegroup = rep(1, nrow(get.edgelist(subgraph))), nodesize = 0)
{
  #V is a call from igraph that creates a vertex sequence of all vertices in a given graph object, so this is a df with the nodes and their community assignment
  nodes <- data.frame(id = V(subgraph)$name, group = vcommunities)
  #set the default font size for the node names
  nodes$font.size<-20
  #This scales the node size
  nodes$value = (nodesize*3000+10)
  #Pull the edge list from the graph object and add color based on edgegroup
  edges <- data.frame(get.edgelist(subgraph), color = edgegroup)
  #Name the columns
  colnames(edges)<-c("from","to", "color")
  plt = visNetwork(nodes, edges, height = "600px")%>%
    visIgraphLayout(layout = "layout_nicely") %>%
    visOptions(selectedBy = "group", 
               highlightNearest = TRUE, 
               nodesIdSelection = TRUE) %>%
    visEdges(color = edgegroup) %>%
    visGroups(groupname = "Root", color = "#F4EF18") %>%
    visGroups(groupname = "Leaf", color = "#E64C2B") %>%
    visGroups(groupname = "Shared", color = "#184FF4") %>%
    visInteraction(keyboard = TRUE,
                   dragNodes = T, 
                   dragView = T, 
                   zoomView = T)
  return(plt)
}

#Custom plotting for difference graphs
#Arguments: graph is the theoretical graph containing all nodes and edges from the 2 graphs compared, subgraph is the difference graph, communities default is length of edge list (assign via diff_plot) and assigns edge colours, nodesize can be PR, v-comm defaults to length of vertex list (assign via high_nodes)
visPlot3 = function(graph, subgraph, communities = rep(1, length(E(graph))), nodesize = 0, v_comm = rep(1, length(V(graph))))
{
  #the node list has the id column with the vertex names of the subgraph and group membership assigned by high_nodes, so that nodes present in one sg are highlighted
  nodes <- data.frame(id = V(graph)$name, group = v_comm)
  nodes$font.size<-20
  #Assign node size (PR value)
  nodes$value = (nodesize)*3000 + 10
  #the edge list taken from the subgraph, but with communities assigned via communities (ie. the diff_plot values)
  edges <- data.frame(from = get.edgelist(graph)[,1], to = get.edgelist(graph)[,2])#, value = comm)
  colnames(edges)<-c("from","to")
  #assign edge color based on the community membership as assigned via diff_plot
  #edges$color <- c('black','red')[edges$value]
  #plot all that data
  plt = visNetwork(nodes, edges, height = "600px")%>%
    visIgraphLayout(layout = "layout_nicely") %>%
    visOptions(highlightNearest = TRUE, 
               nodesIdSelection = TRUE) %>%
    visInteraction(keyboard = TRUE,
                   dragNodes = T, 
                   dragView = T, 
                   zoomView = T)
  return(plt)
}

K_Shell_Dec = function(g)
{
  temp = coreness(g)
  while(min(temp)<4)
  {
    g = delete_vertices(g, names(temp[temp<=3]))
    temp = coreness(g)
  }
  
  return(g)
}