library(DiagrammeR)


grViz("digraph{

      graph[rankdir = LR]
  
      node[shape = rectangle, style = filled]
  
      node[fillcolor = Coral, margin = 0.2]
      A[label = 'Figure 1: Map']
      B[label = 'Figure 2: Metrics']
  
      node[fillcolor = Cyan, margin = 0.2]
      C[label = 'Figures.Rmd']
  
      node[fillcolor = Violet, margin = 0.2]
      D[label = 'Analysis_1.R']
      E[label = 'Analysis_2.R']
  
      subgraph cluster_0 {
        graph[shape = rectangle]
        style = rounded
        bgcolor = Darkorange
    
        label = 'Data Source 1'
        node[shape = rectangle, fillcolor = LemonChiffon, margin = 0.25]
        F[label = 'my_dataframe_1.csv']
        G[label = 'my_dataframe_2.csv']
      }
  
      subgraph cluster_1 {
         graph[shape = rectangle]
         style = rounded
         bgcolor = Gold
    
         label = 'Data Source 2'
         node[shape = rectangle, fillcolor = LemonChiffon, margin = 0.25]
         H[label = 'my_dataframe_3.csv']
         I[label = 'my_dataframe_4.csv']
      }
  
      edge[color = black, arrowhead = vee, arrowsize = 1.25]
      C -> {A B}
      D -> C
      E -> C
      F -> D
      G -> D
      H -> E
      I -> E
      
      }")


grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = box,
        fontname = Helvetica,
        style = filled,  
        color = SeaGreen]
  a [label = 'Hierarchical Cluster Analysis \\n implemented
  \\n Linkage = ward.D \\n Distance = Euclidean']; 
  b [label = 'Fuzzy Correspondence Analysis']; 
  c [label = 'PERMANOVA']; 
  d [label = 'Random Forest']

  node [shape = circle,
        fixedsize = false,
        width = 0.9] // sets as circles
  e [label = 'Gap statistic']; 


  # several 'edge' statements
  # labels used here are for the edges (i.e. arrows)
  a -> e 
  [label = ' Optimal number \\n of groups?']
  e -> d
  [label = ' Most important traits \\n for clustering of taxa?']
  b -> c

}
")

# TODO: Fix borders of tables
# Can we use CSS here?
# See further: https://cyberhelp.sesync.org/blog/visualization-with-diagrammeR.html

grViz("
  digraph D {

  node [shape = plaintext fontname= 'Fira Code' fontsize= '13']

  task_HCA [ label=<
   <table border='1' cellborder='0' cellspacing='1'>
     <tr><td align='left'><b>Hierarchical Cluster <br></br> Analysis</b></td></tr>
     <tr><td align='left'><font color = 'darkgreen'>implemented</font></td></tr>
     <tr><td align='left'>Linkage = ward.D</td></tr>
     <tr><td align='left'>Distance = Euclidean</td></tr>
   </table>>];

  task_GAP [ label=<
   <table border='1' cellborder='0' cellspacing='1'>
     <tr><td align='left'><b>Gap statistic</b></td></tr>
     <tr><td align='left'><font color='darkgreen'>implemented</font></td></tr>
   </table>>];

   task_RF [ label=<
   <table style = 'table-layout: auto; width: 100%;' border='1' cellborder='0' cellspacing='1'>
     <tr><td align='left'><b>Random <br></br> Forest</b></td></tr>
     <tr><td align='left'><font color='darkgreen'>implemented</font></td></tr>
   </table>>];

   task_FCA [ label=<
   <table style = 'table-layout: auto; width: 100%;' border='1' cellborder = '0' cellspacing='1'>
     <tr><td align = 'left'> Fuzzy Correspondence <br></br> Analysis</td></tr>
     <tr><td align='left'><font color='red'> ToDo </font></td></tr>
   </table>>];

   task_PMAV [ label=<
   <table border='1' cellborder='0' cellspacing='1'>
     <tr><td align='left'><b>PERMANOVA </b></td></tr>
     <tr><td align='left'><font color='darkgreen'>implemented</font></td></tr>
   </table>>];

   node [shape = circle,
        fixedsize = false,
         width = 1] // sets as circles
    Q1 [label = 'Testing different \\n linkage methods? \\n Appropriate distance?'];

   task_HCA -> Q1
   task_HCA -> task_GAP [label = ' Optimal number \\n of groups?'] 
   task_GAP -> task_RF [label = ' Most important traits \\n for grouping of taxa?']
   task_FCA -> task_PMAV;
}
")


# can substitue with R expressions 
# ? Can we also put R - output into our flow chart
# grViz("
# digraph a_nice_graph {

# # node definitions with substituted label text
# node [fontname = Helvetica]
# a [label = '@@1']
# b [label = '@@2-1']
# c [label = '@@2-2']
# d [label = '@@2-3']
# e [label = '@@2-4']
# f [label = '@@2-5']
# g [label = '@@2-6']
# h [label = '@@2-7']
# i [label = '@@2-8']
# j [label = '@@2-9']

# # edge definitions with the node IDs
# a -> {b c d e f g h i j}
# }

# [1]: 'top'
# [2]: 10:20
# ")