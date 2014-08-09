import networkx as nx
import matplotlib.pyplot as plt

def makeNetwork(edges):
  G=nx.DiGraph()
  G.add_edges_from(edges)
  return G

def drawNetwork(G):
  #nx.draw(G)
  nx.draw_spring(G)
  plt.show()

def connectivity(G):
  print("create undirected graph: UG=G.to_undirected()")
  UG=G.to_undirected()
  print("test for connectivity: nx.is_connected(UG)\n" + str(nx.is_connected(UG)))
  print("how many connected components? nx.number_connected_components(UG)\n" + str(nx.number_connected_components(UG)))
  print("what are the subgraphs? nx.connected_components(UG)\n" + str(nx.connected_components(UG)))
  

"""
import networkx as nx
G_comp = []
G = net.makeNetwork(net.edges)
UG=G.to_undirected()
for comp in nx.connected_components(UG):
  G_comp.append(G.subgraph(comp))
"""

def getEdge(a, b):
  for edge in edges:
    if edge[0] == a and edge[1] == b:
      print edge
    if edge[0] == b and edge[1] == a:
      print edge

### Remove degree 1 nodes. Not necessarily useful.
"""
deg = UG.degree()
remove = [n for n in deg if deg[n] == 1]
UG.remove_nodes_from(remove)
"""

"""
[(  94   ,  10,{"weight":-0.99906167}), 
(  74   ,  38 ,{"weight":-1.05451978}),
(  25   ,  22 ,{"weight":-0.77492709}),
(  87   ,  97 ,{"weight":-1.16548936}),
(  12   ,  93 ,{"weight":-0.67502484}),
(  51   ,  32 ,{"weight":-1.34700698}),
(  59   ,   6 ,{"weight":-0.52630844}),
(  59   ,   5 ,{"weight":-0.50794252}),
(  23   ,  13 ,{"weight":-1.02944321}),
(  50   ,  33 ,{"weight":-1.03466065}),
(  24   ,  22 ,{"weight":-1.00572381}),
(  58   ,   3 ,{"weight":-0.77431673}),
(  25   ,  62 ,{"weight":-1.09186573}),
(  87   ,  30 ,{"weight":-1.15669341}),
(  76   ,  38 ,{"weight":-0.98950709}),
(  24   ,  76 ,{"weight":-0.70983965}),
(  50   ,  32 ,{"weight":-1.17892331}),
(  36   ,  29 ,{"weight":-0.99012952}),
(  52   ,  30 ,{"weight":-1.39985166}),
(  94   ,  11 ,{"weight":-1.01653161}),
(  87   ,  53 ,{"weight":-0.68632513}),
(  67   ,  41 ,{"weight":-1.39995248}),
(   9   ,  66 ,{"weight":-0.93531419}),
(   2   ,  62 ,{"weight":-1.40158199}),
(  87   ,  67 ,{"weight":-0.88329873}),
(  22   ,  13 ,{"weight":-0.93648341}),
(  46   ,  48 ,{"weight":-0.63045781}),
(  52   ,  32 ,{"weight":-1.24896709}),
(   6   , 100 ,{"weight":-1.37004754}),
(  24   ,  23 ,{"weight":-0.70642538}),
(  23   ,   3 ,{"weight":-0.67746910}),
(  86   ,  21 ,{"weight":-1.32933747}),
(  76   ,  18 ,{"weight":-0.72622137}),
(  83   ,  34 ,{"weight":-1.17575657}),
(   7   ,   5 ,{"weight":-0.70165477}),
(  92   ,   8 ,{"weight":-1.03193445}),
(   9   ,  93 ,{"weight":-0.60056743}),
(  87   , 100 ,{"weight":-1.27971632}),
(  45   ,  90 ,{"weight":-0.75786333}),
(  54   ,  30 ,{"weight":-1.59772492}),
(  47   ,  32 ,{"weight":-1.24287802}),
(  66   ,  26 ,{"weight":-1.07536961}),
(  25   ,  13 ,{"weight":-0.95336102}),
(  23   ,  76 ,{"weight":-0.66299697}),
(  89   ,  19 ,{"weight":-1.37060144}),
(  42   ,  34 ,{"weight":-1.90246392}),
(  24   ,  28 ,{"weight":-0.95442715}),
(  75   ,  38 ,{"weight":-0.99620159}),
(  76   ,  55 ,{"weight":-0.89674595}),
(  48   ,  32 ,{"weight":-1.19765123}),
(  57   ,  19 ,{"weight":-1.30582798}),
(   7   ,  54 ,{"weight":-1.37895332}),
(   7   ,  36 ,{"weight":-0.67291726}),
(   9   ,  65 ,{"weight":-0.91016518}),
(  87   ,  28 ,{"weight":-0.93211341}),
(  23   ,  12 ,{"weight":-0.73371229}),
(  49   ,  96 ,{"weight":-0.84993752}),
(  25   ,  39 ,{"weight":-0.75610110}),
(  89   ,  53 ,{"weight":-0.62761656}),
(  79   ,  40 ,{"weight":-1.24754376}),
(  64   ,  11 ,{"weight":-0.82225326}),
(  45   ,  88 ,{"weight":-0.66060942}),
(  51   ,  52 ,{"weight":-1.13809349}),
(  86   ,  53 ,{"weight":-0.62460934}),
(   2   ,  82 ,{"weight":-0.78019347}),
(  23   ,   2 ,{"weight":-0.67959760}),
(  76   ,  39 ,{"weight":-1.11688569}),
(  58   ,   4 ,{"weight":-0.54442624}),
(  17   ,  20 ,{"weight":-1.58238338}),
(  91   ,  94 ,{"weight":-0.72590126}),
(  82   ,  56 ,{"weight":-0.91958186}),
(  37   ,  29 ,{"weight":-0.86598534}),
(  23   ,  11 ,{"weight":-0.72240525}),
(  92   ,  19 ,{"weight":-1.54004847}),
(  50   ,  96 ,{"weight":-0.80751367}),
(   5   ,  54 ,{"weight":-1.27045656}),
(  16   ,  29 ,{"weight":-1.04180152}),
(  48   ,  96 ,{"weight":-0.85846833}),
(  51   ,  61 ,{"weight":-0.63813459}),
(  80   ,  40 ,{"weight":-1.42671379}),
(  58   ,  10 ,{"weight":-0.82253756}),
(  87   ,  21 ,{"weight":-1.02461679}),
(  75   ,  90 ,{"weight":-0.70117025}),
(   6   ,   5 ,{"weight":-0.90974666}),
(  63   ,  26 ,{"weight":-1.10870797}),
(  92   ,  10 ,{"weight":-1.03644632}),
(  56   ,  34 ,{"weight":-1.27295345}),
(  73   ,  55 ,{"weight":-0.88844131}),
(  24   ,  74 ,{"weight":-0.60146022}),
(  47   ,  96 ,{"weight":-0.83988980}),
(  63   ,  33 ,{"weight":-0.96907083}),
(  81   ,  19 ,{"weight":-1.16246004}),
(  61   ,   7 ,{"weight":-0.55698192}),
(  69   ,  64 ,{"weight":-0.89908872}),
(  63   ,  11 ,{"weight":-0.77892974}),
(  45   ,  68 ,{"weight":-0.54356707}),
(   4   ,  46 ,{"weight":-0.73659845}),
(  56   ,  32 ,{"weight":-1.04486612}),
(  46   ,  52 ,{"weight":-0.58180809}),
(  50   ,  61 ,{"weight":-0.60216580}),
(   9   ,  33 ,{"weight":-0.90288330}),
(  57   ,  34 ,{"weight":-0.93393934}),
(  45   ,  86 ,{"weight":-0.60372771}),
(  63   ,  62 ,{"weight":-1.60245455}),
(  22   ,  69 ,{"weight":-0.55550905}),
(  52   ,  61 ,{"weight":-0.72579846}),
(  88   ,  26 ,{"weight":-1.10415125}),
(  91   ,   8 ,{"weight":-0.82462634}),
(  61   ,  41 ,{"weight":-1.04753955}),
(   6   ,  99 ,{"weight":-1.00247325}),
(  61   ,  39 ,{"weight":-0.84780047}),
(  16   ,  35 ,{"weight":-1.19582862}),
(  28   ,  15 ,{"weight":-0.62593498}),
(  76   ,  40 ,{"weight":-1.40790712}),
(  63   ,  10 ,{"weight":-0.70971354}),
(  23   , 101 ,{"weight": 0.28597114}),
(  10   ,   5 ,{"weight": 0.34465904}),
(  99   ,  93 ,{"weight": 0.67791270}),
(   6   , 101 ,{"weight": 0.19991742}),
(  94   ,   6 ,{"weight": 1.14665061}),
(  63   ,  87 ,{"weight": 0.47717247}),
(  64   , 101 ,{"weight": 0.19779851}),
( 100   ,  24 ,{"weight": 0.32148980}),
(  63   ,  88 ,{"weight": 0.48635535}),
(  55   , 101 ,{"weight": 0.14245890}),
(  69   ,  24 ,{"weight": 0.44337488}),
( 101   ,  50 ,{"weight": 0.53073262}),
(  90   ,   6 ,{"weight": 0.62868719}),
( 100   ,  15 ,{"weight": 0.72928818}),
(  44   ,   6 ,{"weight": 0.24980261}),
(  74   ,  43 ,{"weight": 0.79229448}),
(  46   ,  23 ,{"weight": 0.31378667}),
(  86   ,   6 ,{"weight": 0.32212558}),
(  37   , 101 ,{"weight": 1.07604195}),
(  58   ,  23 ,{"weight": 0.72975539}),
(  65   , 101 ,{"weight": 0.13010747}),
(  29   ,  57 ,{"weight": 1.23968109}),
( 100   ,  23 ,{"weight": 0.31812480}),
(  44   ,   7 ,{"weight": 0.40377123}),
( 101   ,  54 ,{"weight": 0.57212957}),
(  77   ,  44 ,{"weight": 0.38904439}),
(  18   ,  60 ,{"weight": 0.54601566}),
(  91   ,  87 ,{"weight": 0.40290014}),
(  73   ,  61 ,{"weight": 0.52101271}),
(  61   ,  24 ,{"weight": 0.27573143}),
(   4   ,  76 ,{"weight": 0.60411335}),
( 100   ,  78 ,{"weight": 1.75168194}),
(   5   ,  44 ,{"weight": 0.31791081}),
(  20   ,  87 ,{"weight": 0.51135270}),
(  10   ,   6 ,{"weight": 0.18482325}),
(  19   ,  68 ,{"weight": 0.58196358}),
( 100   ,  18 ,{"weight": 0.78586109}),
(  76   ,   6 ,{"weight": 0.21406898}),
(  46   ,   6 ,{"weight": 0.33583645}),
(  31   ,  93 ,{"weight": 0.88476357}),
(  76   ,   5 ,{"weight": 0.30275195}),
(  76   ,  93 ,{"weight": 0.39403515}),
(  69   ,  23 ,{"weight": 0.84354074}),
(  76   ,  61 ,{"weight": 0.54376695}),
(  22   , 101 ,{"weight": 0.13258740}),
(  68   ,  24 ,{"weight": 0.20635145}),
(   2   ,  16 ,{"weight": 0.68148539}),
(  66   ,  86 ,{"weight": 1.11191323}),
(  19   ,  97 ,{"weight": 1.59216157}),
(  87   ,   6 ,{"weight": 0.31633140}),
(   7   ,  75 ,{"weight": 0.39198092}),
(  62   ,  74 ,{"weight": 0.61939764}),
(  45   ,  24 ,{"weight": 0.22205842}),
(  40   ,  65 ,{"weight": 1.54948829}),
(  24   ,  60 ,{"weight": 0.47545924}),
(  17   ,  69 ,{"weight": 0.63531609}),
( 100   ,  14 ,{"weight": 1.35038482}),
(  82   ,  25 ,{"weight": 0.26350582}),
(  17   ,  68 ,{"weight": 0.59918179}),
(  50   ,  44 ,{"weight": 0.66187954}),
(  44   ,  25 ,{"weight": 0.38901889}),
(  17   ,  88 ,{"weight": 0.37108656}),
(  65   ,  57 ,{"weight": 1.34523245}),
(  63   ,  86 ,{"weight": 0.59147248}),
(   8   ,  61 ,{"weight": 0.39559879}),
(  88   ,   6 ,{"weight": 0.38617049}),
( 100   ,  81 ,{"weight": 0.80239428}),
(   3   , 101 ,{"weight": 0.10512044}),
(  63   ,  50 ,{"weight": 1.25820317}),
(  37   ,  93 ,{"weight": 0.45105717}),
(  11   ,   4 ,{"weight": 0.53028448}),
(  33   ,  92 ,{"weight": 0.32874549}),
(  93   ,   6 ,{"weight": 0.27808913}),
(  18   ,  86 ,{"weight": 0.75090568}),
( 100   ,  93 ,{"weight": 0.48563971}),
(  49   ,  44 ,{"weight": 0.53736083}),
(  89   ,   6 ,{"weight": 0.38820098}),
(  25   ,  60 ,{"weight": 0.33484286}),
(  92   , 101 ,{"weight": 0.09140783}),
(  63   , 101 ,{"weight": 0.13865542})]
"""
edges = [(94,10,{"weight":-0.99906167}),(74,38,{"weight":-1.05451978}),(25,22,{"weight":-0.77492709}),(87,97,{"weight":-1.16548936}),(12,93,{"weight":-0.67502484}),(51,32,{"weight":-1.34700698}),(59,6,{"weight":-0.52630844}),(59,5,{"weight":-0.50794252}),(23,13,{"weight":-1.02944321}),(50,33,{"weight":-1.03466065}),(24,22,{"weight":-1.00572381}),(58,3,{"weight":-0.77431673}),(25,62,{"weight":-1.09186573}),(87,30,{"weight":-1.15669341}),(76,38,{"weight":-0.98950709}),(24,76,{"weight":-0.70983965}),(50,32,{"weight":-1.17892331}),(36,29,{"weight":-0.99012952}),(52,30,{"weight":-1.39985166}),(94,11,{"weight":-1.01653161}),(87,53,{"weight":-0.68632513}),(67,41,{"weight":-1.39995248}),(9,66,{"weight":-0.93531419}),(2,62,{"weight":-1.40158199}),(87,67,{"weight":-0.88329873}),(22,13,{"weight":-0.93648341}),(46,48,{"weight":-0.63045781}),(52,32,{"weight":-1.24896709}),(6,100,{"weight":-1.37004754}),(24,23,{"weight":-0.70642538}),(23,3,{"weight":-0.67746910}),(86,21,{"weight":-1.32933747}),(76,18,{"weight":-0.72622137}),(83,34,{"weight":-1.17575657}),(7,5,{"weight":-0.70165477}),(92,8,{"weight":-1.03193445}),(9,93,{"weight":-0.60056743}),(87,100,{"weight":-1.27971632}),(45,90,{"weight":-0.75786333}),(54,30,{"weight":-1.59772492}),(47,32,{"weight":-1.24287802}),(66,26,{"weight":-1.07536961}),(25,13,{"weight":-0.95336102}),(23,76,{"weight":-0.66299697}),(89,19,{"weight":-1.37060144}),(42,34,{"weight":-1.90246392}),(24,28,{"weight":-0.95442715}),(75,38,{"weight":-0.99620159}),(76,55,{"weight":-0.89674595}),(48,32,{"weight":-1.19765123}),(57,19,{"weight":-1.30582798}),(7,54,{"weight":-1.37895332}),(7,36,{"weight":-0.67291726}),(9,65,{"weight":-0.91016518}),(87,28,{"weight":-0.93211341}),(23,12,{"weight":-0.73371229}),(49,96,{"weight":-0.84993752}),(25,39,{"weight":-0.75610110}),(89,53,{"weight":-0.62761656}),(79,40,{"weight":-1.24754376}),(64,11,{"weight":-0.82225326}),(45,88,{"weight":-0.66060942}),(51,52,{"weight":-1.13809349}),(86,53,{"weight":-0.62460934}),(2,82,{"weight":-0.78019347}),(23,2,{"weight":-0.67959760}),(76,39,{"weight":-1.11688569}),(58,4,{"weight":-0.54442624}),(17,20,{"weight":-1.58238338}),(91,94,{"weight":-0.72590126}),(82,56,{"weight":-0.91958186}),(37,29,{"weight":-0.86598534}),(23,11,{"weight":-0.72240525}),(92,19,{"weight":-1.54004847}),(50,96,{"weight":-0.80751367}),(5,54,{"weight":-1.27045656}),(16,29,{"weight":-1.04180152}),(48,96,{"weight":-0.85846833}),(51,61,{"weight":-0.63813459}),(80,40,{"weight":-1.42671379}),(58,10,{"weight":-0.82253756}),(87,21,{"weight":-1.02461679}),(75,90,{"weight":-0.70117025}),(6,5,{"weight":-0.90974666}),(63,26,{"weight":-1.10870797}),(92,10,{"weight":-1.03644632}),(56,34,{"weight":-1.27295345}),(73,55,{"weight":-0.88844131}),(24,74,{"weight":-0.60146022}),(47,96,{"weight":-0.83988980}),(63,33,{"weight":-0.96907083}),(81,19,{"weight":-1.16246004}),(61,7,{"weight":-0.55698192}),(69,64,{"weight":-0.89908872}),(63,11,{"weight":-0.77892974}),(45,68,{"weight":-0.54356707}),(4,46,{"weight":-0.73659845}),(56,32,{"weight":-1.04486612}),(46,52,{"weight":-0.58180809}),(50,61,{"weight":-0.60216580}),(9,33,{"weight":-0.90288330}),(57,34,{"weight":-0.93393934}),(45,86,{"weight":-0.60372771}),(63,62,{"weight":-1.60245455}),(22,69,{"weight":-0.55550905}),(52,61,{"weight":-0.72579846}),(88,26,{"weight":-1.10415125}),(91,8,{"weight":-0.82462634}),(61,41,{"weight":-1.04753955}),(6,99,{"weight":-1.00247325}),(61,39,{"weight":-0.84780047}),(16,35,{"weight":-1.19582862}),(28,15,{"weight":-0.62593498}),(76,40,{"weight":-1.40790712}),(63,10,{"weight":-0.70971354}),(23,101,{"weight":0.28597114}),(10,5,{"weight":0.34465904}),(99,93,{"weight":0.67791270}),(6,101,{"weight":0.19991742}),(94,6,{"weight":1.14665061}),(63,87,{"weight":0.47717247}),(64,101,{"weight":0.19779851}),(100,24,{"weight":0.32148980}),(63,88,{"weight":0.48635535}),(55,101,{"weight":0.14245890}),(69,24,{"weight":0.44337488}),(101,50,{"weight":0.53073262}),(90,6,{"weight":0.62868719}),(100,15,{"weight":0.72928818}),(44,6,{"weight":0.24980261}),(74,43,{"weight":0.79229448}),(46,23,{"weight":0.31378667}),(86,6,{"weight":0.32212558}),(37,101,{"weight":1.07604195}),(58,23,{"weight":0.72975539}),(65,101,{"weight":0.13010747}),(29,57,{"weight":1.23968109}),(100,23,{"weight":0.31812480}),(44,7,{"weight":0.40377123}),(101,54,{"weight":0.57212957}),(77,44,{"weight":0.38904439}),(18,60,{"weight":0.54601566}),(91,87,{"weight":0.40290014}),(73,61,{"weight":0.52101271}),(61,24,{"weight":0.27573143}),(4,76,{"weight":0.60411335}),(100,78,{"weight":1.75168194}),(5,44,{"weight":0.31791081}),(20,87,{"weight":0.51135270}),(10,6,{"weight":0.18482325}),(19,68,{"weight":0.58196358}),(100,18,{"weight":0.78586109}),(76,6,{"weight":0.21406898}),(46,6,{"weight":0.33583645}),(31,93,{"weight":0.88476357}),(76,5,{"weight":0.30275195}),(76,93,{"weight":0.39403515}),(69,23,{"weight":0.84354074}),(76,61,{"weight":0.54376695}),(22,101,{"weight":0.13258740}),(68,24,{"weight":0.20635145}),(2,16,{"weight":0.68148539}),(66,86,{"weight":1.11191323}),(19,97,{"weight":1.59216157}),(87,6,{"weight":0.31633140}),(7,75,{"weight":0.39198092}),(62,74,{"weight":0.61939764}),(45,24,{"weight":0.22205842}),(40,65,{"weight":1.54948829}),(24,60,{"weight":0.47545924}),(17,69,{"weight":0.63531609}),(100,14,{"weight":1.35038482}),(82,25,{"weight":0.26350582}),(17,68,{"weight":0.59918179}),(50,44,{"weight":0.66187954}),(44,25,{"weight":0.38901889}),(17,88,{"weight":0.37108656}),(65,57,{"weight":1.34523245}),(63,86,{"weight":0.59147248}),(8,61,{"weight":0.39559879}),(88,6,{"weight":0.38617049}),(100,81,{"weight":0.80239428}),(3,101,{"weight":0.10512044}),(63,50,{"weight":1.25820317}),(37,93,{"weight":0.45105717}),(11,4,{"weight":0.53028448}),(33,92,{"weight":0.32874549}),(93,6,{"weight":0.27808913}),(18,86,{"weight":0.75090568}),(100,93,{"weight":0.48563971}),(49,44,{"weight":0.53736083}),(89,6,{"weight":0.38820098}),(25,60,{"weight":0.33484286}),(92,101,{"weight":0.09140783}),(63,101,{"weight":0.13865542})]
