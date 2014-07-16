import networkx as nx
import matplotlib.pyplot as plt

def makeNetwork():
  G=nx.DiGraph()
  #edges = [(87 ,97) ,(25 ,22) ,(89 ,45) ,(42 ,45) ,( 5 ,25) ,(68 ,90) ,(80 ,28) ,(68 ,43) ,(89 ,53) ,( 5 , 6) ,(36 , 2) ,(23 ,18) ,(54 ,30) ,( 2 ,62) ,(33 ,97) ,(87 ,98) ,(86 ,98) ,(87 ,00) ,(69 ,66) ,(92 ,19) ,(70 ,43) ,(36 ,21) ,(43 ,58) ,(86 ,21) ,(86 ,97) ,(68 ,89) ,(75 ,83) ,(50 ,87) ,(55 ,33) ,(23 ,76) ,(86 ,99) ,(76 ,83) ,(45 ,88) ,(58 , 2) ,( 7 , 5) ,(14 ,35) ,( 1 ,62) ,( 4 ,61) ,(11 ,21) ,(87 ,21) ,(12 ,24) ,(68 ,66) ,(24 ,11) ,( 3 ,37) ,( 2 ,61)] 
  edges = [(12, 94) ,(12, 93) ,(94, 10) ,(59, 6) ,(87, 97) ,(23, 13) ,(87, 53) ,(50, 32) ,(59, 5) ,(50, 33) ,(51, 32) ,(94, 12) ,(74, 38) ,(24, 76) ,(23, 76) ,(76, 38) ,(25, 22) ,(95, 12) ,(76, 39) ,( 2, 62) ,(58, 3) ,(12, 95) ,(23, 3) ,(86, 21) ,( 9, 93) ,( 6, 5) ,( 6, 100) ,(83, 34) ,(52, 32) ,(46, 48) ,(52, 33) ,(87, 100) ,(25, 39) ,(25, 62) ,(87, 21) ,( 7, 5) ,(22, 13) ,( 6, 99) ,(87, 28) ,( 9, 65) ,(25, 13)]
  G.add_edges_from(edges)
  return G

def drawNetwork(G):
  nx.draw(G)
  plt.show()

def connectivity(G):
  print("create undirected graph: UG=G.to_undirected()")
  UG=G.to_undirected()
  print("test for connectivity: nx.is_connected(UG)\n" + str(nx.is_connected(UG)))
  print("how many connected components? nx.number_connected_components(UG)\n" + str(nx.number_connected_components(UG)))
  print("what are the subgraphs? nx.connected_components(UG)\n" + str(nx.connected_components(UG)))
  

G_comp = []
G = makeNetwork()
UG=G.to_undirected()
for comp in nx.connected_components(UG):
  G_comp.append(G.subgraph(comp))

#[(12,  94) 
#,(12,  93)
#,(94,  10)
#,(59,   6)
#,(87,  97)
#,(23,  13)
#,(87,  53)
#,(50,  32)
#,(59,   5)
#,(50,  33)
#,(51,  32)
#,(94,  12)
#,(74,  38)
#,(24,  76)
#,(23,  76)
#,(76,  38)
#,(25,  22)
#,(95,  12)
#,(76,  39)
#,( 2,  62)
#,(58,   3)
#,(12,  95)
#,(23,   3)
#,(86,  21)
#,( 9,  93)
#,( 6,   5)
#,( 6, 100)
#,(83,  34)
#,(52,  32)
#,(46,  48)
#,(52,  33)
#,(87, 100)
#,(25,  39)
#,(25,  62)
#,(87,  21)
#,( 7,   5)
#,(22,  13)
#,( 6,  99)
#,(87,  28)
#,( 9,  65)
#,(25,  13)]
