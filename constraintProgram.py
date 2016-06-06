# -*- coding:Utf8 -*-

import math
import json
import random
import networkx as nx
from matplotlib import pyplot as plt
from networkx.readwrite import json_graph

dataFile = "newDataSet_tags.txt"
adjacencyDelimiter = ":"
jsonFile = "scalefree.json"


def main( ) :
	g = graphGenerators(100, "scalefree")

def graphGenerators( n, graphName ) :
	if (graphName == "star") :
		return make_reflexive (nx.star_graph (n).to_directed ())
	elif (graphName == "twostars") :
		g = nx.star_graph (n)
		h = nx.star_graph (n)
		newLabels = { }
		for i in range (n + 1) :
			newLabels[i] = i + n + 1
		h = nx.relabel_nodes (h, newLabels)
		gh = nx.union (g, h)
		if (g.degree (0) > 1) :
			gh.add_edge (1, n + 1)
		else :
			gh.add_edge (0, n + 1)
		gh.add_edge (n + 1, n + 2)

		return gh.to_directed ()
	elif (graphName == "complete") :
		return nx.complete_graph (n).to_directed ()
	elif (graphName == "scalefree") :
		return nx.scale_free_graph (n)


def generateRandomWeights(graph):
	avgdeg = 0
	for node in graph.nodes():
		avgdeg+= graph.degree(node)
	avgdeg = avgdeg/graph.order()
	edgeList = []
	for u,v in graph.edges():
		edgeList.append((u,v))
	for u,v in edgeList:
		graph.remove_edge(u,v)
		if not graph.has_edge(u,v):
			#graph.add_edge(u,v,weight=random.normalvariate(avgdeg, 1))
			graph.add_edge(u,v,weight=math.log(len(graph.neighbors(u))+len(graph.neighbors(v))+math.exp(1)))
			#graph.add_edge(u,v,weight=2+math.log(math.exp((len(graph.neighbors(u))+len(graph.neighbors(v)))/2)))

	return graph


def showGraph( graph ) :
	# pos = nx.spring_layout(graph, k=1)
	# nx.draw_networkx_nodes(graph, pos, node_size=4000, node_color="white")
	# nx.draw_networkx_edges(graph, pos, width=3, alpha=0.5, edge_color="black", arrows=True)
	# nx.draw_networkx_labels(graph, pos, font_size=12)
	# nx.draw_networkx_edge_labels(graph, pos, label_pos=0.6, font_size = 10)
	# plt.axis('off')
	nx.draw_spring (graph)
	plt.show ()


# On ne garde pas les noeuds qui n'ont pas de voisins
def createPartialGraphFromDataset( dataFile, size ) :
	g = nx.DiGraph ()
	with open (dataFile) as f :
		try :
			for line in f :
				if g.order() >= size :
					break
				row = line.split (adjacencyDelimiter)
				if len(row)>2:
					x = row[0]
					if len(row)<4:
						raise ValueError("Missing something for node : "+x)
					elif ',' not in row[1] and ',' not in row[2]:
						y = row[1]
						w = row[2]
						t = row[3]
						g.add_edge (x, y, weight=int (w), relation = int(t))
					elif len (row[1].split (',')) == len (row[2].split (',')):
						for j in range (len (row[1].split (','))) :
							if g.order()>=size:
								break
							else:
								y = row[1].split (',')[j]
								w = row[2].split (',')[j]
								t= row[3].split(',')[j]
								g.add_edge (x, y, weight=int (w), relation = int(t))
					else:
						raise ValueError ("Nombre de poids différents du nombre de voisins au noeud " + x)

		except ValueError as er :
			print (er)
			return g

	nodeOutDegrees = {}
	for node in g.nodes():
		nodeOutDegrees[node] = g.out_degree(node)
	nx.set_node_attributes(g, "outdegree", nodeOutDegrees)
	return g


def createGraphFromDataset( dataFile ) :
	g = nx.DiGraph ()
	with open (dataFile) as f :
		for line in f :
			row = line.split (adjacencyDelimiter)
			if (len (row) > 1) :
				x = row[0]
				if (len (row[1].split (",")) == 1) :
					y = row[1]
					g.add_edge (x, y)
				else :
					for y in row[1].split (',') :
						g.add_edge (x, y)
	return g


def findCommunities(graph):
	return graph


def make_reflexive( graph ) :
	for u in graph.nodes () :
		graph.add_edge (u, u)
	return graph


def writeJsonFile( graph, filename ) :
	data = json_graph.node_link_data (graph)
	with open (filename, 'w') as out :
		json.dump (data, out, indent=4)


def readJsonFile(filename):
	g = nx.DiGraph()
	d = json.load(open(filename))
	return json_graph.node_link_graph(d)


def writeGML( graph, fname ) :
	nx.write_gml (graph, fname)


def readGML( fname ) :
	return nx.read_gml (fname)


if __name__ == "__main__" :
	main ()
