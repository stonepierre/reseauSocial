# -*- coding:Utf8 -*-
__author__ = 'Frédéric Simard'

import sys
sys.path.append("/usr/local/lib/python2.7/site-packages")
# import graph_toolCommunity
import networkx as nx
import numpy.random as ran
from matplotlib import pyplot as plt
from graph_tool.all import *

direc = "/output/"

def main():
	h = graph_tool.load_graph("stoneshi0622.gml")
	graph_tool.draw.graph_draw(h, output="stoneshi0622.pdf")


def createGraphsAndCommunities():
	g = nx.scale_free_graph(500, alpha=0.40, beta=0.40, gamma=0.20)
	g1 = nx.powerlaw_cluster_graph(500, 10, 0.2)
	g2 = nx.barabasi_albert_graph(500, 10)
	g3 = nx.newman_watts_strogatz_graph(500, 10, 0.2)
	nx.write_graphml (g, direc+"sfg.graphml")
	nx.write_graphml(g1, direc+"pcg.graphml")
	nx.write_graphml(g2, direc+"bag.graphml")
	nx.write_graphml(g3, direc+"nwsg.graphml")

	graphs = {}
	graphs["sfg"] = graph_tool.load_graph(direc+"sfg.graphml")
	graphs["pcg"] = graph_tool.load_graph(direc+"pcg.graphml")
	graphs["bag"] = graph_tool.load_graph(direc+"bag.graphml")
	graphs["nwsg"] = graph_tool.load_graph(direc+"nwsg.graphml")
	graphs["price"] = graph_tool.generation.price_network(1000)
	
	for i,h in graphs.iteritems():
		s = graph_tool.community.minimize_blockmodel_dl(h)
		b = s.b
		graph_tool.draw.graph_draw(h, vertex_fill_color=b, vertex_shape=b, output=direc+"block"+str(i)+".pdf")
		
		com = graph_tool.community.community_structure(h, 10000, 20)
		graph_tool.draw.graph_draw(h, vertex_fill_color=com, vertex_shape=com, output=direc+"community"+str(i)+".pdf")

		state = graph_tool.community.minimize_nested_blockmodel_dl(h)
		graph_tool.draw.draw_hierarchy(state, output=direc+"nestedblock"+str(i)+".pdf")

		pagerank = graph_tool.centrality.pagerank(h)
		graph_tool.draw.graph_draw(h, vertex_fill_color=pagerank, vertex_size = graph_tool.draw.prop_to_size(pagerank, mi=5, ma=15), vorder=pagerank, output=direc+"pagerank"+str(i)+".pdf")
		h.set_reversed(is_reversed=True)
		pagerank = graph_tool.centrality.pagerank(h)
		graph_tool.draw.graph_draw(h, vertex_fill_color=pagerank, vertex_size = graph_tool.draw.prop_to_size(pagerank, mi=5, ma=15), vorder=pagerank, output=direc+"reversed_pagerank"+str(i)+".pdf")


def createAndSaveBlockModel(graph):
		nx.write_graphml(graph, "graph.graphml")
		h = graph_tool.load_graph("graph.graphml")
		
		s = graph_tool.inference.minimize_blockmodel_dl(h)
		graph_tool.draw.graph_draw(h, vertex_text=h.vertex_properties["name"], vertex_font_size=10, output_size=(700, 700),  vertex_fill_color=s.b, vertex_shape=s.b, output="block.pdf")
		h.properties[("v", "block")] = s.b.copy(value_type="string")
		h.properties[("e", "weight")] = h.ep["weight"].copy(value_type="string")
		h.save("h.graphml")
		k = nx.read_graphml("h.graphml")
		w = nx.get_edge_attributes(k, "weight")
		weights = {}
		for e in k.edges():
			weights[e] = float.fromhex(w[e])
		nx.set_edge_attributes(k, "weight", weights)
		return k


"""	
h = graph_tool.load_graph("graph.graphml")
gr = graph_tool.GraphView(h, vfilt=graph_tool.label_largest_component(h))
kcore = graph_tool.kcore_decomposition(h)
graph_tool.graph_draw(gr, pos=gr.vp["pos"], vertex_fill_color=kcore, vertex_text=kcore, output="netsci-kcore.pdf")
		
def createAndSaveBlockModel(graph):
		nx.write_graphml(graph, "graph.graphml")
		h = graph_tool.load_graph("graph.graphml")

		s = graph_tool.inference.minimize_blockmodel_dl(h)
		graph_tool.draw.graph_draw(h, vertex_text=h.vertex_properties["name"], vertex_font_size=10, output_size=(700, 700),  vertex_fill_color=s.b, vertex_shape=s.b, output="block.pdf")
		#graph_draw(g, vertex_text=v_prop, edge_text=e_prop, edge_pen_width = e_len, vertex_font_size=18, output_size=(800, 800), output="two-nodes.pdf")
		h.properties[("v", "block")] = s.b.copy(value_type="string")
		h.properties[("e", "weight")] = h.ep["weight"].copy(value_type="string")
		h.save("h.graphml")
		k = nx.read_graphml("h.graphml")
		w = nx.get_edge_attributes(k, "weight")
		weights = {}
		for e in k.edges():
			weights[e] = float.fromhex(w[e])
		nx.set_edge_attributes(k, "weight", weights)
		return k		
		"""
if __name__ == "__main__" :
	main ()