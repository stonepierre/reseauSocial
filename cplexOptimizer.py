# -*- coding:Utf8 -*-

from __future__ import division

import scipy.sparse.linalg
import scipy.linalg
import math
import numpy.linalg
import numpy as np
import re
import sys
import networkx as nx
import matplotlib.pyplot as plt
from networkx.readwrite import json_graph
import json
import pickle
from graph_tool.all import *

import os

import constraintProgram as cp
import graphtoolMethods as gt

sys.path.append ('/opt/ibm/ILOG/CPLEX_Studio_Community1263/cplex/python/2.6/x86-64_linux')
import cplex
from cplex.exceptions import CplexError

dataFile = "newDataSet_tags.txt"
dictNodeNumberToId = { }
threshold = 0.00000005
maxexponent = 5
thresholdOfInfluence = 0.5
savedGraphs = { "scalefree": "scalefree.gml" }
pathOfFinalJsonGraph = "/home/franckeinstine/PythonCode/real_code/"
graphdataset = "graphsFiles/base_graphe.graphml"
json_graph_filename = graphdataset.replace(".graphml", ".json")
json_graph = "bds/base_ExtraitBilletterie-SansSalles.json"

def read_json_file(filename):
    with open(filename) as f:
        js_graph = json.load(f)
    return json_graph.node_link_graph(js_graph)

def main( ):
	# graph = cp.readGML(savedGraphs["scalefree"])
	d = dictionnaryFunctions()
	wf = weightFunctions()

	# graph = cp.graphGenerators (500, "scalefree")
	# graph = cp.generateRandomWeights(graph)
	# graph = nx.DiGraph(graph)

	graph = nx.read_graphml(graphdataset)

	graph = gt.createAndSaveBlockModel(graph)
	cp.writeJsonFile(graph, json_graph_filename)

	print("\n\n")
	print(nx.info(graph))
	print("\n\n")

	d.createDictNodeNumberToId(graph)
	d.checkDictionnary()
	w = wf.graphWeightsOnArcs(graph)

	#miProgram(graph, w, "influencersAdjacencyMatrixWithBlocksAndColouringFunction", json_graph_filename)

	for model in ["neighbouringInfluencersWithoutBinaryVariables", "influencersAdjacencyMatrix"]:
		miProgram(graph, w, model, json_graph_filename)

#Dictionnaire qui recuper les noeuds par leurs ID et (i) incrementé
class dictionnaryFunctions:
    def createDictNodeNumberToId( self, graph ):
        i = 0
        for x in graph.nodes():
            if not str(x):
                dictNodeNumberToId[i] = -4
            elif str(x) == '\n':
                dictNodeNumberToId[i] = -3
            elif "notfound" in str(x):
                dictNodeNumberToId[i] = -2
            elif "private" in str(x):
                dictNodeNumberToId[i] = -1                
            elif len(re.findall('\d+', str(x))) > 0:
                dictNodeNumberToId[i] = re.findall('\d+', str(x))[0]
            i += 1

    def getNodeNumberFromID( self, graph, node ):
        for key in dictNodeNumberToId:
            if dictNodeNumberToId[key] == node or dictNodeNumberToId[key] == str(node):
                return key
#Verification de combinaisons
    def checkDictionnary( self ):
        for key, value in dictNodeNumberToId.iteritems():
            if not value:
                print ("Error with key " + str(key))
                break

class ClassList_Influence(object):
	"""docstring for ClassList_Influence"""
	def __init__(self, arg):
		super(ClassList_Influence, self).__init__()
		self.arg = arg
	
	def all_neighbors(graph, node):    
	    if graph.is_directed():
	        values = itertools.chain.from_iterable([graph.predecessors_iter(node), graph.successors_iter(node)])
	    else:
	        values = graph.neighbors_iter(node)

	    return values

class weightFunctions:
    # On suppose que le graphe est dirige
    # w a la structure : {v : {(u,v) : 1}} pour tout u predecesseur de v
    # Créer un dic des poids entre u et v
    def graphIndicatorFunction( self, graph ):
        adjDict = { }
        for v in graph.nodes():
            edgeDict = { }
            for u in graph.predecessors(v):
                edgeDict[str(u) + "," + str(v)] = 1
            if edgeDict:
                adjDict[str(v)] = edgeDict
        return adjDict
    # Recuper uniquement les poid d'un arc
    def graphWeightsOnArcs( self, graph ):
        weights = nx.get_edge_attributes(graph, "weight")
        """
        for v in graph.nodes () :
            edgeDict = { }
            for u in graph.predecessors (v) :
                d = graph.get_edge_attributes(u,v,default=0)
                print(d)
                edgeDict[str (u) + "," + str (v)] = d[0]["weight"]
            if edgeDict :
                weights[str (v)] = edgeDict
                                        """
        w = { }
        for v in graph.nodes():
            eDict = { }
            for edge in weights.keys():
                if edge[1] == v:
                    eDict[str(edge[0]) + "," + str(v)] = weights[edge]
            w[str(v)] = eDict
        return w
    #La somme de toutes de toutes les voisin
    def sumOverNeighborsInWfunction( self, graph, w, u ):
        if u in w:
            edgeDict = w[u] #list de touts les voisin
            wSum = 0   #initialise la var
            for e in edgeDict:
                wSum += edgeDict[e] #la somme de a list de e
            return wSum #retur la somme
        return 0
    #F(x) recuper chaque w(u, v)  
    def getWUV( self, w, u, v ):
        if v in w:
            for edge in w[v]:
                startpoint = edge.split(",")[0] #point initial de recuperation dans la liste
                if str(u) == startpoint or u == startpoint:
                    return w[v][edge] #retour les voisins
        return 0


class modelsOfInfluencers:
	
	def addGlobalConstraint( self, graph, prob, w ):
		d = dictionnaryFunctions()
		wf = weightFunctions()
		rhs = []
		xRowNames = []
		rows = []

		j = 0
		for u in graph.nodes():
			w_u = 0
			w_uv = 0
			for e in w[str(u)]:
				w_uv += w[str(u)][e]
			w_u += w_uv
			rhs.append(w_u / 2)

		for u in graph.nodes():
			constraint = []
			variables = []
			coefficients = []
			varU = str(d.getNodeNumberFromID(graph, u))
			for v in graph.nodes():
				varV = str(d.getNodeNumberFromID(graph, str(v)))
				innerMostSum = 0
				for y in getNeighborhoodIntersectionInDigraph(graph, int(u), int(v)):
					w_uy = wf.getWUV(w, u, str(y))
					w_yv = wf.getWUV(w, str(y), str(v))
					product = w_uy * w_yv
					innerMostSum += product
				variables.append("x" + varU + "_" + varV)
				coefficients.append(innerMostSum)
			constraint.append(variables)
			constraint.append(coefficients)
			rows.append(constraint)
			xRowNames.append("C" + str(j))
			j += 1

		sense = ["G" for i in range(len(rows))]
		if not len(rows) == len(sense) == len(rhs) == len(xRowNames):
			raise ValueError("Lenght of arguments is erroneous")

		prob.linear_constraints.add(lin_expr=rows, senses=sense, rhs=rhs, names=xRowNames)
	
	def addLocalConstraintOnInfluencers( self, graph, prob, w ):
		d, wf = dictionnaryFunctions(), weightFunctions()
		rhs, xRowNames, rows = [], [], []
		j = 0

		for u in graph.nodes(): # REcuper les noeuds on place (nomme) u
			for v in graph.nodes(): # REcuper les noeuds on place (nomme) v
				w_v = wf.sumOverNeighborsInWfunction(graph, w, str(v)) # Appel à la fonction summ pour faire la somme des touts les v
				innersum = 0
				for y in getNeighborhoodIntersectionInDigraph(graph, u, v): #fonction qui fait la somme de tout les voisin de v
					w_uy = wf.getWUV(w, str(u), str(y)) #f(x) qui recuper les u et le v à partir du poit 0
					w_yv = wf.getWUV(w, str(y), str(v)) #f(x) qui recuper les v et le v à partir du poit 0
					prod = w_uy * w_yv #le produi de w(u, y) et w(y, v)
					innersum += prod #la somme de tout
				if w_v > 0: #s'il existe une somme de v > 0
					rhs.append(innersum - 1 / 2 * w_v + threshold) #On addition ds la list rhs la somme -1/2 * 0.00000005 (seuil)
				else:
					rhs.append(-1 / 2 + threshold) #sinon on addition ds la list rhs la somme -1/2 + 0.00000005 (seuil)

		maxRhs = max(rhs) #on recup la val max de rhs
		for i in range(len(rhs)):
			rhs[i] = rhs[i] / maxRhs

		for u in graph.nodes(): #pour les u
			varU = str(d.getNodeNumberFromID(graph, u)) #on recup les id de chaque noeud - donc (u)

			for v in graph.nodes():
				constraint, variables, coefficients = [], [], []

				varV = str(d.getNodeNumberFromID(graph, v))
				variables.append("x" + varU + "_" + varV)
				coefficients.append(1)
				constraint.append(variables)
				constraint.append(coefficients)
				rows.append(constraint)
				xRowNames.append("C" + str(j))
				j += 1

		sense = ["G" for i in range(len(rows))]
		prob.linear_constraints.add(lin_expr=rows, senses=sense, rhs=rhs, names=xRowNames)

	def addBinaryInfluencerVariables( self, graph, prob, w ):
		d = dictionnaryFunctions()
		objFct, xUpBnds, xColNames, xVarType = [], [], [], []

		for u in graph.nodes():
			for v in graph.nodes():
				objFct.append(1)
				xUpBnds.append(1)
				varU = str(d.getNodeNumberFromID(graph, u))
				varV = str(d.getNodeNumberFromID(graph, v))
				if varU not in [-4, -3, -2] and varV not in [-4, -3, -2]:
					xColNames.append("x" + varU + "_" + varV)
				xVarType.append("I")
		prob.variables.add(obj=objFct, ub=xUpBnds, names=xColNames, types=xVarType)

	def addNeighboursModelConstraints( self, graph, prob, w ):
		wf = weightFunctions()
		d = dictionnaryFunctions()
		rhs = []
		xRowNames = []
		rows = []

		j = 0
		for v in graph.nodes():
			wSum = wf.sumOverNeighborsInWfunction(graph, w, str(v))
			rhs.append(wSum / 2)

			constraint = []
			variables = []
			coefficients = []

			for u in graph.predecessors(v):
				varU = str(d.getNodeNumberFromID(graph, u))
				varV = str(d.getNodeNumberFromID(graph, v))

				if varU not in [-4, -3, -2] and varV not in [-4, -3, -2]:
					variables.append("x" + varU + "_" + varV)
					coefficients.append(wf.getWUV(w, str(u), str(v)))

			constraint.append(variables)
			constraint.append(coefficients)

			rows.append(constraint)
			xRowNames.append("C" + str(j))
			j += 1

		sense = ["G" for i in range(len(rows))]
		if not len(rows) == len(sense) == len(rhs) == len(xRowNames):
			raise ValueError("Lenght of arguments is erroneous")

		prob.linear_constraints.add(lin_expr=rows, senses=sense, rhs=rhs, names=xRowNames)

	def addNeighbouringConstraintsWithUnaryVariables( self, graph, prob, w ):
		wf, d = weightFunctions(), dictionnaryFunctions()
		rhs, xRowNames, rows = [], [], []

		j = 0
		for v in graph.nodes():
			wSum = wf.sumOverNeighborsInWfunction(graph, w, str(v))
			rhs.append(thresholdOfInfluence * wSum)

			constraint, variables, coefficients = [], [], []

			for u in graph.predecessors(v):
				varU = str(d.getNodeNumberFromID(graph, u))

				if varU not in [-4, -3, -2]:
					variables.append("x" + varU)
					coefficients.append(wf.getWUV(w, str(u), str(v)))

			constraint.append(variables)
			constraint.append(coefficients)

			rows.append(constraint)
			xRowNames.append("C" + str(j))
			j += 1

		sense = ["G" for i in range(len(rows))]
		if not len(rows) == len(sense) == len(rhs) == len(xRowNames):
			raise ValueError("Lenght of arguments is erroneous")

		prob.linear_constraints.add(lin_expr=rows, senses=sense, rhs=rhs, names=xRowNames)

	def addUnaryInfluencerVariables( self, graph, prob, w ):
		d = dictionnaryFunctions()
		objFunction, xUpperBounds, xColNames, xVarType = [], [], [], []

		for u in graph.nodes():
			objFunction.append(1)
			xUpperBounds.append(1)
			varName = str(d.getNodeNumberFromID(graph, u))
			if varName not in [-4, -3, -2]:
				xColNames.append('x' + varName)
			xVarType.append("I")
		prob.variables.add(obj=objFunction, ub=xUpperBounds, names=xColNames, types=xVarType)

	def constraintWithAdjacencyMatrix( self, graph, prob ):
		d, wf, rhs, xrnames, xr = dictionnaryFunctions(), weightFunctions(), [], [], []
		E = adjMatrixEponential(graph)
		dim = { }  # Les lignes et colonnes de E sont ordonnees selon graph.nodes()
		dim_i = 0
		for v in graph.nodes():
			dim[v] = dim_i
			dim_i += 1

		j = 0
		for v in graph.nodes():
			E_v = E.sum(axis=0).item(dim[v])
			rhs.append(thresholdOfInfluence * E_v)

			constraint, variables, coefficients = [], [], []

			for u in graph.nodes():
				varU = str(d.getNodeNumberFromID(graph, u))

				if varU not in [-4, -3, -2]:
					variables.append("x" + varU)
					coefficients.append(E.item(dim[u], dim[v]))

			if variables:
				constraint.append(variables)
				constraint.append(coefficients)

				xr.append(constraint)
				xrnames.append("C" + str(j))
				j += 1

		sense = ["G" for i in range(len(xr))]
		if not len(xr) == len(sense) == len(rhs) == len(xrnames):
			raise ValueError("Lenght of arguments is invalid")

		prob.linear_constraints.add(lin_expr=xr, senses=sense, rhs=rhs, names=xrnames)

	def constraintWithAdjacencyMatrixAndCommunities( self, graph, prob, colouring, useblockmodel ):
		d, wf, rhs, xrnames, xr = dictionnaryFunctions(), weightFunctions(), [], [], []
		E = adjMatrixEponential(graph)

		block = nx.get_node_attributes(graph, "block")
		blocks = block.values()
		verticesPerBlock = { }
		for b in blocks:
			vertices = []
			for k, v in block.items():
				if v == b:
					vertices.append(k)
			verticesPerBlock[b] = vertices

		colours = { }
		if colouring:
			colours = colournodes(graph, json_graph_filename)

		j = 0
		for v in graph.nodes():
			constraint, variables, coefficients = [], [], []

			if useblockmodel:
				for u in verticesPerBlock[block[v]]:
					self.innersumonneighbours(E, coefficients, colouring, colours, d, graph, rhs, u, v, variables)
			else:
				for u in graph.nodes():
					self.innersumonneighbours(E, coefficients, colouring, colours, d, graph, rhs, u, v, variables)

			if variables:
				constraint.append(variables)
				constraint.append(coefficients)

				xr.append(constraint)
				xrnames.append("C" + str(j))
				j += 1

		sense = ["G" for i in range(len(xr))]
		if not len(xr) == len(sense) == len(rhs) == len(xrnames):
			raise ValueError("Lenght of arguments is invalid")

	def innersumonneighbours( self, E, coefficients, colouring, colours, d, graph, rhs, u, v, variables ):
		E_uv = E.item(d.getNodeNumberFromID(graph, u), d.getNodeNumberFromID(graph, v))
		if colouring:
			rhs.append(thresholdOfInfluence * E_uv * colours[u])

		else:
			rhs.append(thresholdOfInfluence * E_uv)
		varU = str(d.getNodeNumberFromID(graph, u))
		if varU not in [-4, -3, -2]:
			if colouring:
				coefficients.append(E_uv * colours[u])
			else:
				coefficients.append(E_uv)
			if coefficients > 0:
				variables.append("x" + varU)

def createModel( prob, graph, w, constraintName ):
	prob.objective.set_sense(prob.objective.sense.minimize)
	m = modelsOfInfluencers()

	try:
		if constraintName == "localInfluencers":
			m.addBinaryInfluencerVariables(graph, prob, w)
			m.addLocalConstraintOnInfluencers(graph, prob, w)
		elif constraintName == "neighbouringInfluencers":
			m.addBinaryInfluencerVariables(graph, prob, w)
			m.addNeighboursModelConstraints(graph, prob, w)
		elif constraintName == "neighbouringInfluencersWithoutBinaryVariables":
			m.addUnaryInfluencerVariables(graph, prob, w)
			m.addNeighbouringConstraintsWithUnaryVariables(graph, prob, w)
		elif constraintName == "influencersAdjacencyMatrix":
			m.addUnaryInfluencerVariables(graph, prob, w)
			m.constraintWithAdjacencyMatrix(graph, prob)
		elif constraintName == "influencersAdjacencyMatrixWithBlocks":
			m.addUnaryInfluencerVariables(graph, prob, w)
			m.constraintWithAdjacencyMatrixAndCommunities(graph, prob, False, False)
		elif constraintName == "influencersAdjacencyMatrixWithBlocksAndColouringFunction":
			m.addUnaryInfluencerVariables(graph, prob, w)
			m.constraintWithAdjacencyMatrixAndCommunities(graph, prob, True, False)
	except ValueError as e:
		print (e.args)


def miProgram( graph, w, constraintName, jsfilename ):
	try:
		my_prob = cplex.Cplex()
		createModel(my_prob, graph, w, constraintName)
		my_prob.solve()

	except CplexError as exc:
		print (exc)
		return

	print ()
	print ("Solution status = ", my_prob.solution.get_status(), ":")
	print (my_prob.solution.status[my_prob.solution.get_status()])
	print ("Solution value  = ", my_prob.solution.get_objective_value())
	print ("\n Individual solutions : \n")

	writeValuesOfEachVariable(constraintName, graph, my_prob, w)
	my_prob.write("miProgram.lp")
	infDict = extractInfluencersFromSolution(constraintName, graph, my_prob)
	updateJsonWithInfluencers(infDict, graph, jsfilename.replace(".json", "_" + constraintName + ".json"))


def updateJsonWithInfluencers( inflDict, graph, newjsonfile ):
	nx.set_node_attributes(graph, "influence", inflDict)
	if inflDict:
		cp.writeJsonFile(graph, pathOfFinalJsonGraph + newjsonfile)


def extractInfluencersFromSolution( constraintName, graph, my_prob ):
	nodeAndInfluenceSolution = { }
	if constraintName == "neighbouringInfluencersWithoutBinaryVariables" or constraintName == "influencersAdjacencyMatrix" or constraintName == "influencersAdjacencyMatrixWithBlocks" or constraintName == "influencersAdjacencyMatrixWithBlocksAndColouringFunction":
		names = my_prob.variables.get_names()
		influencers = []
		for var in names:
			x = my_prob.solution.get_values(var)
			u = var[1:]
			nodeAndInfluenceSolution[str(dictNodeNumberToId[int(u)])] = x
	return nodeAndInfluenceSolution


def writeValuesOfEachVariable( constraintName, graph, my_prob, w ):
	allZeroes = True
	if constraintName == "localInfluencers" or constraintName == "neighbouringInfluencers":
		names = my_prob.variables.get_names()
		influencers = []
		for var in names:
			x = my_prob.solution.get_values(var)
			if x > 0:
				print ("Variable " + var + " value = %10f" % x)
				u = var[1:].split("_")[0]
				v = var[1:].split("_")[1]
				if not u in influencers and u != v:
					influencers.append(u)
				allZeroes = False
		print ("\n")
		for inf in influencers:
			print (str(dictNodeNumberToId[int(inf)]) + " est un influenceur")

	elif constraintName == "neighbouringInfluencersWithoutBinaryVariables" or constraintName == "influencersAdjacencyMatrix" or constraintName == "influencersAdjacencyMatrixWithBlocks" or constraintName == "influencersAdjacencyMatrixWithBlocksAndColouringFunction":
		names = my_prob.variables.get_names()
		influencers = []
		for var in names:
			x = my_prob.solution.get_values(var)
			if x > 0:
				print ("Variable " + var + " value = %10f" % x)
				u = var[1:]
				if not u in influencers:
					influencers.append(int(u))
				allZeroes = False
		influencers.sort()
		print ("\n")
		for inf in influencers:
			print (str(dictNodeNumberToId[int(inf)]) + " est un influenceur")

	else:
		numcols = my_prob.variables.get_num()
		x = my_prob.solution.get_values()
		for j in range(numcols):
			print ("Column %d:  Value = %10f" % (j, x[j]))

	if (allZeroes):
		print ("All variables are zeroes")

# Recuperation des interctions de entre u et v
def getNeighborhoodIntersectionInDigraph( graph, u, v ):
	n_u = graph.neighbors(u) #Recuperation des voisin de u
	n_u += [u] # somme de tous les voisn de u
	n_v = graph.predecessors(v) ##Recuperation de predecesseurs de v
	if set(n_u) & set(n_v): #Verifie les variables
		return [node for node in n_u if node in n_v] #retour la somme de tout les vois de chaque v 
	return []


def saveJSON( graph, exportFilename ):
	data = json_graph.node_link_data(graph)
	with open(exportFilename, 'w') as out:
		json.dump(data, out, indent=4)


def readJsonFile( filename ):
	d = json.load(open(filename))
	return json_graph.node_link_graph(d)


def updateJson( dic, valname, filename ):
	h = readJsonFile(filename)
	nx.set_node_attributes(h, valname, dic)
	saveJSON(h, filename)


def katzcentrality( digraph, alpha, beta ):
	n = digraph.order()
	A = nx.adjacency_matrix(digraph, weight="weight")
	I = np.identity(n)
	b = np.ones(n)
	inverse = np.linalg.solve(I - alpha * A, np.transpose(b))
	ktz = { }
	sol = beta * inverse
	for u in digraph.nodes():
		ktz[u] = sol[u]
	return ktz


def adjMatrixEponential( graph ):
	A = nx.adjacency_matrix(graph, weight="weight")
	alpha, alpha_v = scipy.sparse.linalg.eigsh(A, 1, which="LM")
	for i in range(2, maxexponent):	#le range(2, maxexponent) limite le comptage que de 2 a la valeur exponentiel déjà declaré plus haut. Donc, on aura 2,3,4 au lieu de 1,2,3,4
		A_i = math.pow(alpha, -i) * numpy.linalg.matrix_power(A, i)
		A += A_i
	return A


def colournodes( graph, jsfile ):
	types = nx.get_node_attributes(graph, "type")
	clients = [node for node in graph.nodes() if types[node] == "Client"]
	others = [node for node in graph.nodes() if node not in clients]
	dic = { }
	for c in clients:
		dic[c] = 1
	for o in others:
		dic[o] = 0
	updateJson(dic, "couleur", jsfile)
	return dic


if __name__ == "__main__":
	main()
