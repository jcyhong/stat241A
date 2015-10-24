# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 00:54:42 2015

@author: pitcany
"""
import numpy as np
import pandas as pd
import networkx as nx
import itertools

obs = pd.read_table('hw4data.data',names=[1,2,3,4,5,6,7])
obs.index=range(1,501,1)

G=nx.Graph()
H=nx.Graph()

G.add_edges_from([(1,5),(1,6),(2,5),(3,4),(3,6),(4,5),(4,7),(6,7)])
H.add_edges_from([(1,2),(1,3),(1,4),(1,7),(2,4),(2,6),(3,4),(3,5),(3,7),
                  (4,6),(5,7)])
#create edges for each graph

J=nx.complete_graph(7)

def shift_one(x):
    return x+1
    
nx.relabel_nodes(J,shift_one,copy=False)

#maximal cliques for G of size 2
#clique potential tables

#g={(i,j):np.ones((2,2)) for i,j in G.edges()}

G15 = np.ones((2,2))
G16 = np.ones((2,2))
G25 = np.ones((2,2))
G34 = np.ones((2,2))
G36 = np.ones((2,2))
G45 = np.ones((2,2))
G47 = np.ones((2,2))
G67 = np.ones((2,2))

#empirical probability distributions
#emp={(i,j):np.ones((2,2)) for i,j in G.edges()}

PG15 = np.ones((2,2))
PG16 = np.ones((2,2))
PG25 = np.ones((2,2))
PG34 = np.ones((2,2))
PG36 = np.ones((2,2))
PG45 = np.ones((2,2))
PG47 = np.ones((2,2))
PG67 = np.ones((2,2))

def fill_2Dtable(X,a,b):
    a = int(a)
    b = int(b)
    for i in range(2):
        for j in range(2):
            X[(i,j)] = len(obs[(obs[a]==i) & (obs[b]==j)])
    X/=500

#fill_2Dtable(e) for e in G.edges()
fill_2Dtable(PG15,1,5)
fill_2Dtable(PG16,1,6)
fill_2Dtable(PG25,2,5)
fill_2Dtable(PG34,3,4)
fill_2Dtable(PG36,3,6)
fill_2Dtable(PG45,4,5)
fill_2Dtable(PG47,4,7)
fill_2Dtable(PG67,6,7)

GJoint=np.ones((2,2,2,2,2,2,2))

#create joint probability table for model marginal probability calculations later
def calculateGJoint():
    binvals=list(itertools.product([0,1],repeat=7))
    for t in binvals:
        GJoint[t]=G15[t[0],t[4]]*G16[t[0],t[5]]*G25[t[1],t[4]]*G34[t[2],t[3]]*G36[t[2],t[5]]*G45[t[3],t[4]]*G47[t[3],t[6]]*G67[t[5],t[6]]

calculateGJoint()
GJoint /= GJoint.sum()

#calculates marginal probabilities fixing an edge
def MarginalProb(X,i,j):
    for i in range(2):
        for j in range(2):
            X[(i,j)]=GJoint[]

#these maximal cliques of H are all of size 3
H124 = np.ones((2,2,2))
H134 = np.ones((2,2,2))
H137 = np.ones((2,2,2))
H357 = np.ones((2,2,2))
H246 = np.ones((2,2,2))

#empirical probability tables
PH124 = np.ones((2,2,2))
PH134 = np.ones((2,2,2))
PH137 = np.ones((2,2,2))
PH357 = np.ones((2,2,2))
PH246 = np.ones((2,2,2))

def fill_3Dtable(X,a,b,c):
    a = int(a)
    b = int(b)
    c = int(c)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                X[(i,j,k)] = len(obs[(obs[a]==i) & (obs[b]==j) & (obs[c]==k)])
    X/=500

fill_3Dtable(PH124,1,2,4)
fill_3Dtable(PH134,1,3,4)
fill_3Dtable(PH137,1,3,7)
fill_3Dtable(PH357,3,5,7)
fill_3Dtable(PH246,2,4,6)

HJoint = np.ones((2,2,2,2,2,2,2))
def calculateHJoint():
    binvals=list(itertools.product([0,1],repeat=7))
    for t in binvals:
        HJoint[t]=H124[t[0],t[1],t[4]]*H134[t[0],t[2],t[3]]*H137[t[0],t[2],t[6]]*H357[t[2],t[4],t[6]]*H246[t[1],t[3],t[5]]

calculateHJoint()
HJoint /= HJoint.sum()

#J
#complete graph...maximal clique is the graph itself
J7D = np.ones((2,2,2,2,2,2,2))
#empirical probability table
PJ7D = np.ones((2,2,2,2,2,2,2))

def fill_7Dtable(X,a,b,c,d,e,f,g):
    a=int(a)
    b=int(b)
    c=int(c)
    d=int(d)
    e=int(e)
    f=int(f)
    g=int(g)
    binvals=list(itertools.product([0,1],repeat=7))
    for t in binvals:
        X[t]=len(obs[(obs[a]==t[0]) & (obs[b]==t[1]) & (obs[c]==t[2]) & 
        (obs[d]==t[3]) & (obs[e]==t[4]) & (obs[f]==t[5]) & (obs[g]==t[6])])
    X/=500

fill_7Dtable(PJ7D,1,2,3,4,5,6,7)
