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

g={(i,j):np.ones((2,2)) for i,j in G.edges()}

#G15 = np.ones((2,2))
#G16 = np.ones((2,2))
#G25 = np.ones((2,2))
#G34 = np.ones((2,2))
#G36 = np.ones((2,2))
#G45 = np.ones((2,2))
#G47 = np.ones((2,2))
#G67 = np.ones((2,2))

#empirical probability distributions
g_emp={(i,j):np.ones((2,2)) for i,j in G.edges()}

#PG15 = np.ones((2,2))
#PG16 = np.ones((2,2))
#PG25 = np.ones((2,2))
#PG34 = np.ones((2,2))
#PG36 = np.ones((2,2))
#PG45 = np.ones((2,2))
#PG47 = np.ones((2,2))
#PG67 = np.ones((2,2))

def fill_2Dtable_g(a,b):
    a = int(a)
    b = int(b)
    for i in range(2):
        for j in range(2):
            g_emp[(a,b)][(i,j)] = len(obs[(obs[a]==i) & (obs[b]==j)])
    g_emp[(a,b)]/=500

#fill_2Dtable(e) for e in G.edges()
for e in G.edges():
    fill_2Dtable_g(e[0],e[1])
#fill_2Dtable(PG15,1,5)
#fill_2Dtable(PG16,1,6)
#fill_2Dtable(PG25,2,5)
#fill_2Dtable(PG34,3,4)
#fill_2Dtable(PG36,3,6)
#fill_2Dtable(PG45,4,5)
#fill_2Dtable(PG47,4,7)
#fill_2Dtable(PG67,6,7)

GJoint=np.ones((2,2,2,2,2,2,2))

#create joint probability table for model marginal probability calculations later
def calculateGJoint():
    binvals=list(itertools.product([0,1],repeat=7))
    for t in binvals:
        GJoint[t]=np.array([g[e][t[e[0]-1],t[e[1]-1]] for e in G.edges()]).prod()

calculateGJoint()
GJoint /= GJoint.sum()


h={(i,j):np.ones((2,2)) for i,j in H.edges()}

h_emp={(i,j):np.ones((2,2)) for i,j in H.edges()}

def fill_2Dtable_h(a,b):
    a = int(a)
    b = int(b)
    for i in range(2):
        for j in range(2):
            h_emp[(a,b)][(i,j)] = len(obs[(obs[a]==i) & (obs[b]==j)])
    h_emp[(a,b)]/=500

#fill_2Dtable(e) for e in G.edges()
for e in H.edges():
    fill_2Dtable_h(e[0],e[1])
    
##these maximal cliques of H are all of size 3
#H124 = np.ones((2,2,2))
#H134 = np.ones((2,2,2))
#H137 = np.ones((2,2,2))
#H357 = np.ones((2,2,2))
#H246 = np.ones((2,2,2))
#
##empirical probability tables
#PH124 = np.ones((2,2,2))
#PH134 = np.ones((2,2,2))
#PH137 = np.ones((2,2,2))
#PH357 = np.ones((2,2,2))
#PH246 = np.ones((2,2,2))

#def fill_3Dtable(X,a,b,c):
#    a = int(a)
#    b = int(b)
#    c = int(c)
#    for i in range(2):
#        for j in range(2):
#            for k in range(2):
#                X[(i,j,k)] = len(obs[(obs[a]==i) & (obs[b]==j) & (obs[c]==k)])
#    X/=500
#
#fill_3Dtable(PH124,1,2,4)
#fill_3Dtable(PH134,1,3,4)
#fill_3Dtable(PH137,1,3,7)
#fill_3Dtable(PH357,3,5,7)
#fill_3Dtable(PH246,2,4,6)

HJoint = np.ones((2,2,2,2,2,2,2))
def calculateHJoint():
    binvals=list(itertools.product([0,1],repeat=7))
    for t in binvals:
        HJoint[t]=np.array([h[e][t[e[0]-1],t[e[1]-1]] for e in H.edges()]).prod()

calculateHJoint()
HJoint /= HJoint.sum()

def getslice(X, axes):
    idx00 = [0 if i in axes else slice(None) for i in range(X.ndim)]    
    idx01 = [0 if i == axes[0] else 1 if i == axes[1] else slice(None) for i in range(X.ndim)]
    idx10 = [1 if i == axes[0] else 0 if i == axes[1] else slice(None) for i in range(X.ndim)]
    idx11 = [1 if i in axes else slice(None) for i in range(X.ndim)]
    d = {(0,0):X[idx00],(0,1):X[idx01],(1,0):X[idx10],(1,1):X[idx11]}
    return d

#calculates marginal probabilities fixing an edge
def MarginalProb(X,a,b):
    Q=np.ones((2,2))
    for i in range(2):
        for j in range(2):
            Q[(i,j)]=getslice(X,(a,b))[(i,j)].sum()
    return Q

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
