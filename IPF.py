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

def getslice(X, axes):
    idx00 = [0 if i in axes else slice(None) for i in range(X.ndim)]    
    idx01 = [0 if i == axes[0] else 1 if i == axes[1] else slice(None) for i in range(X.ndim)]
    idx10 = [1 if i == axes[0] else 0 if i == axes[1] else slice(None) for i in range(X.ndim)]
    idx11 = [1 if i in axes else slice(None) for i in range(X.ndim)]
    d = {(0,0):X[idx00],(0,1):X[idx01],(1,0):X[idx10],(1,1):X[idx11]}
    return d

#calculates marginal probabilities fixing an edge
def MarginalProb(X,e):
    Q=np.ones((2,2))
    for i in range(2):
        for j in range(2):
            Q[(i,j)]=getslice(X,(e[0]-1,e[1]-1))[(i,j)].sum()
    return Q

def normalize(X):
    X /= X.sum()
#maximal cliques for G of size 2
#clique potential tables

g={(i,j):np.ones((2,2)) for i,j in G.edges()}


#empirical probability distributions
g_emp={(i,j):np.ones((2,2)) for i,j in G.edges()}


def fill_2Dtable_g(a,b):
    a = int(a)
    b = int(b)
    for i in range(2):
        for j in range(2):
            g_emp[(a,b)][(i,j)] = len(obs[(obs[a]==i) & (obs[b]==j)])
    normalize(g_emp[(a,b)])

#fill_2Dtable(e) for e in G.edges()
for e in G.edges():
    fill_2Dtable_g(e[0],e[1])

GJoint=np.ones((2,2,2,2,2,2,2))

#create joint probability table for model marginal probability calculations later
    
def calculateGJoint():
    binvals=list(itertools.product([0,1],repeat=7))
    for t in binvals:
        GJoint[t]=np.array([g[e][t[e[0]-1],t[e[1]-1]] for e in G.edges()]).prod()
    return(normalize(GJoint))

calculateGJoint()

h={(i,j):np.ones((2,2)) for i,j in H.edges()}

h_emp={(i,j):np.ones((2,2)) for i,j in H.edges()}

def fill_2Dtable_h(a,b):
    a = int(a)
    b = int(b)
    for i in range(2):
        for j in range(2):
            h_emp[(a,b)][(i,j)] = len(obs[(obs[a]==i) & (obs[b]==j)])
    normalize(h_emp[(a,b)])

#fill empirical probabilities for each edge in H.edges()
for e in H.edges():
    fill_2Dtable_h(e[0],e[1])


HJoint = np.ones((2,2,2,2,2,2,2))
def calculateHJoint():
    binvals=list(itertools.product([0,1],repeat=7))
    for t in binvals:
        HJoint[t]=np.array([h[e][t[e[0]-1],t[e[1]-1]] for e in H.edges()]).prod()
    return normalize(HJoint)
    
calculateHJoint()

def IPF_G_step():
    for e in G.edges():
        #update potentials
        if np.allclose(g_emp[e]/MarginalProb(GJoint,e), np.ones((2,2)), .01) == False:
            g[e] *= g_emp[e]/MarginalProb(GJoint,e)
            #update joint probability
            calculateGJoint()

def IPF_G():
    #create a list telling us if the ratios for empirical to model are close to 1 for each edge
    #ratios={e:g_emp[e]/MarginalProb(GJoint,e) for e in G.edges()}
    ratios=[g_emp[e]/MarginalProb(GJoint,e) for e in G.edges()]
    while(all(np.allclose(ratio,np.ones((2,2))) for ratio in ratios) == False):
        IPF_G_step()
    

def IPF_H_step():
    for e in H.edges():
        #update potentials
        if np.allclose(h_emp[e]/MarginalProb(HJoint,e), np.ones((2,2)), .01) == False:
            h[e] *= h_emp[e]/MarginalProb(HJoint,e)
            #update joint probability
            calculateHJoint()

def IPF_H():
    #create a list telling us if the ratios for empirical to model are close to 1 for each edge
    #ratios={e:g_emp[e]/MarginalProb(GJoint,e) for e in G.edges()}
    ratios=[h_emp[e]/MarginalProb(HJoint,e) for e in H.edges()]
    while(all(np.allclose(ratio,np.ones((2,2))) for ratio in ratios) == False):
        IPF_H_step()


##J
##complete graph...maximal clique is the graph itself
#J7D = np.ones((2,2,2,2,2,2,2))
##empirical probability table
#PJ7D = np.ones((2,2,2,2,2,2,2))
#
#def fill_7Dtable(X,a,b,c,d,e,f,g):
#    a=int(a)
#    b=int(b)
#    c=int(c)
#    d=int(d)
#    e=int(e)
#    f=int(f)
#    g=int(g)
#    binvals=list(itertools.product([0,1],repeat=7))
#    for t in binvals:
#        X[t]=len(obs[(obs[a]==t[0]) & (obs[b]==t[1]) & (obs[c]==t[2]) & 
#        (obs[d]==t[3]) & (obs[e]==t[4]) & (obs[f]==t[5]) & (obs[g]==t[6])])
#    X/=500
#
#fill_7Dtable(PJ7D,1,2,3,4,5,6,7)

