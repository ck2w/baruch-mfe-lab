import math
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from structure import *
from black_model import black_impv

class RoughBergomiModel:
    
    def __init__(self, H, rho, nu, fvc, r=0):
        
        self.H = H
        self.rho = rho
        self.nu = nu
        self.fvc = fvc
        self.r = r
        
        # binomial map:
        # 0: xi = +1, zeta = +1
        # 1: xi = +1, zeta = -1
        # 2: xi = -1, zeta = +1
        # 3: xi = -1, zeta = -1
        self.bmap = [(1,1),(1,-1),(-1,1), (-1,-1)]
        
        # The following are derived intermedate parameters. DO NOT CHANGE.
        self.rho2 = math.sqrt(1-rho*rho)
        self.alpha = H - 0.5
        self.twoH = 2*H
        self.CH = math.sqrt((2*H*math.gamma(1.5-H))/(math.gamma(0.5+H)*math.gamma(2-2*H)))
        self.twoNuCH = 2*nu*self.CH
        self.NuCH2H = nu*nu*self.CH*self.CH/H
        
        print("Hurst parameter H:", self.H)
        print("Correlation coefficient 1:", round(self.rho,3))
        print("Correlation coefficient 2:", round(self.rho2,3))
        print("Volatility of volatility:", self.nu)
        print("Forward variance curve (flat):", self.fvc(0))
        print("α: ", self.alpha)
        print("2H:", self.twoH)
        print("C_H:", round(self.CH,6))
        print("2*ν*C_H:", round(self.twoNuCH,6))
        print("ν^2*C_H^2/H:", round(self.NuCH2H,6))
        
        # Stock prices
        self.X = None
        
    def draw_fbm_tree(self, N, T):
        
        NN = total_num_nodes(N)
        dt = T/N
        rdt = math.sqrt(dt)
        
        print("Number of time steps:", N)
        print("Total number of nodes:", NN)
        print("Time horizon:", T)
        print("Time step width:", round(dt,3))
        print("Square root of time step:", round(rdt, 3))
        print("-------------------------------")
        
        # binomial map:
        # 0: xi = +1, zeta = +1
        # 1: xi = +1, zeta = -1
        # 2: xi = -1, zeta = +1
        # 3: xi = -1, zeta = -1
        bmap = [(1,1),(1,-1),(-1,1), (-1,-1)]
        
        G = np.zeros(NN)
        G[0] = 0.0

        for nn in range(1, NN):
            path = path_from_index(nn)
            n = len(path)
            GW = 0.0
            for k in range(n):
                move = path[k]
                xi, zeta = bmap[move][0], bmap[move][1]
                GW += xi*(((n-k)*dt)**self.alpha)
            GW *= rdt
            G[nn] = GW

        for n in range(N):
            nn1 = total_num_nodes(n)
            nn2 = total_num_nodes(n+1)
            print("Tree span at step", n+1, ":", round(G[nn1:nn2].min(),3), round(G[nn1:nn2].max(),3))
            
        # Draw volatility tree
        g = nx.Graph()

        # Draw nodes
        for nn in range(NN):
            path = path_from_index(nn)
            n = len(path)
            t = n*dt
            g.add_node(nn,pos=(t,G[nn]))

        # Draw edges
        for nn in range(NN):
            next_nn = 4*nn+1
            if next_nn < NN:
                g.add_edge(nn, next_nn+0, weight=0.1)
                g.add_edge(nn, next_nn+1, weight=0.1)
                g.add_edge(nn, next_nn+2, weight=0.1)
                g.add_edge(nn, next_nn+3, weight=0.1)

        pos=nx.get_node_attributes(g,'pos')
                
        plt.ylabel("$\mathcal{V}_t$", fontsize=25)
        plt.xlabel("Time", fontsize=25)
        plt.title("H = " + str(self.H), fontsize=25)
        
        nx.draw_networkx(g,pos,arrows=False, with_labels=False, node_size=50,node_color="skyblue",alpha=0.2, width=0.1)
        plt.grid('on')
        
        plt.rcParams['xtick.labelsize']=25
        plt.rcParams['ytick.labelsize']=25
        plt.rcParams["figure.figsize"] = (15,10)
        plt.show()
        print("H =", self.H)

    def create_price_tree(self, N, T, S0, EVOLVE_METHOD, DEBUG):
        
        NN = total_num_nodes(N)
        dt = T/N
        rdt = math.sqrt(dt)
        X0 = math.log(S0)
        assert EVOLVE_METHOD in {'dS', 'dLogS'}
        DLOGS = True if EVOLVE_METHOD == 'dLogS' else False
        
        print("Rough Bergomi - Number of time steps:", N)
        print("Rough Bergomi - Total number of nodes:", NN)
        print("Rough Bergomi - Time horizon:", T)
        #print("Rough Bergomi - Time step width:", round(dt,3))
        #print("Rough Bergomi - Square root of time step:", round(rdt, 3))

        # binomial map:
        # 0: xi = +1, zeta = +1
        # 1: xi = +1, zeta = -1
        # 2: xi = -1, zeta = +1
        # 3: xi = -1, zeta = -1
        bmap = [(1,1),(1,-1),(-1,1), (-1,-1)]
        
        # Allocation
        self.X = np.zeros(NN)

        # Initialization
        if DLOGS: self.X[0] = X0
        else: self.X[0] = S0

        # Fill X at time step #1
        nn1 = total_num_nodes(0)
        nn2 = total_num_nodes(1)
        for nn in range(nn1,nn2):
            ########################################################################
            path = path_from_index(nn)
            move = path[0]
            xi, zeta = bmap[move][0], bmap[move][1]
            Phi = self.fvc(0)
            PhiDt = Phi*dt
            rPhiDt = math.sqrt(PhiDt)
            if DLOGS: self.X[nn] = self.X[0]-0.5*PhiDt+rPhiDt*(self.rho*xi+self.rho2*zeta)
            else: self.X[nn] = self.X[0]*(1+rPhiDt*(self.rho*xi+self.rho2*zeta))
            ########################################################################
            if DEBUG: 
                logstr = ( "Time 1: (" + str(xi) + "," + str(zeta) + ")" + 
                           ", X[prev_n] =" + str(round(self.X[0], 3)) + 
                           ", σ×dB = " + str(round(-0.5*PhiDt+rPhiDt*(self.rho*xi+self.rho2*zeta),3)) +
                           ", X[n] = " + str(round(self.X[0]-0.5*PhiDt+rPhiDt*(self.rho*xi+self.rho2*zeta), 3)) )
                print(logstr)

        # Fill X at later time steps
        for nn in range(nn2, NN):
            ########################################################################
            path = path_from_index(nn)
            n = len(path)
            t = n*dt
            prev_t = (n-1)*dt
            t2H = prev_t**self.twoH
            ########################################################################
            GW = 0.0
            for k in range(n-1):
                move = path[k]
                xi, zeta = bmap[move][0], bmap[move][1]
                GW += xi*(((n-1-k)*dt)**self.alpha)
            GW *= rdt
            Phi = self.fvc(prev_t)*math.exp(self.twoNuCH*GW-self.NuCH2H*t2H)
            PhiDt = Phi*dt
            rPhiDt = math.sqrt(PhiDt)
            ########################################################################
            move = path[-1]
            prev_nn = (nn-move-1)//4
            xi, zeta = bmap[move][0], bmap[move][1]
            if DLOGS: self.X[nn] = self.X[prev_nn]-0.5*PhiDt+rPhiDt*(self.rho*xi+self.rho2*zeta)
            else: self.X[nn] = self.X[prev_nn]*(1+rPhiDt*(self.rho*xi+self.rho2*zeta))
            ########################################################################

            if DEBUG:
                logstr = "Time " + str(n) + ": "
                ####################################################
                prev_path = path_from_index(prev_nn)
                for k in range(len(prev_path)):
                    move = prev_path[k]
                    xi, zeta = bmap[move][0], bmap[move][1]
                    logstr += "(" + str(xi) + "," + str(zeta) + ")"
                ####################################################
                if len(prev_path)>0: logstr += " -> "
                for k in range(n):
                    move = path[k]
                    xi, zeta = bmap[move][0], bmap[move][1]
                    logstr += "(" + str(xi) + "," + str(zeta) + "),"

                logstr += ( " X[prev_n] =" + str(round(self.X[prev_nn], 3)) + 
                            ", σ×dB = " + str(round(-0.5*PhiDt+rPhiDt*(self.rho*xi+self.rho2*zeta),3)) + 
                            ", X[n] = " + str(round(self.X[prev_nn]-0.5*PhiDt+rPhiDt*(self.rho*xi+self.rho2*zeta), 3)) + 
                            ", GW =" + str(round(GW, 3)) )
                ####################################################
                print(logstr)

        if DLOGS: self.X = np.exp(self.X)
        
        # Check martingality at terminal
        EX = 0.0
        NI = start_index_at_step(N)
        for nn in range(NI, NN):
            EX += self.X[nn]
        EX /= (NN-NI)
        #print("E[X(T)] =", EX)
        
        # Print tree span
        for n in range(N):
            nn1 = total_num_nodes(n)
            nn2 = total_num_nodes(n+1)
            Xmin = self.X[nn1:nn2].min()
            Xmax = self.X[nn1:nn2].max()
            
        print("Rough Bergomi - Tree span at step", n+1, ":", round(Xmin,3), round(Xmax,3))
            
        return Xmin, Xmax
        
    def draw_price_tree(self, N, T):
        
        NN = total_num_nodes(N)
        dt = T/N
        
        g = nx.Graph()

        # Draw nodes
        for nn in range(NN):
            path = path_from_index(nn)
            n = len(path)
            t = n*dt
            g.add_node(nn,pos=(t,self.X[nn]))

        # Draw edges
        for nn in range(NN):
            next_nn = 4*nn+1
            if next_nn < NN:
                g.add_edge(nn, next_nn+0, weight=0.1)
                g.add_edge(nn, next_nn+1, weight=0.1)
                g.add_edge(nn, next_nn+2, weight=0.1)
                g.add_edge(nn, next_nn+3, weight=0.1)
            
        pos=nx.get_node_attributes(g,'pos')
        
        plt.ylabel("$\mathcal{S}_t$", fontsize=25)
        plt.xlabel("Time", fontsize=25)
        plt.title("H = " + str(self.H), fontsize=25)
        
        nx.draw_networkx(g,pos,arrows=False,with_labels=False,node_size=50,node_color="skyblue",alpha=0.5)
        plt.grid('off')
        
        plt.rcParams['xtick.labelsize']=25
        plt.rcParams['ytick.labelsize']=25
        plt.rcParams["figure.figsize"] = (15,10)
        plt.show()
        print("H =", self.H)
        
    def price_option(self, N, T, strikes, opttype, earlyexercise, r):
        
        NN = total_num_nodes(N)
        assert len(self.X) == NN
        assert opttype in {1,-1}
        
        NI = start_index_at_step(N)
        Xmin = self.X[NI:NN].min()
        Xmax = self.X[NI:NN].max()

        p = 0.25
        dt = T/N
        disc = math.exp(-r*dt)
                
        Klist = np.atleast_1d(strikes)
        Vlist = []
        EarlyExerciseBoundaryList = []
        for K in Klist:
            assert (opttype == 1 and K <= Xmax) or (opttype == -1 and K >= Xmin)
            # Initialization
            V = np.zeros(NN)

            # Terminal payoff
            NI = start_index_at_step(N)
            for nn in range(NI, NN):
                if opttype == 1: 
                    v = self.X[nn]-K
                else: 
                    v = K-self.X[nn]
                if v>0: V[nn] = v

            # Backward valuation
            for n in range(N-1,-1,-1):
                NI = start_index_at_step(n)
                NF = start_index_at_step(n+1)
                for nn in range(NI, NF):
                    next_nn = 4*nn+1
                    V[nn] = disc*p*(V[next_nn+0]+V[next_nn+1]+V[next_nn+2]+V[next_nn+3])
                    if earlyexercise:
                        if opttype == 1: v = self.X[nn]-K
                        else: v = K-self.X[nn]
                        if v > V[nn]: V[nn] = v

            Vlist.append(V[0])
            
            # Early exercise boundary
            #if earlyexercise:
            #    boundary = []
            #    for n in range(N+1):
            #        NI = start_index_at_step(n)
            #        NF = start_index_at_step(n+1)
            #        
            #        if opttype == 1: stoppingX = self.X[NI]
            #        else: stoppingX = stoppingX = self.X[NF-1]
            #        
            #        for nn in range(NI, NF):
            #            if opttype == 1: v = self.X[nn]-K
            #            else: v = K-self.X[nn]
            #            
            #            if v > V[nn]:
            #                if opttype == 1:
            #                    if self.X[nn] < stoppingX:
            #                        stoppingX = self.X[nn]
            #                else:
            #                    if self.X[nn] > stoppingX:
            #                        stoppingX = self.X[nn]
            #                        
            #        boundary.append(stoppingX)
            #        
            #    EarlyExerciseBoundaryList.append(boundary)
            
        #sigmaList = None
        #if getvol and not earlyexercise:
        #    sigmaList= []
        #    for k, v in zip(Klist, Vlist): 
        #        sigmaList.append(black_impv(k, T, self.X[0], v, 1))
                
        return Vlist, EarlyExerciseBoundaryList
    