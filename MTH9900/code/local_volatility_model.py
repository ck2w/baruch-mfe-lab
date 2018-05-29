import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

class LocalVolatilityModel:
    
    def __init__(self, r=0):
        
        self.r = r
        self.S = None # Implied stock tree
        self.p = None # Transition probability tree
        self.lam = None # Arrow-Debreu price tree
        
    def create_tree(self, N, T, S0, rBergomi, rBergomiN, DEBUG):
        
        #############################################
        # The storage is implemented as an array of #
        # size M=(N+1)(N+2)/2 where N is the number #
        # of time steps. The i'th node at the k'th  #
        # step is indexed at k(k+1)/2 + i, i<=k.    #
        #                                           #
        #                      15                   #
        #                   10                      #
        #                 6    16                   #
        #               3   11                      #
        #             1   7    17                   #
        #           0   4   12                      #
        #             2   8    18                   #
        #               5   13                      #
        #                 9    19                   #
        #                   14                      #
        #                      20                   #
        #                                           #
        #############################################
        
        NN = (N+1)*(N+2)//2;
        dt = T/N
        rdt = math.sqrt(dt)
        
        print("Number of time steps:", N)
        print("Total number of nodes:", NN)
        print("Time horizon:", T)
        print("Time step width:", round(dt,3))
        print("Square root of time step:", round(rdt, 3))
        print("=======================================================================")
        
        # Allocation
        self.S = np.zeros(NN)
        self.p = np.zeros(NN)
        self.lam = np.zeros(NN)
                
        # Initialization
        self.S[0] = S0
        self.lam[0] = 1.0
        
        # Compounding and discounting factor
        comp = math.exp(self.r*dt)
        disc = 1.0/comp
        
        # Evolve forward in time
        for k in range(1,N+1):
            
            #if k>4: break
            
            # Current time
            t = k*dt
            
            # Index root of this time step
            index_root = (k*(k+1))//2
            prev_index_root = ((k-1)*k)//2
            pivot = k//2
            
            if DEBUG: print("Current time step:", k, 
                            ", Current time:", t, 
                            ", pivot:", pivot, 
                            ", pivot index:", index_root + pivot)
            
            # Create rough Bergomi tree
            Smin, Smax = rBergomi.create_price_tree(rBergomiN, t, S0, EVOLVE_METHOD="dLogS", DEBUG=False)
            
            # Odd number of nodes at this time step
            if k%2 == 0:
                pivot_index = index_root + pivot # In this case, k is an even number, k/2 does not have remainder.
                self.S[pivot_index] = S0
                
                # Use call options
                if DEBUG: print("Decreasing from", pivot-1, "to 0 to use call options.")
                for i in range(pivot, 0, -1):
                    curr_index = index_root + i
                    prev_index = prev_index_root + (i-1)
                    if DEBUG: print("hey~", i-1, ", curr_index:", curr_index, ", prev_index:", prev_index)
                    si = self.S[prev_index]
                    Si = self.S[curr_index]
                    Fi = comp*si
                    Call = rBergomi.price_option(rBergomiN, t, si, opttype=1, earlyexercise=False, r=self.r)[0][0]
                    lam = self.lam[prev_index]
                    Sigma = 0
                    for j in range(i-1):
                        if DEBUG: print("Calculating Sigma, adding ", prev_index_root+j)
                        Sigma += self.lam[prev_index_root+j]*(comp*self.S[prev_index_root+j]-si)
                        
                    if DEBUG: print("Si =", Si, 
                                    ", si =", si, 
                                    ", Fi =", Fi, 
                                    ", Call =", Call, 
                                    ", lambda:", lam, 
                                    ", Sigma =", Sigma)
                        
                    if DEBUG: print("Double check, numerator:", Si*(comp*Call-Sigma)-lam*si*(Fi-Si), 
                                    ", denominator:", (comp*Call-Sigma-lam*(Fi-Si)))
                        
                    if DEBUG: print("Fill", curr_index-1, 
                                    "using", curr_index, 
                                    "and", prev_index)
                        
                    self.S[curr_index-1] = (Si*(comp*Call-Sigma)-lam*si*(Fi-Si))/(comp*Call-Sigma-lam*(Fi-Si))
                    
                    # Mannual overrides
                    if self.S[curr_index-1] <= comp*si or ( (i-2) >= 0 and self.S[curr_index-1] >= comp*self.S[prev_index-1] ):
                        print("Mannual override calls (A.1):", curr_index-1,curr_index, prev_index, prev_index+1)
                        self.S[curr_index-1] = self.S[curr_index]*self.S[prev_index]/self.S[prev_index+1]
                    if self.S[curr_index-1] > Smax: 
                        print("Mannual override calls (A.2):", curr_index-1, Smax)
                        self.S[curr_index-1] = Smax
                    
                # Use put options
                if DEBUG: print("Increasing from", pivot+1, "to", k, "to use put options.")
                    
                for i in range(pivot, k, 1):
                    curr_index = index_root + i
                    prev_index = prev_index_root + i
                    
                    if DEBUG: print("hey~", i+1, 
                                    ", curr_index:", curr_index, 
                                    ", prev_index:", prev_index)
                    
                    si = self.S[prev_index]
                    Si1 = self.S[curr_index]
                    Fi = comp*si
                    Put = rBergomi.price_option(rBergomiN, t, si, opttype=-1, earlyexercise=False, r=self.r)[0][0]
                    lam = self.lam[prev_index]
                    Sigma = 0
                    for j in range(i+1, k):
                        if DEBUG: print("Calculating Sigma, adding ", prev_index_root+j)
                        Sigma += self.lam[prev_index_root+j]*(si-comp*self.S[prev_index_root+j])
                    
                    if DEBUG: print("S(i+1) =", Si1, 
                                    ", si =", si, 
                                    ", Fi =", Fi, 
                                    ", Put =", Put, 
                                    ", lambda:", lam, 
                                    ", Sigma =", Sigma)
                        
                    if DEBUG: print("Double check, numerator:", (Si1*(comp*Put-Sigma)+lam*si*(Fi-Si1)), 
                                    ", denominator:", (comp*Put-Sigma+lam*(Fi-Si1)))
                        
                    if DEBUG: print("Fill", curr_index+1, "using", curr_index, "and", prev_index)
                        
                    self.S[curr_index+1] = (Si1*(comp*Put-Sigma)+lam*si*(Fi-Si1))/(comp*Put-Sigma+lam*(Fi-Si1))
                    
                    # Mannual overrides
                    if ( self.S[curr_index+1] >= comp*si or 
                         ( (i+1) <= (k-1) and self.S[curr_index+1] <= comp*self.S[prev_index+1] ) ):
                        print("Mannual override puts (A.1):", curr_index+1, curr_index, prev_index, prev_index-1)
                        self.S[curr_index+1] = self.S[curr_index]*self.S[prev_index]/self.S[prev_index-1]
                    if self.S[curr_index+1] < Smin: 
                        print("Mannual override puts (A.2):", curr_index+1, Smin)
                        self.S[curr_index+1] = Smin
                    
            # Even number of nodes at this time step
            else:
                # Current time step index range
                pivot_index = index_root + pivot # In this case, k is an odd number, k/2 has a remainder 1.
                start_index = index_root
                end_index = index_root + k
                
                if DEBUG: print("-----------------------------------------------------------------------")
                if DEBUG: print("start_index:", start_index, ", end_index:", end_index, ", pivot_index:", pivot_index)
                
                # Previous time step index range
                prev_pivot_index = prev_index_root + pivot
                prev_start_index = prev_index_root
                prev_end_index = prev_index_root + (k-1) # Within previous time step #(k-1), index from 0,...,(k-1)
                
                if DEBUG: print("prev_start_index:", prev_start_index, 
                                ", prev_end_index:", prev_end_index, 
                                ", prev_pivot_index:", prev_pivot_index)
                
                si = self.S[prev_pivot_index]
                if DEBUG: print("si =", si)
                assert si == S0
                
                Fi = comp*si
                if DEBUG: print("Fi = ", Fi)
                    
                Call = rBergomi.price_option(rBergomiN, t, si, opttype=1, earlyexercise=False, r=self.r)[0][0]
                if DEBUG: print("C(S,t') =", Call)
                
                lam = self.lam[prev_pivot_index]
                if DEBUG: print("lambda_i =", lam)
                
                Sigma = 0
                for j in range(pivot): #range(prev_pivot_index+1, prev_end_index+1):
                    if DEBUG: print("Calculating Sigma, adding ", prev_index_root+j)
                    Sigma += self.lam[prev_index_root+j]*(comp*self.S[prev_index_root+j]-si)
                if DEBUG: print("Sigma =", Sigma)
                
                # Setup pivot
                self.S[pivot_index] = si*(comp*Call + lam*si - Sigma)/(lam*Fi - comp*Call + Sigma)
                self.S[pivot_index+1] = si*si/self.S[pivot_index]
                
                # Use call options
                if DEBUG: print("Decreasing from", pivot-1, "to 0 to use call options.")
                for i in range(pivot, 0, -1):
                    curr_index = index_root + i
                    prev_index = prev_index_root + i - 1
                    si = self.S[prev_index]
                    Si = self.S[curr_index]
                    Fi = comp*si
                    Call = rBergomi.price_option(rBergomiN, t, si, opttype=1, earlyexercise=False, r=self.r)[0][0]
                    lam = self.lam[prev_index]
                    Sigma = 0
                    for j in range(i-1):
                        if DEBUG: print("Calculating Sigma, adding ", prev_index_root+j, ", intra-index:", j)
                        Sigma += self.lam[prev_index_root+j]*(comp*self.S[prev_index_root+j]-si)
                    
                    if DEBUG: print("i:", i, 
                                    ", Si =", Si, 
                                    ", si =", si, 
                                    ", Fi =", Fi, 
                                    ", Call =", Call, 
                                    ", lambda =", lam, 
                                    ", Sigma =", Sigma)
                        
                    if DEBUG: print("Fill", curr_index-1, "using", curr_index, "and", prev_index)
                    
                    self.S[curr_index-1] = (Si*(comp*Call-Sigma)-lam*si*(Fi-Si))/(comp*Call-Sigma-lam*(Fi-Si))
                    
                    # Mannual overrides
                    if self.S[curr_index-1] <= comp*si or ( (i-2) >= 0 and self.S[curr_index-1] >= comp*self.S[prev_index-1] ):
                        print("Mannual override calls (B.1):", curr_index-1,curr_index, prev_index, prev_index+1)
                        self.S[curr_index-1] = self.S[curr_index]*self.S[prev_index]/self.S[prev_index+1]
                    if self.S[curr_index-1] > Smax: 
                        print("Mannual override calls (B.2):", curr_index-1, Smax)
                        self.S[curr_index-1] = Smax
                    
                # Use put options
                if DEBUG: print("Increasing from", pivot+2, "to", k, "to use put options.")
                for i in range(pivot+1, k, 1):
                    curr_index = index_root + i
                    prev_index = prev_index_root + i
                    si = self.S[prev_index]
                    Si1 = self.S[curr_index]
                    Fi = comp*si
                    Put = rBergomi.price_option(rBergomiN, t, si, opttype=-1, earlyexercise=False, r=self.r)[0][0]
                    lam = self.lam[prev_index]
                    Sigma = 0
                    for j in range(i+1, k):
                        if DEBUG: print("Calculating Sigma, adding ", prev_index_root+j, ", intra-index:", j)
                        Sigma += self.lam[prev_index_root+j]*(si-comp*self.S[prev_index_root+j])
                    
                    if DEBUG: print("i:", i, 
                                    ", S(i+1) =", Si1, 
                                    ", si =", si, 
                                    ", Fi =", Fi, 
                                    ", Put =", Put, 
                                    ", lambda =", lam, 
                                    ", Sigma =", Sigma)
                        
                    if DEBUG: print("Fill", curr_index+1, "using", curr_index, "and", prev_index)
                    self.S[curr_index+1] = (Si1*(comp*Put-Sigma)+lam*si*(Fi-Si1))/(comp*Put-Sigma+lam*(Fi-Si1))
                    
                    # Mannual overrides
                    if ( self.S[curr_index+1] >= comp*si or 
                         ( (i+1) <= (k-1) and self.S[curr_index+1] <= comp*self.S[prev_index+1] ) ):
                        print("Mannual override puts (B.1):", curr_index+1, curr_index, prev_index, prev_index-1)
                        self.S[curr_index+1] = self.S[curr_index]*self.S[prev_index]/self.S[prev_index-1]
                    if self.S[curr_index+1] < Smin: 
                        print("Mannual override puts (B.2):", curr_index+1, Smin)
                        self.S[curr_index+1] = Smin
            
            if DEBUG: print("------------------------ Calculate p's --------------------------------")
            for i in range(k):
                prev_index = prev_index_root + i
                curr_index = index_root + i
                
                if DEBUG: print("Calculate p at index:", prev_index, 
                                ", numerator:", (comp*self.S[prev_index]-self.S[curr_index+1]), 
                                ", denominator:", (self.S[curr_index]-self.S[curr_index+1]))
                
                self.p[prev_index] = (comp*self.S[prev_index]-self.S[curr_index+1])/(self.S[curr_index]-self.S[curr_index+1])
            
            if DEBUG: print("------------------------ Calculate lambda's ---------------------------")
            
            if DEBUG: print("Top lambda at", index_root, " calculated using previous index at", prev_index_root)
            self.lam[index_root] = disc*self.p[prev_index_root]*self.lam[prev_index_root]
            
            if DEBUG: print("Bottom lambda at", index_root+k, " calculated using previous index at", prev_index_root+k-1)
            
            self.lam[index_root+k] = disc*(1-self.p[prev_index_root+k-1])*self.lam[prev_index_root+k-1]
            for i in range(1,k):
                if DEBUG: print("Middle lambda at", index_root+i, 
                                "calculated using previous index at", prev_index_root+i, 
                                "and", prev_index_root+i-1)
                self.lam[index_root+i] = disc*( self.p[prev_index_root+i]*self.lam[prev_index_root+i] +
                                                (1-self.p[prev_index_root+i-1])*self.lam[prev_index_root+i-1])
            
            if DEBUG: print("=======================================================================")
            
            if DEBUG: 
                print("Stock prices:")
                print(localVol.S)
                print("p:")
                print(localVol.p)
                print("Lamda:")
                print(localVol.lam)
            
            if DEBUG: print("**************************************************************************************")
            
    def draw_price_tree(self, N, T):
        
        NN = ((N+1)*(N+2))//2
        dt = T/N
        
        g = nx.Graph()

        # Draw nodes
        for k in range(N+1):
            t = k*dt
            index_root = (k*(k+1))//2
            for i in range(k+1):
                nn = index_root + i
                g.add_node(nn,pos=(t,self.S[nn]))

        # Draw edges
        for k in range(N):
            index_root = (k*(k+1))//2
            next_index_root = ((k+1)*(k+2))//2
            for i in range(k+1):
                nn = index_root + i
                g.add_edge(nn, next_index_root+i, weight=0.1)
                g.add_edge(nn, next_index_root+i+1, weight=0.1)
                
        pos=nx.get_node_attributes(g,'pos')
        nx.draw(g,pos,node_size=50,node_color="skyblue",alpha=0.5)

        plt.rcParams["figure.figsize"] = (10,6)
        plt.show()
        
    def price_option(self, N, T, strikes, opttype, earlyexercise, r):
        
        NN = (N+1)*(N+2)//2
        assert len(self.S) == NN
        assert opttype in {1,-1}
        
        NI = (N*(N+1))//2
        Smin = self.S[NI:NN].min()
        Smax = self.S[NI:NN].max()

        dt = T/N
        disc = math.exp(-r*dt)
                
        Klist = np.atleast_1d(strikes)
        Vlist = []
        EarlyExerciseBoundaryList = []
        for K in Klist:
            #assert (opttype == 1 and K <= Xmax) or (opttype == -1 and K >= Xmin)
            # Initialization
            V = np.zeros(NN)
            
            # Terminal payoff
            NI = (N*(N+1))//2
            for nn in range(NI, NN):
                if opttype == 1: 
                    v = self.S[nn]-K
                else: 
                    v = K-self.S[nn]
                if v>0: V[nn] = v
                    
            # Backward valuation
            for n in range(N-1,-1,-1):
                NI = n*(n+1)//2
                NF = (n+1)*(n+2)//2
                for nn in range(NI, NF):
                    next_up = NF + nn - NI
                    next_down = NF + nn - NI + 1
                    V[nn] = disc*(self.p[nn]*V[next_up]+(1-self.p[nn])*V[next_down])
                    if earlyexercise:
                        if opttype == 1: 
                            v = self.S[nn]-K
                        else:
                            v = K-self.S[nn]
                        if v > V[nn]: 
                            V[nn] = v

            Vlist.append(V[0])
            
            # Early exercise boundary
            #if earlyexercise:
            #    boundary = []
            #    for n in range(N+1):
            #        NI = n*(n+1)//2
            #        NF = (n+1)*(n+2)//2
            #        
            #        if opttype == 1: stoppingS = self.S[NI]
            #        else: stoppingX = stoppingS = self.S[NF-1]
            #        
            #        for nn in range(NI, NF):
            #            if opttype == 1: v = self.S[nn]-K
            #            else: v = K-self.S[nn]
            #            
            #            if v > V[nn]:
            #                if opttype == 1:
            #                    if self.S[nn] < stoppingS:
            #                        stoppingS = self.S[nn]
            #                else:
            #                    if self.S[nn] > stoppingS:
            #                        stoppingS = self.S[nn]
            #                        
            #        boundary.append(stoppingS)
            #        
            #    EarlyExerciseBoundaryList.append(boundary)
                
        return Vlist, EarlyExerciseBoundaryList
            