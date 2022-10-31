import numpy as np
import pandas as pd
import math
import time

start_time = time.time()

class Optimize:
    def __init__(self,P_max,P_min,def_lim):
        self.P_max=P_max
        self.P_min = P_min
        self.def_lim=def_lim
    def Generate_allowed_piles(self,B_min,B_max,L_min,L_max,dx,dy,d_edge,N,L):
        nx = (((B_max-B_min)/2-d_edge-dx/2)/dx+1)//1
        ny = (((L_max-L_min)/2 - d_edge-dy/2) / dy+1)//1
        pile_x = []
        pile_y = []
        all_piles = []
        for i in range(int(nx)):
            pile_x.append(B_min+d_edge+dx*i)
        for i in range(int(nx)):
            pile_x.append(B_max-d_edge-dx*i)
        if  (B_max-B_min)/2-d_edge-(nx-1)*dx>dx:
            pile_x.append((B_max+B_min)/2)



        for i in range(int(ny)):
            pile_y.append(L_min+d_edge+dy*i)
        for i in range(int(ny)):
            pile_y.append(L_max-d_edge-dy*i)
        if  (L_max-L_min)/2-d_edge-(ny-1)*dy>dy:
            pile_y.append((L_max+L_min)/2)


        pile_x = np.sort(pile_x)
        pile_y = np.sort(pile_y)


        k = 0
        for x in pile_x:
            for y in pile_y:
                all_piles.append({'Nr': k, 'x':x,'y':y,'z':0,'Type':1,'N':N,'L':L,'redu':1})
                k += 1

        self.all_piles = pd.DataFrame(all_piles, index=None)




    def Generate_piles(self,genom):

        piles = self.all_piles
        piles['Angle'] = 0
        piles['Acti'] = 1

        for item, row in piles.iterrows():
            piles.loc[piles['Nr']== item,'Acti']=genom[2*item]
            piles.loc[piles['Nr'] == item, 'Angle'] = genom[2 * item+1]*90

        self.piles = piles[["Nr",'Type','Angle','x','y','z','N','L','redu','Acti']]







if __name__ == "__main__":
    Optimiser = Optimize()
    Loads = pd.read_csv('Loads.csv', delimiter=";")
    Optimiser.Generate_allowed_piles(B_min=-10,B_max=10,L_min=-4,L_max=4,dx=2,dy=2,d_edge=0.6,N=4,L=30)
    Optimiser.Generate_piles()
    print("--- %s seconds ---" % (time.time() - start_time))


