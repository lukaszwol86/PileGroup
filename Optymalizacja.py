import numpy as np
import pandas as pd
import math
import time

start_time = time.time()

class Optimize:
    def __init__(self):
        pass
    def Generate_allowed_piles(self,B_min,B_max,L_min,L_max,dx,dy,d_edge):
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


        self.pile_x = np.sort(pile_x)
        self.pile_y = np.sort(pile_y)

        x_m = sum(pile_x)/len(pile_x)
        y_m = sum(pile_y) / len(pile_y)
        k = 0
        for x in self.pile_x:
            for y in self.pile_y:
                r = ((x-x_m)**2+(y-y_m)**2)**0.5
                all_piles.append({'Ind': k, 'x':x,'y':y,'r':r})
                k += 1

        self.all_piles = pd.DataFrame(all_piles, index=None)
        prob = []
        prob_x = []
        prob_y = []
        r_2=0
        x_2=0
        y_2 = 0
        for index, row in self.all_piles.iterrows():
            r_2 += (row.r)**2
            x_2 += (row.x) ** 2
            y_2 += (row.y) ** 2
        for index, row in self.all_piles.iterrows():
            prob.append(row.r/r_2)
            prob_x.append(abs(row.x) / x_2)
            prob_y.append(abs(row.y) / y_2)
        self.all_piles['prob']=prob
        self.all_piles['prob_x'] = prob_x
        self.all_piles['prob_y'] = prob_y

    def Initial_number_piles(self,Loads,P_cap):
        self.n_p = math.ceil(max(Loads.Fz)/P_cap)
        n_y = 0

        for index, row in Loads.iterrows():
            n_y = max(n_y,abs(row.Fy)/(abs(row.Fy)+abs(row.Fx)))
        self.n_y = n_y
        self.n_x = 1-n_y

    def Generate_piles(self,N,L):
        n_piles = int(min((2*self.n_p)//4*4,len(self.all_piles)))
        ny_piles = max(int((self.n_y*n_piles)//4*4),4)
        print(self.n_y)

        nx_piles = n_piles-ny_piles
        piles = self.all_piles.nlargest(n=n_piles, columns='prob')
        piles['Nr']= np.arange(0, n_piles, 1, dtype=int)
        piles['z'] = 0
        piles['Type'] = 1
        piles['N'] = N
        piles['L'] = L
        piles['redu'] = 1
        piles['Angle'] = 0
        piles['Acti'] = np.random.randint(2, size=n_piles)
        piles_y = piles.nlargest(n=ny_piles, columns='prob_y')
        for item, row in piles_y.iterrows():
            if row.y>0:
                piles.loc[piles['Nr']== row.Nr,'Angle']=90
            else:
                piles.loc[piles['Nr']== row.Nr,'Angle']= 270

        for item, row in piles.iterrows():
            if row.Angle == 0:
                if row.x>0:
                    piles.loc[piles['Nr'] == row.Nr, 'Angle'] = 0
                else:
                    piles.loc[piles['Nr'] == row.Nr, 'Angle'] = 180

        self.piles = piles[["Nr",'Type','Angle','x','y','z','N','L','redu','Acti']]

    def Update_piles(self, action):
        pile_no = action//4
        modif = action%4
        if modif ==0:
            self.piles.loc[self.piles['Nr'] == pile_no, 'Acti'] = 0
        elif modif ==1:
            self.piles.loc[self.piles['Nr'] == pile_no, 'Acti'] = 1
        elif modif ==2:
            self.piles.loc[self.piles['Nr'] == pile_no, 'Angle'] -= 90
        elif modif ==3:
            self.piles.loc[self.piles['Nr'] == pile_no, 'Angle'] += 90



Optimiser = Optimize()
Loads = pd.read_csv('Loads.csv', delimiter=";")

Optimiser.Generate_allowed_piles(B_min=-10,B_max=10,L_min=-4,L_max=4,dx=2,dy=2,d_edge=0.6)
Optimiser.Initial_number_piles(Loads=Loads,P_cap=2000)



Optimiser.Generate_piles(N=4,L=30)

print("Used piles:")
print(Optimiser.piles)

Optimiser.Update_piles(0)
print(Optimiser.piles)
Optimiser.Update_piles(6)
print(Optimiser.piles)
print("--- %s seconds ---" % (time.time() - start_time))


