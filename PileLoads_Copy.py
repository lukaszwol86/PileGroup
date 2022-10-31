import numpy as np
import pandas as pd
import math

np.set_printoptions(formatter={'float': '{: 0.2f}'.format})

Ep = 33*10**9    # [Pa]  Pile stiffness
Gp = 12*10**9       # [Pa]  Pile rotation stiffness
Ap = 109378*10**-6   # [m2] pile area
Jp = 58791*10**(4-4*3)  # [m4] pile moment of inertia
#L = 15          # [m] pile length

soil_type = "cohesive" # "no" / "cohesive" / "friction"

m = 0           # type of connection 0- hinge connection, 1 - fix connection

kd = 4615*1000        #[Pa] soil stiffness
nh = 3*10^6         # [N]


piles = pd.read_csv('piles.csv', delimiter=";")
Loads = pd.read_csv('Loads.csv', delimiter=";")

R = []
for index, row in Loads.iterrows():
    r = np.array([row.Fx*1000, row.Fy*1000, row.Fz*1000, row.Mx*1000, row.My*1000, row.Mz*1000]) #.T
    R.append(r)





Kp = []
for index, row in piles.iterrows():
    L=row.L
    k=np.zeros((6,6))
    if soil_type == "no":
        k[0,0] =  m*3*Ep*Jp/L**3
        k[1, 1] = k[0,0]
        k[2, 2] = Ap * Ep / L
        k[3, 3] = m*3*Ep*Jp/L
        k[4, 4] = k[3, 3]
        k[5, 5] = m *  Gp * Jp / L
    elif soil_type == "cohesive":
        Le = (4*Ep*Jp/kd)**(1/4)
        k[0,0] =  (m+1)*2*Ep*Jp/Le**3
        k[1, 1] = k[0,0]
        k[2, 2] = Ap * Ep / L
        k[3, 3] = m*2*Ep*Jp/Le
        k[4, 4] = k[3, 3]
        k[5, 5] = m *  Gp * Jp / L
        k[0, 4] = m * 2 * Ep * Jp / Le ** 2
        k[4, 0] =k[0, 4]
        k[1, 3] = -k[0, 4]
        k[3, 1] = -k[0, 4]
        #print(f"Pile {index+1} Le= {Le: 0.2f}")

    elif soil_type == "friction":
        Li = 1.8*(Ep*Jp/nh)**(1/5)
        k[0,0] =  (3*m+1)*3*Ep*Jp/Li**3
        k[1, 1] = k[0,0]
        k[2, 2] = Ap * Ep / L
        k[3, 3] = m*4*Ep*Jp/Li
        k[4, 4] = k[3, 3]
        k[5, 5] = m *  Gp * Jp / L
        k[0, 4] = m * 6 * Ep * Jp / Li ** 2
        k[4, 0] =k[0, 4]
        k[1, 3] = -k[0, 4]
        k[3, 1] = -k[0, 4]

    Kp.append(k)






Dp = []
Cp = []
Api = []
S = np.zeros((6,6))

for index, row in piles.iterrows():
    cp = np.array([[1, 0, 0, 0, 0, 0],
                  [0, 1, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0],
                  [0, -row.z, row.y, m, 0, 0],
                  [row.z, 0, -row.x, 0, m, 0],
                  [-row.y, row.x, 0, 0, 0, m]])
    beta = math.atan(1/row.N)
    alfa = math.radians(row.Angle)

 #   ap = np.array([[math.cos(beta)*math.cos(alfa), -math.sin(alfa), math.sin(beta)*math.cos(alfa), 0, 0, 0],
 #                 [math.cos(beta)*math.sin(alfa),  math.cos(alfa), math.sin(beta)*math.sin(alfa), 0, 0, 0],
 #                 [-math.sin(beta)*math.cos(alfa), 0, math.cos(beta), 0, 0, 0],
 #                 [0, 0, 0, math.cos(beta)*math.cos(alfa), -math.sin(alfa), math.sin(beta)*math.cos(alfa)],
 #                 [0, 0, 0, math.cos(beta)*math.sin(alfa),  math.cos(alfa), math.sin(beta)*math.sin(alfa)],
 #                 [0, 0, 0, -math.sin(beta)*math.cos(alfa), 0, math.cos(beta)]])


    ap = np.array([[math.cos(beta)*math.cos(alfa), -math.sin(alfa), math.sin(beta)*math.cos(alfa), 0, 0, 0],
                  [math.cos(beta)*math.sin(alfa),  math.cos(alfa), math.sin(beta)*math.sin(alfa), 0, 0, 0],
                  [-math.sin(beta), 0, math.cos(beta), 0, 0, 0],
                  [0, 0, 0, math.cos(beta)*math.cos(alfa), -math.sin(alfa), math.sin(beta)*math.cos(alfa)],
                  [0, 0, 0, math.cos(beta)*math.sin(alfa),  math.cos(alfa), math.sin(beta)*math.sin(alfa)],
                  [0, 0, 0, -math.sin(beta), 0, math.cos(beta)]])

    dp = np.matmul(cp,ap)
    sp = np.matmul(np.matmul(dp,Kp[index]),dp.T)

    Dp.append(dp)
    S=np.add(S,sp)

    Api.append(ap)
    Cp.append(cp)

U = []
DispP = np.zeros((len(Loads),len(piles),6))
ForcsP = np.zeros((len(Loads),len(piles),6))
for i, reaction in enumerate(R):
    u = np.matmul(np.linalg.inv(S),reaction)
    U.append(u)
    for n, dp in enumerate(Dp):
        disp = np.matmul(dp.T, u)
        force = np.matmul(Kp[n], disp)/1000
        #print(disp)
        DispP[i,n]= disp
        ForcsP[i,n] = force

#print(ForcsP)

N_P_LC = np.zeros((len(Loads),len(piles)))

for pile in range(len(piles)):
    for lc in range(len(Loads)):
        N_P_LC[lc,pile] = ForcsP[lc,pile,2]

print(N_P_LC)

#N_P_LC.tofile('Results.csv', sep = ',')