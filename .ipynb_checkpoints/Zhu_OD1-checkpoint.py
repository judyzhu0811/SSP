import numpy as np

def angular_momentum(pos, vel):
    m = 58.3760923664
    L = m * np.cross(pos, vel) #angular momentum 
    return L
def main():
    Zhu_OD1 = np.loadtxt('Zhu_OD1Input.txt') #extract array
    pos = Zhu_OD1[0] #first row of array
    vel = Zhu_OD1[1] #second row of array
    print(angular_momentum(pos, vel))