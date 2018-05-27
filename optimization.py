'''
from recordtype import recordtype

Point = recordtype('Point', 'x y')
p = Point(1, 3)
p.x = 2
#assert p.x == 2
#assert p.y == 3
print(p.x)
'''


import numpy as np

class Individuo(object):
    def __init__ (self,name,val):
        self.name = name
        self.val = val

class Pop(object):
    def __init__(self,n,name):
        self.individuos = []
        for i in range (n):
            self.individuos.append(Individuo(name[i],(i+1)*10))
        
        
    def imprime(self,n):
        for i in range(n):    
            print("Nome: {}".format(self.individuos[i].name))
            print("Val: {}".format(self.individuos[i].val))

    def troca(self,p):
        self.individuos[1] = p.individuos[1]
        
if __name__ == '__main__':
    name1 = ["Pedro","Carol"]
    name2 = ["Lucas", "Larissa"]
    p1 = Pop(2,name1)
    p2 = Pop(2,name2)
    #print("IMpressaona main {}".format(p2.individuos[].name))

    p2.individuos[0].val = 33
    p2.individuos[1].val = 66
    p1.imprime(2)
    p2.imprime(2)
    p1.troca(p2)
    p1.imprime(2)
    p2.imprime(2)
    
    for i in range(10):
        print(np.random.uniform(20,10))
    
