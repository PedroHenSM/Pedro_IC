'''
from recordtype import recordtype

Point = recordtype('Point', 'x y')
p = Point(1, 3)
p.x = 2
#assert p.x == 2
#assert p.y == 3
print(p.x)
'''

import random
import string
import numpy as np
#from operator import attrgetter
import operator


class Individuo(object):
    def __init__ (self,name,val):
        self.name = name
        self.val = val
    
    def __repr__(self):
        return str(self.__dict__)

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
     
    def __repr__(self):
        return str(self.__dict__)
        
class Produto(object):
   
  def __init__(self, nome, valor):
    self.__nome = nome
    self.__valor = valor
     
  def __repr__(self):
    return "nome:%s valor:%s" % (self.__nome, self.__valor)
 
  def get_nome(self):
    return self.__nome
 
  def get_valor(self):
    return self.__val
        

class Aluno:
    def __init__(self, nome, matricula):
        self.nome = nome
        self.matricula = matricula
 
    def __str__(self):
        return "%s - %s" % (str(self.nome), str(self.matricula))
    
    
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
        
    
    alunos = [Aluno("".join(random.sample(string.ascii_letters, 5)), random.randint(0, 10)) for i in range(10)]
    for aluno in alunos:
        print (aluno)
    alunos.sort(key=operator.attrgetter("matricula","nome"))
    print("Ordenados !!!")
    for aluno in alunos:
        print (aluno)

    print(p1)
    
    
            
    
