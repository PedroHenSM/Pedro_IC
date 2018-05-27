#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 20 17:53:50 2018

@author: pedrohen
"""

import numpy as np
import sys
from functions import Functions
from random import uniform

'''
typedef struct Individuo {
    float x[TAM_X_MAXIMO];
    float sigma[TAM_X_MAXIMO]; // Strategy Parameters ( Parametros de Controle )
    float funcaoObjetivo[TAM_X_MAXIMO];
    float g[NUM_RESTRICOES];
    float h[NUM_RESTRICOES];
    //float penaltyCoefficients[NUM_RESTRICOES]; // For APM
    float fitness; // Fitness for APM
    //float avgObjFunc; // For APM
    float v[2*NUM_RESTRICOES]; // Vetor de violacoes
    float violacao; // Soma das violacoes // Soma das violacoes
}Individuo;
'''

'''
n = search space dimension
sigma = ES parameter
objectiveFunction = Objective function to be minimized
g = Number of inequalities constraints
h = Number of equalities constraints
violations = Array of violations
violatonSum = Sum of violations
fitness = Fitness of each individual (for APM)


'''



class Individual(object):
    #population = []
    
    
    def __init__(self,n,objectiveFunction = [], g = [], h = [],violations = [], violationSum = 0):
        self.n = n
        self.objectiveFunction = objectiveFunction
        self.g = g
        self.h = h
        self.violations = violations
        self.violationSum = violationSum
    
    '''
    def __init__(self,n,objectiveFunction,g,h):
        #self.n.append(n)
        self.n = n
        self.objectiveFunction = objectiveFunction
        self.g = g
        self.h = h 
    '''
    
    
    '''
    def __init__(self,n,sigma,objectiveFunction,g,h,violations,violationSum,fitness):
        Individual.n = n
        Individual.sigma = sigma
        Individual.objectiveFunction = objectiveFunction
        Individual.g = g
        Individual.h = h
        Individual.violations = violations
        Individual.violationSum = violationSum
        Individual.fitness = fitness
    '''
    
    '''
    def __init__(self,name,vec):
        self.name = name
        self.vec = vec
    '''
        
        
'''
    def sumViolations(self,popSize,g,h):
        for i in range(popSize):
            for j in range (g+h):
                if(j < g):
                    #self.individuals[i].violations[j]
                    #self.individuals[i].violations.
'''


'''
Recebe vetor de inviduos e inicializa pop
Ou pode chamar um for na main, passando individo e incializando 
separadamente a populacao

'''


class Population(object):
    def __init__(self,popSize,nSize,function): 
        self.individuals = []
        for i in range (popSize): #
            values = []
            for j in range (nSize): # Dimension
                if (function == 1):
                    values.append(np.random.uniform(0,10))
                elif (function == 2):
                    values.append(np.random.uniform(-5,5))
                elif (function == 3 or function == 12 or function == 14 or function == 15):
                    values.append(np.random.uniform(-1000,1000))
                elif (function == 4 or function == 18):
                    values.append(np.random.uniform(-50,50))
                elif (function == 5 or function == 6):
                    values.append(np.random.uniform(-600,600))
                elif (function == 7 or function == 8):
                    values.append(np.random.uniform(-140,140))
                elif (function == 9 or function == 10 or function == 13):
                    values.append(np.random.uniform(-500,500))
                elif (function == 11):
                    values.append(np.random.uniform(-100,100))
                elif (function == 16 or function == 17):
                    values.append(np.random.uniform(-10,10))
                else:
                    print("Function not encountered")
                    sys.exit("Function not encountered")
                #self.individuals[i].append(Individual(value))
                #self.individuals.append(Individual(vecInd[i].name,vecInd[i].vec))
            list1 = []
            list1.append(1)
            list1.append(2)
            list1.append(3)
            #print("Lista")
            #print(list1)
            #print("Values")
            #print(values)
            self.individuals.append(Individual(values))
          
    def printPopulation(self,popSize,n):
        for i in range(popSize):
            print("Individual: {}".format(i))
            for j in range(n):
                print("N: {}".format(self.individuals[i].n[j]))
                
    
    def evaluate(self,popSize,function,nSize,gSize,hSize,fe): # verificar parametro fe 
        for i in range(popSize):
            fe = fe + 1
            if (function == 1):
                Functions.C01(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 2):
                Functions.C02(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 3):
                Functions.C03(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 4):
                Functions.C04(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 5):
                Functions.C05(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 6):
                Functions.C06(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 7):
                Functions.C07(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 8):
                Functions.C08(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 9):
                Functions.C09(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 10):
                Functions.C10(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 11):
                Functions.C11(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 12):
                Functions.C12(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 13):
                Functions.C13(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 14):
                Functions.C14(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 15):
                Functions.C15(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 16):
                Functions.C16(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 17):
                Functions.C17(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 18):
                Functions.C18(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
            elif (function == 99):
                print("AvaliouFuncaoTeste")
            else:
                print("Function not encountered")
                sys.exit("Function not encountered")
        return fe
    
    
    def sumViolations(self,popSize,gSize,hSize):
        for i in range(popSize):
            idxG = 0
            idxH = 0
            for j in range(gSize + hSize):
                if (j < gSize):
                    self.individuals[i].violations[j] = self.individuals[i].g[idxG]
                    idxG = idxG + 1
                else:
                    self.individuals[i].violations[j] = self.individuals[i].h[idxH] - 0.0001 # Converts equality on inequality
                    idxH = idxH + 1
            self.individuals[i].violationSum = np.sum(self.individuals[i].violations)
        
    
    def selection(self,parents,parentsSize,nSize,offspringsSize): # Compares violations and objective function
        for i in range(offspringsSize):
            idx1 = np.random.randint(parentsSize)
            idx2 = np.random.randint(parentsSize)
            if (parents.individuals[idx1].violationSum < parents.individuals[idx2].violationSum):
                self.individuals[i] = parents.individuals[idx1]
            elif (parents.individuals[idx1].violationSum > parents.individuals[idx2].violationSum):
                self.individuals[i] = parents.individuals[idx2]
            elif (parents.individuals[idx1].objectiveFunction[0] < parents.individuals[idx2].objectiveFunction[0]):
                self.individuals[i] = parents.individuals[idx1]
            else:
                self.individuals[i] = parents.individuals[idx2]
              
                
    def standardCrossover(self,nSize,offspringsSize): # Randomizes uniform value between n1 and n2   
        for i in range(offspringsSize):
            for j in range(nSize):
                n1 = self.individuals[i].n[j]
                n2 = self.individuals[i+1].n[j]
                self.individuals[i].n[j] = np.random.uniform(n1,n2)
                self.individuals[i+1].n[j] = np.random.uniform(n1,n2)
    
    
    
    """
    eta : Crowding degree of the crossover. A high eta will produce
    children resembling to their parents, while a small eta will
    produce solutions much more different., usually 2 to 5
    """
    def sbCrossover(self,eta,nSize,offspringsSize): # SBX (Simulated Binary Crossover)
        for i in range(0,offspringsSize,2):
            for j in range(nSize):
                rand = np.random.rand() # Generates random value between (0,1) 
                if (rand < 0.5):
                    beta = 2. * rand
                else:
                    beta = 1. / (2. * (1. - rand))
                
                #beta** = 1. /(eta + 1)
                beta = np.power(beta, 1. /(eta+1))
                x1 = self.individuals[i].n[j]
                x2 = self.individuals[i+1].n[j]
                self.individuals[i].n[j] = 0.5 * (((1 + beta) * x1) + ((1 - beta) * x2))
                self.individuals[i+1].n[j] = 0.5 * (((1 - beta) * x1) + ((1 + beta) * x2))
                
                
                #ind1[i] = 0.5 * (((1 + beta) * x1) + ((1 - beta) * x2))
                #ind2[i] = 0.5 * (((1 - beta) * x1) + ((1 + beta) * x2))
        
    def lol(self):
        imprimeTeste()
        
        
        
    def testeFor(self,tamanho,soma):
        print("Tamanho: {}".format(tamanho))
        for j in range(tamanho):
            soma = soma + 1
            print("j: {}".format(j))
        return soma
    
    
    
    def listaMuda(self,lista):
        print("LIsta em lista Muda",lista)
        lista.append(99)
        for i in range(len(lista)):
            lista[i] = i+1


def crossoverProbability(crossoverProb):
    if (np.random.randint(100) < crossoverProb):
        return 1
    else:
        return 0


def imprimeTeste():
    print("Funciona")



def initializeConstraints(function):
    if(function == 1):
        g = 2
        h = 0
        return (g,h,g+h)
    elif (function == 2 or function == 17):
        g = 2
        h = 1
        return (g,h,g+h)
    elif (function == 3 or function == 9 or function == 10 or function == 11):
        g = 0
        h = 1
        return (g,h,g+h)
    elif (function == 4):
        g = 0
        h = 4
        return (g,h,g+h)
    elif (function == 5 or function == 6):
        g = 0
        h = 2
        return (g,h,g+h)
    elif (function == 7 or function == 8):
        g = 1
        h = 0
        return (g,h,g+h)
    elif (function == 12 or function == 18):
        g = 1
        h = 1
        return (g,h,g+h)
    elif (function == 13 or function == 14 or function == 15):
        g = 3
        h = 0
        return (g,h,g+h)
    elif (function == 16):
        g = 2
        h = 2
        return (g,h,g+h)
    else:
        print("Function not encountered")
        sys.exit("Function not encountered")


'''
function = function to be minimized  --- TIPO FUNCAO
seed = seed to be used --- SEED
penaltyMethod = Standard or Adaptive Penalty Method (APM)
popSize = Size of population
generatedOffspring = Number of offsprings generated by each parent
maxFE = Max number of function evaluations
probCrossover = Crossover probability (0 - 100)
esType = Type of evolution strategies ( + or , )
globalSigma = Sigma parameter for ES. If sigmaGlobal equals 1 , each individual has 1 sigma. Else,
each individual has an array of sigma.

'''


def GA(function,seed,penaltyMethod,parentsSize,nSize,generatedOffspring,maxFE,crossoverProb,esType,globalSigma): # Genetic Algorithm
    np.random.seed(seed)
    crossoverType = 1
    functionEvaluations = 0
    offspringsSize = parentsSize * generatedOffspring
    np.random.seed(seed) 
    gSize,hSize,numConstraintsSize = initializeConstraints(function)
    parents = Population(parentsSize,nSize,function) # Initialize parents population
    offsprings = Population(offspringsSize,nSize,function)
    parents.evaluate(parentsSize,function,nSize,gSize,hSize,functionEvaluations)
    if (penaltyMethod == 1):
        parents.sumViolations(parentsSize,gSize,hSize)
    elif (penaltyMethod == 2):
        print("Adaptar codigo heder deopis")
        sys.exit(1)
    else:
        print("Penalthy method not encountered")
        sys.exit("Penalty method not encountered")
    while (functionEvaluations < maxFE):
        if (penaltyMethod == 1):
            offsprings.selection(parents,parentsSize,nSize,offspringsSize)
        elif (penaltyMethod == 2):
            print("Adaptar codigo heder deopis")
            sys.exit(1)
        if(crossoverProbability(crossoverProb)):
            if(crossoverType == 0):
                offsprings.standardCrossover(nSize,offspringsSize)
            elif(crossoverType == 1):
                offsprings.sbCrossover(2,nSize,offspringsSize)
            else:
                print("Crossover type not encountered")
                sys.exit("Crossover type not encountered")
    
    
    
    
    

if __name__ == '__main__':
    np.random.seed(1)
    pSize = 5
    nSize = 10
    l = []
    l.append(11)
    l.append(22)
    l.append(33)
    population = Population(pSize,nSize,1)
    #populations = Population(10,10,1)
    population.printPopulation(pSize,nSize)
    #populations.printPopulation(10,10)
    fes = 0
    population.lol()
    #population.evaluate(popSize,function,nSize,gSize,hSize,functionEvaluations)
    fes = population.evaluate(pSize,99,nSize,1,1,fes)
    fes = population.testeFor(pSize,fes)
    print(fes)
    print(l)
    population.listaMuda(l)
    print(l)
    print(np.sum(l))
    print("randons")
    for j in range(50):
        print(np.random.rand())
    '''
    population = Population(2,vecInd)
    print("opa")
    print(population)
    print(population.individuals[0].name)
    print(population.individuals[1].name)
    population.printPopulation(2)
    '''
    
    
    