#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 20 17:53:50 2018

@author: pedrohen
"""

import numpy as np
import sys
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
numG = Number of inequalities constraints
numH = Number of equalities constraints
violations = Array of violations
violatonSum = Sum of violations
fitness = Fitness of each individual (for APM)


'''



class Individual(object):
    population = []
    
    def __init__(self,n):
        #self.n.append(n)
        self.n = n
    '''def __init__(self,n,sigma,objectiveFunction,numG,numH,violations,violationSum,fitness):
        Individual.n = n
        Individual.sigma = sigma
        Individual.objectiveFunction = objectiveFunction
        Individual.numG = numG
        Individual.numH = numH
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
    def sumViolations(self,popSize,numG,numH):
        for i in range(popSize):
            for j in range (numG+numH):
                if(j < numG):
                    #self.individuals[i].violations[j]
                    #self.individuals[i].violations.
'''


'''
Recebe vetor de inviduos e inicializa pop
Ou pode chamar um for na main, passando individo e incializando 
separadamente a populacao

'''


class Population(object):
    def __init__(self,popSize,n,function): 
        self.individuals = []
        for i in range (popSize): #
            values = []
            for j in range (n): # Dimension
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
            print("Lista")
            print(list1)
            print("Values")
            print(values)
            self.individuals.append(Individual(values))
          
    def printPopulation(self,popSize,n):
        for i in range(popSize):
            print("Individual: {}".format(i))
            for j in range(n):
                print("N: {}".format(self.individuals[i].n[j]))
    






def initializeConstraints(function):
    if(function == 1):
        numG = 2
        numH = 0
        return (numG,numH,numG+numH)
    elif (function == 2 or function == 17):
        numG = 2
        numH = 1
        return (numG,numH,numG+numH)
    elif (function == 3 or function == 9 or function == 10 or function == 11):
        numG = 0
        numH = 1
        return (numG,numH,numG+numH)
    elif (function == 4):
        numG = 0
        numH = 4
        return (numG,numH,numG+numH)
    elif (function == 5 or function == 6):
        numG = 0
        numH = 2
        return (numG,numH,numG+numH)
    elif (function == 7 or function == 8):
        numG = 1
        numH = 0
        return (numG,numH,numG+numH)
    elif (function == 12 or function == 18):
        numG = 1
        numH = 1
        return (numG,numH,numG+numH)
    elif (function == 13 or function == 14 or function == 15):
        numG = 3
        numH = 0
        return (numG,numH,numG+numH)
    elif (function == 16):
        numG = 2
        numH = 2
        return (numG,numH,numG+numH)
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
def GA(function,seed,penaltyMethod,popSize,n,generatedOffspring,maxFE,probCrossover,esType,globalSigma): # Genetic Algorithm
    functionEvaluations = 0
    offspringPopSize = popSize * generatedOffspring
    np.random.seed(seed) 
    numG,numH,numConstraints = initializeConstraints(function)
    population = Population(popSize,n,function)
    population.printPopulation(popSize,n)
    
    
    
    

if __name__ == '__main__':
    np.random.seed(1)
    population = Population(2,2,1)
    population.printPopulation(2,2)
    print("Funcionando")
    #print(population.individuals[1].n[0])
    print("teste")
    for i in range(5):
        print(np.random.uniform(-5,5))
    print("Rodando AG")
    GA(1,1,1,40,10,1,20000,1,1,1)
    
    
    '''
    population = Population(2,vecInd)
    print("opa")
    print(population)
    print(population.individuals[0].name)
    print(population.individuals[1].name)
    population.printPopulation(2)
    '''
    
    
    