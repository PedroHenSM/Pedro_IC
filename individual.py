#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 20 17:53:50 2018

@author: pedrohen
"""

import numpy as np
import sys
import operator as op
import ctypes
from functions import Functions


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
    def __init__(self,n = None, objectiveFunction = None, g = None, h =None,violations = None ,sigma = None, violationSum = None, fitness = None):
        self.n = [-1 for i in range (30)] if n is None else n
        self.objectiveFunction = [-1 for i in range(1)] if objectiveFunction is None else objectiveFunction
        self.g = [-1 for i in range (5)] if g is None else g
        self.h = [-1 for i in range (5)] if h is None else h
        self.sigma = [-1 for i in range (30)] if sigma is None else sigma
        self.violations = [-1 for i in range (10) ] if violations is None else violations
        self.violationSum = -1 if violationSum is None else violationSum
        self.fitness = -1 if fitness is None else fitness

    def copyIndividual(self, individual, nSize, objectiveFunctionSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod):
        for j in range(nSize):  # Copies n
            self.n[j] = individual.n[j]
        for j in range(objectiveFunctionSize):
            self.objectiveFunction[j] = individual.objectiveFunction[j]
        for j in range(gSize):
            self.g[j] = individual.g[j]
        for j in range(hSize):
            self.h[j] = individual.h[j]
        if penaltyMethod == 1:  # Standard
            if globalSigma == 1:  # 1 sigma for each individual, utilizes only the first position of sigma array
                self.sigma[0] = individual.sigma[0]
            elif globalSigma == 2:  # 1 sigma for each n of each individual
                for j in range(nSize):
                    self.sigma[j] = individual.sigma[j]
            for j in range(constraintsSize): # gSize + hSize (violations)
                self.violations[j] = individual.violations[j]
            self.violationSum = individual.violationSum
        elif penaltyMethod == 2:  # APM
            self.fitness = individual.fitness

    def __repr__(self):
        return str(self.__dict__)

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
                    #value = np.random.uniform(0,10)
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
                    #sys.exit("Function not encountered")
                """   
                #self.individuals[i].append(Individual(value))
                if ( j == 0):
                    self.individuals.append()
                else: # 
                    self.individuals[i].n.append(value)
                self.individuals.insert(0,).n.append()
                """

            self.individuals.append(Individual(values))

    """
    def __repr__(self):
        return str(self.__dict__)
    """
    # self,n = None, objectiveFunction = None, g = None, h =None,violations = None ,sigma = None, violationSum = None, fitness = None
    def copyIndividual(self, idxDest, idxToBeCopy, population, nSize, objectiveFunctionSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod):
        for j in range(nSize):  # Copies n
            self.individuals[idxDest].n[j] = population.individuals[idxToBeCopy].n[j]
        for j in range(objectiveFunctionSize):
            self.individuals[idxDest].objectiveFunction[j] = population.individuals[idxToBeCopy].objectiveFunction[j]
        for j in range(gSize):
            self.individuals[idxDest].g[j] = population.individuals[idxToBeCopy].g[j]
        for j in range(hSize):
            self.individuals[idxDest].h[j] = population.individuals[idxToBeCopy].h[j]
        if penaltyMethod == 1:  # Standard
            if globalSigma == 1:  # 1 sigma for each individual, utilizes only the first position of sigma array
                self.individuals[idxDest].sigma[0] = population.individuals[idxToBeCopy].sigma[0]
            elif globalSigma == 2:  # 1 sigma for each n of each individual
                for j in range(nSize):
                    self.individuals[idxDest].sigma[j] = population.individuals[idxToBeCopy].sigma[j]
            for j in range(constraintsSize): # gSize + hSize (violations)
                self.individuals[idxDest].violations[j] = population.individuals[idxToBeCopy].violations[j]
            self.individuals[idxDest].violationSum = population.individuals[idxToBeCopy].violationSum
        elif penaltyMethod == 2:  # APM
            self.individuals[idxDest].fitness = population.individuals[idxToBeCopy].fitness

    def evaluate(self,popSize,function,nSize,gSize,hSize,fe):
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
        #print(self.individuals)
        return fe

    def sumViolations(self,popSize,gSize,hSize):
        for i in range(popSize):
            idxG = 0
            idxH = 0
            for j in range(gSize + hSize):
                if (j < gSize):
                    self.individuals[i].violations[j] = self.individuals[i].g[idxG]
                    #self.individuals[i].violations.append(self.individuals[i].g[idxG])
                    idxG = idxG + 1
                else:
                    self.individuals[i].violations[j] = self.individuals[i].h[idxH] - 0.0001 # Converts equality on inequality
                    #self.individuals[i].violations.append(self.individuals[i].h[idxH] - 0.0001)
                    idxH = idxH + 1
            # Only sums positives values. Negatives violations are not summed
            self.individuals[i].violationSum = sum(l for l in self.individuals[i].violations if l > 0)

    def selection(self,parents, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod):
        if (penaltyMethod == 1): # Not apm, 'standard' method | Compares violations and objective function
            for i in range(offspringsSize):
                idx1 = np.random.randint(0,parentsSize)
                idx2 = np.random.randint(0,parentsSize)
                #print("\n\nIndividuo Pai 1")
                #print(parents.individuals[idx1])
                #print("\n\nIndividuo Pai 2")
                #print(parents.individuals[idx2])
                if (parents.individuals[idx1].violationSum < parents.individuals[idx2].violationSum):
                    self.copyIndividual(i, idx1, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
                    # self.individuals[i] = parents.individuals[idx1]
                elif (parents.individuals[idx1].violationSum > parents.individuals[idx2].violationSum):
                    self.copyIndividual(i, idx2, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma,penaltyMethod)
                    # self.individuals[i] = parents.individuals[idx2]
                elif (parents.individuals[idx1].objectiveFunction[0] < parents.individuals[idx2].objectiveFunction[0]):
                    self.copyIndividual(i, idx1, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
                    # self.individuals[i] = parents.individuals[idx1]
                else:
                    self.copyIndividual(i, idx2, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
                    # self.individuals[i] = parents.individuals[idx2]
                #print("\n\nIndividuo Filho Escolhido")
                #print(self.individuals[i])
        elif (penaltyMethod == 2): # APM | Compares fitness of each individual
            for i in range(offspringsSize):
                idx1 = np.random.randint(0,parentsSize)
                idx2 = np.random.randint(0,parentsSize)
                if (parents.individuals[idx1].fitness < parents.individuals[idx2].fitness):
                    self.individuals[i] = parents.individuals[idx1]
                elif (parents.individuals[idx1].fitness >= parents.individuals[idx2].fitness):
                    self.individuals[i] = parents.individuals[idx2]


    def standardCrossover(self,nSize,offspringsSize):
        # Randomizes uniform value between n1 and n2
        for i in range(0, offspringsSize, 2):
            for j in range(nSize):
                n1 = self.individuals[i].n[j]
                n2 = self.individuals[i+1].n[j]
                if n1 < n2: # Just certificates that n1 is greater than n2
                    aux = n1
                    n1 = n2
                    n2 = aux
                self.individuals[i].n[j] = np.random.uniform(n1, n2)
                self.individuals[i+1].n[j] = np.random.uniform(n1, n2)

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


    def mutation(self,nSize,offspringsSize):
        MEAN = 0
        STD = 1
        for i in range(offspringsSize):
            for j in range(nSize):
                self.individuals[i].n[j] = self.individuals[i].n[j] + np.random.normal(MEAN,STD)

    def bounding(self,nSize,function,popSize):
        nMin = 0
        nMax = 0
        if (function == 1):
            nMin = 0
            nMax = 10
        elif (function == 2):
            nMin = -5.12
            nMax = 5.12
        elif (function == 3 or function == 12 or function == 14 or function == 15):
            nMin = -1000
            nMax = 1000
        elif (function == 4 or function == 18):
            nMin = -50
            nMax = 50
        elif (function == 5 or function == 6):
            nMin = -600
            nMax = 600
        elif (function == 7 or function == 8):
            nMin = -140
            nMax = 140
        elif (function == 9 or function == 10 or function == 13):
            nMin = -500
            nMax = 500
        elif (function == 11):
            nMin = -100
            nMax = 100
        elif (function == 16 or function == 17):
            nMin = -10
            nMax = 10
        else:
            print("Function not encountered")
            sys.exit("Function not encountered")

        for i in range(popSize):
            for j in range(nSize):
                if (self.individuals[i].n[j] > nMax):
                    self.individuals[i].n[j] = nMax
                elif (self.individuals[i].n[j] < nMin):
                    self.individuals[i].n[j] = nMin


    def sort(self,offsprings,penaltyMethod):
        if (penaltyMethod == 1): # Sort based on violatiom and then objective function
            self.individuals.sort(key = op.attrgetter ("violationSum",
                                                       "objectiveFunction"))
            offsprings.individuals.sort( key = op.attrgetter (
                "violationSum",  "objectiveFunction"))
        elif (penaltyMethod == 2): # Sort based on fitness and then objective function
            self.individuals.sort(key = op.attrgetter ("fitness",
                                                       'objectiveFunction'))
            offsprings.individuals.sort(key = op.attrgetter("fitness",
                                                            "objectiveFunction"))

    def elitism(self,offsprings,parentsSize, nSize, gSize, hSize,
                constraintsSize, globalSigma, penaltyMethod):
        if (penaltyMethod == 1): # Not apm
            copyStart = 5
            i = 0
            for j in range(copyStart,parentsSize): # J iterates on parents
                self.copyIndividual(j, i, offsprings, nSize, 1, gSize, hSize,
                                    constraintsSize, globalSigma,
                                    penaltyMethod)
                # self.individuals[j] = offsprings.individuals[i]
                i = i + 1
        elif (penaltyMethod == 2): # APM | Chooses the best individual beetwen parent-offspring
            for i in range(parentsSize):
                if (offsprings.individuals[i].fitness < self.individuals[i].fitness):
                    self.individuals[i] = offsprings.individuals[i]
                elif (offsprings.individuals[i].fitness == self.individuals[i].fitness):
                    if (offsprings.individuals[i].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]):
                        self.individuals[i] = offsprings.individuals[i]


    def printBest(self,parentsSize,penaltyMethod):
        best = bestIndividual(self,parentsSize,penaltyMethod)
        if (penaltyMethod == 1): # not apm
            print("Violation\t{:e}\tObjectiveFunction\t{:e}\n".format(
                best.violationSum,best.objectiveFunction[0]))
        elif (penaltyMethod == 2): # APM
            print("Fitness\t{:e}\tObjectiveFunction\t{:e}\n".format(
               best.fitness,best.objectiveFunction[0]))
            if best.fitness == best.objectiveFunction[0]:
                print("Fitness == objectiveFunction")



    def DESelection(self,offsprings,parentsSize, nSize, gSize, hSize,
                    constraintsSize, penaltyMethod):
        if penaltyMethod == 1:
            for i in range(parentsSize):
                if (offsprings.individuals[i].violationSum < self.individuals[i].violationSum): # If offspring is better than parent
                    self.copyIndividual(i, i, offsprings, nSize, 1, gSize,
                                        hSize, constraintsSize, -1, penaltyMethod)
                    #self.individuals[i] = offsprings.individuals[i] # Offspring
                    #  "becomes" parent
                elif (offsprings.individuals[i].violationSum == self.individuals[i].violationSum): # Compares if violatiomSun is equal for offspring and parents
                    if (offsprings.individuals[i].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]): # Offspring better than parent
                        self.copyIndividual(i, i, offsprings, nSize, 1, gSize,
                                            hSize, constraintsSize, -1,
                                            penaltyMethod)
                        # self.individuals[i] = offsprings.individuals[i]
        elif penaltyMethod == 2:
            print("Au")
            for i in range(parentsSize):
                if offsprings.individuals[i].fitness < self.individuals[i].fitness:
                    self.copyIndividual(i, i, offsprings, nSize, 1, gSize,
                                        hSize, constraintsSize, -1,
                                        penaltyMethod)
                elif offsprings.individuals[i].fitness == self.individuals[i].fitness:
                    if (offsprings.individuals[i].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]): # Offspring better than parent
                        self.copyIndividual(i, i, offsprings, nSize, 1, gSize,
                                            hSize, constraintsSize, -1,
                                            penaltyMethod)

    def initializeEvolutionStrategy(self,offsprings,nSize,parentsSize,
                                    offspringsSize,globalSigma):
        STD = 1
        for i in range(parentsSize): # Initialize sigma parameter for parents
            if (globalSigma == 1):
                self.individuals[i].sigma[0] = STD
            else:
                for j in range(nSize):
                    self.individuals[i].sigma[j] = STD

        for i in range(offspringsSize): # Initialize sigma parameter for offsprings
            if (globalSigma == 1):
                offsprings.individuals[i].sigma[0] = STD
            else:
                for j in range(nSize):
                    offsprings.individuals[i].sigma[j] = STD

    def sigmaSelfAdaptation(self,offsprings,nSize,parentsSize,generatedOffspring,globalSigma):
        MEAN = 0
        STD = 1
        l = 0  # iterates over offsprings
        for i in range(parentsSize):
            for j in range(generatedOffspring): # Each parents generate X offsprings
                tau = 1 / np.sqrt(nSize)
                epsilon = tau * np.random.normal(MEAN,STD) # Normal distribution
                if (globalSigma == 1): # 1 sigma for each individual, utilizes only the first position of sigma array
                    offsprings.individuals[l].sigma[0] = self.individuals[i].sigma[0] * np.exp(epsilon)
                for k in range(nSize):
                    if (globalSigma == 1): # 1 sigma for each individual, utilizes only the first position of sigma array
                        offsprings.individuals[l].n[k] = self.individuals[i].n[k] + offsprings.individuals[l].sigma[0] * np.random.normal(MEAN,STD)
                    else:
                        offsprings.individuals[l].sigma[k] = self.individuals[i].sigma[k] * np.exp(epsilon) * np.exp(epsilon) # NOTE: Is the double produ   ct necessary?
                        offsprings.individuals[l].n[k] = self.individuals[i].n[k] + offsprings.individuals[l].sigma[k] * np.random.normal(MEAN,STD)
                l = l + 1

    def elitismES(self,offsprings,parentsSize,offspringsSize,nSize,
                  gSize, hSize, constraintsSize, globalSigma, esType,
                  generatedOffspring, penaltyMethod):
        if (esType == 0): # Es + | Pick bests individuals among parents and offsprings
            # parents = Population(parentsSize,nSize,function) # Initialize
            # parents population
            aux = Population(parentsSize + offspringsSize, nSize, 1)
            k = 0
            for i in range(parentsSize + offspringsSize):
                if (i < parentsSize):
                    aux.copyIndividual(i, i, self, nSize, 1, gSize, hSize,
                                       constraintsSize, globalSigma,
                                       penaltyMethod)
                    # aux.individuals.append(self.individuals[i])
                else:
                    aux.copyIndividual(i, k, offsprings, nSize, 1, gSize,
                                       hSize, constraintsSize, globalSigma,
                                       penaltyMethod)
                    # aux.individuals.append(offsprings.individuals[k])
                    k = k + 1
            if (penaltyMethod == 1):
                aux.individuals.sort(key = op.attrgetter ("violationSum","objectiveFunction"))
            elif (penaltyMethod == 2):
                aux.individuals.sort(key = op.attrgetter ("fitness","objectiveFunction"))

            for i in range(parentsSize):
                self.copyIndividual(i, i, aux, nSize, 1, gSize, hSize,
                                    constraintsSize, globalSigma,
                                    penaltyMethod)
                # self.individuals[i] = aux.individuals[i]
        elif (esType == 1): # Es , | Each offspring only "exists" for 1 generation. Pick best offsprings of each parent
            j = 0
            for i in range(parentsSize):
                best = Individual()
                best.copyIndividual(offsprings.individuals[j], nSize, 1,
                                    gSize, hSize, constraintsSize,
                                    globalSigma, penaltyMethod)
                # best = offsprings.individuals[j]
                while (j < generatedOffspring*(i+1)): # Goes through each offspring of each parent
                    if (penaltyMethod == 1): # Not apm
                        if (offsprings.individuals[j].violationSum < best.violationSum):
                            best.copyIndividual(offsprings.individuals[j],
                                                nSize, 1, gSize, hSize,
                                                constraintsSize,
                                                globalSigma, penaltyMethod)
                            # best = offsprings.individuals[j]
                        elif (offsprings.individuals[j].violationSum == best.violationSum):
                            if (offsprings.individuals[j].objectiveFunction[0] < best.objectiveFunction[0]):
                                best.copyIndividual(offsprings.individuals[j],
                                                    nSize, 1, gSize, hSize,
                                                    constraintsSize,
                                                    globalSigma, penaltyMethod)
                                # best = offsprings.individuals[j]
                    elif (penaltyMethod == 2): # APM
                        if (offsprings.individuals[j].fitness < best.fitness):
                            best.copyIndividual(offsprings.individuals[j],
                                                nSize, 1, gSize, hSize,
                                                constraintsSize,
                                                globalSigma, penaltyMethod)
                            # best = offsprings.individuals[j]
                        elif (offsprings.individuals[j].fitness == best.fitness):
                            if (offsprings.individuals[j].objectiveFunction[0] < best.objectiveFunction[0]):
                                best.copyIndividual(offsprings.individuals[j],
                                                    nSize, 1, gSize, hSize,
                                                    constraintsSize,
                                                    globalSigma, penaltyMethod)
                                # best = offsprings.individuals[j]
                    j = j + 1
                # self.individuals[i] = best
                self.individuals[i].copyIndividual(best, nSize, 1, gSize,
                                                   hSize, constraintsSize,
                                                   globalSigma, penaltyMethod)
            if (penaltyMethod == 1):
                self.individuals.sort(key = op.attrgetter ("violationSum","objectiveFunction"))
            elif (penaltyMethod == 2):
                self.individuals.sort(key = op.attrgetter ("fitness","objectiveFunction"))
        else:
            print("Es type not encountered")
            sys.exit("Es type not encountered")

    def uniteConstraints(self,parentsSize,gSize,hSize):
        for i in range(parentsSize):
            idxG = 0
            idxH = 0
            for j in range(gSize+hSize):
                if (j < gSize):
                    self.individuals[i].violations[j] = self.individuals[i].g[idxG]
                    idxG = idxG + 1
                else:
                    self.individuals[i].violations[j] = np.abs(self.individuals[i].g[idxH]) - 0.0001
                    idxH = idxH + 1

    def calculatePenaltyCoefficients(self,popSize,numberOfConstraints,
                                     penaltyCoefficients,
                                     averageObjectiveFunctionValues):
        sumObjectiveFunction = 0

        #foreach candidate solution
        for i in range(popSize):
            sumObjectiveFunction = sumObjectiveFunction + self.individuals[i].objectiveFunction[0]
        # the absolute value of sumObjectiveFunction
        if (sumObjectiveFunction < 0):
            sumObjectiveFunction = sumObjectiveFunction * -1
        # the average of the objective function values
        averageObjectiveFunctionValues = sumObjectiveFunction / popSize
        #the denominator of the equation of the penalty coefficients
        denominator = 0
        #the sum of the constraint violation values
        #these values are recorded to be used in the next situation
        sumViolation = []
        for l in range(numberOfConstraints):
            sumViolation.append(0)
            for i in range(popSize):
                if (self.individuals[i].violations[l] > 0): # TODO: inser 'fancy' IF here
                    sumViolation[l] = sumViolation[l] + self.individuals[i].violations[l]
                """
                else:
                    sumViolation[l] = sumViolation[l] + 0 # TODO pode
                    # colocar um pass aqui, acho
                """
            denominator = denominator + sumViolation[l] * sumViolation[l]
        # the penalty coefficients  are calculated
        for j in range(numberOfConstraints): # TODO: insert 'fancy' IF here
            if (denominator == 0):
                penaltyCoefficients[j] = 0
                # penaltyCoefficients.append(0)
            else:
                penaltyCoefficients[j] = (sumObjectiveFunction /
                                          denominator) * sumViolation[j]
                # penaltyCoefficients.append((sumObjectiveFunction /
                # denominator) * sumViolation[j])

        return averageObjectiveFunctionValues  # returns list and value


    def calculateAllFitness(self,popSize,numberOfConstraints,penaltyCoefficients,averageObjectiveFunctionValues):
        for i in range(popSize):
            # indicates if the candidate solution is feasible
            infeasible = 0 # boolean
            # the penalty value
            penalty = 0
            for j in range(numberOfConstraints):
                if (self.individuals[i].violations[j] > 0):
                    # the candidate solution is infeasible if some constraint is violated
                    infeasible = 1
                    # the panalty value is updated
                    penalty = penalty + penaltyCoefficients[j] * self.individuals[i].violations[j]


            #fitness is the sum of the objective function and penalty values
            #the candidate solution is infeasible and just the objective function value,
            #otherwise
            #TODO: fancy if here
            if (infeasible):
                if (self.individuals[i].objectiveFunction[0] > averageObjectiveFunctionValues):
                    self.individuals[i].fitness = self.individuals[i].objectiveFunction[0] + penalty
                else:
                    self.individuals[i].fitness = averageObjectiveFunctionValues + penalty
            else:
                self.individuals[i].fitness = self.individuals[i].objectiveFunction[0]

    def printPopulation(self,popSize):
        for i in range(popSize):
            print(self.individuals[i])

def bestIndividual(parents,parentsSize,penaltyMethod):
    best = parents.individuals[0]
    if (penaltyMethod == 1): # not apm
        for i in range(1,parentsSize):
            if (parents.individuals[i].violationSum < best.violationSum):
                best = parents.individuals[i]
            elif (parents.individuals[i].violationSum == best.violationSum):
                if(parents.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]):
                    best = parents.individuals[i]
    elif (penaltyMethod == 2):
        for i in range(1,parentsSize):
            if (parents.individuals[i].fitness < best.fitness):
                best = parents.individuals[i]
            elif (parents.individuals[i].fitness == best.fitness):
                if (parents.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]):
                    best = parents.individuals[i]
    return best

def crossoverProbability(crossoverProb):
    if (np.random.randint(0,100) < crossoverProb):
        return 1
    else:
        return 0

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

def populationPick(solution, flags, parentsSize):
    while (1):
        contains = 0
        idx = np.random.randint(0,parentsSize)
        if (idx != solution ):
            for l in range(3):
                if (idx == flags[l]):
                    contains = 1
                    break
            if (contains == 0):
                break
    return idx

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


def GA(function, seed, penaltyMethod, parentsSize, nSize, generatedOffspring,
       maxFE, crossoverProb, esType, globalSigma): # Genetic
    # Algorithm
    esType = globalSigma = -1
    np.random.seed(seed)
    crossoverType = 1
    functionEvaluations = 0
    offspringsSize = parentsSize * generatedOffspring
    gSize,hSize,constraintsSize = initializeConstraints(function) # constraintsSize is the sum of gSize and hSize
    penaltyCoefficients = [-1 for i in range(constraintsSize)]
    avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
    parents = Population(parentsSize,nSize,function) # Initialize parents population
    offsprings = Population(offspringsSize,nSize,function)
    functionEvaluations = parents.evaluate(parentsSize,function,nSize,gSize,hSize,functionEvaluations)
    if (penaltyMethod == 1):
        parents.sumViolations(parentsSize,gSize,hSize)
    elif (penaltyMethod == 2):
        parents.uniteConstraints(parentsSize,gSize,hSize)
        avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
        parents.calculateAllFitness(parentsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
    else:
        print("Penalthy method not encountered")
        sys.exit("Penalty method not encountered")
    while (functionEvaluations < maxFE):
        offsprings.selection(parents, parentsSize, offspringsSize, nSize,
                             gSize, hSize, constraintsSize, globalSigma,
                             penaltyMethod) # Selection
        if(crossoverProbability(crossoverProb)):
            if(crossoverType == 0):
                offsprings.standardCrossover(nSize,offspringsSize) # Standard crossover
            elif(crossoverType == 1):
                offsprings.sbCrossover(2,nSize,offspringsSize) # SBX Crossover
            else:
                print("Crossover type not encountered")
                sys.exit("Crossover type not encountered")
        offsprings.mutation(nSize,offspringsSize)
        offsprings.bounding(nSize,function,offspringsSize)
        functionEvaluations = offsprings.evaluate(offspringsSize,function,nSize,gSize,hSize,functionEvaluations)
        if (penaltyMethod == 1):
            offsprings.sumViolations(offspringsSize,gSize,hSize)
        elif (penaltyMethod == 2):
            offsprings.uniteConstraints(offspringsSize,gSize,hSize)
            avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
            offsprings.calculateAllFitness(offspringsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
        parents.sort(offsprings, penaltyMethod)
        parents.elitism(offsprings, parentsSize, nSize, gSize, hSize,
                        constraintsSize, globalSigma, penaltyMethod)
        parents.sort(offsprings, penaltyMethod)
        print("FE: {}".format(functionEvaluations))
        if (penaltyMethod == 1):
            parents.printBest(parentsSize,penaltyMethod)
        elif (penaltyMethod == 2):
            print("implementar apm depois")

def DE(function,seed,penaltyMethod,parentsSize, nSize, generatedOffspring,
       maxFE, crossoverProb, esType, globalSigma): # Differential
    # Evolution
    crossoverProb = esType = globalSigma = -1
    np.random.seed(seed)
    CR = 0.9
    F = 0.5
    functionEvaluations = 0
    generatedOffspring = 1 # One parents only generates one offsprings
    offspringsSize = parentsSize * generatedOffspring
    gSize,hSize,constraintsSize = initializeConstraints(function) # Initialize constraints
    penaltyCoefficients = [-1 for i in range(constraintsSize)]
    avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
    parents = Population(parentsSize,nSize,function) # Initialize parents population
    offsprings = Population(offspringsSize,nSize,function) # Initialize offsprings population
    functionEvaluations = parents.evaluate(parentsSize,function,nSize,gSize,hSize,functionEvaluations) # Evaluate parents
    if(penaltyMethod == 1): # Padrao?  (not apm)
        parents.sumViolations(parentsSize,gSize,hSize)
    elif(penaltyMethod == 2): # // Adaptive Penalty Method ( APM )
        parents.uniteConstraints(parentsSize,gSize,hSize)
        avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
        parents.calculateAllFitness(parentsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
    else:
        print("Penalthy method not encountered")
        sys.exit("Penalty method not encountered")

    while (functionEvaluations < maxFE):
        offsprings.selection(parents, parentsSize, offspringsSize, nSize,
                             gSize, hSize, constraintsSize, globalSigma,
                             penaltyMethod)
        flags = [-1,-1,-1]
        for i in range(parentsSize):
            for l in range(len(flags)):
                flags[l] = populationPick(i,flags,parentsSize)
            R = np.random.randint(0,nSize) # Random index
            for j in range(nSize):
                Ri = np.random.rand() # Generates random number between (0,1)
                if (Ri < CR or j == R):
                    offsprings.individuals[i].n[j] = parents.individuals[flags[0]].n[j] + F * (parents.individuals[flags[1]].n[j] - parents.individuals[flags[2]].n[j])
                else:
                    offsprings.individuals[i].n[j] = parents.individuals[i].n[j]
        offsprings.bounding(nSize,function,offspringsSize)
        functionEvaluations = offsprings.evaluate(offspringsSize,function,nSize,gSize,hSize,functionEvaluations)
        if (penaltyMethod == 1): # Not apm
            offsprings.sumViolations(offspringsSize,gSize,hSize)
        elif (penaltyMethod == 2):
            offsprings.uniteConstraints(offspringsSize,gSize,hSize)
            avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
            offsprings.calculateAllFitness(offspringsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
        parents.DESelection(offsprings, parentsSize, nSize, gSize, hSize,
                            constraintsSize, penaltyMethod)
        parents.printBest(parentsSize,penaltyMethod)

def ES(function, seed, penaltyMethod, parentsSize, nSize, generatedOffspring,
       maxFE, crossoverProb, esType, globalSigma): # Evolution Strategy
    crossoverProb = -1
    np.random.seed(seed)
    functionEvaluations = 0
    offspringsSize = parentsSize * generatedOffspring
    gSize,hSize,constraintsSize = initializeConstraints(function) # Initialize constraints
    penaltyCoefficients = [-1 for i in range(constraintsSize)]
    avgObjFunc = 0  # will be subscribed on 'calculatePenaltyCoefficients'
    parents = Population(parentsSize,nSize,function) # Initialize parents population
    offsprings = Population(offspringsSize,nSize,function) # Initialize offsprings population
    functionEvaluations = parents.evaluate(parentsSize,function,nSize,gSize,hSize,functionEvaluations) # Evaluate parents
    if (penaltyMethod == 1): # Padrao?  (not apm)
        parents.sumViolations(parentsSize,gSize,hSize)
    elif (penaltyMethod == 2): # // Adaptive Penalty Method ( APM )
        parents.uniteConstraints(parentsSize,gSize,hSize)
        avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
        parents.calculateAllFitness(parentsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
    else:
        print("Penalthy method not encountered")
        sys.exit("Penalty method not encountered")
    parents.initializeEvolutionStrategy(offsprings,nSize,parentsSize,offspringsSize,globalSigma)
    while (functionEvaluations < maxFE):
        parents.sigmaSelfAdaptation(offsprings,nSize,parentsSize,generatedOffspring,globalSigma)
        offsprings.bounding(nSize,function,offspringsSize)
        functionEvaluations = offsprings.evaluate(offspringsSize,function,nSize,gSize,hSize,functionEvaluations)
        if (penaltyMethod == 1):
            offsprings.sumViolations(offspringsSize,gSize,hSize)
        else:
            offsprings.uniteConstraints(offspringsSize,gSize,hSize)
            avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
            offsprings.calculateAllFitness(offspringsSize,constraintsSize,penaltyCoefficients,avgObjFunc)
        parents.sort(offsprings,penaltyMethod)
        parents.elitismES(offsprings, parentsSize, offspringsSize, nSize,
                          gSize, hSize, constraintsSize, globalSigma,
                          esType, generatedOffspring, penaltyMethod)
        parents.printBest(parentsSize,penaltyMethod)


if __name__ == '__main__':
    # function,seed,penaltyMethod,parentsSize,nSize,generatedOffspring,
    # maxFE,crossoverProb,esType,globalSigma
    function = 1
    seed = 1
    penaltyMethod = 1
    parentsSize = 50
    nSize = 10
    generatedOffspring = 10
    maxFE = 20000
    crossoverProb = 60
    esType = 1
    globalSigma = 1

    """
    DE(function, seed, penaltyMethod, parentsSize, nSize, generatedOffspring,
       maxFE, crossoverProb, esType, globalSigma)

    GA(function, seed, penaltyMethod, parentsSize, nSize, generatedOffspring, 
    maxFE, crossoverProb, esType, globalSigma)
    
    ES(function, seed, penaltyMethod, parentsSize, nSize, generatedOffspring,
       maxFE, crossoverProb, esType, globalSigma)
    """

    # gcc -shared -o libCFunctions.so -fPIC -Wall -g cFunctions.c  # gerar .so
    cLib = ctypes.cdll.LoadLibrary("./libCFunctions.so")
    cLib.helloWorld()
    print(cLib.soma(5,2))

    # arrays
    cLib.change_array.restype = None # void
    # Seta 1 param( ponteiro para inteiro) e segundo param (inteiro)
    cLib.change_array.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_int)
    x = [2,3,4]
    y = []
    param = (ctypes.c_int * len(x))(*x)
    cLib.change_array(param, len(x))
    print(list(param))
    print(x)
    for i in range (len(param)):
        x[i] = param[i]
    print(x)


    # print(res)