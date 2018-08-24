#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 19:19:34 2018

@author: pedrohen
"""

import numpy as np
from functions import Functions

#ISSO NAO FUNCIONA ?!
"""
class Individuo(object):
    def __init__(self,x = [-1,-1,-1], y = [-1,-1,-1]):
        self.x = x
        self.y = y
    def __repr__(self): # Apenas para impressão
        return str(self.__dict__)


"""
import sys



class Individual(object): # ESSE FUNCIONA CACHOEIRA
    def __init__(self,n = None, objectiveFunction = None, g = None, h = None, violations = None, violationSum = None):
        self.n = [-1,-1] if n is None else n
        self.objectiveFunction = [-1 for i in range(1)] if objectiveFunction is None else objectiveFunction
        self.g = [-1 for i in range(2)] if g is None else g
        self.h = [-1] if h is None else h
        self.violations = [-1,-1] if violations is None else violations
        self.violationSum = [-1] if violationSum is None else violationSum
    def __repr__(self): # Apenas para impressão
        return str(self.__dict__)

    
class Population(object):
    """
    def __init__(self, popSize, nSize):
        # self.popSize = 50 if popSize is None else popSize
        # self.nSize = 2 if popSize is None else nSize
        self.individuals = []
        for i in range (popSize):
            values = []
            for j in range(nSize):
                values.append(np.random.uniform(0,10))
            self.individuals.append(Individual(values))
    """

    def __init__(self, alist):
        self.individuals = alist

    # @staticmethod
    def evaluate(self,popSize, nSize, gSize, hSize, fe):
        for i in range(popSize):
            fe = fe + 1
            Functions.C01(self.individuals[i].n,self.individuals[i].objectiveFunction,self.individuals[i].g,self.individuals[i].h,nSize,1,gSize,hSize)
        return fe

    def sumViolations(self,popSize,gSize):
        for i in range(popSize):
            for j in range(gSize):
                self.individuals[i].violations[j] = self.individuals[i].g[j]
            self.individuals[i].violationSum = sum(l for l in self.individuals[i].violations if l > 0)

    def selection(self,parents,popSize):
        """
        print("Offsprings(self) Impresso dentro de selection (INICIO).")
        self.printPopulation(popSize)
        print("Parents Impresso dentro de selection (INICIO).")
        parents.printPopulation(popSize)
        """
        for i in range(popSize):
            idx1 = np.random.randint(0,popSize)
            idx2 = np.random.randint(0,popSize)
            if parents.individuals[idx1].violationSum < parents.individuals[idx2].violationSum:
                self.individuals[i] = parents.individuals[idx1]
            elif parents.individuals[idx1].violationSum > parents.individuals[idx2].violationSum:
                self.individuals[i] = parents.individuals[idx2]
            elif parents.individuals[idx1].objectiveFunction[0] < parents.individuals[idx2].objectiveFunction[0]:
                self.individuals[i] = parents.individuals[idx1]
            else:
                self.individuals[i] = parents.individuals[idx2]
        """
        print("Offsprings(self) Impresso dentro de selection (FIM).")
        self.printPopulation(popSize)
        print("Parents Impresso dentro de selection (FIM).")
        parents.printPopulation(popSize)
        """

    def printPopulation(self, popSize):
        for i in range(popSize):
            print(self.individuals[i])




    def printBest(self,popSize):
        best = bestIndividual(self, popSize)
        print("Violation: {:e}\tObjectiveFunction: {:e}\n".format(best.violationSum, best.objectiveFunction[0]))

    def bounding(self,popSize,nSize):
        nMax = 10
        nMin = 0
        for i in range(popSize):
            for j in range(nSize):
                if self.individuals[i].n[j] > nMax:
                    self.individuals[i].n[j] = nMax
                elif self.individuals[i].n[j] < nMin:
                    self.individuals[i].n[j] = nMin

    def selectionDE(self,offsprings, popSize):
        for i in range(popSize):
            if offsprings.individuals[i].violationSum < self.individuals[i].violationSum:
                self.individuals[i] = offsprings.individuals[i]
            elif offsprings.individuals[i].violationSum == self.individuals[i].violationSum:
                if offsprings.individuals[i].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
                    self.individuals[i] = offsprings.individuals[i]
                    
                    
    def atribuiValorN(self,idxPop,idxN,val):
        self.individuals[idxPop].n[idxN] = val


def populationPick(solution, flags, popSize):
    while 1:
        contains = 0
        idx = np.random.randint(0,popSize)
        if idx != solution:
            for l in range(3):
                if idx == flags[l]:
                    contains = 1
                    break
            if contains == 0:
                break
    return idx

def bestIndividual(parents, popSize):
    # best = Individual()
    best = parents.individuals[0]
    #print("Imprimindo best: {}".format(best))
    #print(type(best))
    #print("UAAAAAAAAAAA")
    #sys.exit("imprimiu o tipo")
    for i in range(1, popSize):
        if (parents.individuals[i].violationSum < best.violationSum):
            best = parents.individuals[i]
        elif (parents.individuals[i].violationSum == best.violationSum):
            if (parents.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]):
                best = parents.individuals[i]
    return best


if __name__ == '__main__':
    np.random.seed(1)
    CR = 0.9
    F = 0.7
    maxFE = 200
    fe = 0
    nSize = 2
    popSize = 5
    gSize = 2
    hSize = 0
    indiP = []
    indiO = []
    for i in range(popSize):
        values = []
        for j in range(nSize):
            values.append(np.random.uniform(0, 10))
        indiP.append(Individual(values))
        
    for i in range(popSize):
        values = []
        for j in range(nSize):
            values.append(np.random.uniform(0, 10))
        indiO.append(Individual(values))

    parents = Population(indiP)
    offsprings = Population(indiO)
    #parents = Population(popSize, nSize)
    #offsprings = Population(popSize, nSize)
    fe = parents.evaluate(popSize, nSize, gSize, hSize, fe)
    fe = offsprings.evaluate(popSize,nSize,gSize,hSize, 0)
    parents.sumViolations(popSize, gSize)
    offsprings.sumViolations(popSize,gSize)



    
    print("Parents Antes do selection do offsprings. FE = {}".format(fe))
    parents.printPopulation(popSize)
    parents.printBest(popSize)
    print("Offsprings Antes do selection do offsprings. FE = {}".format(fe))
    offsprings.printPopulation(popSize)
    offsprings.printBest(popSize)
    
    while (fe < maxFE):
        offsprings.selection(parents,popSize)
        
        """
        O individuo de parents que não sofre "selecao", ou seja
        não é sorteado na funcao offsprings.selection(), é mantido corretamente
        após os cálculos realizados abaixo, em offsprings.
        """
        
        print("Parents DEPOIS do selection do offsprings. FE = {}".format(fe))
        parents.printPopulation(popSize)
        parents.printBest(popSize)
        print("Offsprings DEPOIS do selection do offsprings. FE = {}".format(fe))
        offsprings.printPopulation(popSize)
        offsprings.printBest(popSize)
                
        flags = [-1,-1,-1]
        for i in range(popSize):
            # print("entrou no for")
            for l in range(len(flags)):
                flags[l] = populationPick(i,flags,popSize)
            R = np.random.randint(0,nSize)
            # print(flags)
            for j in range(nSize):
                
                ## Conta de teste
                """
                x = parents.individuals[i].n[0]
                y = parents.individuals[i].n[1]
                z = (x + y) / 2
                offsprings.atribuiValorN(i,j,z)
                """
                
                print("Offspring {} Antes de realizar conta para n[{}]:".format(i,j))
                #print(offsprings.individuals[i])
                offsprings.printPopulation(popSize)
                print("Parents {} Antes de realizar conta para n[{}]:".format(i,j))
                #print(parents.individuals[i])
                parents.printPopulation(popSize)
                
                offsprings.individuals[i].n[j] = (parents.individuals[i].n[0] + parents.individuals[i].n[1]) / 2
                
                print("Offspring {} Após realizar conta para n[{}]:".format(i,j))
                #print(offsprings.individuals[i])
                offsprings.printPopulation(popSize)
                print("Parents {} Após realizar conta para n[{}]:".format(i,j))
                #print(parents.individuals[i])
                parents.printPopulation(popSize)
                
                """
                Ri = np.random.rand()
                if Ri < 0.9 or j == R:
                    # print("J pelo if:{}".format(j))
                    idx0 = flags[0]
                    idx1 = flags[1]
                    idx2 = flags[2]
                    # print("{} {} {}".format(idx0,idx1,idx2))
                    #individuoTemp = parents.individuals[idx0].n[j] + 0.7 * (parents.individuals[idx1].n[j] - parents.individuals[idx2].n[j])
                    #offsprings.individuals[i].n[j] = individuoTemp
                    # offsprings.individuals[i].n[j] = parents.individuals[flags[0]].n[j] + F * (parents.individuals[flags[1]].n[j] - parents.individuals[flags[2]].n[j])
                else:
                    # print("J pelo else:{}".format(j))
                    #individuoTemp1 = parents.individuals[i].n[j]
                    offsprings.individuals[i].n[j] = parents.individuals[i].n[j]
                """

        print("Parents DEPOIS do calculo do DE do offsprings. FE = {}".format(fe))
        parents.printPopulation(popSize)
        parents.printBest(popSize)
        print("Offsprings DEPOIS do calculo do DE do offsprings. FE = {}".format(fe))
        offsprings.printPopulation(popSize)
        offsprings.printBest(popSize)
        sys.exit("i'm giving up")
        
        offsprings.bounding(popSize,nSize)
        fe = offsprings.evaluate(popSize,nSize,gSize,hSize,fe)
        offsprings.sumViolations(popSize,gSize)
        parents.selectionDE(offsprings,popSize)
        print("Parents no FIM do while FE: {}".format(fe))
        parents.printPopulation(popSize)
        parents.printBest(popSize)