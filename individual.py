#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 20 17:53:50 2018

@author: pedrohen
"""
import sys
sys.path.insert(0, "eureka")
import eureka
import argparse
import numpy as np
import operator as op
import random
from scipy.stats import truncnorm
# import ctypes
from functions import Functions

"""
n = search space dimension
sigma = ES parameter
objectiveFunction = Objective function to be minimized
g = Number of inequalities constraints
h = Number of equalities constraints
violations = Array of violations
violatonSum = Sum of violations
fitness = Fitness of each individual (for APM)
"""


class Individual(object):
    # noinspection PyUnusedLocal
    def __init__(self, n=None, objectiveFunction=None, g=None, h=None, violations=None, sigma=None, violationSum=None, fitness=None):
        self.n = [-1 for i in range(1000)] if n is None else n
        self.objectiveFunction = [-1 for i in range(1)] if objectiveFunction is None else objectiveFunction
        self.g = [-1 for i in range(1000)] if g is None else g
        self.h = [-1 for i in range(1000)] if h is None else h
        self.sigma = [-1 for i in range(1000)] if sigma is None else sigma
        self.violations = [-1 for i in range(1000)] if violations is None else violations
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
            for j in range(constraintsSize):  # gSize + hSize (violations)
                self.violations[j] = individual.violations[j]
            self.violationSum = individual.violationSum
        elif penaltyMethod == 2:  # APM
            self.fitness = individual.fitness

    def __repr__(self):
        return str(self.__dict__)


class Population(object):
    def __init__(self, popSize, nSize, function, lowerBound=None, upperBound=None):  # TODO Remover parametro truss?
        if lowerBound is not None and upperBound is not None:
            mu = (lowerBound + upperBound) / 2  # midpoint of interval
            sigma = (upperBound - lowerBound) / 4  # std which is 1/4 of the diameter
        strFunction = str(function)
        self.individuals = []
        for i in range(popSize):
            values = []
            for j in range(nSize):  # Dimension
                if strFunction[0] == "1":  # "Normal" problems (functions to be minimized)
                    if function == 11:
                        values.append(np.random.uniform(0, 10))  # value = np.random.uniform(0,10)
                    elif function == 12:
                        values.append(np.random.uniform(-5, 5))
                    elif function == 13 or function == 112 or function == 114 or function == 115:
                        values.append(np.random.uniform(-1000, 1000))
                    elif function == 14 or function == 118:
                        values.append(np.random.uniform(-50, 50))
                    elif function == 15 or function == 16:
                        values.append(np.random.uniform(-600, 600))
                    elif function == 17 or function == 18:
                        values.append(np.random.uniform(-140, 140))
                    elif function == 19 or function == 110 or function == 113:
                        values.append(np.random.uniform(-500, 500))
                    elif function == 111:
                        values.append(np.random.uniform(-100, 100))
                    elif function == 116 or function == 117:
                        values.append(np.random.uniform(-10, 10))
                elif strFunction[0] == "2":  # Truss problems
                    """
                    s = get_truncated_normal(mu, sigma, lowerBound, upperBound).rvs(nSize)
                    s = s.tolist()
                    for it in range(nSize):
                        values.append(s[it])
                    """
                    """
                    flag = True
                    while flag:
                        # val = np.random.normal(mu, sigma)
                        val = random.gauss(mu, sigma)
                        if val >= lowerBound and val <= upperBound:  # lowerBound <= val <= upperBound
                            values.append(val)
                            flag = False
                    """
                    # print(values)
                    # values.append(np.random.uniform(lowerBound, upperBound))
                else:
                    print("Function not encountered")  # sys.exit("Function not encountered")
            s = get_truncated_normal(mu, sigma, lowerBound, upperBound).rvs(nSize)
            s = s.tolist()
            for it in range(nSize):
                values.append(s[it])
            self.individuals.append(Individual(values))

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
            for j in range(constraintsSize):  # gSize + hSize (violations)
                self.individuals[idxDest].violations[j] = population.individuals[idxToBeCopy].violations[j]
            self.individuals[idxDest].violationSum = population.individuals[idxToBeCopy].violationSum
        elif penaltyMethod == 2:  # APM
            self.individuals[idxDest].fitness = population.individuals[idxToBeCopy].fitness

    def evaluate(self, popSize, function, nSize, gSize, hSize, fe, truss=None):
        strFunction = str(function)
        for i in range(popSize):
            fe = fe + 1
            if strFunction[0] == "1":  # "Normal" problems (functions to be minimized)
                if function == 11:
                    Functions.C01(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                    """
                    cLib.C01.restype = None  # void
                    # Seta 4 param( ponteiros para float) e 4 ultimos param (inteiros)
                    cLib.C01.argtypes = (ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER( ctypes.c_float),ctypes.POINTER(ctypes.c_float),  ctypes.c_int,ctypes.c_int, ctypes.c_int,ctypes.c_int)
                    nParam = (ctypes.c_float * nSize) (*(self.individuals[i].n))
                    objFuncParam = (ctypes.c_float * nSize)(*(self.individuals[i].objectiveFunction))
                    gParam = (ctypes.c_float * nSize)(*(self.individuals[i].g))
                    hParam = (ctypes.c_float * nSize)(*(self.individuals[i].h))
                    cLib.C01(nParam, objFuncParam, gParam, hParam, nSize, 1, gSize, hSize)
                    self.individuals[i].n = [x for x in nParam]  # copia array
                    self.individuals[i].objectiveFunction = [x for x in objFuncParam]
                    self.individuals[i].g = [x for x in gParam]  # copia array
                    self.individuals[i].h = [x for x in hParam]  # copia array
                    """
                elif function == 12:
                    Functions.C02(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 13:
                    Functions.C03(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 14:
                    Functions.C04(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 15:
                    Functions.C05(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 16:
                    Functions.C06(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 17:
                    Functions.C07(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 18:
                    Functions.C08(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 19:
                    Functions.C09(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 110:
                    Functions.C10(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 111:
                    Functions.C11(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 112:
                    Functions.C12(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 113:
                    Functions.C13(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 114:
                    Functions.C14(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 115:
                    Functions.C15(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 116:
                    Functions.C16(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 117:
                    Functions.C17(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
                elif function == 118:
                    Functions.C18(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
            elif strFunction[0] == "2":  # evaluate method on c++ : evalute(*vector, *values) where value is the dimension of the problem and values is the objFunction + constraints
                valuesArraySize = truss.getNumberObjectives() + truss.getNumberConstraints()  # the length will be objFunction (1) + gSize
                dimensionArray = eureka.new_doubleArray(truss.getDimension())  # creates an array
                valuesArray = eureka.new_doubleArray(valuesArraySize)  # the length will be objFunct(1) + gSize
                build_array(dimensionArray, self.individuals[i].n, truss.getDimension())  # transfers values to C++ array
                valuesList = self.individuals[i].objectiveFunction + self.individuals[i].g  # concatenates the two lists
                build_array(valuesArray, valuesList, valuesArraySize)
                truss.evaluation(dimensionArray, valuesArray)
                # truss.evaluate(dimensionArray, valuesArray)
                build_list(self.individuals[i].n, dimensionArray, 0, truss.getDimension())  # transfers values to python list
                self.individuals[i].objectiveFunction[0] = eureka.doubleArray_getitem(valuesArray, 0)
                build_list(self.individuals[i].g, valuesArray, 1, valuesArraySize)
                eureka.delete_doubleArray(dimensionArray)
                eureka.delete_doubleArray(valuesArray)
            else:
                print("Function not encountered")
                sys.exit("Function not encountered")
        # print(self.individuals)
        return fe

    def sumViolations(self, popSize, gSize, hSize):
        for i in range(popSize):
            idxG = 0
            idxH = 0
            for j in range(gSize + hSize):
                if j < gSize:
                    self.individuals[i].violations[j] = self.individuals[i].g[idxG]
                    # self.individuals[i].violations.append(self.individuals[i].g[idxG])
                    idxG = idxG + 1
                else:
                    self.individuals[i].violations[j] = self.individuals[i].h[idxH] - 0.0001  # Converts equality on inequality
                    # self.individuals[i].violations.append(self.individuals[i].h[idxH] - 0.0001)
                    idxH = idxH + 1
            # Only sums positives values. Negatives violations are not summed
            self.individuals[i].violationSum = sum(l for l in self.individuals[i].violations if l > 0)

    def selection(self, parents, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod):
        if penaltyMethod == 1:  # Not apm, 'standard' method | Compares violations and objective function
            for i in range(offspringsSize):
                idx1 = np.random.randint(0, parentsSize)
                idx2 = np.random.randint(0, parentsSize)
                # print("\n\nIndividuo Pai 1")
                # print(parents.individuals[idx1])
                # print("\n\nIndividuo Pai 2")
                # print(parents.individuals[idx2])
                if parents.individuals[idx1].violationSum < parents.individuals[idx2].violationSum:
                    self.copyIndividual(i, idx1, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)  # self.individuals[i] = parents.individuals[idx1]
                elif parents.individuals[idx1].violationSum > parents.individuals[idx2].violationSum:
                    self.copyIndividual(i, idx2, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)  # self.individuals[i] = parents.individuals[idx2]
                elif parents.individuals[idx1].objectiveFunction[0] < parents.individuals[idx2].objectiveFunction[0]:
                    self.copyIndividual(i, idx1, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)  # self.individuals[i] = parents.individuals[idx1]
                else:
                    self.copyIndividual(i, idx2, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)  # self.individuals[i] = parents.individuals[idx2]  # print("\n\nIndividuo Filho Escolhido")  # print(self.individuals[i])
        elif penaltyMethod == 2:  # APM | Compares fitness of each individual
            for i in range(offspringsSize):
                idx1 = np.random.randint(0, parentsSize)
                idx2 = np.random.randint(0, parentsSize)
                if parents.individuals[idx1].fitness < parents.individuals[idx2].fitness:
                    self.copyIndividual(i, idx1, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
                elif parents.individuals[idx1].fitness >= parents.individuals[idx2].fitness:
                    self.copyIndividual(i, idx2, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)

    def standardCrossover(self, nSize, offspringsSize):
        # Randomizes uniform value between n1 and n2
        for i in range(0, offspringsSize, 2):
            for j in range(nSize):
                n1 = self.individuals[i].n[j]
                n2 = self.individuals[i + 1].n[j]
                if n1 < n2:  # Just certificates that n1 is greater than n2
                    aux = n1
                    n1 = n2
                    n2 = aux
                self.individuals[i].n[j] = np.random.uniform(n1, n2)
                self.individuals[i + 1].n[j] = np.random.uniform(n1, n2)

    """
    eta : Crowding degree of the crossover. A high eta will produce
    children resembling to their parents, while a small eta will
    produce solutions much more different., usually 2 to 5
    """

    def sbCrossover(self, eta, nSize, offspringsSize):  # SBX (Simulated Binary Crossover)
        for i in range(0, offspringsSize, 2):
            for j in range(nSize):
                rand = np.random.rand()  # Generates random value between (0,1)
                if rand < 0.5:
                    beta = 2. * rand
                else:
                    beta = 1. / (2. * (1. - rand))

                # beta** = 1. /(eta + 1)
                beta = np.power(beta, 1. / (eta + 1))
                x1 = self.individuals[i].n[j]
                x2 = self.individuals[i + 1].n[j]
                self.individuals[i].n[j] = 0.5 * (((1 + beta) * x1) + ((1 - beta) * x2))
                self.individuals[i + 1].n[j] = 0.5 * (((1 - beta) * x1) + ((1 + beta) * x2))

                # ind1[i] = 0.5 * (((1 + beta) * x1) + ((1 - beta) * x2))  # ind2[i] = 0.5 * (((1 - beta) * x1) + ((1 + beta) * x2))

    def mutation(self, nSize, offspringsSize):
        MEAN = 0
        STD = 1
        for i in range(offspringsSize):
            for j in range(nSize):
                self.individuals[i].n[j] = self.individuals[i].n[j] + np.random.normal(MEAN, STD)

    def bounding(self, nSize, function, popSize, lowerBound=None, upperBound=None):
        strFunction = str(function)
        nMin = nMax = 0
        if strFunction[0] == "1":
            if function == 11:
                nMin = 0
                nMax = 10
            elif function == 12:
                nMin = -5.12
                nMax = 5.12
            elif function == 13 or function == 112 or function == 114 or function == 115:
                nMin = -1000
                nMax = 1000
            elif function == 14 or function == 118:
                nMin = -50
                nMax = 50
            elif function == 15 or function == 16:
                nMin = -600
                nMax = 600
            elif function == 17 or function == 18:
                nMin = -140
                nMax = 140
            elif function == 19 or function == 110 or function == 113:
                nMin = -500
                nMax = 500
            elif function == 111:
                nMin = -100
                nMax = 100
            elif function == 116 or function == 117:
                nMin = -10
                nMax = 10
        elif strFunction[0] == "2":  # Truss problems
            nMin = lowerBound
            nMax = upperBound
        else:
            print("Function not encountered")
            sys.exit("Function not encountered")

        for i in range(popSize):
            for j in range(nSize):
                if self.individuals[i].n[j] > nMax:
                    self.individuals[i].n[j] = nMax
                elif self.individuals[i].n[j] < nMin:
                    self.individuals[i].n[j] = nMin

    def sort(self, offsprings, penaltyMethod):
        if penaltyMethod == 1:  # Sort based on violatiom and then objective function
            self.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
            offsprings.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
        elif penaltyMethod == 2:  # Sort based on fitness and then objective function
            self.individuals.sort(key=op.attrgetter("fitness", 'objectiveFunction'))
            offsprings.individuals.sort(key=op.attrgetter("fitness", "objectiveFunction"))

    def elitism(self, offsprings, parentsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod):
        copyStart = 5
        i = 0
        for j in range(copyStart, parentsSize):  # J iterates on parents
            self.copyIndividual(j, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
            # self.individuals[j] = offsprings.individuals[i]
            i = i + 1

    def printBest(self, nSize, parentsSize, penaltyMethod):
        best = bestIndividual(self, parentsSize, penaltyMethod)
        if penaltyMethod == 1:  # not apm
            # print("Violation\t{:e}\tObjectiveFunction\t{:e}\t".format(best.violationSum, best.objectiveFunction[0]), end = " ")
            print("{:}\t{:}\t".format(best.violationSum, best.objectiveFunction[0]), end=" ")
            for i in range(nSize):
                print("{}\t".format(best.n[i]), end=" ")
            print("")
        elif penaltyMethod == 2:  # APM
            # print("Fitness\t{:e}\tObjectiveFunction\t{:e}\t".format(best.fitness, best.objectiveFunction[0]), end = " ")
            if best.fitness == best.objectiveFunction[0]:
                print("{:}\t{:}\t".format(best.fitness, best.objectiveFunction[0]), end=" ")
                for i in range(nSize):
                    print("{}\t".format(best.n[i]), end=" ")
                print("")

            """
            if best.fitness == best.objectiveFunction[0]:
                print("Fitness == objectiveFunction")
            """

    def printBestFO(self, parentsSize, penaltyMethod):
        best = bestIndividual(self, parentsSize, penaltyMethod)
        if penaltyMethod == 1:  # not apm
            print("Violation\t{:e}\tObjectiveFunction\t{:e}\n".format(best.violationSum, best.objectiveFunction[0]))
        elif penaltyMethod == 2:  # APM
            if best.fitness == best.objectiveFunction[0]:
                print("Fitness\t{:e}\tObjectiveFunction\t{:e}\n".format(best.fitness, best.objectiveFunction[0]))
                # if best.fitness == best.objectiveFunction[0]:
                #  print("Fitness == objectiveFunction")

    def DESelection(self, offsprings, generatedOffspring, parentsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod):
        if penaltyMethod == 1:
            j = 0
            for i in range(parentsSize):
                bestIdx = j
                while j < generatedOffspring * (i + 1):  # walks through every n offsprings of each parent
                    # get the best individual among the offsprings
                    if offsprings.individuals[j].violationSum < offsprings.individuals[bestIdx].violationSum:
                        bestIdx = j
                    elif offsprings.individuals[j].violationSum == offsprings.individuals[bestIdx].violationSum:
                        if offsprings.individuals[j].objectiveFunction[0] < offsprings.individuals[bestIdx].violationSum:
                            bestIdx = j
                    j = j + 1
                # get the best individual among the parent and the best offspring
                if offsprings.individuals[bestIdx].violationSum < self.individuals[i].violationSum:
                    self.copyIndividual(i, bestIdx, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
                elif offsprings.individuals[bestIdx].violationSum == self.individuals[i].violationSum:
                    if offsprings.individuals[bestIdx].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
                        self.copyIndividual(i, bestIdx, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
        elif penaltyMethod == 2:
            j = 0
            for i in range(parentsSize):
                bestIdx = j
                while j < generatedOffspring * (i + 1):  # walks through every n offspring of each parent
                    if offsprings.individuals[j].fitness < offsprings.individuals[bestIdx].fitness:
                        bestIdx = j  # self.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
                    elif offsprings.individuals[j].fitness == offsprings.individuals[bestIdx].fitness:
                        if offsprings.individuals[j].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:  # Offspring better than parent
                            bestIdx = j  # self.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
                    j = j + 1
                # get the best individual among the parent and the best offspring
                if offsprings.individuals[bestIdx].fitness < self.individuals[i].fitness:
                    self.copyIndividual(i, bestIdx, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
                elif offsprings.individuals[bestIdx].fitness == self.individuals[i].fitness:
                    if offsprings.individuals[bestIdx].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
                        self.copyIndividual(i, bestIdx, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)

    def DESelectionWorking(self, offsprings, parentsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod):
        if penaltyMethod == 1:
            for i in range(parentsSize):
                if offsprings.individuals[i].violationSum < self.individuals[i].violationSum:  # If offspring is better than parent
                    self.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)  # self.individuals[i] = offsprings.individuals[i] # Offspring  #  "becomes" parent
                elif offsprings.individuals[i].violationSum == self.individuals[i].violationSum:  # Compares if violationSun is equal for offspring and parents
                    if offsprings.individuals[i].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:  # Offspring better than parent
                        self.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)  # self.individuals[i] = offsprings.individuals[i]
        elif penaltyMethod == 2:
            for i in range(parentsSize):
                if offsprings.individuals[i].fitness < self.individuals[i].fitness:
                    self.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
                elif offsprings.individuals[i].fitness == self.individuals[i].fitness:
                    if offsprings.individuals[i].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:  # Offspring better than parent
                        self.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)

    def initializeEvolutionStrategy(self, offsprings, nSize, parentsSize, offspringsSize, globalSigma):
        STD = 1
        for i in range(parentsSize):  # Initialize sigma parameter for parents
            if globalSigma == 1:
                self.individuals[i].sigma[0] = STD
            else:
                for j in range(nSize):
                    self.individuals[i].sigma[j] = STD

        for i in range(offspringsSize):  # Initialize sigma parameter for offsprings
            if globalSigma == 1:
                offsprings.individuals[i].sigma[0] = STD
            else:
                for j in range(nSize):
                    offsprings.individuals[i].sigma[j] = STD

    def sigmaSelfAdaptation(self, offsprings, nSize, parentsSize, generatedOffspring, globalSigma):
        #  λ ≥ 5n, µ ≈ λ/4
        MEAN = 0
        STD = 1
        m = 0  # iterates over offsprings
        for i in range(parentsSize):
            for j in range(generatedOffspring):  # Each parents generate X offsprings
                tau = 1 / np.sqrt(nSize)
                epsilon = tau * np.random.normal(MEAN, STD)  # Normal distribution
                if globalSigma == 1:  # 1 sigma for each individual, utilizes only the first position of sigma array
                    offsprings.individuals[m].sigma[0] = self.individuals[i].sigma[0] * np.exp(epsilon) * np.exp(epsilon)
                for k in range(nSize):
                    if globalSigma == 1:  # 1 sigma for each individual, utilizes only the first position of sigma array
                        offsprings.individuals[m].n[k] = self.individuals[i].n[k] + offsprings.individuals[m].sigma[0] * np.random.normal(MEAN, STD)
                    else:
                        offsprings.individuals[m].sigma[k] = self.individuals[i].sigma[k] * np.exp(epsilon) * np.exp(epsilon)  # NOTE: Is the double product necessary?
                        offsprings.individuals[m].n[k] = self.individuals[i].n[k] + offsprings.individuals[m].sigma[k] * np.random.normal(MEAN, STD)
                m = m + 1

    def elitismES(self, offsprings, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, esType, generatedOffspring, penaltyMethod, strFunction, truss, lowerBound, upperBound):
        if esType == 0:  # Es + | Pick bests individuals among parents and offsprings
            # parents = Population(parentsSize,nSize,function) # Initialize
            # parents population
            if strFunction[0] == "2":  # truss
                aux = Population(parentsSize + offspringsSize, nSize, 2, lowerBound, upperBound)
            elif strFunction[0] == "1":
                aux = Population(parentsSize + offspringsSize, nSize, 11)

            k = 0
            for i in range(parentsSize + offspringsSize):
                if i < parentsSize:
                    aux.copyIndividual(i, i, self, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)  # aux.individuals.append(self.individuals[i])
                else:
                    aux.copyIndividual(i, k, offsprings, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
                    # aux.individuals.append(offsprings.individuals[k])
                    k = k + 1
            if penaltyMethod == 1:
                aux.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
            elif penaltyMethod == 2:
                aux.individuals.sort(key=op.attrgetter("fitness", "objectiveFunction"))

            for i in range(parentsSize):
                self.copyIndividual(i, i, aux, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)  # self.individuals[i] = aux.individuals[i]
        elif esType == 1:  # Es , | Each offspring only "exists" for 1 generation. Pick best offsprings of each parent
            j = 0
            for i in range(parentsSize):
                bestIdx = j
                # best = offsprings.individuals[j]
                while j < generatedOffspring * (i + 1):  # Goes through each offspring of each parent
                    # get the best individual among the offsprings
                    if penaltyMethod == 1:  # Not apm
                        if offsprings.individuals[j].violationSum < offsprings.individuals[bestIdx].violationSum:
                            bestIdx = j
                        elif offsprings.individuals[j].violationSum == offsprings.individuals[bestIdx].violationSum:
                            if offsprings.individuals[j].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:
                                bestIdx = j
                    elif penaltyMethod == 2:  # APM
                        if offsprings.individuals[j].fitness < offsprings.individuals[bestIdx].fitness:
                            bestIdx = j
                        elif offsprings.individuals[j].fitness == offsprings.individuals[bestIdx].fitness:
                            if offsprings.individuals[j].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:
                                bestIdx = j
                    j = j + 1
                # self.individuals[i] = best
                self.copyIndividual(i, bestIdx, offsprings, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
            if penaltyMethod == 1:
                self.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
            elif penaltyMethod == 2:
                self.individuals.sort(key=op.attrgetter("fitness", "objectiveFunction"))
        else:
            print("Es type not encountered")
            sys.exit("Es type not encountered")

    def uniteConstraints(self, parentsSize, gSize, hSize):
        for i in range(parentsSize):
            idxG = 0
            idxH = 0
            for j in range(gSize + hSize):
                if j < gSize:
                    self.individuals[i].violations[j] = self.individuals[i].g[idxG]
                    idxG = idxG + 1
                else:
                    self.individuals[i].violations[j] = np.abs(self.individuals[i].h[idxH]) - 0.0001
                    idxH = idxH + 1

    # noinspection PyUnusedLocal
    def calculatePenaltyCoefficients(self, popSize, numberOfConstraints, penaltyCoefficients, averageObjectiveFunctionValues):
        sumObjectiveFunction = 0

        # foreach candidate solution
        for i in range(popSize):
            sumObjectiveFunction = sumObjectiveFunction + self.individuals[i].objectiveFunction[0]
        # the absolute value of sumObjectiveFunction
        if sumObjectiveFunction < 0:
            sumObjectiveFunction = sumObjectiveFunction * -1
        # the average of the objective function values
        averageObjectiveFunctionValues = sumObjectiveFunction / popSize
        # the denominator of the equation of the penalty coefficients
        denominator = 0
        # the sum of the constraint violation values
        # these values are recorded to be used in the next situation

        sumViolation = []
        for l in range(numberOfConstraints):
            sumViolation.append(0)
            for i in range(popSize):
                if self.individuals[i].violations[l] > 0:  # TODO: insert 'fancy' IF here
                    sumViolation[l] = sumViolation[l] + self.individuals[i].violations[l]
                """
                else:
                    sumViolation[l] = sumViolation[l] + 0 # TODO pode
                    # colocar um pass aqui, acho
                """
            denominator = denominator + sumViolation[l] * sumViolation[l]
        # the penalty coefficients  are calculated
        for j in range(numberOfConstraints):  # TODO: insert 'fancy' IF here
            if denominator == 0:
                penaltyCoefficients[j] = 0  # penaltyCoefficients.append(0)
            else:
                penaltyCoefficients[j] = (sumObjectiveFunction / denominator) * sumViolation[j]  # penaltyCoefficients.append((sumObjectiveFunction /  # denominator) * sumViolation[j])

        return averageObjectiveFunctionValues  # returns list and value

    def calculateAllFitness(self, popSize, numberOfConstraints, penaltyCoefficients, averageObjectiveFunctionValues):
        for i in range(popSize):
            # indicates if the candidate solution is feasible
            infeasible = 0  # boolean
            # the penalty value
            penalty = 0
            for j in range(numberOfConstraints):
                if self.individuals[i].violations[j] > 0:
                    # the candidate solution is infeasible if some constraint is violated
                    infeasible = 1
                    # the penalty value is updated
                    penalty = penalty + penaltyCoefficients[j] * self.individuals[i].violations[j]

            # fitness is the sum of the objective function and penalty values
            # the candidate solution is infeasible and just the objective function value,
            # otherwise
            # TODO: fancy if here
            if infeasible:
                if self.individuals[i].objectiveFunction[0] > averageObjectiveFunctionValues:
                    self.individuals[i].fitness = self.individuals[i].objectiveFunction[0] + penalty
                else:
                    self.individuals[i].fitness = averageObjectiveFunctionValues + penalty
            else:
                self.individuals[i].fitness = self.individuals[i].objectiveFunction[0]

    def printPopulation(self, popSize):
        for i in range(popSize):
            print(self.individuals[i])

    def printDimensionsAndViolationPopulation(self, popSize, nSize):
        for i in range(popSize):
            for j in range(nSize):
                print(self.individuals[i].n[j], end="\t")
            print(self.individuals[i].violationSum)

    def calculateTrussWeight(self, parentsSize, penaltyMethod, bars):
        PropertyE = 0.1
        weight = 0
        best = bestIndividual(self, parentsSize, penaltyMethod)
        for i in range(len(bars)):  # len(bars) == len(n) (nSize)
            weight = weight + best.n[i] * PropertyE * bars[i][-1]
        return weight

    def calculateTrussWeightGrouping(self, parentsSize, penaltyMethod, bars, grouping, function):
        PropertyE = 0.1
        weight = 0
        j = 0
        best = bestIndividual(self, parentsSize, penaltyMethod)
        if function == 210:  # function 210
            for i in range(len(bars)):  # len(bars) == len(n) (nSize)
                weight = weight + best.n[i] * PropertyE * bars[i][-1]
        elif function == 260:
            for i in range(len(bars)):  # len(bars) == len(n) (nSize)
                weight = weight + best.n[grouping[i]] * PropertyE * bars[i][-1]
        elif function == 2942:
            for i in range(len(grouping)):  # len(bars) == len(n) (nSize)
                integer = int(best.n[i])
                if integer > 200:
                    integer = 200
                while j < grouping[i]:
                    weight = weight + integer * PropertyE * bars[j][-1]
                    j = j + 1
        else:
            for i in range(len(grouping)):  # len(bars) == len(n) (nSize)
                while j < grouping[i]:
                    weight = weight + best.n[i] * PropertyE * bars[j][-1]
                    j = j + 1
        return weight


def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)


def copyOneIndividual(individual, idxToBeCopy, population, nSize, objectiveFunctionSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod):
    for j in range(nSize):  # Copies n
        individual.n[j] = population.individuals[idxToBeCopy].n[j]  # self.individuals[idxDest].n[j] = population.individuals[idxToBeCopy].n[j]
    for j in range(objectiveFunctionSize):
        individual.objectiveFunction[j] = population.individuals[idxToBeCopy].objectiveFunction[j]  # self.individuals[idxDest].objectiveFunction[j] = population.individuals[idxToBeCopy].objectiveFunction[j]
    for j in range(gSize):
        individual.g[j] = population.individuals[idxToBeCopy].g[j]  # self.individuals[idxDest].g[j] = population.individuals[idxToBeCopy].g[j]
    for j in range(hSize):
        individual.h[j] = population.individuals[idxToBeCopy].h[j]  # self.individuals[idxDest].h[j] = population.individuals[idxToBeCopy].h[j]
    if penaltyMethod == 1:  # Standard
        if globalSigma == 1:  # 1 sigma for each individual, utilizes only the first position of sigma array
            individual.sigma[0] = population.individuals[idxToBeCopy].sigma[0]  # self.individuals[idxDest].sigma[0] = population.individuals[idxToBeCopy].sigma[0]
        elif globalSigma == 2:  # 1 sigma for each n of each individual
            for j in range(nSize):
                individual.sigma[j] = population.individuals[idxToBeCopy].sigma[j]  # self.individuals[idxDest].sigma[j] = population.individuals[idxToBeCopy].sigma[j]
        for j in range(constraintsSize):  # gSize + hSize (violations)
            individual.violations[j] = population.individuals[idxToBeCopy].violations[j]  # self.individuals[idxDest].violations[j] = population.individuals[idxToBeCopy].violations[j]
        individual.violationSum = population.individuals[idxToBeCopy].violationSum  # self.individuals[idxDest].violationSum = population.individuals[idxToBeCopy].violationSum
    elif penaltyMethod == 2:  # APM
        individual.fitness = population.individuals[idxToBeCopy].fitness  # self.individuals[idxDest].fitness = population.individuals[idxToBeCopy].fitness


def bestIndividual(parents, parentsSize, penaltyMethod):
    best = parents.individuals[0]
    if penaltyMethod == 1:  # not apm
        for i in range(1, parentsSize):
            if parents.individuals[i].violationSum < best.violationSum:
                best = parents.individuals[i]
            elif parents.individuals[i].violationSum == best.violationSum:
                if parents.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]:
                    best = parents.individuals[i]
    elif penaltyMethod == 2:
        for i in range(1, parentsSize):
            if parents.individuals[i].fitness < best.fitness:
                best = parents.individuals[i]
            elif parents.individuals[i].fitness == best.fitness:
                if parents.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]:
                    best = parents.individuals[i]
    return best


def crossoverProbability(crossoverProb):
    if np.random.randint(0, 100) < crossoverProb:
        return 1
    else:
        return 0


def initializeConstraints(function):
    if function == 11:
        g = 2
        h = 0
        return g, h, g + h
    elif function == 12 or function == 117:
        g = 2
        h = 1
        return g, h, g + h
    elif function == 13 or function == 19 or function == 110 or function == 111:
        g = 0
        h = 1
        return g, h, g + h
    elif function == 14:
        g = 0
        h = 4
        return g, h, g + h
    elif function == 15 or function == 16:
        g = 0
        h = 2
        return g, h, g + h
    elif function == 17 or function == 18:
        g = 1
        h = 0
        return g, h, g + h
    elif function == 112 or function == 118:
        g = 1
        h = 1
        return g, h, g + h
    elif function == 113 or function == 114 or function == 115:
        g = 3
        h = 0
        return g, h, g + h
    elif function == 116:
        g = 2
        h = 2
        return g, h, g + h
    else:
        print("Function not encountered")
        sys.exit("Function not encountered")


def initializeTruss(function):  # Initializes truss and bounds
    if function == 210:  # Truss 10 bars
        truss = eureka.F101Truss10Bar()
    elif function == 225:  # Truss 25 bars
        truss = eureka.F103Truss25Bar()
    elif function == 260:  # Truss 60 bars
        truss = eureka.F105Truss60Bar()
    elif function == 272:  # Truss 72 bars
        truss = eureka.F107Truss72Bar()
    elif function == 2942:  # Truss 942 bars
        truss = eureka.F109Truss942Bar()
    else:
        print("Function not encountered")
        sys.exit("Function not encountered")
    bounds = eureka.new_doubleddArray(truss.getDimension())
    bounds = eureka.castToDouble(truss.getBounds())
    lowerBound = eureka.doubleddArray_getitem(bounds, 0, 0)
    upperBound = eureka.doubleddArray_getitem(bounds, 0, 1)
    return truss, lowerBound, upperBound


def initializeConstraintsTrusses(truss):
    return truss.getNumberConstraints(), 0, truss.getNumberConstraints() + 0


def populationPick(solution, flags, parentsSize):
    while 1:
        contains = 0
        idx = np.random.randint(0, parentsSize)
        if idx != solution:
            for l in range(3):
                if idx == flags[l]:
                    contains = 1
                    break
            if contains == 0:
                break
    return idx


# Python function to pass values of a list to a C++ array
def build_array(a, l, size):
    for i in range(size):
        eureka.doubleArray_setitem(a, i, l[i])  # sets on array "a" at idx "i" the value of "l[i]"


# Python function to pass values of a C++ array to a list
# startIdx is from where (on array) the values will be copied
def build_list(l, a, startIdx, size):
    for i in range(size):
        l[i] = eureka.doubleArray_getitem(a, startIdx)  # Gets item of array "a" at idx "idxStart"
        startIdx = startIdx + 1


'''
function = function to be minimized  --- TIPO FUNCAO
seed = seed to be used --- SEED
penaltyMethod = Standard or Adaptive Penalty Method (APM)
popSize = Size of population
offspringsSize = λ is number of offsprings, offsprings population size
maxFE = Max number of function evaluations
probCrossover = Crossover probability (0 - 100)
esType = Type of evolution strategies ( + or , )
globalSigma = Sigma parameter for ES. If sigmaGlobal equals 1 , each individual has 1 sigma. Else,
each individual has an array of sigma.

'''
# nodes = [idxNode, x, y, z]
# bars = [idxBar, nodeIdx1, nodeIdx2]


def calculateBarLength(nodes, bars):
    for i in range(len(bars)):
        idx1 = bars[i][1] - 1  # 5
        idx2 = bars[i][2] - 1  # 3
        x = np.power(nodes[idx1][1] - nodes[idx2][1], 2)  # (x1 - x2)²
        y = np.power(nodes[idx1][2] - nodes[idx2][2], 2)
        z = np.power(nodes[idx1][3] - nodes[idx2][3], 2)
        bars[i].append(np.sqrt(x+y+z))  # appends length of each bar
    # return np.sqrt(x+y+z)


# noinspection PyTypeChecker,PyTypeChecker,PyShadowingNames
def readTrussInput(function):
    nodes = []
    bars = []
    groupSize = -1
    if function == 210:  # working
        groupSize = -1
        file = open("input10.dat", "r")
    elif function == 225:  # working
        groupSize = 8
        file = open("input25.dat", "r")
    elif function == 260:  # not working
        groupSize = 60
        file = open("input60.dat", "r")
    elif function == 272:  # working
        groupSize = 16
        file = open("input72.dat", "r")
    elif function == 2942:
        groupSize = 59
        file = open("input942.dat", "r")

    file.readline()  # Ignores first line (name of truss)
    numNodes = int(file.readline().split()[0])
    # print(numNodes)
    for i in range(numNodes):
        buffer = file.readline().split()
        for j in range(3):  # removes useless information
            buffer.pop(1)
        buffer[0] = int(buffer[0])
        for index in range(1, len(buffer)):  # saves nodes's idx and coordenates x, y and z
            buffer[index] = float(buffer[index])
        nodes.append(buffer)  # Insert the line in to the list
    ignore = int(file.readline().split()[-1])
    if function == 260:  # ignore da leitura, fora de padrao?!
        ignore = 10
    elif function == 272:
        ignore = 8
    for i in range(ignore):  # ignores the next lines
        file.readline()
    numBars = int(file.readline().split()[1])
    # print(numBars)
    # sys.exit()
    PropertyE = 10000000.
    for i in range(numBars):  # bars properties | ignores (always 10000000)
        file.readline()
    for i in range(numBars):  # which bar is connected to who
        buffer = file.readline().split()
        buffer = buffer[:-2]
        for index in range(len(buffer)):  # converts to int
            buffer[index] = int(buffer[index])
        bars.append(buffer)
    # print(*bars, sep="\n")
    # print(len(bars))
    # sys.exit("uaai")
    calculateBarLength(nodes, bars)
    # print(*bars, sep="\n")
    if function == 210:
        grouping = -1
    elif function == 225:
        grouping = [0 for i in range(groupSize)]
        grouping[0] = 1
        grouping[1] = 5
        grouping[2] = 9
        grouping[3] = 11
        grouping[4] = 13
        grouping[5] = 17
        grouping[6] = 21
        grouping[7] = 25
    elif function == 260:
        grouping = [0 for i in range(groupSize)]
        grouping[0] = 1
        grouping[1] = 2
        grouping[2] = 3
        grouping[3] = 4
        grouping[4] = 5
        grouping[5] = 6
        grouping[6] = 7
        grouping[7] = 8
        grouping[8] = 9
        grouping[9] = 10
        grouping[10] = 11
        grouping[11] = 12
        grouping[12] = 1
        grouping[13] = 2
        grouping[14] = 3
        grouping[15] = 4
        grouping[16] = 5
        grouping[17] = 6
        grouping[18] = 7
        grouping[19] = 8
        grouping[20] = 9
        grouping[21] = 10
        grouping[22] = 11
        grouping[23] = 12
        grouping[24] = 13
        grouping[25] = 14
        grouping[26] = 15
        grouping[27] = 16
        grouping[28] = 17
        grouping[29] = 18
        grouping[30] = 19
        grouping[31] = 20
        grouping[32] = 21
        grouping[33] = 22
        grouping[34] = 23
        grouping[35] = 24
        grouping[36] = 13
        grouping[37] = 14
        grouping[38] = 15
        grouping[39] = 16
        grouping[40] = 17
        grouping[41] = 18
        grouping[42] = 19
        grouping[43] = 20
        grouping[44] = 21
        grouping[45] = 22
        grouping[46] = 23
        grouping[47] = 24
        grouping[48] = 0
        grouping[49] = 0
        grouping[50] = 0
        grouping[51] = 0
        grouping[52] = 0
        grouping[53] = 0
        grouping[54] = 0
        grouping[55] = 0
        grouping[56] = 0
        grouping[57] = 0
        grouping[58] = 0
        grouping[59] = 0
    elif function == 272:
        grouping = [0 for i in range(groupSize)]
        grouping[0] = 4
        grouping[1] = 12
        grouping[2] = 16
        grouping[3] = 18
        grouping[4] = 22
        grouping[5] = 30
        grouping[6] = 34
        grouping[7] = 36
        grouping[8] = 40
        grouping[9] = 48
        grouping[10] = 52
        grouping[11] = 54
        grouping[12] = 58
        grouping[13] = 66
        grouping[14] = 70
        grouping[15] = 72
    elif function == 2942:
        grouping = [0 for i in range(groupSize)]
        grouping[0] = 2
        grouping[1] = 10
        grouping[2] = 18
        grouping[3] = 34
        grouping[4] = 46
        grouping[5] = 58
        grouping[6] = 82
        grouping[7] = 86
        grouping[8] = 90
        grouping[9] = 98
        grouping[10] = 106
        grouping[11] = 122
        grouping[12] = 130
        grouping[13] = 162
        grouping[14] = 170
        grouping[15] = 186
        grouping[16] = 194
        grouping[17] = 226
        grouping[18] = 234
        grouping[19] = 258
        grouping[20] = 270
        grouping[21] = 318
        grouping[22] = 330
        grouping[23] = 338
        grouping[24] = 342
        grouping[25] = 350
        grouping[26] = 358
        grouping[27] = 366
        grouping[28] = 382
        grouping[29] = 390
        grouping[30] = 398
        grouping[31] = 430
        grouping[32] = 446
        grouping[33] = 462
        grouping[34] = 486
        grouping[35] = 498
        grouping[36] = 510
        grouping[37] = 558
        grouping[38] = 582
        grouping[39] = 606
        grouping[40] = 630
        grouping[41] = 642
        grouping[42] = 654
        grouping[43] = 702
        grouping[44] = 726
        grouping[45] = 750
        grouping[46] = 774
        grouping[47] = 786
        grouping[48] = 798
        grouping[49] = 846
        grouping[50] = 870
        grouping[51] = 894
        grouping[52] = 902
        grouping[53] = 906
        grouping[54] = 910
        grouping[55] = 918
        grouping[56] = 926
        grouping[57] = 934
        grouping[58] = 942

    # print(grouping)
    # print(len(grouping))
    # print(*bars, sep="\n")
    # sys.exit("saiu na barra")
    return bars, grouping


def GA(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):  # Genetic Algorithm
    strFunction = str(function)
    esType = globalSigma = -1
    np.random.seed(seed)
    crossoverType = 1
    functionEvaluations = 0
    generatedOffspring = int(offspringsSize / parentsSize)
    lowerBound = upperBound = truss = 0
    if strFunction[0] == "2":  # solving trusses
        truss, lowerBound, upperBound = initializeTruss(function)
        nSize = truss.getDimension()
        gSize, hSize, constraintsSize = initializeConstraintsTrusses(truss)
        penaltyCoefficients = [-1 for i in range(constraintsSize)]
        avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
        parents = Population(parentsSize, nSize, function, lowerBound, upperBound)
        offsprings = Population(offspringsSize, nSize, function, lowerBound, upperBound)
        functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
    else:
        gSize, hSize, constraintsSize = initializeConstraints(function)  # constraintsSize is the sum of gSize and hSize
        penaltyCoefficients = [-1 for i in range(constraintsSize)]
        avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
        parents = Population(parentsSize, nSize, function)  # Initialize parents population
        offsprings = Population(offspringsSize, nSize, function)
        functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations)
    if penaltyMethod == 1:
        parents.sumViolations(parentsSize, gSize, hSize)
    elif penaltyMethod == 2:
        parents.uniteConstraints(parentsSize, gSize, hSize)
        avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
        parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
    else:
        print("Penalthy method not encountered")
        sys.exit("Penalty method not encountered")
    while functionEvaluations < maxFE:
        offsprings.selection(parents, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)  # Selection
        if crossoverProbability(crossoverProb):
            if crossoverType == 0:
                offsprings.standardCrossover(nSize, offspringsSize)  # Standard crossover
            elif crossoverType == 1:
                offsprings.sbCrossover(2, nSize, offspringsSize)  # SBX Crossover
            else:
                print("Crossover type not encountered")
                sys.exit("Crossover type not encountered")
        offsprings.mutation(nSize, offspringsSize)
        if strFunction[0] == "2":
            offsprings.bounding(nSize, function, offspringsSize, lowerBound, upperBound)
            functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
        else:
            offsprings.bounding(nSize, function, offspringsSize)
            functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
        if penaltyMethod == 1:
            offsprings.sumViolations(offspringsSize, gSize, hSize)
        elif penaltyMethod == 2:
            offsprings.uniteConstraints(offspringsSize, gSize, hSize)
            avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
            offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
        parents.sort(offsprings, penaltyMethod)
        parents.elitism(offsprings, parentsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
        parents.sort(offsprings, penaltyMethod)
        # print("FE: {}".format(functionEvaluations))
        parents.printBest(nSize, parentsSize, penaltyMethod)


def DE(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):  # Differential Evolution
    strFunction = str(function)
    crossoverProb = esType = globalSigma = -1
    np.random.seed(seed)
    CR = 0.9
    F = 0.5
    functionEvaluations = 0
    generatedOffspring = int(offspringsSize / parentsSize)
    lowerBound = upperBound = truss = 0
    if strFunction[0] == "2":  # solving trusses
        bars, grouping = readTrussInput(function)
        truss, lowerBound, upperBound = initializeTruss(function)
        nSize = truss.getDimension()
        gSize, hSize, constraintsSize = initializeConstraintsTrusses(truss)  # TODO: Juntar com o initializeConstraints?!
        penaltyCoefficients = [-1 for i in range(constraintsSize)]
        avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
        parents = Population(parentsSize, nSize, function, lowerBound, upperBound)
        offsprings = Population(offspringsSize, nSize, function, lowerBound, upperBound)
        functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
        # lol = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)  # TODO: LINHA CODIGO ROBSON
    else:  # Solving 'normal' functions
        gSize, hSize, constraintsSize = initializeConstraints(function)  # Initialize constraints
        penaltyCoefficients = [-1 for i in range(constraintsSize)]
        avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
        parents = Population(parentsSize, nSize, function)  # Initialize parents population
        offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
        functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations)  # Evaluate parents
    if penaltyMethod == 1:  # Padrao?  (not apm)
        parents.sumViolations(parentsSize, gSize, hSize)
        offsprings.sumViolations(offspringsSize, gSize, hSize)  # TODO: LINHA CODIGO ROBSON
        offsprings.printDimensionsAndViolationPopulation(offspringsSize, nSize)
        # sys.exit("saiu antes ainda")
    elif penaltyMethod == 2:  # // Adaptive Penalty Method ( APM )
        parents.uniteConstraints(parentsSize, gSize, hSize)
        avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
        parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
    else:
        print("Penalthy method not encountered")
        sys.exit("Penalty method not encountered")
    firstIteration = True
    while functionEvaluations < maxFE:
        offsprings.selection(parents, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
        if firstIteration:  # gets first set of offsprings
            # print("wtf")
            offsprings.printDimensionsAndViolationPopulation(offspringsSize, nSize)
            firstIteration = False
            # sys.exit("hue")
        flags = [-1, -1, -1]
        offspringIdx = 0
        for i in range(parentsSize):
            for k in range(generatedOffspring):
                for l in range(len(flags)):
                    flags[l] = populationPick(i, flags, parentsSize)
                R = np.random.randint(0, nSize)  # Random index
                for j in range(nSize):
                    Ri = np.random.rand()  # Generates random number between (0,1)
                    # print("offspringIdx: {}".format(offspringIdx))
                    if Ri < CR or j == R:
                        offsprings.individuals[offspringIdx].n[j] = parents.individuals[flags[0]].n[j] + F * (parents.individuals[flags[1]].n[j] - parents.individuals[flags[2]].n[j])
                    else:
                        offsprings.individuals[offspringIdx].n[j] = parents.individuals[i].n[j]
                offspringIdx = offspringIdx + 1
        if strFunction[0] == "2":
            # Modelo Robson. Primeira vez, passar x(n) e violationSum(cSum). Retornará x(n) e e violationSum(csum). Selecionar os melhores e avaliar no simulador.
            offsprings.bounding(nSize, function, offspringsSize, lowerBound, upperBound)
            functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
        else:
            offsprings.bounding(nSize, function, offspringsSize)
            functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
        if penaltyMethod == 1:  # Not apm
            offsprings.sumViolations(offspringsSize, gSize, hSize)
        elif penaltyMethod == 2:
            offsprings.uniteConstraints(offspringsSize, gSize, hSize)
            avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
            offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
        parents.DESelection(offsprings, generatedOffspring, parentsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod)
        # parents.printBest(nSize, parentsSize, penaltyMethod)
        parents.printBestFO(parentsSize, penaltyMethod)
        # weight = parents.calculateTrussWeight(parentsSize, penaltyMethod, bars)
        weight = parents.calculateTrussWeightGrouping(parentsSize, penaltyMethod, bars, grouping, function)
        print("Weigth: {:e}".format(weight))


def DERobson(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):  # Differential Evolution
    strFunction = str(function)
    crossoverProb = esType = globalSigma = -1
    np.random.seed(seed)
    CR = 0.9
    F = 0.5
    functionEvaluations = 0
    generatedOffspring = int(offspringsSize / parentsSize)
    lowerBound = upperBound = truss = 0
    if strFunction[0] == "2":  # solving trusses
        bars, grouping = readTrussInput(function)
        truss, lowerBound, upperBound = initializeTruss(function)
        nSize = truss.getDimension()
        gSize, hSize, constraintsSize = initializeConstraintsTrusses(truss)  # TODO: Juntar com o initializeConstraints?!
        penaltyCoefficients = [-1 for i in range(constraintsSize)]
        avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
        parents = Population(parentsSize, nSize, function, lowerBound, upperBound)
        offsprings = Population(offspringsSize, nSize, function, lowerBound, upperBound)
        functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
        functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
        # lol = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)  # TODO: LINHA CODIGO ROBSON
    else:  # Solving 'normal' functions
        gSize, hSize, constraintsSize = initializeConstraints(function)  # Initialize constraints
        penaltyCoefficients = [-1 for i in range(constraintsSize)]
        avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
        parents = Population(parentsSize, nSize, function)  # Initialize parents population
        offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
        functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations)  # Evaluate parents
    if penaltyMethod == 1:  # Padrao?  (not apm)
        parents.sumViolations(parentsSize, gSize, hSize)
        offsprings.sumViolations(offspringsSize, gSize, hSize)  # TODO: LINHA CODIGO ROBSON
        # TODO: TREINAR MODELO ROBSON
        # offsprings.printDimensionsAndViolationPopulation(offspringsSize, nSize)  # TODO: LINHA CODIGO ROBSON
        # sys.exit("saiu antes ainda")
    elif penaltyMethod == 2:  # // Adaptive Penalty Method ( APM )
        parents.uniteConstraints(parentsSize, gSize, hSize)
        avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
        parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
    else:
        print("Penalthy method not encountered")
        sys.exit("Penalty method not encountered")
    # print("funcEvaluation: {}".format(functionEvaluations))
    # sys.exit()
    while functionEvaluations < maxFE:
        if functionEvaluations != 250:  # Not first iteration
            offsprings.selection(parents, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
        flags = [-1, -1, -1]
        offspringIdx = 0
        for i in range(parentsSize):
            for k in range(generatedOffspring):
                for l in range(len(flags)):
                    flags[l] = populationPick(i, flags, parentsSize)
                R = np.random.randint(0, nSize)  # Random index
                for j in range(nSize):
                    Ri = np.random.rand()  # Generates random number between (0,1)
                    # print("offspringIdx: {}".format(offspringIdx))
                    if Ri < CR or j == R:
                        offsprings.individuals[offspringIdx].n[j] = parents.individuals[flags[0]].n[j] + F * (parents.individuals[flags[1]].n[j] - parents.individuals[flags[2]].n[j])
                    else:
                        offsprings.individuals[offspringIdx].n[j] = parents.individuals[i].n[j]
                offspringIdx = offspringIdx + 1
        if strFunction[0] == "2":
            # Modelo Robson. Primeira vez, passar x(n) e violationSum(cSum). Retornará x(n) e e violationSum(csum). Selecionar os melhores e avaliar no simulador.
            offsprings.bounding(nSize, function, offspringsSize, lowerBound, upperBound)
            functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
        else:
            offsprings.bounding(nSize, function, offspringsSize)
            functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
        if penaltyMethod == 1:  # Not apm
            offsprings.sumViolations(offspringsSize, gSize, hSize)
        elif penaltyMethod == 2:
            offsprings.uniteConstraints(offspringsSize, gSize, hSize)
            avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
            offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
        parents.DESelection(offsprings, generatedOffspring, parentsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod)
        # parents.printBest(nSize, parentsSize, penaltyMethod)
        parents.printBestFO(parentsSize, penaltyMethod)
        # weight = parents.calculateTrussWeight(parentsSize, penaltyMethod, bars)
        weight = parents.calculateTrussWeightGrouping(parentsSize, penaltyMethod, bars, grouping, function)
        print("Weigth: {:e}".format(weight))


def ES(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):  # Evolution Strategy
    strFunction = str(function)
    crossoverProb = -1
    np.random.seed(seed)
    functionEvaluations = 0
    generatedOffspring = int(offspringsSize / parentsSize)
    lowerBound = upperBound = truss = 0
    if strFunction[0] == "2":  # solving trusses
        truss, lowerBound, upperBound = initializeTruss(function)
        nSize = truss.getDimension()
        gSize, hSize, constraintsSize = initializeConstraintsTrusses(truss)
        penaltyCoefficients = [-1 for i in range(constraintsSize)]
        avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
        parents = Population(parentsSize, nSize, function, lowerBound, upperBound)
        offsprings = Population(offspringsSize, nSize, function, lowerBound, upperBound)
        functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
    else:
        gSize, hSize, constraintsSize = initializeConstraints(function)  # Initialize constraints
        penaltyCoefficients = [-1 for i in range(constraintsSize)]
        avgObjFunc = 0  # will be subscribed on 'calculatePenaltyCoefficients'
        parents = Population(parentsSize, nSize, function)  # Initialize parents population
        offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
        functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations)  # Evaluate parents
    if penaltyMethod == 1:  # Padrao?  (not apm)
        parents.sumViolations(parentsSize, gSize, hSize)
    elif penaltyMethod == 2:  # // Adaptive Penalty Method ( APM )
        parents.uniteConstraints(parentsSize, gSize, hSize)
        avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
        parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
    else:
        print("Penalthy method not encountered")
        sys.exit("Penalty method not encountered")
    parents.initializeEvolutionStrategy(offsprings, nSize, parentsSize, offspringsSize, globalSigma)
    while functionEvaluations < maxFE:
        parents.sigmaSelfAdaptation(offsprings, nSize, parentsSize, generatedOffspring, globalSigma)
        if strFunction[0] == "2":
            offsprings.bounding(nSize, function, offspringsSize, lowerBound, upperBound)
            functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
        else:
            offsprings.bounding(nSize, function, offspringsSize)
            functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
        if penaltyMethod == 1:
            offsprings.sumViolations(offspringsSize, gSize, hSize)
        else:
            offsprings.uniteConstraints(offspringsSize, gSize, hSize)
            avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
            offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
        parents.elitismES(offsprings, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, esType, generatedOffspring, penaltyMethod, strFunction, truss, lowerBound, upperBound)
        parents.printBest(nSize, parentsSize, penaltyMethod)


# noinspection PyShadowingNames
def algorithm(algorithm, function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):
    if algorithm == "GA":  # Genetic Algorithm
        GA(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma)
    elif algorithm == "DE":  # Differential Evolution
        #DE(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma)
        DERobson(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma)
    elif algorithm == "ES":  # Evolution Strategy
        ES(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma)
    else:
        print("Algorithm not encountered")
        sys.exit("Algorithm not encountered")


def main():
    # function,seed,penaltyMethod,parentsSize,nSize,offspringsSize,maxFE,crossoverProb,esType,globalSigma
    # ES µ ≈ λ/4
    parser = argparse.ArgumentParser(description="Evolutionary Algorithms")
    parser.add_argument("--algorithm", "-a", type=str, default="DE", help="Algorithm to be used (GA, ES or DE)")
    parser.add_argument("--function", "-f", type=int, default=210, help="Truss to be solved (10, 25, 60, 72 or 942 bars). "
                        "For the truss problem, the first digit must be 2, followed by the number of the bars in the problem. "
                        "Example: 225, is for the truss of 25 bars")
    parser.add_argument("--seed", "-s", type=int, default=1, help="Seed to be used")
    parser.add_argument("--penaltyMethod", "-p", type=int, default=1, help="Penalty method to be used. 1 for Deb Penalty or 2 for APM")
    parser.add_argument("--parentsSize", "-u", type=int, default=50, help="µ is the parental population size")  # u from µ (mi) | µ ≈ λ/4
    parser.add_argument("--nSize", "-n", type=int, default=10, help="Search space dimension")
    parser.add_argument("--offspringsSize", "-l", type=int, default=200, help="λ is number of offsprings, offsprings population size")  # l from λ (lambda) | µ ≈ λ/4
    parser.add_argument("--maxFE", "-m", type=int, default=20000, help="The max number of functions evaluations")
    parser.add_argument("--crossoverProb", "-c", type=int, default=100, help="The crossover probability [0,100]")
    parser.add_argument("--esType", "-e", type=int, default=0, help="The type of ES. 0 for ES(µ + λ) or 1 for ES(µ , λ)")
    parser.add_argument("--globalSigma", "-g", type=int, default=0, help="If the σ parameter is global or not. 1 for global σ or 0 if not")
    args = parser.parse_args()
    """
    args.algorithm = "DE"
    args.function = 11
    args.seed = 1
    args.penaltyMethod = 1
    args.parentsSize = 50  # µ | µ ≈ λ/4
    args.nSize = 10
    args.offspringsSize = 50  # λ | µ ≈ λ/4
    args.maxFE = 20000
    args.crossoverProb = 100
    args.esType = 0  # 0 Es + and 1 Es ,
    args.globalSigma = 1
    """
    # args.algorithm = "ES"
    # args.globalSigma = 0
    # args.maxFE = 2000
    # args.esType = 1
    # args.penaltyMethod = 2
    # args.function = 2942
    # args.offspringsSize = args.parentsSize
    # 10, 72 - 10 000
    # 942 - 150 000
    # Demais - 15 000
    # args.seed = 2
    # readTrussInput()
    algorithm(args.algorithm, args.function, args.seed, args.penaltyMethod, args.parentsSize, args.nSize, args.offspringsSize, args.maxFE, args.crossoverProb, args.esType, args.globalSigma)
    # print(args)


if __name__ == '__main__':
    main()
