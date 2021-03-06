#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 20 17:53:50 2018

@author: pedrohen
"""
from scipy.linalg import eig
import sys
# sys.path.append("eureka")
sys.path.append("cec20/NoConstraints")
sys.path.append("cec20/NoConstraints/C version")
# sys.path.append("cec17NoConstraints")
# sys.path.append("..")
import eureka
import cec17NoConstraints
import cec20NoConstraints
import argparse
import numpy as np
import operator as op
import collections
from scipy.stats import truncnorm
from timeit import default_timer as timer
# import ctypes
from functions import Functions
import modelo_GPR_log
# import modelo_GPR_log_corrigido as modelo_GPR_log  # modelo novo
# import modelo_GPR_log_corrigido
import pandas as pd

"""
n = search space dimension
sigma = ES parameter
objectiveFunction = Objective function to be minimized
g = Number of inequalities constraints
h = Number of equalities constraints
violations = Array of violations
violationSum = Sum of violations
fitness = Fitness of each individual (for APM)
"""


class Individual(object):
	# noinspection PyUnusedLocal
	def __init__(self, n=None, objectiveFunction=None, g=None, h=None, violations=None, sigma=None, violationSum=None, fitness=None, arz=None):
		self.n = [-1 for i in range(1000)] if n is None else n
		self.objectiveFunction = [-1 for i in range(1)] if objectiveFunction is None else objectiveFunction
		self.g = [-1 for i in range(1000)] if g is None else g
		self.h = [-1 for i in range(1000)] if h is None else h
		self.sigma = [-1 for i in range(1000)] if sigma is None else sigma
		self.violations = [-1 for i in range(1000)] if violations is None else violations
		self.violationSum = -1 if violationSum is None else violationSum
		self.fitness = -1 if fitness is None else fitness
		self.arz = [-1 for i in range(1000)] if arz is None else arz

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

	def isObjectiveFunctionEqual(self, individual):
		if self.objectiveFunction[0] != individual.objectiveFunction[0]:
			return False
		return True

	def printIndividual(self, nSize, parentsSize, penaltyMethod, hasConstraints):
		if hasConstraints:
			if penaltyMethod == 1:  # not apm
				# print("Violation\t{:e}\tObjectiveFunction\t{:e}\t".format(best.violationSum, best.objectiveFunction[0]), end = " ")
				print("{:}\t{:}\t".format(self.violationSum, self.objectiveFunction[0]), end=" ")
				for i in range(nSize):
					print("{}\t".format(self.n[i]), end=" ")
				print("")
			elif penaltyMethod == 2:  # APM
				# print("Fitness\t{:e}\tObjectiveFunction\t{:e}\t".format(best.fitness, best.objectiveFunction[0]), end = " ")
				if self.fitness == self.objectiveFunction[0]:
					print("{:}\t{:}\t".format(self.fitness, self.objectiveFunction[0]), end=" ")
					for i in range(nSize):
						print("{}\t".format(self.n[i]), end=" ")
					print("")

				"""
				if best.fitness == best.objectiveFunction[0]:
					print("Fitness == objectiveFunction")
				"""
		else:
			print("{:}\t".format(self.objectiveFunction[0]), end=" ")
			for i in range(nSize):
				print("{}\t".format(self.n[i]), end=" ")
			print("")

	def printFO(self, parentsSize, penaltyMethod, hasConstraints):
		if hasConstraints:
			if penaltyMethod == 1:  # not apm
				# print("Violation\t{:e}\tObjectiveFunction\t{:e}\n".format(best.violationSum, best.objectiveFunction[0]))
				print("{}\t{}".format(self.violationSum, self.objectiveFunction[0]), end="\n")
			elif penaltyMethod == 2:  # APM
				if self.fitness == self.objectiveFunction[0]:
					print("Fitness\t{:e}\tObjectiveFunction\t{:e}\n".format(self.fitness, self.objectiveFunction[0]))
					# if best.fitness == best.objectiveFunction[0]:
					#  print("Fitness == objectiveFunction")
		else:
			print("{}".format(self.objectiveFunction[0]), end="\n")

	def __repr__(self):
		return str(self.__dict__)

	def __hash__(self):
		# return hash((self.n, self.objectiveFunction, self.g, self.h, self.sigma, self.violations, self.violationSum, self.fitness))
		return hash((self.objectiveFunction[0], self.violationSum))

	def __eq__(self, other):
		# return (self.n, self.objectiveFunction, self.g, self.h, self.sigma, self.violations, self.violationSum, self.fitness) == (other.n, other.objectiveFunction, other.g, other.h, other.sigma, other.violations, other.violationSum, other.fitness)
		return self.objectiveFunction[0], self.violationSum == other.objectiveFunction[0], other.violationSum

	# def __iter__(self):
	#    return (self.n, self.objectiveFunction)


class Population(object):
	def __init__(self, popSize, nSize, function, initalizeValues=True, lowerBound=None, upperBound=None):
		initalizeValues = False  # TODO initalizeValues=False nos parametros não tá funcionando?
		mu = False
		sigma = False
		if lowerBound is not None and upperBound is not None:
			# print("IF lower bound")
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
					values.append(np.random.uniform(lowerBound, upperBound))
				elif strFunction[0] == "8" or strFunction[0] == "7":  # cec17NoConstraints and cec20NoConstraints
					values.append(np.random.uniform(-100,100))  # search range is [-100,100] for all functions
				elif strFunction[0] == "9": # Functions with no constraints
					if function == 91:  # Rosenbrock Function
						# print("Initializing Rosenbrock Function")
						values.append(np.random.uniform(-5, 10))
				else:
					sys.exit("Function not encountered, exiting")
					# print("Function not encountered")  
			if initalizeValues:  # intialize values(n) with gaussian distribution
				print("Entrou aqui")
				s = get_truncated_normal(mu, sigma, lowerBound, upperBound).rvs(nSize)
				s = s.tolist()
				for it in range(nSize):
					values.append(s[it])
			# else:  # not initializing values
			# 	for it in range(nSize):
			# 		print("AH FILHA DA PUTA")
			# 		values.append(-1)
			self.individuals.append(Individual(values))

	def __iter__(self):
		#  Generator function
		for individual in self.individuals:
			yield individual

	# self,n = None, objectiveFunction = None, g = None, h =None,violations = None ,sigma = None, violationSum = None, fitness = None

	def copyIndividual(self, idxDest, idxToBeCopy, population, nSize, objectiveFunctionSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod):
		for j in range(nSize):  # Copies n and ARZ
			self.individuals[idxDest].n[j] = population.individuals[idxToBeCopy].n[j]
			self.individuals[idxDest].arz[j] = population.individuals[idxToBeCopy].arz[j]
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
			if strFunction[0] == "1":  # "Normal constrained" problems (functions to be minimized)
				if function == 11:
					Functions.C01(self.individuals[i].n, self.individuals[i].objectiveFunction, self.individuals[i].g, self.individuals[i].h, nSize, 1, gSize, hSize)
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
				else:
					print("Function not encountered\n")
					exit(1)

			elif strFunction[0] == "2":  # evaluate method on c++ : evalute(*vector, *values) where value is the dimension of the problem and values is the objFunction + constraints
				valuesArraySize = truss.getNumberObjectives() + truss.getNumberConstraints()  # the length will be objFunction (1) + gSize
				dimensionArray = eureka.new_doubleArray(truss.getDimension())  # creates an array
				valuesArray = eureka.new_doubleArray(valuesArraySize)  # the length will be objFunct(1) + gSize
				build_array(dimensionArray, self.individuals[i].n, 0, truss.getDimension(), strFunction)  # transfers values to C++ array
				valuesList = self.individuals[i].objectiveFunction + self.individuals[i].g  # concatenates the two lists
				build_array(valuesArray, valuesList, 0, valuesArraySize, strFunction)
				truss.evaluation(dimensionArray, valuesArray)
				# truss.evaluate(dimensionArray, valuesArray)
				build_list(self.individuals[i].n, dimensionArray, 0, truss.getDimension(), strFunction)  # transfers values to python list
				self.individuals[i].objectiveFunction[0] = eureka.doubleArray_getitem(valuesArray, 0)
				build_list(self.individuals[i].g, valuesArray, 1, valuesArraySize, strFunction)
				eureka.delete_doubleArray(dimensionArray)
				eureka.delete_doubleArray(valuesArray)
			
			elif strFunction[0] == "7":	# no constraints cec2020 competition functions
				funcNum = int(strFunction[1:])	# gets all numbers except the first, in this case, the 7 is ignored
				# 	creates C arrays
				xArray = cec20NoConstraints.new_doubleArray(nSize)	# creates an array with nSize Dimension
				objFuncArray = cec20NoConstraints.new_doubleArray(1)	# creates an array with 1 dimension, for(objective function)
				# print(self.individuals[i].n)
				build_array(xArray, self.individuals[i].n, 0, nSize, strFunction)	# transfers values form list to C array
				build_array(objFuncArray, self.individuals[i].objectiveFunction, 0, 1, strFunction)
				# test_func(x, f, dimension, population_size,func_num);
				cec20NoConstraints.cec20_test_func(xArray, objFuncArray, nSize, 1, funcNum)
				# transfers values back to individuals
				build_list(self.individuals[i].n, xArray, 0, nSize, strFunction)  # transfers project variables
				self.individuals[i].objectiveFunction[0] = cec20NoConstraints.doubleArray_getitem(objFuncArray, 0)  # transfers objective function
				# clean mess
				cec20NoConstraints.delete_doubleArray(xArray)
				cec20NoConstraints.delete_doubleArray(objFuncArray)
			elif strFunction[0] == "8":  # no constraints cec2017 competition functions
				funcNum = int(strFunction[1:])  # gets all numbers except the first, in this case , the 8 is ignored
				#  creates C++ arrays
				xArray = cec17NoConstraints.new_doubleArray(nSize)  # creates an array with nSize dimension
				objFuncArray = cec17NoConstraints.new_doubleArray(1)  # creates an array with 1 dimension (for objective function)
				# tranfers values from individuals to C++ array for functions be evaluated
				# print(self.individuals[i].n)
				
				build_array(xArray, self.individuals[i].n, 0, nSize, strFunction)  # transfers values from list to C++ array
				build_array(objFuncArray, self.individuals[i].objectiveFunction, 0, 1, strFunction)
				# cec17NoConstraints.doubleArray_setitem(objFuncArray, 0, self.individuals[i].objectiveFunction[0])

				"""
				print("individuals[{}].n before evaluating: ".format(i))
				print(self.individuals[i].n)
				print("individuals[{}].objectiveFunction[0] before evaluating: ".format(i))
				print(self.individuals[i].objectiveFunction[0])

				print("Array with values of individuals before evaluating(C structure):")
				cec17NoConstraints.printDoubleArray(xArray, nSize)
				print("Array with objectiveFunction of individuals before evaluating(C structure):")
				cec17NoConstraints.printDoubleArray(objFuncArray, 1)
				"""

				# EVALUATE FUNCTION
				# parameters of cec17_test_func(x, objFunc, dimension, population_size, number of function)

				cec17NoConstraints.cec17_test_func(xArray, objFuncArray, nSize, 1, funcNum)

				"""
				print("Array with values of individuals after evaluating(C structure):")
				cec17NoConstraints.printDoubleArray(xArray, nSize)
				print("Array with objectiveFunction of individuals after evaluating(C structure):")
				cec17NoConstraints.printDoubleArray(objFuncArray, 1)
				"""

				# transfers values back to individuals
				build_list(self.individuals[i].n, xArray, 0, nSize, strFunction)  # transfers project variables
				self.individuals[i].objectiveFunction[0] = cec17NoConstraints.doubleArray_getitem(objFuncArray, 0)  # transfers objective function
				"""
				print("individuals[{}].n after evaluating: ".format(i))
				print(self.individuals[i].n)
				print("individuals[{}].objectiveFunction[0] after evaluating: ".format(i))
				print(self.individuals[i].objectiveFunction[0])


				print("Array with values of individuals after evaluating(C structure):")
				cec17NoConstraints.printDoubleArray(xArray, nSize)
				print("Array with objectiveFunction of individuals after evaluating(C structure):")
				cec17NoConstraints.printDoubleArray(objFuncArray, 1)
				"""

				# clean mess
				cec17NoConstraints.delete_doubleArray(xArray)
				cec17NoConstraints.delete_doubleArray(objFuncArray)
				# print("Executing 2017 cec competitions functions (noConstraints)")

				"""
				funcNum = int(strFunction[1:])  # gets all numbers exepect the first, in this case , the 8 is ignored
				#  creates C++ arrays
				xArray = cec17NoConstraints.new_longDoubleArray(nSize)  # creates an array with nSize dimension
				objFuncArray = cec17NoConstraints.new_longDoubleArray(1)  # creates an array with 1 dimension (for objective function)
				# tranfers values from individuals to C++ array for functions be evaluated
				build_array(xArray, self.individuals[i].n, 0, nSize, strFunction)  # transfers values from list to C++ array
				build_array(objFuncArray, self.individuals[i].objectiveFunction, 0, 1, strFunction)
				# cec17NoConstraints.doubleArray_setitem(objFuncArray, 0, self.individuals[i].objectiveFunction[0])
				print("individuals[{}].n before evaluating: ".format(i))
				print(self.individuals[i].n)
				print("individuals[{}].objectiveFunction[0] before evaluating: ".format(i))
				print(self.individuals[i].objectiveFunction[0])

				print("Array with values of individuals before evaluating(C structure):")
				cec17NoConstraints.printLongDoubleArray(xArray, nSize)
				print("Array with objectiveFunction of individuals before evaluating(C structure):")
				cec17NoConstraints.printLongDoubleArray(objFuncArray, 1)

				# EVALUATE FUNCTION
				# parameters of cec17_test_func(x, objFunc, dimension, population_size, number of function)

				cec17NoConstraints.cec17_test_func(xArray, objFuncArray, nSize, 1, funcNum)

				# transfers values back to individuals
				build_list(self.individuals[i].n, xArray, 0, nSize, strFunction)  # transfers project variables
				self.individuals[i].objectiveFunction[0] = cec17NoConstraints.longDoubleArray_getitem(objFuncArray, 0)  # transfers objective function
				print("individuals[{}].n after evaluating: ".format(i))
				print(self.individuals[i].n)
				print("individuals[{}].objectiveFunction[0] after evaluating: ".format(i))
				print(self.individuals[i].objectiveFunction[0])

				print("Array with values of individuals after evaluating(C structure):")
				cec17NoConstraints.printLongDoubleArray(xArray, nSize)
				print("Array with objectiveFunction of individuals after evaluating(C structure):")
				cec17NoConstraints.printLongDoubleArray(objFuncArray, 1)

				# clean mess
				cec17NoConstraints.delete_longDoubleArray(xArray)
				cec17NoConstraints.delete_longDoubleArray(objFuncArray)
				# print("Executing 2017 cec competitions functions (noConstraints)")
				"""

			elif strFunction[0] == "9":  # functions without constraints
				if function == 91:  # RosenbrockFunction
					sumRosen = 0
					for k in range(nSize - 1):
						sumRosen = sumRosen + 100*np.power((self.individuals[i].n[k+1] - np.power(self.individuals[i].n[k], 2)), 2) + np.power((1 - self.individuals[i].n[k]), 2)
					self.individuals[i].objectiveFunction[0] = sumRosen
			else:
				print("Function not encountered")
				sys.exit("Function not encountered")
		# print(self.individuals)
		"""
		# Function that evaluates all population in one take
		if strFunction[0] == "8":
			print("Evaluating all population for once")
			funcNum = int(strFunction[1:])  # gets all numbers exepect the first, in this case , the 8 is ignored
			xArray = cec17NoConstraints.new_doubleArray(nSize*popSize)  # creates an array with nSize*popSize dimension
			objFuncArray = cec17NoConstraints.new_doubleArray(popSize)  # creates an array with popSize dimension (for objectives functions)
			for j in range(popSize):  # xArray and objFuncArray has values of the entire population
				build_array(xArray, self.individuals[j].n, j*nSize, nSize, strFunction)
				build_array(objFuncArray, self.individuals[j].objectiveFunction, j, 1, strFunction)

			# Evaluates all population
			cec17NoConstraints.cec17_test_func(xArray, objFuncArray, nSize, popSize, funcNum)

			for j in range(popSize):
				self.individuals[j].objectiveFunction[0] = cec17NoConstraints.doubleArray_getitem(objFuncArray, j)
				# build_list(self.individuals[j].objectiveFunction)

			cec17NoConstraints.delete_doubleArray(xArray)
			cec17NoConstraints.delete_doubleArray(objFuncArray)
		"""
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

	def selection(self, parents, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod, hasConstraints):
		if hasConstraints:
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
		else:
			for i in range(offspringsSize):
				idx1 = np.random.randint(0, parentsSize)
				idx2 = np.random.randint(0, parentsSize)
				if parents.individuals[idx1].objectiveFunction[0] < parents.individuals[idx2].objectiveFunction[0]:
					self.copyIndividual(i, idx1, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)  # self.individuals[i] = parents.individuals[idx1]
				else:
					self.copyIndividual(i, idx2, parents, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)  # self.individuals[i] = parents.individuals[idx2]  # print("\n\nIndividuo Filho

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
		elif strFunction[0] == "8" or strFunction[0] == "7":  # cec2017NoConstraints or cec2020NoConstraints problems
			nMin = -100
			nMax = 100
		elif strFunction[0] == "9":
			if function == 91:  # Rosenbrock Function
				nMin = -5
				nMax = 10
		else:
			print("Function not encountered")
			sys.exit("Function not encountered")

		for i in range(popSize):
			for j in range(nSize):
				if self.individuals[i].n[j] > nMax:
					self.individuals[i].n[j] = nMax
				elif self.individuals[i].n[j] < nMin:
					self.individuals[i].n[j] = nMin

	def sort(self, offsprings, penaltyMethod, hasConstraints):
		if hasConstraints:
			if penaltyMethod == 1:  # Sort based on violatiom and then objective function
				self.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
				offsprings.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
			elif penaltyMethod == 2:  # Sort based on fitness and then objective function
				self.individuals.sort(key=op.attrgetter("fitness", 'objectiveFunction'))
				offsprings.individuals.sort(key=op.attrgetter("fitness", "objectiveFunction"))
		else:
			self.individuals.sort(key=op.attrgetter("objectiveFunction"))
			offsprings.individuals.sort(key=op.attrgetter("objectiveFunction"))

	def sortPopulation(self, penaltyMethod):
		if penaltyMethod == 1:  # Sort based on violatiom and then objective function
			self.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
		elif penaltyMethod == 2:  # Sort based on fitness and then objective function
			self.individuals.sort(key=op.attrgetter("fitness", 'objectiveFunction'))

	def elitism(self, offsprings, parentsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod):
		copyStart = 5
		i = 0
		for j in range(copyStart, parentsSize):  # J iterates on parents
			self.copyIndividual(j, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
			# self.individuals[j] = offsprings.individuals[i]
			i = i + 1

	def printBest(self, nSize, parentsSize, penaltyMethod, hasConstraints):
		best = bestIndividual(self, parentsSize, penaltyMethod, hasConstraints)
		if hasConstraints:
			if penaltyMethod == 1:  # DEB
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
		else:
			print("{:}\t".format(best.objectiveFunction[0]), end=" ")
			for i in range(nSize):
				print("{}\t".format(best.n[i]), end=" ")
			print("")

	def printBestFO(self, parentsSize, penaltyMethod, hasConstraints):
		best = bestIndividual(self, parentsSize, penaltyMethod, hasConstraints)
		if hasConstraints:
			if penaltyMethod == 1:  # not apm
				# print("Violation\t{:e}\tObjectiveFunction\t{:e}\n".format(best.violationSum, best.objectiveFunction[0]))
				print("{}\t{}".format(best.violationSum, best.objectiveFunction[0]), end="\n")
			elif penaltyMethod == 2:  # APM
				if best.fitness == best.objectiveFunction[0]:
					print("Fitness\t{:e}\tObjectiveFunction\t{:e}\n".format(best.fitness, best.objectiveFunction[0]))
					# if best.fitness == best.objectiveFunction[0]:
					#  print("Fitness == objectiveFunction")
		else:
			print("{}".format(best.objectiveFunction[0]), end="\t")

	def DESelection(self, offsprings, generatedOffspring, parentsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod, hasConstraints):
		if hasConstraints:
			if penaltyMethod == 1:
				j = 0
				for i in range(parentsSize):
					bestIdx = j
					while j < generatedOffspring * (i + 1):  # walks through every n offsprings of each parent
						# get the best individual among the offsprings
						if offsprings.individuals[j].violationSum < offsprings.individuals[bestIdx].violationSum:
							bestIdx = j
						elif offsprings.individuals[j].violationSum == offsprings.individuals[bestIdx].violationSum:
							if offsprings.individuals[j].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:
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
		else:  # no constraints
			j = 0
			for i in range(parentsSize):
				bestIdx = j
				while j < generatedOffspring * (i + 1):  # walks through every n offsprings of each parent
					# get the best individual among the offsprings
					if offsprings.individuals[j].violationSum < offsprings.individuals[bestIdx].violationSum:
						bestIdx = j
					elif offsprings.individuals[j].violationSum == offsprings.individuals[bestIdx].violationSum:
						if offsprings.individuals[j].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:
							bestIdx = j
					j = j + 1
				if offsprings.individuals[bestIdx].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
					self.copyIndividual(i, bestIdx, offsprings, nSize, 1, -1, -1, -1, -1, penaltyMethod)

	def DESelectionModelGPR(self, bestsIndividuals, parentsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod):
		if penaltyMethod == 1:
			for i in range(parentsSize):
				# get the best individual among the parent and the best offspring
				if bestsIndividuals.individuals[i].violationSum < self.individuals[i].violationSum:
					self.copyIndividual(i, i, bestsIndividuals, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
				elif bestsIndividuals.individuals[i].violationSum == self.individuals[i].violationSum:
					if bestsIndividuals.individuals[i].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
						self.copyIndividual(i, i, bestsIndividuals, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
		elif penaltyMethod == 2:
			for i in range(parentsSize):
				# get the best individual among the parent and the best offspring
				if bestsIndividuals.individuals[i].fitness < self.individuals[i].fitness:
					self.copyIndividual(i, i, bestsIndividuals, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
				elif bestsIndividuals.individuals[i].fitness == self.individuals[i].fitness:
					if bestsIndividuals.individuals[i].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
						self.copyIndividual(i, i, bestsIndividuals, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)

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

	# noinspection PyUnboundLocalVariable
	def elitismES(self, offsprings, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, esType, generatedOffspring, penaltyMethod, strFunction, truss, lowerBound, upperBound, hasConstraints):
		if hasConstraints:
			if esType == 0:  # Es + | Pick bests individuals among parents and offsprings
				# parents = Population(parentsSize,nSize,function) # Initialize
				# parents population
				# print("LB{} and UP{}".format(lowerBound, upperBound))
				if strFunction[0] == "2":  # truss
					aux = Population(parentsSize + offspringsSize, nSize, 2, True, lowerBound, upperBound)
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
						if penaltyMethod == 1:  # Deb
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
			elif esType == 2: # CMAES, temporary?!
				for i in range(parentsSize):  # copies all individuals from offsprings to parents
					self.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)

				if penaltyMethod == 1:
					self.individuals.sort(key=op.attrgetter("objectiveFunction"))
				elif penaltyMethod == 2:
					self.individuals.sort(key=op.attrgetter("objectiveFunction"))
			else:
				print("Es type not encountered")
				sys.exit("Es type not encountered")

		else:  # Function has no constraints
			# Elitism for ES Type 1 ( ES ,), will work just for CMAES, for testing
			print("ta caidnoa qui neh")
			
			for i in range(parentsSize):  # copies all individuals from offsprings to parents
				self.copyIndividual(i, i, offsprings, nSize, 1, 0, 0, 0, 0, 0)

			self.individuals.sort(key=op.attrgetter("objectiveFunction"))

	def elitismNoConstraints(self, offsprings, parentsSize, offspringsSize, nSize):
		for i in range(parentsSize):  # copies all individuals from offsprings to parents
			self.copyIndividual(i, i, offsprings, nSize, 1, 0, 0, 0, 0, 0)

		self.individuals.sort(key=op.attrgetter("objectiveFunction"))

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

	def printViolationSum(self, popSize):
		for i in range(popSize):
			print(self.individuals[i].violationSum, end="\t")

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

	def calculateTrussWeightGroupingBest(self, parentsSize, penaltyMethod, bars, grouping, function):
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

	def calculateTrussWeightGrouping(self, popSize, bars, grouping, function):
		PropertyE = 0.1
		# weight = 0
		# j = 0
		for it in range(popSize):
			weight = 0
			j = 0
			if function == 210:  # function 210
				for i in range(len(bars)):  # len(bars) == len(n) (nSize)
					weight = weight + self.individuals[it].n[i] * PropertyE * bars[i][-1]
			elif function == 260:
				for i in range(len(bars)):  # len(bars) == len(n) (nSize)
					weight = weight + self.individuals[it].n[grouping[i]] * PropertyE * bars[i][-1]
			elif function == 2942:
				for i in range(len(grouping)):  # len(bars) == len(n) (nSize)
					integer = int(self.individuals[it].n[i])
					if integer > 200:
						integer = 200
					while j < grouping[i]:
						weight = weight + integer * PropertyE * bars[j][-1]
						j = j + 1
			else:  # 25 and 72 bars trusses
				for i in range(len(grouping)):  # len(bars) == len(n) (nSize)
					while j < grouping[i]:
						weight = weight + self.individuals[it].n[i] * PropertyE * bars[j][-1]
						j = j + 1
			self.individuals[it].objectiveFunction[0] = weight

	def makePopulationList(self, popSize, nSize, booleanViolSum, penaltyMethod):  # big list of all individuals and vilationSum for model_gpr
		populationList = []
		for i in range(popSize):
			for j in range(nSize):
				populationList.append(self.individuals[i].n[j])  # copies n values
			if booleanViolSum:
				if penaltyMethod == 1:  # deb
					populationList.append(self.individuals[i].violationSum)  # copies csum
				elif penaltyMethod == 2:  # apm
					populationList.append(self.individuals[i].fitness)  # copies fitness
		# print("population list {}".format(populationList))
		return populationList

	def evaluateOnGPRModel(self, popSize, nSize, modelGPR, bars, grouping, function, penaltyMethod):
		populationList = self.makePopulationList(popSize, nSize, False, penaltyMethod)  # not passing csum to modelGpr
		df = pd.DataFrame(np.array(populationList).reshape(popSize, nSize))  # converts do dataframe
		# print(type(df))
		csums = modelGPR.predict(np.log(df))  # predict cSums
		# print(type(csums))
		# sys.exit("printando pop list")
		csums = csums.tolist()  # converts np.array to list
		# print(type(csums))
		# print(len(csums))
		# print(csums)
		for i in range(popSize):  # copies csums predicted values to offsprings individuals
			if penaltyMethod == 1:  # deb
				self.individuals[i].violationSum = csums[i]
			elif penaltyMethod == 2:  # apm
				self.individuals[i].fitness = csums[i]
		# if penaltyMethod == 1:  # calculates objective function only if deb penalization
		self.calculateTrussWeightGrouping(popSize, bars, grouping, function)  # calculates weight(objectiveFunction) for each individual

	def evaluateOnlyBestsFromModelGPR(self, bestsIndividuals, generatedOffspring, offspringsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod, function, functionEvaluations, truss):
		parentsSize = int(offspringsSize / generatedOffspring)
		# print(parentsSize)
		# bestsIndividuals = Population(parentsSize, nSize, 2, False, -1, -1)  # last two parameters is lowerBound and upperBound, will be subscribed
		auxIdxBestsIndividuals = []  # saves the indexes of bests individuals of offsprings
		# self.individuals[j].violationSum
		if penaltyMethod == 1:
			j = 0
			for i in range(parentsSize):  # offspringsSIze
				bestIdx = j  # gets idx of best offspring of each parent
				while j < generatedOffspring * (i + 1):  # walks through every n offsprings of each parent
					# get the best individual among the offsprings
					# print(j)
					if self.individuals[j].violationSum < self.individuals[bestIdx].violationSum:
						bestIdx = j
					elif self.individuals[j].violationSum == self.individuals[bestIdx].violationSum:
						if self.individuals[j].objectiveFunction[0] < self.individuals[bestIdx].objectiveFunction[0]:
							bestIdx = j
					j = j + 1
				# aux.copyIndividual(i, i, self, nSize, 1, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)  # aux.individuals.append(self.individuals[i])
				# auxIdxBestsIndividuals.append(bestIdx)
				bestsIndividuals.copyIndividual(i, bestIdx, self, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)  # adds best offsprings to be evaluate on simulator
			functionEvaluations = bestsIndividuals.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)  # evaluate on simulator
			"""
			for i in range(parentsSize):  # copies best individuals to offsprings
				idxDest = auxIdxBestsIndividuals[i]  # gets the index to where the evaluate offspring will be copied
				self.copyIndividual(idxDest, i, bestsIndividuals, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
			return functionEvaluations
			"""
			return bestsIndividuals, functionEvaluations

		elif penaltyMethod == 2:
			j = 0
			for i in range(parentsSize):
				bestIdx = j
				while j < generatedOffspring * (i + 1):  # walks through every n offspring of each parent
					if self.individuals[j].fitness < self.individuals[bestIdx].fitness:
						bestIdx = j  # self.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
					elif self.individuals[j].fitness == self.individuals[bestIdx].fitness:
						if self.individuals[j].objectiveFunction[0] < self.individuals[bestIdx].objectiveFunction[0]:  # Offspring better than parent
							bestIdx = j  # self.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
					j = j + 1
				bestsIndividuals.copyIndividual(i, bestIdx, self, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)  # adds best offsprings to be evaluate on simulator
			functionEvaluations = bestsIndividuals.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)  # evaluate on simulator
			return bestsIndividuals, functionEvaluations


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


def bestIndividual(parents, parentsSize, penaltyMethod, hasConstraints):
	best = parents.individuals[0]
	if hasConstraints:
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
	else:
		for i in range(1, parentsSize):
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
def build_array(a, l, startIdx, size, strFunction):
	if strFunction[0] == "2":  # truss problems, utilizes eureka
		for i in range(size):
			# eureka.doubleArray_setitem(a, i, l[i])  # sets on array "a" at idx "i" the value of "l[i]"
			eureka.doubleArray_setitem(a, startIdx, l[i])  # sets on array "a" at idx "i" the value of "l[i]"
			startIdx = startIdx + 1
	elif strFunction[0] == "7":  # cec2020NoConstraints problems, utilizes cec20NoConstraints
		for i in range(size):
			cec20NoConstraints.doubleArray_setitem(a, startIdx, l[i])
			startIdx = startIdx + 1
	elif strFunction[0] == "8":  # cec2017NoConstraints problems, utilizes cec17NoConstraints
		for i in range(size):
			# cec17NoConstraints.doubleArray_setitem(a, i, l[i])  # sets on array "a" at idx "i" the value of "l[i]"
			# print("transfering element l[{}] = {}".format(i, l[i]))
			cec17NoConstraints.doubleArray_setitem(a, startIdx, l[i])
			startIdx = startIdx + 1


# Python function to pass values of a C++ array to a list
# startIdx is from where (on array) the values will be copied
def build_list(l, a, startIdx, size, strFunction):
	if strFunction[0] == "2":  # truss problems, utilizes eureka
		for i in range(size):
			l[i] = eureka.doubleArray_getitem(a, startIdx)  # Gets item of array "a" at idx "idxStart"
			startIdx = startIdx + 1
	elif strFunction[0] == "7":# cec2020NoConstraints problems, utilizes cec20NoConstraints
		for i in range(size):
			l[i] = cec20NoConstraints.doubleArray_getitem(a, startIdx)
			startIdx = startIdx + 1
	elif strFunction[0] == "8":  # cec2017NoConstraints problems, utilizes cec17NoConstraints
		for i in range(size):
			l[i] = cec17NoConstraints.doubleArray_getitem(a, startIdx)
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
	file.close()
	return bars, grouping


def compareLists(s, t):
	return collections.Counter(s) == collections.Counter(t)


def adjustGPRModel(offsprings, offspringsSize, nSize, penaltyMethod, functionEvaluations):
	populationList = offsprings.makePopulationList(offspringsSize, nSize, True, penaltyMethod)  # true stands for adding violationSum to populationList
	# print("offspringsSize {}".format(offspringsSize))
	# print("population list {}".format(len(populationList)))
	# print("Imprimindo offspringsSize(popTrainingIdx) {} ".format(offspringsSize))
	modelGPR, crossVMean, crossVStd = modelo_GPR_log.surGPR_training(populationList, offspringsSize, nSize + 1, functionEvaluations)  # TODO: TREINAR MODELO ROBSON
	# modelGPR, crossVMean, crossVStd = modelo_GPR_log.surGPR_training(populationList, off)
	# modelGPR, crossVMean, crossVstd = modelo_GPR_log_corrigido.surGPR_training(populationList, offspringsSize, nSize + 1)
	# print("Modelo atualizado")
	return modelGPR, crossVMean, crossVStd


def GA(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):  # Genetic Algorithm
	strFunction = str(function)
	esType = globalSigma = -1
	np.random.seed(seed)
	hasConstraints = False
	if strFunction[0] == "1" or strFunction[0] == "2":
		hasConstraints = True
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
		parents = Population(parentsSize, nSize, function, True, lowerBound, upperBound)
		offsprings = Population(offspringsSize, nSize, function, True, lowerBound, upperBound)
		functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
	else:
		if hasConstraints:
			gSize, hSize, constraintsSize = initializeConstraints(function)  # constraintsSize is the sum of gSize and hSize
			penaltyCoefficients = [-1 for i in range(constraintsSize)]
			avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
			parents = Population(parentsSize, nSize, function)  # Initialize parents population
			offsprings = Population(offspringsSize, nSize, function)
			functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations)
		else:
			parents = Population(parentsSize, nSize, function)  # Initialize parents population
			offsprings = Population(offspringsSize, nSize, function)
			functionEvaluations = parents.evaluate(parentsSize, function, nSize, -1, -1, functionEvaluations)
	if hasConstraints:
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
		if hasConstraints:
			offsprings.selection(parents, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod, hasConstraints)  # Selection
		else:
			offsprings.selection(parents, parentsSize, offspringsSize, nSize, -1, -1, -1, globalSigma, penaltyMethod, hasConstraints)  # Selection
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
			if hasConstraints:
				offsprings.bounding(nSize, function, offspringsSize)
				functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
			else:
				offsprings.bounding(nSize, function, offspringsSize)
				functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, -1, -1, functionEvaluations)
		if hasConstraints:
			if penaltyMethod == 1:
				offsprings.sumViolations(offspringsSize, gSize, hSize)
			elif penaltyMethod == 2:
				offsprings.uniteConstraints(offspringsSize, gSize, hSize)
				avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
				offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
		parents.sort(offsprings, penaltyMethod, hasConstraints)
		if hasConstraints:
			parents.elitism(offsprings, parentsSize, nSize, gSize, hSize, constraintsSize, globalSigma, penaltyMethod)
		else:
			parents.elitism(offsprings, parentsSize, nSize, -1, -1, -1, globalSigma, penaltyMethod)
		parents.sort(offsprings, penaltyMethod, hasConstraints)
		# print("FE: {}".format(functionEvaluations))
		# parents.printBest(nSize, parentsSize, penaltyMsethod, hasConstraints)
		# parents.printBestFO(parentsSize, penaltyMethod, hasConstraints)
		parents.printBestFO(parentsSize, penaltyMethod, hasConstraints)


def DE(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):  # Differential Evolution
	strFunction = str(function)
	crossoverProb = esType = globalSigma = -1
	np.random.seed(seed)
	hasConstraints = False
	if strFunction[0] == "1" or strFunction[0] == "2":
		hasConstraints = True
	if strFunction[0] == "7":	# cec2020 no constraints competition
		if(nSize == 5):	
			maxFE = 50000
		elif (nSize == 10):
			maxFE = 1000000
		elif (nSize == 15):
			maxFE = 3000000
		elif (nSize == 20):
			maxFE = 10000000
		else:
			print("Dimension not defined")
	if strFunction[0] == "8":
		maxFE = nSize * 10000
	CR = 0.9
	F = 0.6
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
		parents = Population(parentsSize, nSize, function, True, lowerBound, upperBound)
		offsprings = Population(offspringsSize, nSize, function, True, lowerBound, upperBound)
		functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)  # TODO Verificar
		functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
	else:  # Solving 'normal' functions
		if hasConstraints:
			gSize, hSize, constraintsSize = initializeConstraints(function)  # Initialize constraints
			penaltyCoefficients = [-1 for i in range(constraintsSize)]
			avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
			parents = Population(parentsSize, nSize, function)  # Initialize parents population
			offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
			functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations)  # Evaluate parents
			# functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
		else:
			parents = Population(parentsSize, nSize, function)  # Initialize parents population
			offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
			functionEvaluations = parents.evaluate(parentsSize, function, nSize, -1, -1, functionEvaluations)  # Evaluate parents
			functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, -1, -1, functionEvaluations)

	if hasConstraints:
		if penaltyMethod == 1:  # Padrao?  (not apm)
			# parents.sumViolations(parentsSize, gSize, hSize)
			offsprings.sumViolations(offspringsSize, gSize, hSize)
			offsprings.sortPopulation(penaltyMethod)
			for i in range(parentsSize):
				parents.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
			# offsprings.printDimensionsAndViolationPopulation(offspringsSize, nSize)
			# sys.exit("saiu antes ainda")
		elif penaltyMethod == 2:  # // Adaptive Penalty Method ( APM )
			parents.uniteConstraints(parentsSize, gSize, hSize)
			avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
			parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
		else:
			print("Penalthy method not encountered")
			sys.exit("Penalty method not encountered")
	while functionEvaluations < maxFE:
		flags = [-1, -1, -1]
		offspringIdx = 0
		for i in range(parentsSize):
			for k in range(generatedOffspring):
				flags = [-1, -1, -1]
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
			offsprings.bounding(nSize, function, offspringsSize, lowerBound, upperBound)
			functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
		else:
			if hasConstraints:
				offsprings.bounding(nSize, function, offspringsSize)
				functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
			else:
				offsprings.bounding(nSize, function, offspringsSize, hasConstraints)
				functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, -1, -1, functionEvaluations)
		if hasConstraints:
			if penaltyMethod == 1:  # Deb method
				offsprings.sumViolations(offspringsSize, gSize, hSize)
			elif penaltyMethod == 2:	# APM
				offsprings.uniteConstraints(offspringsSize, gSize, hSize)
				avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
				offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
				parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
		if hasConstraints:
			parents.DESelection(offsprings, generatedOffspring, parentsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod, hasConstraints)
		else:
			parents.DESelection(offsprings, generatedOffspring, parentsSize, nSize, -1, -1, -1, penaltyMethod, hasConstraints)

		# parents.printBest(nSize, parentsSize, penaltyMethod, hasConstraints)

		# parents.printBestFO(parentsSize, penaltyMethod, hasConstraints)
		# print("Ta rodando o dE")
		# weight = parents.calculateTrussWeight(parentsSize, penaltyMethod, bars)
		# weight = parents.calculateTrussWeightGroupingBest(parentsSize, penaltyMethod, bars, grouping, function)
		# print("Weigth: {:e}".format(weight))
	# parents.printBestFO(parentsSize, penaltyMethod, hasConstraints)
	if strFunction[0] == "8":
		best = bestIndividual(parents, parentsSize, penaltyMethod, hasConstraints)
		# parents.printBestFO(parentsSize, penaltyMethod, hasConstraints)
		best.objectiveFunction[0] = best.objectiveFunction[0] - 100 * int(strFunction[1:])
	parents.printBestFO(parentsSize, penaltyMethod, hasConstraints)


def DERobson(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma, windowSize):  # Differential Evolution
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
		avgObjFuncPar = -1
		file = open('crossVDatas.txt', 'w')
		parents = Population(parentsSize, nSize, function, True, lowerBound, upperBound)
		offsprings = Population(offspringsSize, nSize, function, True, lowerBound, upperBound)
		# functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
		functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
	else:  # Solving 'normal' functions
		gSize, hSize, constraintsSize = initializeConstraints(function)  # Initialize constraints
		penaltyCoefficients = [-1 for i in range(constraintsSize)]
		avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
		parents = Population(parentsSize, nSize, function)  # Initialize parents population
		offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
		functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations)  # Evaluate parents
	if penaltyMethod == 1:  # Padrao?  (not apm)
		# parents.sumViolations(parentsSize, gSize, hSize)
		offsprings.sumViolations(offspringsSize, gSize, hSize)  # TODO: LINHA CODIGO ROBSON
		# parents.sort(offsprings, penaltyMethod)
		offsprings.sortPopulation(penaltyMethod)
		for i in range(parentsSize):
			parents.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
		# populationList = offsprings.makePopulationList(offspringsSize, nSize, True)  # true stands for adding violationSum to populationList
		# print(len(populationList))
		# sys.exit("aaa")
		# file = open('crossVDatas.txt', 'w')
		# modelGPR, crossVMean, crossVStd = modelo_GPR_log.surGPR_training(populationList, offspringsSize, nSize + 1)  # TODO: TREINAR MODELO ROBSON
		modelGPR, crossVMean, crossVStd = adjustGPRModel(offsprings, offspringsSize, nSize, penaltyMethod, functionEvaluations)
		# print(crossVMean, crossVStd)
		file.write("{}\t{}\n".format(crossVMean, crossVStd))
		# offsprings.printDimensionsAndViolationPopulation(offspringsSize, nSize)  # TODO: LINHA CODIGO ROBSON
		# sys.exit("saiu antes ainda")
	elif penaltyMethod == 2:  # // Adaptive Penalty Method ( APM )
		# parents.uniteConstraints(parentsSize, gSize, hSize)
		offsprings.uniteConstraints(offspringsSize, gSize, hSize)
		# avgObjFuncPar = parents.calculatePenaltyCoefficients(parentsSize, constraintsSize, penaltyCoefficients, avgObjFuncPar)
		avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
		# parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFuncPar)
		offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
		offsprings.sortPopulation(penaltyMethod)
		for i in range(parentsSize):
			parents.copyIndividual(i, i, offsprings, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)
		modelGPR, crossVMean, crossVStd = adjustGPRModel(offsprings, offspringsSize, nSize, penaltyMethod, functionEvaluations)
		# print(crossVMean, crossVStd)
	else:
		print("Penalthy method not encountered")
		sys.exit("Penalty method not encountered")
	# print("funcEvaluation: {}".format(functionEvaluations))
	# sys.exit()
	cont = 0
	popForTraining = Population(750, nSize, 2, False, -1, -1)
	popForTrainingIdx = 0
	while functionEvaluations < maxFE:
		if functionEvaluations != 200:  # Not first iteration
			if crossVMean > 0.6:  # modelo "bom"
				cont = cont + 1
				if cont == windowSize:  # Treina modelo pela janela
					# modelGPR, crossVMean, crossVStd = adjustGPRModel(offsprings, offspringsSize, nSize, penaltyMethod, functionEvaluations)
					# print(setForTraining)
					modelGPR, crossVMean, crossVStd = adjustGPRModel(popForTraining, popForTrainingIdx, nSize, penaltyMethod, functionEvaluations)

					# print(crossVMean, crossVStd)
					file.write("{}\t{}\n".format(crossVMean, crossVStd))

					popForTraining = Population(750, nSize, 2, False, -1, -1)
					popForTrainingIdx = 0
					cont = 0
			else:  # Treina modelo pela qualidade
				# modelGPR, crossVMean, crossVStd = adjustGPRModel(offsprings, offspringsSize, nSize, penaltyMethod, functionEvaluations)
				# print(crossVMean, crossVStd)
				modelGPR, crossVMean, crossVStd = adjustGPRModel(popForTraining, popForTrainingIdx, nSize, penaltyMethod, functionEvaluations)
				file.write("{}\t{}\n".format(crossVMean, crossVStd))
				popForTraining = Population(750, nSize, 2, False, -1, -1)
				popForTrainingIdx = 0
		flags = [-1, -1, -1]
		offspringIdx = 0
		for i in range(parentsSize):
			for k in range(generatedOffspring):
				flags = [-1, -1, -1]
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
			offsprings.evaluateOnGPRModel(offspringsSize, nSize, modelGPR, bars, grouping, function, penaltyMethod)
			# bestIndividuals contains the best individuals(offsprings) evaluated on modelGPR and then they are evaluated on simulator
			bestsIndividuals = Population(parentsSize, nSize, 2, False, -1, -1)  # last two parameters is lowerBound and upperBound, will be subscribed
			bestsIndividuals, functionEvaluations = offsprings.evaluateOnlyBestsFromModelGPR(bestsIndividuals, generatedOffspring, offspringsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod, function, functionEvaluations, truss)
			"""
			for i in range(parentsSize):  # len(bestIndividuals)
				adds = 0
				for j in range(len(popForTraining.individuals)):  # returns true if lists are equal
					if not compareLists(popForTraining.individuals[j].n, bestsIndividuals.individuals[i].n):  # and compareLists(popForTraining.individuals[i].objectiveFunction, bestsIndividuals.individuals[i].objectiveFunction) and compareLists()
						adds = adds + 1
				if adds == len(popForTraining.individuals):  # if individuals are different
					popForTraining.copyIndividual(popForTrainingIdx, i, bestsIndividuals, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)  # copies bests to training pop
					popForTrainingIdx = popForTrainingIdx + 1
			"""
			for i in range(parentsSize):  # len(bestIndividuals)
				popForTraining.copyIndividual(popForTrainingIdx, i, bestsIndividuals, nSize, 1, gSize, hSize, constraintsSize, -1, penaltyMethod)  # copies bests to training pop
				popForTrainingIdx = popForTrainingIdx + 1
			# print(bestsIndividuals)
			# sys.exit("Printando bestsIndividuals")
			# def evaluateOnlyBestsFromModelGPR(self, bestsIndividuals, generatedOffspring, offspringsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod, function, functionEvaluations, truss):
			# offsprings.evaluateOnlyBestsFromModelGPR(bestINd)
			# sys.exit("acima csums")
			# functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
		else:
			offsprings.bounding(nSize, function, offspringsSize)
			functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
		if penaltyMethod == 1:  # Not apm
			# offsprings.sumViolations(offspringsSize, gSize, hSize)
			bestsIndividuals.sumViolations(parentsSize, gSize, hSize)
		elif penaltyMethod == 2:
			# calculate fitness biased on simulator
			bestsIndividuals.uniteConstraints(parentsSize, gSize, hSize)
			avgObjFunc = bestsIndividuals.calculatePenaltyCoefficients(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
			bestsIndividuals.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
			parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
		# parents.DESelection(offsprings, generatedOffspring, parentsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod)
		# print("parents")
		# parents.printDimensionsAndViolationPopulation(parentsSize, nSize)
		parents.DESelectionModelGPR(bestsIndividuals, parentsSize, nSize, gSize, hSize, constraintsSize, penaltyMethod)
		# print("Best individuals")
		# bestsIndividuals.printDimensionsAndViolationPopulation(parentsSize, nSize)
		# print("parents")
		# parents.printDimensionsAndViolationPopulation(parentsSize, nSize)
		# parents.printBest(nSize, parentsSize, penaltyMethod)
		# parents.printBestFO(parentsSize, penaltyMethod)
		# weight = parents.calculateTrussWeight(parentsSize, penaltyMethod, bars)
		# weight = parents.calculateTrussWeightGroupingBest(parentsSize, penaltyMethod, bars, grouping, function)
		# print("Weight: {:e}".format(weight))
		# parents.printBest(nSize, parentsSize, penaltyMethod)
	print("{}\t{}".format(crossVMean, crossVStd), end="\t")
	parents.printBestFO(parentsSize, penaltyMethod)


def ES(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):  # Evolution Strategy
	strFunction = str(function)
	crossoverProb = -1
	np.random.seed(seed)
	hasConstraints = False
	if strFunction[0] == "1" or strFunction[0] == "2":
		hasConstraints = True
	functionEvaluations = 0
	generatedOffspring = int(offspringsSize / parentsSize)
	lowerBound = upperBound = truss = 0
	if strFunction[0] == "2":  # solving trusses
		truss, lowerBound, upperBound = initializeTruss(function)
		nSize = truss.getDimension()
		gSize, hSize, constraintsSize = initializeConstraintsTrusses(truss)
		penaltyCoefficients = [-1 for i in range(constraintsSize)]
		avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
		parents = Population(parentsSize, nSize, function, True, lowerBound, upperBound)
		offsprings = Population(offspringsSize, nSize, function,True, lowerBound, upperBound)
		functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
	else:
		if hasConstraints:
			gSize, hSize, constraintsSize = initializeConstraints(function)  # Initialize constraints
			penaltyCoefficients = [-1 for i in range(constraintsSize)]
			avgObjFunc = 0  # will be subscribed on 'calculatePenaltyCoefficients'
			parents = Population(parentsSize, nSize, function)  # Initialize parents population
			offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
			functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations)  # Evaluate parents
		else:
			parents = Population(parentsSize, nSize, function)  # Initialize parents population
			offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
			functionEvaluations = parents.evaluate(parentsSize, function, -1, -1, -1, functionEvaluations)  # Evaluate parents
	if hasConstraints:
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
			if hasConstraints:
				offsprings.bounding(nSize, function, offspringsSize)
				functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
			else:
				offsprings.bounding(nSize, function, offspringsSize)
				functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, -1, -1, functionEvaluations)
		if hasConstraints:
			if penaltyMethod == 1:
				offsprings.sumViolations(offspringsSize, gSize, hSize)
			else:
				offsprings.uniteConstraints(offspringsSize, gSize, hSize)
				avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
				offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
		if hasConstraints:
			parents.elitismES(offsprings, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, esType, generatedOffspring, penaltyMethod, strFunction, truss, lowerBound, upperBound, hasConstraints)
		else:
			parents.elitismES(offsprings, parentsSize, offspringsSize, nSize, -1, -1, -1, globalSigma, esType, generatedOffspring, penaltyMethod, strFunction, truss, lowerBound, upperBound, hasConstraints)
		# parents.printBest(nSize, parentsSize, penaltyMethod, hasConstraints)
	parents.printBestFO(parentsSize, penaltyMethod, hasConstraints)

"""
lambda = parentsSize
mi = offspringsSize
w = weights
sigma = steo size
B ∈ R n, an orthogonal matrix. Columns of B are eigenvectors of C with unit length and correspond to the diagonal elements of D.
C(g) ∈ R n×n, covariance matrix at generation g. cii, diagonal elements of C.
cc ≤ 1, learning rate for cumulation for the rank-one update of the covariance matrix, see (24) and (45), and Table 1.
c1 ≤ 1 − cµ, learning rate for the rank-one update of the covariance matrix update, see (28), (30), and (47), and Table 1.
cµ ≤ 1 − c1, learning rate for the rank-µ update of the covariance matrix update, see (16), (30), and (47), and Table 1.
cσ < 1, learning rate for the cumulation for the step-size control, see (31) and (43), and  Table 1
D ∈ Rn, a diagonal matrix. The diagonal elements of D are square roots of eigenvalues ofC and correspond to the respective columns of B.
di > 0, diagonal elements of diagonal matrix D, d²i are eigenvalues of C.
dσ ≈ 1, damping parameter for step-size update, see (32), (37), and (44).
"""


def ESCMAColunaTests(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):  # Evolution Strategy
	strFunction = str(function)
	esType = 1
	crossoverProb = -1
	np.random.seed(seed)
	hasConstraints = False
	if strFunction[0] == "1" or strFunction[0] == "2":
		hasConstraints = True
	functionEvaluations = 0
	# User defined parameters\

	"""
	λ ≥ 2, population size, sample size, number of offspring, see (5).
	µ ≤ λ parent number, number of (positively) selected search points in the population, number
	of strictly positive recombination weights, see (6).
	sigma = 0.5
	xmean = np.random.randn(nSize)  # np.random.randn(nSize, 1)
	# maxFE = 30

	# Strategy parameters setting: Selection
	parentsSize = 4 + np.floor(3 * np.log(nSize))  # parentsSize is biased on nSize
	parentsSize = int(parentsSize)
	mu = parentsSize / 2  # mu is NOT offspringsSize.
	muList = [i + 1 for i in range(int(mu))]
	weights = np.log(mu+1/2)-np.log(muList)  # muXone recombination weights
	mu = np.floor(mu)
	# mu = int(mu)
	weights = weights/np.sum(weights)
	mueff = np.power(sum(weights), 2) / np.sum(np.power(weights, 2))

	# Strategy parameter setting: Adaptation
	cc = (4+mueff / nSize) / (nSize+4 + 2*mueff/nSize) # time constant for cumulation for C
	cs = (mueff+2) / (nSize+mueff+5)  # t-const for cumulation for sigma control
	c1 = 2 / (np.power(nSize + 1.3, 2) + mueff)  # learning rate for rank one update of C
	cmu = np.minimum(1 - c1, 2 * (mueff - 2 + 1/mueff)/(np.power(nSize+2, 2) + 2*mueff/2))  # and for rank-mu update
	damps = 1 + 2*np.maximum(0, np.sqrt((mueff -1) / (nSize + 1)) -1 ) + cs  # damping for sigma

	# Initiliaze dynamic (internal) strategy parameters  and constants
	pc = np.zeros(nSize) # evolutions paths for C
	ps = np.zeros(nSize)  # evolutions paths for sigma
	B = np.eye(nSize)  # B defines de coordinate system
	D = np.eye(nSize)  # diagonal matrix D defines the scaling
	C = np.matmul(np.matmul(B, D), np.matmul(B, D).transpose())  #  covariance matrix
	eigenval = 0  # B and D update at counteval == 0
	chinN = nSize**(0.5) * (1-1/(4*nSize)+1 / (21*np.power(nSize, 2)))  # expectation of ||N(0,I)|| == norm(randn(N,1))
	"""

	# CODIGO ANTIGO

	# generatedOffspring = int(offspringsSize / parentsSize)  # TODO Verifiar isso depois
	generatedOffspring = 1
	offspringsSize = parentsSize
	lowerBound = upperBound = truss = 0
	if strFunction[0] == "2":  # solving trusses
		truss, lowerBound, upperBound = initializeTruss(function)
		nSize = truss.getDimension()
		print("nSize{}".format(nSize))
		gSize, hSize, constraintsSize = initializeConstraintsTrusses(truss)
		penaltyCoefficients = [-1 for i in range(constraintsSize)]
		avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
		parents = Population(parentsSize, nSize, function, True, lowerBound, upperBound)
		offsprings = Population(offspringsSize, nSize, function, True, lowerBound, upperBound)
		functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
	else:
		if hasConstraints:
			gSize, hSize, constraintsSize = initializeConstraints(function)  # Initialize constraints
			penaltyCoefficients = [-1 for i in range(constraintsSize)]
			avgObjFunc = 0  # will be subscribed on 'calculatePenaltyCoefficients'
			parents = Population(parentsSize, nSize, function)  # Initialize parents population
			offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
			functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations)  # Evaluate parents
		else:  # Functions with no constraints
			parents = Population(parentsSize, nSize, function)  # Initialize parents population
			offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
			functionEvaluations = parents.evaluate(parentsSize, function, nSize, -1, -1, functionEvaluations)
	# Choose method of penalization
	if hasConstraints:
		if penaltyMethod == 1:  # Padrao?  (not apm)
			parents.sumViolations(parentsSize, gSize, hSize)
		elif penaltyMethod == 2:  # // Adaptive Penalty Method ( APM )
			parents.uniteConstraints(parentsSize, gSize, hSize)
			avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
			parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
		else:
			print("Penalthy method not encountered")
			sys.exit("Penalty method not encountered")


	# User defined parameters
	sigma = 0.5
	xmean = np.random.randn(nSize)  # np.random.randn(nSize, 1)
	# maxFE = 30
	"""
	λ ≥ 2, population size, sample size, number of offspring, see (5).
	µ ≤ λ parent number, number of (positively) selected search points in the population, number
	of strictly positive recombination weights, see (6).
	"""
	# Strategy parameters setting: Selection
	parentsSize = 4 + np.floor(3 * np.log(nSize))  # parentsSize is biased on nSize
	parentsSize = int(parentsSize)
	offspringsSize = parentsSize  # temporay?
	mu = parentsSize / 2  # mu is NOT offspringsSize.
	muList = [i + 1 for i in range(int(mu))]
	weights = np.log(mu+1/2)-np.log(muList).conj().T  # muXone recombination weights
	mu = np.floor(mu)
	# mu = int(mu)
	weights = weights/np.sum(weights)
	mueff = np.power(sum(weights), 2) / np.sum(np.power(weights, 2))

	# Strategy parameter setting: Adaptation
	cc = (4+mueff / nSize) / (nSize+4 + 2*mueff/nSize) # time constant for cumulation for C
	cs = (mueff+2) / (nSize+mueff+5)  # t-const for cumulation for sigma control
	c1 = 2 / (np.power(nSize + 1.3, 2) + mueff)  # learning rate for rank one update of C
	cmu = np.minimum(1 - c1, 2 * (mueff - 2 + 1/mueff)/(np.power(nSize+2, 2) + 2*mueff/2))  # and for rank-mu update
	damps = 1 + 2*np.maximum(0, np.sqrt((mueff -1) / (nSize + 1)) -1 ) + cs  # damping for sigma

	# Initiliaze dynamic (internal) strategy parameters  and constants
	pc = np.zeros(nSize) # evolutions paths for C
	ps = np.zeros(nSize)  # evolutions paths for sigma
	B = np.eye(nSize)  # B defines de coordinate system
	D = np.eye(nSize)  # diagonal matrix D defines the scaling
	C = B @ D @ (B@D).conj().T
	# C = np.matmul(np.matmul(B, D), np.matmul(B, D).transpose())  #  covariance matrix
	# print(C)
	# exit("sai viado")
	eigenval = 0  # B and D update at counteval == 0
	chinN = nSize**(0.5) * (1-1/(4*nSize)+1 / (21*np.power(nSize, 2)))  # expectation of ||N(0,I)|| == norm(randn(N,1))



	# parents.initializeEvolutionStrategy(offsprings, nSize, parentsSize, offspringsSize, globalSigma)
	# FIM CODIGO ANTIGO

	counteval = 0
	while functionEvaluations < maxFE:
		# Generate and evaluate lambda offspring
		arzAuxList = []
		arxAuxList = []
		"""
		print("B")
		print(B)
		print("D")
		print(D)
		"""
		for i in range(parentsSize):
			arz = np.random.randn(nSize)  # standard normally distributed vector  (38)
			# arx = xmean + sigma * (np.dot(np.matmul(B, D), arz))  # add mutation N(m,sigma²C) (40)
			arx = xmean + sigma * (B @ D @ arz)  # add mutation N(m,sigma²C) (40)
			if functionEvaluations > 10:
				"""
				print("Conta (np.dot(B*D, arz))")
				print((np.dot(B*D, arz)))
				print("ArzLoop")
				print(arz)
				print("ArxLoop")
				print(arx)
				"""
			# arx = xmean + sigma * (np.dot(np.matmul(B, D), arz))
			arzAuxList.append(arz)
			arxAuxList.append(arx)
			# arx é  individuals[i].n
			for j in range(nSize):  # Copies arx to individual[].n TODO This can be done as offsprings.individuals[k].n = [i for i in arx] ?!
				offsprings.individuals[i].arz[j] = arz[j]
				offsprings.individuals[i].n[j] = arx[j]
			counteval = counteval + 1
		arz = np.vstack(arzAuxList)  # matrix nd.array with all values calculated above
		arx = np.vstack(arxAuxList)  # matrix nd.array with all values calculated above
		if functionEvaluations > 10:
			"""
			print("sigma")
			print(sigma)
			print("xmean")
			print(xmean)
			print("arz")
			print(arz)
			print("arx")
			print(arx)
			"""
		# Evaluate function
		if strFunction[0] == "2":
			offsprings.bounding(nSize, function, offspringsSize, lowerBound, upperBound)
			functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
		else:
			if hasConstraints:
				offsprings.bounding(nSize, function, offspringsSize) # Commented for rosenbrock
				functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
			else:
				offsprings.bounding(nSize, function, offspringsSize) # Commented for rosenbrock
				functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, -1, -1, functionEvaluations)

		if hasConstraints:
			# Choose method of penalization
			if penaltyMethod == 1:  # Deb penalization
				offsprings.sumViolations(offspringsSize, gSize, hSize)
			else:  # APM penalization
				offsprings.uniteConstraints(offspringsSize, gSize, hSize)
				avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
				offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)

		# Elitismo apenas ordena os filhos. esType == 1(pai não é levado em conta (ES ,))
		# Sort individuals
		if hasConstraints:
			parents.elitismES(offsprings, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, esType, generatedOffspring, penaltyMethod, strFunction, truss, lowerBound, upperBound, hasConstraints)
		else:
			parents.elitismES(offsprings, parentsSize, offspringsSize, nSize, -1, -1, -1, globalSigma, esType, generatedOffspring, penaltyMethod, strFunction, truss, lowerBound, upperBound, hasConstraints)

		# parents.elitismForRosenbrockCMAES(offsprings, parentsSize, offspringsSize, nSize)

		# Individuals are sorted
		# print(arx)
		"""
		for i in range(parentsSize):  # TODO Individuos armazenados em LINHA, não em COLUNAS, como no octave !! ATENCAO
			for j in range(nSize):
				arx[i][j] = parents.individuals[i].n[j]
				arz[i][j] = parents.individuals[i].arz[j]  #
		"""

		# arxB4 = np.copy(arx)
		arz = arz.T
		arx = arx.T
		for j in range(parentsSize):  # TODO Individuos armazenados em COLUNA, não em LINHAS, como no octave !! ATENCAO
			for i in range(nSize):
				arx[i][j] = parents.individuals[j].n[i]
				arz[i][j] = parents.individuals[j].arz[i]  #
			# del offsprings.individuals[i].n[nSize:]

		# print(arx)
		# print(arz)
		# xmeans = []
		# zmeans = []
		# weights = list(weights)
		# del offsprings.individuals[0].n[nSize:]  # delete the unused part of an list TODO seria bom fazer isso para todos os atributos de individuos, após a inicialização da população ser feita
		"""
		### INICIO Primeiro jeito de fazer, manualmente
		for j in range(nSize):
			somaX = 0
			somaZ = 0
			for i in range(int(mu)):
				# print(np.dot(np.asarray(offsprings.individuals[i].n), weights)
				somaX = offsprings.individuals[i].n[j] * weights[i] + somaX
				somaZ = offsprings.individuals[i].arz[j] * weights[i] + somaZ
				# xmeans = np.dot(np.asarray(offsprings.individuals[i].n), weights, xmeans)
			xmeans.append(somaX)
			zmeans.append(somaZ)
		### FIM Primeiro jeito de fazer, manualmente
		"""
		# print("Xmean x zmean por individuos\n")
		#print("{}\n{}".format(xmeans, zmeans))

		### INICIO Segundo jeito de fazer, com numpy
		# arx = np.asarray(arx)
		# arz = np.asarray(arz)
		"""
		# If individuals are stored in lines, muBestX and Z has to be tranposed
		muBestX = np.delete(arx, np.s_[int(mu):], 0)  # remove as linhas  de mu em diante da matrix arx
		muBestZ = np.delete(arz, np.s_[int(mu):], 0)  # remove as linhas  de mu em diante da matrix arx
		xmean = np.dot(muBestX.transpose(), weights)  # xmean is one array with nSize positions
		zmean = np.dot(muBestZ.transpose(), weights)  # zmeanis one array with nSize positions
		"""
		# Individuals already stored in columns, no need to tranpose
		muBestX = np.delete(arx, np.s_[int(mu):], 1)  # remove as colunas  de mu em diante da matrix arx
		muBestZ = np.delete(arz, np.s_[int(mu):], 1)  # remove as colunas  de mu em diante da matrix arx
		xmean = np.matmul(muBestX, weights)  # xmean is one array with nSize positions
		zmean = np.matmul(muBestZ, weights)  # zmeanis one array with nSize positions
		# xmeanTest = np.dot(arx[0:mu].transpose(), weights)  # works too, without needing to slice the arx vec
		"""
		print("zmean")
		print(zmean)
		print("xmean")
		print(xmean)
		"""
		### FIM Segundo jeito de fazer, com numpy
		# print("{}\n{}".format(xmean, zmean))
		# TODO Testar essa conta, com os mesmos valores do octave. TESTADO, CALCUCLO ABAIXO FUNCIONANDO PARA OS VALORES DA 1 ITERACAO DO OCTAVE
		ps = (1-cs)*ps + (np.sqrt(cs*(2-cs)*mueff)) * np.matmul(B, zmean)  # Eq. 43
		hsig = True if np.linalg.norm(ps) / np.sqrt(1-np.power((1-cs), (2*counteval/nSize)))/chinN < 1.4 + 2/(nSize + 1) else False
		# TODO Testar essa conta, com os mesmos valores do octave. TESTADO, CALCUCLO ABAIXO FUNCIONANDO PARA OS VALORES DA 1 ITERACAO DO OCTAVE
		pc = (1-cc)*pc + hsig * np.sqrt(cc*(2-cc)*mueff) * np.matmul(np.matmul(B, D), zmean)  # Eq. 45



		# Individuals already stored in columns, no need to tranpose
		# print("C before recalculing it")
		# print(C)
		# TODO Testar essa conta, com os mesmos valores do octave. TESTADO, CALCUCLO ABAIXO FUNCIONANDO PARA OS VALORES DA 1 ITERACAO DO OCTAVE
		"""
		Caux1 = (1-c1-cmu) * C
		Caux2 = np.outer(pc, pc) + ((1-hsig) * cc * (2-cc) * C)
		BDARGZ = np.matmul(np.matmul(B, D), muBestZ)
		Caux3 = c1*Caux2
		Caux4 = (cmu * BDARGZ)
		Caux5 = np.matmul(np.diag(weights), BDARGZ.transpose())
		Caux6 = np.matmul(Caux4, Caux5)
		C = Caux1 + Caux3 + Caux6
		"""
		# **************************************** NEW CALCULUS OF C (IM GONNA CRY)
		# TODO: CONTA TESTADA E "FUNCIONANDO" PARA A SEGUNDA ITERACAO COM AS MATRIZES DO OCTAVE
		C = (1-c1-cmu) * C + c1 * (np.outer(pc, pc) + (1-hsig) * cc*(2-cc) * C) + cmu * (B@D@muBestZ) @ np.diag(weights) @ (B@D@muBestZ).conj().T
		"""
		C = ((1-c1-cmu) * C #
			 + c1 * (np.outer(pc, pc)
			+ (1-hsig) * cc * (2 - cc) * C)  # THIS LAST MULTIPLICATION SHOULD BE @, MATLAB IS CRAP
			 + cmu
			 * (B@D@muBestZ)
			 @ np.diag(weights) @ (B@D@muBestZ).T
			)

		 """
		# print("FunctionEvaluations: {}".format(functionEvaluations))
		# print("C after recalculing it")
		# print(C)
		#  Adapt step-size sigma
		sigma = sigma * np.exp((cs/damps) * (np.linalg.norm(ps)/chinN - 1))  # Adapt sigma step-size Eq. 44

		# np.set_printoptions(suppress=True)

		# Update B and D from C
		if counteval - eigenval > parentsSize / (c1 + cmu) / nSize / 10:  # to achieve 0(N^2)
			eigenval = counteval
			# print("C before enforce symmetry")
			# print(C)
			C = np.triu(C) + np.triu(C, 1).conj().T  # enforce symmetry
			"""
			print("C DENTRO DE UPDATE B AND D")
			print("Printando C")
			print(C)
			"""
			# print("C after enforce symmetry")
			# print(C)
			D, B = np.linalg.eig(C)  # eigen decomposition, B == normalized eigenvector  # B está dando diferente do MATLAB, deve ser por ser autovetor
			# D = np.diag(D) # function above returns D as a diagonal matrix, (on octave), and on python return just a array. Transforms array on diag matrix
			D = np.diag(np.sqrt(D))
			# print("a")
			# print(D)
		"""
		print("B")
		print(B)
		print("D")
		print(D)
		"""
		# TODO can implement flatfitness if later
		parents.printBest(nSize, parentsSize, penaltyMethod, hasConstraints)


def initilaizeCMAESParemeters(cc, cs, c1, cmu, damps, pc, ps, B, D, C, eigenval, chinN):
	# Strategy parameter setting: Adaptation
	cc = (4+mueff / nSize) / (nSize+4 + 2*mueff/nSize) # time constant for cumulation for C
	cs = (mueff+2) / (nSize+mueff+5)  # t-const for cumulation for sigma control
	c1 = 2 / (np.power(nSize + 1.3, 2) + mueff)  # learning rate for rank one update of C
	cmu = np.minimum(1 - c1, 2 * (mueff - 2 + 1/mueff)/(np.power(nSize+2, 2) + 2*mueff/2))  # and for rank-mu update
	damps = 1 + 2*np.maximum(0, np.sqrt((mueff -1) / (nSize + 1)) -1 ) + cs  # damping for sigma

	# Initiliaze dynamic (internal) strategy parameters and constants
	pc = np.zeros(nSize) # evolutions paths for C
	ps = np.zeros(nSize)  # evolutions paths for sigma
	B = np.eye(nSize)  # B defines de coordinate system
	D = np.eye(nSize)  # diagonal matrix D defines the scaling
	C = B @ D @ (B@D).conj().T # covariance matrix
	eigenval = 0  # B and D update at counteval == 0
	chinN = nSize**(0.5) * (1-1/(4*nSize)+1 / (21*np.power(nSize, 2)))  # expectation of ||N(0,I)|| == norm(randn(N,1))

def ESCMAColuna(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):  # Evolution Strategykj
	strFunction = str(function)
	esType = 1
	crossoverProb = -1
	np.random.seed(seed)
	hasConstraints = False
	generatedOffspring = 1
	if strFunction[0] == "1" or strFunction[0] == "2":
		hasConstraints = True
	if strFunction[0] == "7":	# cec2020 no constraints competition
		if(nSize == 5):	
			maxFE = 50000
		elif (nSize == 10):
			maxFE = 1000000
		elif (nSize == 15):
			maxFE = 3000000
		elif (nSize == 20):
			maxFE = 10000000
		else:
			print("Dimension not defined")
	if strFunction[0] == "8":
		maxFE = nSize * 10000
	functionEvaluations = 0
	maxFE = 110000
	# User defined parameters
	sigma = 0.5 # coordinate wise standard deviation (step-size)
	xmean = np.random.randn(nSize)  # objective variables initial point

	"""
	λ ≥ 2, population size, sample size, number of offspring, see (5).
	µ ≤ λ parent number, number of (positively) selected search points in the population, number
	of strictly positive recombination weights, see (6).
	"""

	# Strategy parameters setting: Selection
	# parentsSize = 4 + np.floor(3 * np.log(nSize))  # population size (λ), offsprings number (same as parents)
	# parentsSize = int(parentsSize)
	parentsSize = nSize * 4 # population size, offsprings number (same as parents)
	offspringsSize = parentsSize 
	mu = parentsSize / 2   # number of parents selected (selected search points in the population)
	muList = [i + 1 for i in range(int(mu))]
	weights = np.log(mu+1/2)-np.log(muList).conj().T  # muXone recombination weights
	mu = np.floor(mu) # number of parents selected (µ)(selected search points in the population)
	weights = weights/np.sum(weights) # normalize recombination weights array
	mueff = np.power(np.sum(weights), 2) / np.sum(np.power(weights, 2)) # variance-effective size of mu

	# Strategy parameter setting: Adaptation
	cc = (4+mueff / nSize) / (nSize+4 + 2*mueff/nSize) # time constant for cumulation for C
	cs = (mueff+2) / (nSize+mueff+5)  # t-const for cumulation for sigma control
	c1 = 2 / (np.power(nSize + 1.3, 2) + mueff)  # learning rate for rank one update of C
	cmu = np.minimum(1 - c1, 2 * (mueff - 2 + 1/mueff)/(np.power(nSize+2, 2) + 2*mueff/2))  # and for rank-mu update
	damps = 1 + 2*np.maximum(0, np.sqrt((mueff -1) / (nSize + 1)) -1 ) + cs  # damping for sigma

	# Initiliaze dynamic (internal) strategy parameters and constants
	pc = np.zeros(nSize) # evolutions paths for C
	ps = np.zeros(nSize)  # evolutions paths for sigma
	B = np.eye(nSize)  # B defines de coordinate system
	D = np.eye(nSize)  # diagonal matrix D defines the scaling
	C = B @ D @ (B@D).conj().T # covariance matrix
	eigenval = 0  # B and D update at counteval == 0
	chinN = nSize**(0.5) * (1-1/(4*nSize)+1 / (21*np.power(nSize, 2)))  # expectation of ||N(0,I)|| == norm(randn(N,1))

	# Initialize Individuals
	lowerBound = upperBound = truss = 0
	if strFunction[0] == "2":  # solving trusses
		truss, lowerBound, upperBound = initializeTruss(function)
		nSize = truss.getDimension()
		gSize, hSize, constraintsSize = initializeConstraintsTrusses(truss)
		penaltyCoefficients = [-1 for i in range(constraintsSize)]
		avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
		parents = Population(parentsSize, nSize, function, True, lowerBound, upperBound)
		offsprings = Population(offspringsSize, nSize, function, True, lowerBound, upperBound)
		functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
	else:
		if hasConstraints:
			gSize, hSize, constraintsSize = initializeConstraints(function)  # Initialize constraints
			penaltyCoefficients = [-1 for i in range(constraintsSize)]
			avgObjFunc = 0  # will be subscribed on 'calculatePenaltyCoefficients'
			parents = Population(parentsSize, nSize, function)  # Initialize parents population
			offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
			functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations)  # Evaluate parents
		else:  # Functions with no constraints
			parents = Population(parentsSize, nSize, function)  # Initialize parents population
			offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
			# functionEvaluations = parents.evaluate(parentsSize, function, nSize, -1, -1, functionEvaluations)
	# Choose method of penalization
	if hasConstraints:
		if penaltyMethod == 1:  # Padrao?  (not apm)
			parents.sumViolations(parentsSize, gSize, hSize)
		elif penaltyMethod == 2:  # // Adaptive Penalty Method ( APM )
			parents.uniteConstraints(parentsSize, gSize, hSize)
			avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
			parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
		else:
			print("Penalthy method not encountered")
			sys.exit("Penalty method not encountered")

	tempBestInd = Individual()  # in case of covariance matrix degenerates
	# bestOfEachPopulation = 0  # selects the best individual among all populations
	bestOfEachPopulation = Individual()  # selects the best individual among all populations
	deg = 0
	while functionEvaluations < maxFE:
		arzAuxList = []
		arxAuxList = []
		complexIdx = 0
		# Generate and evaluate lambda offspring
		for i in range(parentsSize):
			arz = np.random.randn(nSize)  # standard normally distributed vector  (38)
			# print("Arz: {}".format(arz))
			arx = xmean + sigma * (B @ D @ arz)  # add mutation N(m,sigma²C) (40)
			# print("B@D@arz: {}".format(B @ D @ arz))
			# print("Arx: {}".format(arx))
			# if np.any(arx  > 1000000):
			# 	print(arx)
			# 	print("deg: {}".format(deg))
			# 	sys.exit("explodiu")
			arzAuxList.append(arz)
			arxAuxList.append(arx)
			# TODO Use the vector on individuals used on others ES, instead of using individuals.arz
			for j in range(nSize):  # Copies arx to individual[].n TODO: This can be done as offsprings.individuals[k].n = [i for i in arx] ?!
				offsprings.individuals[i].arz[j] = arz[j]
				offsprings.individuals[i].n[j] = arx[j]
				complexxx = np.iscomplexobj(offsprings.individuals[i].n)
				# if complexxx:
				# 	complexIdx = i
					# print(offsprings.individuals[complexIdx].n)
			# counteval = counteval + 1 # TODO: Pode ser removido (equivalente ao functionsEvaluations)

		if complexxx:
			# print("offsprings.[{}]:".format(complexIdx))
			
			sys.exit("Complex number in individuals")
		arz = np.vstack(arzAuxList)  # matrix nd.array with all values calculated above
		arx = np.vstack(arxAuxList)  # matrix nd.array with all values calculated above

		# Evaluate function
		if strFunction[0] == "2":
			offsprings.bounding(nSize, function, offspringsSize, lowerBound, upperBound)
			functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
		else:
			if hasConstraints:
				offsprings.bounding(nSize, function, offspringsSize) # Commented for rosenbrock
				functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize, functionEvaluations)
			else:
				offsprings.bounding(nSize, function, offspringsSize) # Commented for rosenbrock
				functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, -1, -1, functionEvaluations)

		if hasConstraints:
			# Choose method of penalization
			if penaltyMethod == 1:  # Deb penalization
				offsprings.sumViolations(offspringsSize, gSize, hSize)
			else:  # APM penalization
				offsprings.uniteConstraints(offspringsSize, gSize, hSize)
				avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
				offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)

		# if(bestOfEachPopulation != 0):
		print("The value of bestOfEachPopulation.objFunction ONDE VOCE MUDA CARA 1: {}".format(bestOfEachPopulation.objectiveFunction[0]))

		# Sort offsprings and put then (all) on parents 
		if hasConstraints:
			parents.elitismES(offsprings, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma, esType, generatedOffspring, penaltyMethod, strFunction, truss, lowerBound, upperBound, hasConstraints)
		else:
			parents.elitismES(offsprings, parentsSize, offspringsSize, nSize, -1, -1, -1, globalSigma, esType, generatedOffspring, penaltyMethod, strFunction, truss, lowerBound, upperBound, hasConstraints)
		# if(bestOfEachPopulation != 0):
		print("The value of bestOfEachPopulation.objFunction ONDE VOCE MUDA CARA 2: {}".format(bestOfEachPopulation.objectiveFunction[0]))

		#  Stored individuals , just like matlab code
		arz = arz.T
		arx = arx.T
		for j in range(parentsSize):
			for i in range(nSize):
				arx[i][j] = parents.individuals[j].n[i]
				arz[i][j] = parents.individuals[j].arz[i]  #
			# del offsprings.individuals[i].n[nSize:]

		# del offsprings.individuals[0].n[nSize:]  # delete the unused part of an list TODO seria bom fazer isso para todos os atributos de individuos, após a inicialização da população ser feita

		# Individuals already stored in columns, no need to tranpose
		muBestX = np.delete(arx, np.s_[int(mu):], 1)  # Remove columns after mu index (Select mu individuals)
		muBestZ = np.delete(arz, np.s_[int(mu):], 1)  # Remove columns after mu index (Select mu individuals)
		xmean = np.matmul(muBestX, weights)  # xmean is one array with nSize positions. Recombination Eq. 42
		zmean = np.matmul(muBestZ, weights)  # zmeanis one array with nSize positions. == Dˆ-1*B’*(xmean-xold)/sigma

		# Cumulation: Updatte evolution paths
		ps = (1-cs)*ps + (np.sqrt(cs*(2-cs)*mueff)) * B@zmean  # Eq. 43
		hsig = True if np.linalg.norm(ps) / np.sqrt(1-np.power((1-cs), (2*functionEvaluations/parentsSize)))/chinN < 1.4 + 2/(nSize + 1) else False
		pc = (1-cc)*pc + hsig * np.sqrt(cc*(2-cc)*mueff) * B@D@zmean # Eq. 45

		# print("C before recalculing it")
		# print(C)

		# Adapt covariance matrix C. Eq. 47
		# C = (1-c1-cmu) * C + c1 * (np.outer(pc, pc) + (1-hsig) * cc*(2-cc) * C) + cmu * (B@D@muBestZ) @ np.diag(weights) @ (B@D@muBestZ).conj().T # Eq. 47
		C = ((1-c1-cmu) * C + # regard old matrix
			c1 * (np.outer(pc, pc) + # plus rank on update
			(1-hsig) * cc*(2-cc) * C) + # minor correction
			cmu *  # plus rank mu update
			(B@D@muBestZ) @ np.diag(weights) @ (B@D@muBestZ).conj().T)

		# cIsComplex = np.iscomplexobj(C)
		# cIsNan= np.isnan(C)
		# cIsInf = np.isinf(C)
		# print("cIsComplex: {}".format(cIsComplex))
		# print("cIsNan: {}".format(cIsNan))
		# print("cIsInf: {}".format(cIsInf))
		
		# print("C after recalculing it")
		# print(C)

		#  Adapt step-size sigma
		sigma = sigma * np.exp((cs/damps) * (np.linalg.norm(ps)/chinN - 1)) # Eq. 44

		# First population
		if(bestOfEachPopulation.objectiveFunction[0] == -1):
		# if(bestOfEachPopulation == 0):
			print("tá entrandoa qui o porra")
			bestOfEachPopulation = bestIndividual(parents, parentsSize, penaltyMethod, hasConstraints)
		else:
			# print("The value of bestOfEachPopulation.objFunction before calculating the new presentBest: {}".format(bestOfEachPopulation.objectiveFunction[0]))
			presentBest = bestIndividual(parents, parentsSize, penaltyMethod, hasConstraints)
			# print("The value of bestOfEachPopulation.objFunction just after calculating the new presentBest: {}".format(bestOfEachPopulation.objectiveFunction[0]))
			if(presentBest.objectiveFunction[0] < bestOfEachPopulation.objectiveFunction[0]):
				print("The value of bestOfEachPopulation.objFunction: {}".format(bestOfEachPopulation.objectiveFunction[0]))
				print("The value of presentBest.objFunction: {}".format(presentBest.objectiveFunction[0]))
				if hasConstraints:
					print("has constraints, need to implement")
				else:
					print("ah bao")
					# bestOfEachPopulation.copyIndividual(presentBest, nSize, 1, -1, -1, -1, -1, -1)
					bestOfEachPopulation.copyIndividual(presentBest, nSize, 1, 0, 0, 0, 0, 0)
					print("The value of bestOfEachPopulation.objFunction after being exchange with presentBest: {}".format(bestOfEachPopulation.objectiveFunction[0]))
					print("finalizou uma geracao")

		# Update B and D from C
		if functionEvaluations - eigenval > parentsSize / (c1 + cmu) / nSize / 10:  # to achieve 0(N^2)
			eigenval = functionEvaluations
			# print("C before enforce symmetry")
			# print(C)
			C = np.triu(C) + np.triu(C, 1).conj().T  # enforce symmetry

			# print("C after enforce symmetry")
			# print(C)
			# # # # # print("D before recalculating")
			# # # # # print(D)
			# np.linalg returns eigenvalues (D) as an array
			D, B = np.linalg.eig(C)  # eigen decomposition, B == normalized eigenvector 
			if(np.any(D<=0)): # degenerates
				deg = deg + 1
				print("Degenerou {} vezes".format(deg))
				# print("D has negative values in it")
				# print(D)
				# print(arx)
				currentBest = bestIndividual(parents, parentsSize, penaltyMethod, hasConstraints)
				# tempBestInd already exists 
				if(tempBestInd.objectiveFunction[0] != -1):
					print("The value of tempBestInd.objFunction: {}".format(tempBestInd.objectiveFunction[0]))
					print("The value of currentBest.objFunction: {}".format(currentBest.objectiveFunction[0]))
					# print("The values of tempBestInd.n: {}".format(tempBestInd.n))
					# current individual is better than tempBestInd, it becames tempBestInd
					if(currentBest.objectiveFunction[0] <= tempBestInd.objectiveFunction[0]):
						if hasConstraints:
							print("has constraints, need to implement")
						else:
							print("New better than older")
							tempBestInd.copyIndividual(currentBest, nSize, 1, -1, -1, -1, -1, -1)
				else: # First times it degenaretes current best is the 
					tempBestInd.copyIndividual(currentBest, nSize, 1, -1, -1, -1, -1, -1)
					
				# parents.printBest(nSize, parentsSize, penaltyMethod, hasConstraints)
				# tempBestInd.printFO(parentsSize, penaltyMethod, hasConstraints)
				
				
				# User defined parameters
				sigma = 0.5 # coordinate wise standard deviation (step-size)
				xmean = np.random.randn(nSize)  # objective variables initial point

				# Strategy parameter setting: Adaptation
				cc = (4+mueff / nSize) / (nSize+4 + 2*mueff/nSize) # time constant for cumulation for C
				cs = (mueff+2) / (nSize+mueff+5)  # t-const for cumulation for sigma control
				c1 = 2 / (np.power(nSize + 1.3, 2) + mueff)  # learning rate for rank one update of C
				cmu = np.minimum(1 - c1, 2 * (mueff - 2 + 1/mueff)/(np.power(nSize+2, 2) + 2*mueff/2))  # and for rank-mu update
				damps = 1 + 2*np.maximum(0, np.sqrt((mueff -1) / (nSize + 1)) -1 ) + cs  # damping for sigma

				# Initiliaze dynamic (internal) strategy parameters and constants
				pc = np.zeros(nSize) # evolutions paths for C
				ps = np.zeros(nSize)  # evolutions paths for sigma
				B = np.eye(nSize)  # B defines de coordinate system
				D = np.eye(nSize)  # diagonal matrix D defines the scaling
				C = B @ D @ (B@D).conj().T # covariance matrix
				# eigenval = 0  # B and D update at counteval == 0 #TODO: Verify is eigenval need to be 0 again.
				chinN = nSize**(0.5) * (1-1/(4*nSize)+1 / (21*np.power(nSize, 2)))  # expectation of ||N(0,I)|| == norm(randn(N,1))
				# print("D after restartring it(after sqrt, final D)")
				# print(D)
				# sys.exit("Existe valores negativos em D")
			else:
				D = np.diag(np.sqrt(D)) # D containts standard deviations now
			# DdiagSqrt = np.diag(np.sqrt(D)) # D containts standard deviations now
			# Ddiag = np.diag(D)  # function above returns D as a diagonal matrix, (on octave), and on python return just a array. Transforms array on diag matrix
			# if(deg > 0):
			# 	print(D)
			# 	print(DdiagSqrt)
			# 	print(Ddiag)
			# 	sys.exit("pronto")
			# print("D after recalculating(after sqrt, final D)")
			
			# print(D)

		# TODO: can implement flatfitness if later



		# print(functionEvaluations)
		# parents.printBest(nSize, parentsSize, penaltyMethod, hasConstraints)

	bestFromCurrentPopulation = bestIndividual(parents, parentsSize, penaltyMethod, hasConstraints)
	# best individual from current is better than the best individuals of previous populations
	# print("The value of tempBestInd.objFunction: {}".format(tempBestInd.objectiveFunction[0]))
	# print("The value of bestFromCurrentPopulation.objFunction: {}".format(bestFromCurrentPopulation.objectiveFunction[0]))
	if (bestFromCurrentPopulation.objectiveFunction[0] < tempBestInd.objectiveFunction[0]):
		if hasConstraints:
			print("Need to implement, todo")
		else:
			tempBestInd.copyIndividual(bestFromCurrentPopulation, nSize, 1, -1, -1, -1, -1, -1)
	# Prints objective function of the best individual
	print("Best individual (based on actual and previous populations)")
	tempBestInd.printFO(parentsSize, penaltyMethod, hasConstraints)
	# print("Best individual among all populations (bestOfEachPopulation)")
	# bestOfEachPopulation.printFO(parentsSize, penaltyMethod, hasConstraints)

		
	
# noinspection PyShadowingNames
def algorithm(algorithm, function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma, windowSize):
	if algorithm == "GA":  # Genetic Algorithm
		GA(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma)
	elif algorithm == "DE":  # Differential Evolution
		"""
		start = timer()
		DE(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma)
		end = timer()
		print(end - start)
		"""
		start = timer()
		DE(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma)
		# DERobson(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma, windowSize)
		end = timer()
		# print("{} seconds".format(end - start), end="")
	elif algorithm == "ES":  # Evolution Strategy
		ES(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma)
	elif algorithm == "ESCMA":
		ESCMAColuna(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma)
	elif algorithm == "ESCMATest":
		ESCMAColunaTests(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma )
	else:
		print("Algorithm not encountered")
		sys.exit("Algorithm not encountered")


def main():
	# function,seed,penaltyMethod,parentsSize,nSize,offspringsSize,maxFE,crossoverProb,esType,globalSigma
	# ES µ ≈ λ/4
	parser = argparse.ArgumentParser(description="Evolutionary Algorithms")
	parser.add_argument("--algorithm", "-a", type=str, default="ESCMA", help="Algorithm to be used (GA, ES, DE or ESCMA)")
	parser.add_argument("--function", "-f", type=int, default=71, help="Truss to be solved (10, 25, 60, 72 or 942 bars). "
						"For the truss problem, the first digit must be 2, followed by the number of the bars in the problem. "
						"Example: 225, is for the truss of 25 bars")
	parser.add_argument("--seed", "-s", type=int, default=1, help="Seed to be used")
	parser.add_argument("--penaltyMethod", "-p", type=int, default=1, help="Penalty method to be used. 1 for Deb Penalty or 2 for APM")
	parser.add_argument("--parentsSize", "-u", type=int, default=50, help="µ is the parental population size")  # u from µ (mi) | µ ≈ λ/4
	parser.add_argument("--nSize", "-n", type=int, default=5, help="Search space dimension")
	parser.add_argument("--offspringsSize", "-l", type=int, default=50, help="λ is number of offsprings, offsprings population size")  # l from λ (lambda) | µ ≈ λ/4
	parser.add_argument("--maxFE", "-m", type=int, default=10000, help="The max number of functions evaluations")
	parser.add_argument("--crossoverProb", "-c", type=int, default=100, help="The crossover probability [0,100]")
	parser.add_argument("--esType", "-e", type=int, default=0, help="The type of ES. 0 for ES(µ + λ) or 1 for ES(µ , λ)")
	parser.add_argument("--globalSigma", "-g", type=int, default=0, help="If the σ parameter is global or not. 1 for global σ or 0 if not")
	parser.add_argument("--windowSize", "-w", type=int, default=5, help="Size of the window for updating gaussian model")
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
	# args.function = 225
	# args.algorithm = "ES"
	# args.globalSigma = 0
	# args.maxFE = 15000
	# args.esType = 1
	# args.penaltyMethod = 2
	# args.function = 272

	# args.offspringsSize = args.parentsSize
	# 10, 72 - 10 000
	# 942 - 150 000
	# Demais - 15 000
	# args.seed = 2
	# readTrussInput()
	# F1-F5 and F8-F10: D = 5, 10, 15,20
	# F6 and F7: D = 10, 15, 20
	algorithm(args.algorithm, args.function, args.seed, args.penaltyMethod, args.parentsSize, args.nSize, args.offspringsSize, args.maxFE, args.crossoverProb, args.esType, args.globalSigma, args.windowSize)
	print(args)
	# algorithms = ["DE", "ESCMA"]
	# print(algorithms)s
	# if args.seed == 1:
	#     print("DE\tESCMA")
	# for i in algorithms:
	#     algorithm(i, args.function, args.seed, args.penaltyMethod, args.parentsSize, args.nSize, args.offspringsSize, args.maxFE, args.crossoverProb, args.esType, args.globalSigma, args.windowSize)

if __name__ == '__main__':
	main()
