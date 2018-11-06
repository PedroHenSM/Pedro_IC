import numpy as np
# noinspection PyUnresolvedReferences
import sys


# noinspection PyShadowingNames
def defineParameters():
    totalAlgorithms = 3
    algorithms = ['DE', 'GA', 'ES']
    a = 0

    totalFunctions = 1  # Funcoes
    functions = [210, 225, 260, 272, 2942]
    func = 0

    totalSeeds = 3  # Seeds
    seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
    s = 0

    totalPenaltyMethods = 2  # Penalty Methods
    penaltyMethods = [1, 2]
    pm = 0

    totalPopulation = 1  # Tamanho Populacao
    populations = [50, 100]
    p = 0

    """
    totalDimensions = 1  # TamX
    # dimensions=(10 30)
    dimensions = [10, 8, 25, 16, 59] # 10, 25, 60 , 72 and 942 bars and respctives dimensions
    d = 0
    """
    totaloffsprings = 2
    # offsprings = [50, 100, 200, 400]
    offsprings = [50, 200]
    l = 0

    totalFE = 1  # Max funtions evaluations (inicalmente 1 - apenas 20000)
    functionEvaluations = [300, 20000, 100000, 500000]
    fe = 0

    totalProbCrossover = 1  # Probabilidade Crossover
    probCrossovers = [80, 100]
    c = 0

    totalTipoES = 2  # Tipos de ES (+ ,)
    tipoES = [0, 1]  # 0 + | 1 ,
    es = 0

    totalSigma = 2  # Sigma global
    sigmas = [1, 2]  # 1 sigmaGlobal | (!=1)sigmaIndividual
    sig = 0

    return algorithms, totalAlgorithms, totalFunctions, functions, totalSeeds, seeds, totalPenaltyMethods, penaltyMethods, totalPopulation, populations, totaloffsprings, offsprings, totalFE, functionEvaluations, totalProbCrossover, probCrossovers, totalTipoES, tipoES, totalSigma, sigmas


# noinspection PyShadowingNames
def cabecalho(functions, func):
    print("C{:02}\tMean\t\tStd\t\tBest\t\tWorst".format(functions[func]))


# noinspection PyShadowingNames
def headerLatexObjFunc(functions, func):
    print("\\begin{table}[h]")  # h de Here
    print("\centering")
    print("\caption{{Function C{:02}}}".format(functions[func]))
    print("\\vspace{0.5cm}")
    print("\\begin{tabular}{@{} l | r r r r @{}}")  # r Divide primeira coluna das outras e rl alinha a direita
    print("\hline")
    print("Algorithm & Mean & Std & Best & Worst \\\\")
    print("\hline")


# noinspection PyShadowingNames,SpellCheckingInspection
def headerLatexObjFunc1(functions, func, algorithmsStr):
    print("\\begin{table}[h]")  # h de Here
    print("\centering")
    print("\caption{{Function C{:02}}}".format(functions[func]))
    print("\\vspace{0.5cm}")
    print("\\resizebox{\\textwidth}{!} {%")
    print("\\begin{tabular}{@{} l | r r r r r r@{}}")  # r Divide primeira coluna das outras e rl alinha a direita
    print("\hline")
    print("Algorithm & {} & {} & {} & {} & {} & {} \\\\".format(algorithmsStr[0], algorithmsStr[1], algorithmsStr[2], algorithmsStr[3], algorithmsStr[4], algorithmsStr[5]))
    print("\hline")


# noinspection PyShadowingNames,SpellCheckingInspection
def headerLatexVariableProjects(functions, func, algorithmsStr):
    print("\\begin{table}[h]")  # h de Here
    print("\centering")
    print("\caption{{Project Variables for Function C{:02}}}".format(functions[func]))
    print("\\vspace{0.5cm}")
    print("\\resizebox{\\textwidth}{!} {%")
    print("\\begin{tabular}{@{} c | r r r r r r@{}}")  # r Divide primeira coluna das outras e rl alinha a direita
    print("\hline")
    print("Variables & {} & {} & {} & {} & {} & {} \\\\".format(algorithmsStr[0], algorithmsStr[1], algorithmsStr[2], algorithmsStr[3], algorithmsStr[4], algorithmsStr[5]))
    print("\hline")


# noinspection PyShadowingNames
def footerLatex(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm):
    print("\end{tabular}")
    print("\\textbf{{FE}}: {}\t\\textbf{{Population}}: {}\t\\textbf{{Offsprings}}: {}\t\\textbf{{Crossover}}: {}\%\t\\textbf{{PenalthyMethod}}: {} \\\\ ".format(functionEvaluations[fe], populations[p], offsprings[l], probCrossovers[c], penaltyMethodsStr[pm]))
    print("\end{table}")
    print("\n")


# noinspection PyShadowingNames
def footerLatex1(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm):
    print("\end{tabular}%")
    print("}")
    print("\\textbf{{FE}}: {}\t\\textbf{{Population}}: {}\t\\textbf{{Offsprings}}: {}\t\\textbf{{Crossover}}: {}\%\t\\textbf{{PenalthyMethod}}: {} \\\\ ".format(functionEvaluations[fe], populations[p], offsprings[l], probCrossovers[c], penaltyMethodsStr[pm]))
    print("\end{table}")
    print("\n")


# noinspection PyShadowingNames
def footerLatexProjectVariables(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm):
    print("\end{tabular}%")
    print("}")
    print("\\textbf{{FE}}: {}\t\\textbf{{Population}}: {}\t\\textbf{{Offsprings}}: {}\t\\textbf{{Crossover}}: {}\%\t\\textbf{{PenalthyMethod}}: {} \\\\ ".format(functionEvaluations[fe], populations[p], offsprings[l], probCrossovers[c], penaltyMethodsStr[pm]))
    print("\end{table}")
    print("\n")


# noinspection PyShadowingNames
def latexModel(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm):
    headerLatexObjFunc(functions, func)
    footerLatex(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm)


# noinspection PyShadowingNames
def imprime(de, es01, es02, es11, es12, ag):
    print("DE\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(de), np.std(de), np.amin(de), np.amax(de)))
    print("AG\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(ag), np.std(ag), np.amin(ag), np.amax(ag)))
    print("ES + G\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es01), np.std(es01), np.amin(es01), np.amax(es01)))
    print("ES + I\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es02), np.std(es02), np.amin(es02), np.amax(es02)))
    print("ES , G\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es11), np.std(es11), np.amin(es11), np.amax(es11)))
    print("ES , I\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es12), np.std(es12), np.amin(es12), np.amax(es12)))


# def printLatex(de,es01, es02, es11, es12, ga, dimensions, functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, d, fe, p, l, c, pm, totalSeeds):
# noinspection PyShadowingNames
def objectiveFunctionAnalysis(algorithmsVariations, totalAlgorithmsVariations, algorithmsStr, statisticalAnalysisStr, totalSeeds):
    totalStatisticalAnalysis = 4  # Mean, std, max and min
    meansObjectiveFunction = []
    stdsObjectiveFunction = []
    bestsObjectiveFunction = []
    worstsObjectiveFunction = []
    # objectiveFunctions = [[]] * totalAlgorithmsVariations # doesnt work, same id for lists
    objectiveFunctions = [[] for i in range(totalAlgorithmsVariations)]
    # print("objectionFunctions")
    # print(objectiveFunctions)
    # print(*algorithmsVariations, sep="\n")
    for a in range(totalAlgorithmsVariations):
        for i in range(totalSeeds):
            # print("algorithmsVariations[{}][{}][0]: {}".format(a, i, algorithmsVariations[a][i][0]))
            objectiveFunctions[a].append(algorithmsVariations[a][i][0])
        meansObjectiveFunction.append(np.mean(objectiveFunctions[a]))
        stdsObjectiveFunction.append(np.std(objectiveFunctions[a]))
        bestsObjectiveFunction.append(np.amin(objectiveFunctions[a]))
        worstsObjectiveFunction.append(np.amax(objectiveFunctions[a]))

    """
    print(*algorithmsVariations, sep="\n")
    print("objectiveFunction[a]")
    print(objectiveFunctions[0])
    print("objectiveFunction")
    print(objectiveFunctions, sep="\n")
    """

    statisticalValues = ([])
    statisticalValues.append(meansObjectiveFunction)
    statisticalValues.append(stdsObjectiveFunction)
    statisticalValues.append(bestsObjectiveFunction)
    statisticalValues.append(worstsObjectiveFunction)
    bestStatisticalValues = [statisticalValues[0][0], statisticalValues[1][0], statisticalValues[2][0], statisticalValues[3][0]]  # means, stds, bests, worsts | first idx is the analysis(means,std,etc)and second is the algorithm(de,ga,etc)
    # TODO: Note: Could be improved using numpy (gets the minimum value of each dimension, no need to do the double loop
    for a in range(1, totalAlgorithmsVariations):  # Pick best means | stds | bests | worsts
        if np.mean(objectiveFunctions[a]) < bestStatisticalValues[0]:  # Mean
            bestStatisticalValues[0] = np.mean(objectiveFunctions[a])
        if np.std(objectiveFunctions[a]) < bestStatisticalValues[1]:  # Std
            bestStatisticalValues[1] = np.std(objectiveFunctions[a])
        if np.amin(objectiveFunctions[a]) < bestStatisticalValues[2]:  # Best
            bestStatisticalValues[2] = np.amin(objectiveFunctions[a])
        if np.amax(objectiveFunctions[a]) < bestStatisticalValues[3]:  # Worst
            bestStatisticalValues[3] = np.amax(objectiveFunctions[a])

    return bestStatisticalValues, statisticalValues, totalStatisticalAnalysis, bestsObjectiveFunction


# noinspection PyShadowingNames
def printObjectiveFunctionsAnalysis(bestStatisticalValues, statisticalValues, totalStatisticalAnalysis, algorithmsStr, statisticalAnalysisStr, penaltyMethodsStr, algorithmsOnColumns, totalAlgorithmsVariations, functions, func, functionEvaluations, populations, offsprings, probCrossovers, fe, p, l, c, pm):
    if algorithmsOnColumns:
        headerLatexObjFunc1(functions, func, algorithmsStr)
        for s in range(totalStatisticalAnalysis):
            print(statisticalAnalysisStr[s] + " &", end=" ")
            # print(algorithmsStr[a] + " &", end=' ')
            for j in range(totalAlgorithmsVariations):
                if bestStatisticalValues[s] == statisticalValues[s][j]:  # if actual algorithm is the best, put bold
                    if j == totalAlgorithmsVariations - 1:  # just for formatting
                        print("\\textbf{{{:e}}}".format(statisticalValues[s][j]), end=' ')
                    else:
                        print("\\textbf{{{:e}}} &".format(statisticalValues[s][j]), end=' ')
                else:
                    if j == totalAlgorithmsVariations - 1:  # just for formatting
                        print("{:e}".format(statisticalValues[s][j]), end=' ')
                    else:
                        print("{:e} &".format(statisticalValues[s][j]), end=' ')
            print("\\\\")
        footerLatex1(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm)
    else:
        headerLatexObjFunc(functions, func)
        for a in range(totalAlgorithmsVariations):
            print(algorithmsStr[a] + " &", end=' ')
            for j in range(totalStatisticalAnalysis):
                if bestStatisticalValues[j] == statisticalValues[j][a]:  # if actual algorithm is the best, put bold
                    if j == totalStatisticalAnalysis - 1:  # just for formatting
                        print("\\textbf{{{:e}}}".format(statisticalValues[j][a]), end=' ')
                    else:
                        print("\\textbf{{{:e}}} &".format(statisticalValues[j][a]), end=' ')
                else:
                    if j == totalStatisticalAnalysis - 1:  # just for formatting
                        print("{:e}".format(statisticalValues[j][a]), end=' ')
                    else:
                        print("{:e} &".format(statisticalValues[j][a]), end=' ')
            print("\\\\")
        footerLatex(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm)


# noinspection PyShadowingNames
def printProjectVariables(bestOfEachAlgorithm, bestsObjectiveFunction, totalAlgorithmsVariations, functions, func, algorithmsStr, functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm):
    totalProjectVariables = len(bestOfEachAlgorithm[0][0])  # best of each algorithm has the objective function and project variable (all the same dimension)
    # print(bestOfEachAlgorithm[0][0])
    # print("bestOfEachAlgorithm")
    # print(*bestOfEachAlgorithm, sep="\n")
    # print(bestsObjectiveFunction)
    headerLatexVariableProjects(functions, func, algorithmsStr)

    for v in range(totalProjectVariables):
        # print(statisticalAnalysisStr[s] + " &", end=" ")
        if v != totalProjectVariables - 1:  # print project variables
            if func == 0:  # If func 0, the functions is 210 (10 bar truss), prints "v"
                print("{} &".format(v+1), end=" ")
            else:  # prints A
                print("A_{} &".format(v+1), end=" ")
        else:  # print weight
            print("W &", end=" ")
        for a in range(totalAlgorithmsVariations):  # the idx 0 is the objective function, will be printed after the project variables
            if v != totalProjectVariables - 1:  # print project variables
                if a == totalAlgorithmsVariations - 1:  # just for formatting
                    print("{:.4f}".format(bestOfEachAlgorithm[a][0][v+1]), end=' ')
                else:
                    print("{:.4f} &".format(bestOfEachAlgorithm[a][0][v+1]), end=' ')
            else:  # last iteration, print weight (objectiveFunction)
                if np.amin(bestsObjectiveFunction) == bestOfEachAlgorithm[a][0][0]:  # if actual algorithm is the best, put bold
                    if a == totalAlgorithmsVariations - 1:  # just for formatting
                        print("\\textbf{{{:e}}}".format(bestOfEachAlgorithm[a][0][0]), end=' ')
                    else:
                        print("\\textbf{{{:e}}} &".format(bestOfEachAlgorithm[a][0][0]), end=' ')
                else:  # just print
                    if a == totalAlgorithmsVariations - 1:  # just for formatting
                        print("{:e}".format(bestOfEachAlgorithm[a][0][0]), end=' ')
                    else:
                        print("{:e} &".format(bestOfEachAlgorithm[a][0][0]), end=' ')
        print("\\\\")
    footerLatexProjectVariables(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm)


# noinspection PyShadowingNames
def printLatex(bestOfEachAlgorithm, algorithmsVariations, de, functionEvaluations, populations, offsprings, probCrossovers, fe, p, l, c, pm, totalSeeds, totalAlgorithmsVariations, functions, func):
    algorithmsStr = ["DE", "GA", "ES + G", "ES + I", "ES , G", "ES , I"]
    statisticalAnalysisStr = ["Mean", "Std", "Best", "Worst"]
    penaltyMethodsStr = ["Deb", "APM"]
    algorithmsOnColumns = False
    bestStatisticalValues, statisticalValues, totalStatisticalAnalysis, bestsObjectiveFunction = objectiveFunctionAnalysis(algorithmsVariations, totalAlgorithmsVariations, algorithmsStr, statisticalAnalysisStr, totalSeeds)
    printObjectiveFunctionsAnalysis(bestStatisticalValues, statisticalValues, totalStatisticalAnalysis, algorithmsStr, statisticalAnalysisStr, penaltyMethodsStr, algorithmsOnColumns, totalAlgorithmsVariations, functions, func, functionEvaluations, populations, offsprings, probCrossovers, fe, p, l, c, pm)
    printProjectVariables(bestOfEachAlgorithm, bestsObjectiveFunction, totalAlgorithmsVariations, functions, func, algorithmsStr, functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm)


if __name__ == '__main__':
    algorithms, totalAlgorithms, totalFunctions, functions, totalSeeds, seeds, totalPenaltyMethods, penaltyMethods, totalPopulation, populations, totaloffsprings, offsprings, totalFE, functionEvaluations, totalProbCrossover, probCrossovers, totalTipoES, tipoES, totalSigma, sigmas = defineParameters()
    fe = c = 0
    # print(totalSeeds)
    # print(algorithms)
    # print(fe)
    for func in range(totalFunctions):
        for l in range(totaloffsprings):
            for p in range(totalPopulation):
                for pm in range(totalPenaltyMethods):
                    totalAlgorithmsVariations = 6
                    bestOfEachAlgorithm = ([[], [], [], [], [], []])
                    # bestOfEachAlgorithm = [[] for i in range(totalAlgorithmsVariations)]  # de, ga, es01, es02, es11, es12

                    de = []  # DE
                    es_0_1 = []  # ES 0 sigma global (+ 1)
                    es_0_2 = []  # ES 0 sigma indivudal (+ 2)
                    es_1_1 = []  # ES 1 sigma global (, 1)
                    es_1_2 = []  # ES 1 sigma individual (, 2)
                    ga = []  # AG
                    algorithmsVariations = ([])
                    algorithmsVariations.append(de)
                    algorithmsVariations.append(ga)
                    algorithmsVariations.append(es_0_1)
                    algorithmsVariations.append(es_0_2)
                    algorithmsVariations.append(es_1_1)
                    algorithmsVariations.append(es_1_2)
                    # print(algorithmsVariations)
                    # sys.exit("ui")
                    for a in range(totalAlgorithms):
                        for es in range(totalTipoES):
                            for sig in range(totalSigma):
                                for s in range(totalSeeds):
                                    if algorithms[a] == 'DE' and es == 0 and sig == 0:
                                        idxAlgorithmVariation = 0
                                        file = open("/home/pedrohen/Documentos/PedroIC/Resultados/trusses/DE/DE_truss_{}_{}_{}_{}_{}_{}.txt".format(functions[func], seeds[s], penaltyMethods[pm], populations[p], offsprings[l], functionEvaluations[fe]), "r")  # de.append(bestLine)
                                    elif algorithms[a] == 'GA' and es == 0 and sig == 0:
                                        idxAlgorithmVariation = 1
                                        file = open("/home/pedrohen/Documentos/PedroIC/Resultados/trusses/GA/GA_truss_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func], seeds[s], penaltyMethods[pm], populations[p], offsprings[l], functionEvaluations[fe], probCrossovers[c]), "r")
                                    elif algorithms[a] == 'ES':
                                        if tipoES[es] == 0:  # Es0 +
                                            if sigmas[sig] == 1:  # Sigma Global
                                                idxAlgorithmVariation = 2
                                                file = open("/home/pedrohen/Documentos/PedroIC/Resultados/trusses/ES/ES0/sigmaGlobal/ES_truss_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func], seeds[s], penaltyMethods[pm], populations[p], offsprings[l], functionEvaluations[fe], tipoES[es], sigmas[sig]), "r")
                                            elif sigmas[sig] == 2:  # Sigma Individual
                                                idxAlgorithmVariation = 3
                                                file = open("/home/pedrohen/Documentos/PedroIC/Resultados/trusses/ES/ES0/sigmaIndividual/ES_truss_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func], seeds[s], penaltyMethods[pm], populations[p], offsprings[l], functionEvaluations[fe], tipoES[es], sigmas[sig]), "r")
                                        elif tipoES[es] == 1:  # Es1 ,
                                            if sigmas[sig] == 1:  # Sigma Global
                                                idxAlgorithmVariation = 4
                                                file = open("/home/pedrohen/Documentos/PedroIC/Resultados/trusses/ES/ES1/sigmaGlobal/ES_truss_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func], seeds[s], penaltyMethods[pm], populations[p], offsprings[l], functionEvaluations[fe], tipoES[es], sigmas[sig]), "r")
                                            elif sigmas[sig] == 2:  # Sigma Individual
                                                idxAlgorithmVariation = 5
                                                file = open("/home/pedrohen/Documentos/PedroIC/Resultados/trusses/ES/ES1/sigmaIndividual/ES_truss_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func], seeds[s], penaltyMethods[pm], populations[p], offsprings[l], functionEvaluations[fe], tipoES[es], sigmas[sig]), "r")
                                    splittedLine = file.readline().split("\t")  # | Line Model is Violation\t ObjectiveFunction \t dimensions
                                    bestLine = splittedLine
                                    for line in file:  # continues reading file
                                        splittedLine = line.split("\t")
                                        if penaltyMethods[pm] == 1:  # Deb Penalty
                                            if float(splittedLine[0]) == 0 and splittedLine[1] < bestLine[1]:  # Tests if violation is 0 and then compares objecive Function
                                                bestLine = splittedLine
                                        elif penaltyMethods[pm] == 2:  # APM
                                            if float(splittedLine[0]) == float(splittedLine[1]) and splittedLine[1] < bestLine[1]:  # Test if fitness is equal to objectiveFunction and compares objectiveFunction
                                                bestLine = splittedLine
                                    bestLine = bestLine[1:-1]  # removes violation and '\n' from list, now the first idx is the ObjFunction
                                    bestLine = [float(i) for i in bestLine]  # converts all str values to float
                                    # print("Fora do printLatex")
                                    # print(*algorithmsVariations, sep="\n")
                                    if idxAlgorithmVariation != -1:
                                        algorithmsVariations[idxAlgorithmVariation].append(bestLine)
                                        if s == 0:  # If is the first iteration (seed 0)
                                            bestOfEachAlgorithm[idxAlgorithmVariation].append(bestLine)
                                        else:
                                            if bestLine[0] < bestOfEachAlgorithm[idxAlgorithmVariation][0][0]:  # compares objective function
                                                bestOfEachAlgorithm[idxAlgorithmVariation][0] = bestLine
                                        idxAlgorithmVariation = -1
                    if populations[p] == 50 and offsprings[l] == 50 or populations[p] == 50 and offsprings[l] == 200:
                        """
                        print("algorithmsVariation")
                        print(algorithmsVariations)
                        print("bestOfEachAlgorithm")
                        print(*bestOfEachAlgorithm, sep="\n")
                        """
                        # sys.exit("wtf")
                        printLatex(bestOfEachAlgorithm, algorithmsVariations, de, functionEvaluations, populations, offsprings, probCrossovers, fe, p, l, c, pm, totalSeeds, totalAlgorithmsVariations, functions, func)  # printLatex(de, dimensions, functionEvaluations, populations, offsprings, probCrossovers, d, fe, p, l, c, pm, totalSeeds, totalAlgorithmsVariations)
                        # print("Best of each algorithm")
                        # print(*bestOfEachAlgorithm, sep="\n")
