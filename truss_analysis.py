import numpy as np
# noinspection PyUnresolvedReferences
import sys


# noinspection PyShadowingNames
def defineParameters():
    totalAlgorithms = 3
    algorithms = ['DE', 'ES', 'GA']
    a = 0

    totalFunctions = 1  # Funcoes
    functions = [210, 225, 260, 272, 2942]
    func = 0

    totalSeeds = 3  # Seeds
    seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
    s = 0

    totalPenaltyMethods = 1  # Penalty Methods
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
def headerLatex(functions, func):
    print("\\begin{table}[h]")  # h de Here
    print("\centering")
    print("\caption{{Function C{:02}}}".format(functions[func]))
    # print("\\label{Label}")
    print("\\vspace{0.5cm}")
    print("\\begin{tabular}{@{} l | r r r r @{}}")  # r Divide primeira coluna das outras e rl alinha a direita
    print("\hline")
    print("Algorithm & Mean & Std & Best & Worst \\\\")
    print("\hline")


# noinspection PyShadowingNames
def footerLatex(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm):
    print("\end{tabular}")
    print("\end{table}")
    print("\\textbf{{FE}}: {}\t\\textbf{{Population}}: {}\t\\textbf{{Offsprings.}}: {}\t\\textbf{{Crossover}}: {}\%\t\\textbf{{PenalthyMethod}}: {} \\\\ ".format(functionEvaluations[fe], populations[p], offsprings[l], probCrossovers[c], penaltyMethodsStr[pm]))
    print("\n\n")


# noinspection PyShadowingNames
def latexModel(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm):
    headerLatex(functions, func)
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
def printLatex(algorithmsVariations, de, functionEvaluations, populations, offsprings, probCrossovers, fe, p, l, c, pm, totalSeeds, totalAlgorithmsVariations):
    algorithmsStr = ["DE", "GA", "ES + G", "ES + I", "ES , G", "ES , I"]
    penaltyMethodsStr = ["Deb", "APM"]
    means = []
    stds = []
    bests = []
    worsts = []
    # objectiveFunctions = [[]] * totalAlgorithmsVariations # doesnt work
    objectiveFunctions = [[] for i in range(totalAlgorithmsVariations)]
    # print("objectionFunctions")
    print(objectiveFunctions)
    # print(*algorithmsVariations, sep="\n")
    for a in range(totalAlgorithmsVariations):
        for i in range(totalSeeds):
            # print("algorithmsVariations[{}][{}][0]: {}".format(a, i, algorithmsVariations[a][i][0]))
            objectiveFunctions[a].append(algorithmsVariations[a][i][0])
        means.append(np.mean(objectiveFunctions[a]))
        stds.append(np.std(objectiveFunctions[a]))
        bests.append(np.amin(objectiveFunctions[a]))
        worsts.append(np.amax(objectiveFunctions[a]))
    """
    print(*algorithmsVariations, sep="\n")
    print("objectiveFunction[a]")
    print(objectiveFunctions[0])
    print("objectiveFunction")
    print(objectiveFunctions, sep="\n")
    """
    sys.exit("tchaa")
    values = ([])
    values.append(means)
    values.append(stds)
    values.append(bests)
    values.append(worsts)
    bestValues = [values[0][0], values[1][0], values[2][0], values[3][0]]  # means, stds, bests, worsts | first idx is the analysis(means,std,etc)and second is the algorithm(de,ga,etc)

    for a in range(1, totalAlgorithmsVariations):  # Pick best means | stds | bests | worsts
        if np.mean(objectiveFunctions[a]) < bestValues[0]:  # Mean
            bestValues[0] = np.mean(objectiveFunctions[a])
        if np.std(objectiveFunctions[a]) < bestValues[1]:  # Std
            bestValues[1] = np.std(objectiveFunctions[a])
        if np.amin(objectiveFunctions[a]) < bestValues[2]:  # Best
            bestValues[2] = np.amin(objectiveFunctions[a])
        if np.amax(objectiveFunctions[a]) < bestValues[3]:  # Worst
            bestValues[3] = np.amax(objectiveFunctions[a])
    headerLatex(functions, func)
    for a in range(totalAlgorithmsVariations):
        print(algorithmsStr[a] + " &", end=' ')
        for j in range(4):
            if bestValues[j] == values[j][a]:  # if actual algorithm is the best, put bold
                if j == 3:  # just for formatting
                    print("\\textbf{{{:e}}}".format(values[j][a]), end=' ')
                else:
                    print("\\textbf{{{:e}}} &".format(values[j][a]), end=' ')
            else:
                if j == 3:  # just for formatting
                    print("{:e}".format(values[j][a]), end=' ')
                else:
                    print("{:e} &".format(values[j][a]), end=' ')
        print("\\\\")

    footerLatex(functionEvaluations, populations, offsprings, probCrossovers, penaltyMethodsStr, fe, p, l, c, pm)


if __name__ == '__main__':
    algorithms, totalAlgorithms, totalFunctions, functions, totalSeeds, seeds, totalPenaltyMethods, penaltyMethods, totalPopulation, populations, totaloffsprings, offsprings, totalFE, functionEvaluations, totalProbCrossover, probCrossovers, totalTipoES, tipoES, totalSigma, sigmas = defineParameters()
    fe = c = 0
    print(totalSeeds)
    print(algorithms)
    print(fe)
    for func in range(totalFunctions):
        for l in range(totaloffsprings):
            for p in range(totalPopulation):
                for pm in range(totalPenaltyMethods):
                    totalAlgorithmsVariations = 6
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
                                        if float(splittedLine[0]) == 0 and splittedLine[1] < bestLine[1]:  # Tests if violation is 0 and then compares objecive Function
                                            bestLine = splittedLine
                                    bestLine = bestLine[1:-1]  # removes violation and '\n' from list, now the first idx is the ObjFunction
                                    bestLine = [float(i) for i in bestLine]  # converts all str values to float
                                    # print("Fora do printLatex")
                                    # print(*algorithmsVariations, sep="\n")
                                    if idxAlgorithmVariation != -1:
                                        algorithmsVariations[idxAlgorithmVariation].append(bestLine)
                                        idxAlgorithmVariation = -1
                    if populations[p] == 50 and offsprings[l] == 50 or populations[p] == 50 and offsprings[l] == 200:
                        printLatex(algorithmsVariations, de, functionEvaluations, populations, offsprings, probCrossovers, fe, p, l, c, pm, totalSeeds, totalAlgorithmsVariations)  # printLatex(de, dimensions, functionEvaluations, populations, offsprings, probCrossovers, d, fe, p, l, c, pm, totalSeeds, totalAlgorithmsVariations)
