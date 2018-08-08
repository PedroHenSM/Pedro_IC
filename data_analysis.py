# -*- coding: utf-8 -*-
#import sys
#print (sys.version_info) # Versão python
import  numpy as np

def cabecalho(functions,func):
    print("C{:02}\tMean\t\tStd\t\tBest\t\tWorst".format(functions[func]))
    
def headerLatex(functions,func):
    print("\\begin{table}[h]") # h de Here
    print("\centering")
    print("\caption{{Function C{:02}}}".format(functions[func]))
    #print("\\label{Label}")
    print("\\vspace{0.5cm}")
    print("\\begin{tabular}{@{} l | r r r r @{}}") # r Divide primeira coluna das outras e rl alinha a direita
    print("\hline")
    print("Algorithm & Mean & Std & Best & Worst \\\\")
    print("\hline")
    
def footerLatex():
    print("\end{tabular}")
    print("\end{table}")
    print("\\textbf{{N}}: {}\t\\textbf{{FE}}: {}\t\\textbf{{Population}}: {}\t\\textbf{{Children/Gen.}}: {}\t\\textbf{{Crossover}}: {}\%\t\\textbf{{PenalthyMethod}}: {} \\\\ ".format(dimensions[d],functionEvaluations[fe],populations[p],filhosGerados[fi],probCrossovers[c],penaltyMethodsStr[pm]))
    print("\n\n")

def latexModel():
    headerLatex(functions,func)
    footerLatex()

def imprime(de,es01,es02,es11,es12,ag):
    print("DE\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(de),np.std(de),np.amin(de),np.amax(de)))
    print("AG\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(ag),np.std(ag),np.amin(ag),np.amax(ag)))
    print("ES + G\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es01),np.std(es01),np.amin(es01),np.amax(es01)))
    print("ES + I\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es02),np.std(es02),np.amin(es02),np.amax(es02)))
    print("ES , G\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es11),np.std(es11),np.amin(es11),np.amax(es11)))
    print("ES , I\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es12),np.std(es12),np.amin(es12),np.amax(es12)))
    
def printLatex(de,es01,es02,es11,es12,ag):
    algorithmsStr = ["DE","AG", "ES + G","ES + I","ES , G","ES , I"]
    
    algorithms = []
    algorithms.append(de)
    algorithms.append(ag)
    algorithms.append(es01)
    algorithms.append(es02)
    algorithms.append(es11)
    algorithms.append(es12)

    
    #print(algorithms)

    means = []
    stds = []
    bests = []
    worsts = []
    for i in range (6):
        means.append(np.mean(algorithms[i]))
        stds.append(np.std(algorithms[i]))
        bests.append(np.amin(algorithms[i]))
        worsts.append(np.amax(algorithms[i]))
        
    values = []
    values.append(means)
    values.append(stds)
    values.append(bests)
    values.append(worsts)

    #print(values)
    #print(values[1][0])
    bestValues = [np.mean(de),np.std(de),np.amin(de),np.amax(de)]
    #print(bestValues)
    for m in range (1,6): # Pick best means | stds | bests | worsts
        if np.mean(algorithms[m]) < bestValues[0]: # Mean
            bestValues[0] = np.mean(algorithms[m])
        if np.std(algorithms[m]) < bestValues[1]: # Std
            bestValues[1] = np.std(algorithms[m])
        if np.amin(algorithms[m]) < bestValues[2]: # Best
            bestValues[2] = np.amin(algorithms[m])
        if np.amax(algorithms[m]) < bestValues[3]: # Worst
            bestValues[3] = np.amax(algorithms[m])

    headerLatex(functions,func)
    
    for m in range(6):
        print(algorithmsStr[m]+" &", end = ' ')
        for p in range(4):
            if bestValues[p] == values[p][m]:
                if p == 3:
                    print("\\textbf{{{:e}}}".format(values[p][m]),end = ' ')
                #print("\\text{{{:e}}}\t".format(values[p][m]), end = ' ')
                else:
                    print("\\textbf{{{:e}}} &".format(values[p][m]),end = ' ')
            #else:
                #print("{:e}\t".format(values[p][m]),end = ' ')
            else:
                if p == 3:
                    print("{:e}".format(values[p][m]),end = ' ')
                else:
                    print("{:e} &".format(values[p][m]),end = ' ')    
        print("\\\\")

    footerLatex()
    
    
totalMethods = 3
methods = ['DE','ES','AG']
m = 0

totalFunctions = 18 # Funcoes
functions = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
func = 0


totalSeeds = 30 # Seeds
seeds = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30] 
s = 0

totalPenaltyMethods = 2 # Penalty Methods
penaltyMethods = [1,2]
penaltyMethodsStr = ["Padrao", "APM"]
pm = 0

totalPopulation = 2 # Tamanho Populacao
populations = [50,25]
p = 0

totalDimensions = 1 # TamX
#dimensions=(10 30)
dimensions=[10,30]
d = 0

totalFilhosGerados = 2
#filhosGerados=(1 100)
filhosGerados = [1,100]
fi = 0

totalFE = 1 # Max funtions evaluations (inicalmente 1 - apenas 20000)
functionEvaluations = [20000,100000,500000]
fe = 0

totalProbCrossover = 1 # Probabilidade Crossover
#probCrossovers=(100 80)
probCrossovers=[100,80]
c = 0

totalTipoES = 2 # Tipos de ES (+ ,)
tipoES = [0,1] # 0 + | 1 ,
es = 0

totalSigma = 2 # Sigma global
sigmas = [1,2] # 1 sigmaGlobal | (!=1)sigmaIndividual
sig = 0


for func in range (totalFunctions):
    for fi in range (totalFilhosGerados):
        for p in range(totalPopulation):
            for pm in range(totalPenaltyMethods):
                de = [] # DE
                es_0_1 = [] # ES 0 sigma global (+ 1)
                es_0_2 = [] # ES 0 sigma indivudal (+ 2)
                es_1_1 = [] # ES 1 sigma global (, 1)
                es_1_2 = [] # ES 1 sigma individual (, 2)
                ag = [] # AG
                for m in range(totalMethods):
                    for es in range(totalTipoES):
                        for sig in range(totalSigma):
                            for s in range (totalSeeds):
                                if methods[m] == 'DE' and es == 0 and sig == 0:
                                    file = open("/home/pedrohen/Documentos/PedroIC/Resultados/DE/DE_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],penaltyMethods[pm],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe]),"r")
                                    text = file.readlines() # Le arquivo txt
                                    lastLine = text[-1] # Ultima linha do txt
                                    values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                                    de.append(float(values[-1])) # Adiciona FO ao vetor
                                elif methods[m] == 'ES': #ES_FUNC _ SEED _ POP _ X _ FILHOS _ FE _ TIPOES _ SIGMAGLOBAL
                                    if tipoES[es] == 0: # Es0 +
                                        if sigmas[sig] == 1: # Sigma Global
                                            file = open("/home/pedrohen/Documentos/PedroIC/Resultados/ES/ES0/sigmaGlobal/ES_{}_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],penaltyMethods[pm],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe],tipoES[es],sigmas[sig]),"r")
                                            text = file.readlines() # Le arquivo txt
                                            lastLine = text[-1] # Ultima linha do txt
                                            values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                                            es_0_1.append(float(values[-1])) # Adiciona FO ao vetor
                                        elif sigmas[sig] == 2: # Sigma Individual
                                            file = open("/home/pedrohen/Documentos/PedroIC/Resultados/ES/ES0/sigmaIndividual/ES_{}_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],penaltyMethods[pm],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe],tipoES[es],sigmas[sig]),"r")
                                            text = file.readlines() # Le arquivo txt
                                            lastLine = text[-1] # Ultima linha do txt
                                            values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                                            es_0_2.append(float(values[-1])) # Adiciona FO ao vetor
                                    elif tipoES[es] == 1: # Es1 ,
                                        if sigmas[sig] == 1: # Sigma Global
                                            file = open("/home/pedrohen/Documentos/PedroIC/Resultados/ES/ES1/sigmaGlobal/ES_{}_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],penaltyMethods[pm],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe],tipoES[es],sigmas[sig]),"r")
                                            text = file.readlines() # Le arquivo txt
                                            lastLine = text[-1] # Ultima linha do txt
                                            values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                                            es_1_1.append(float(values[-1])) # Adiciona FO ao vetor
                                        elif sigmas[sig] == 2: # Sigma Individual
                                            file = open("/home/pedrohen/Documentos/PedroIC/Resultados/ES/ES1/sigmaIndividual/ES_{}_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],penaltyMethods[pm],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe],tipoES[es],sigmas[sig]),"r")
                                            text = file.readlines() # Le arquivo txt
                                            lastLine = text[-1] # Ultima linha do txt
                                            values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                                            es_1_2.append(float(values[-1])) # Adiciona FO ao vetor
                                elif methods[m] == 'AG' and es == 0 and sig == 0: #AG_FUNC _ SEED _ POP _ X _ FILHOS _ FE _ PROBCROSSOVER
                                    file = open("/home/pedrohen/Documentos/PedroIC/Resultados/AG/AG_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],penaltyMethods[pm],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe],probCrossovers[c]),"r")
                                    text = file.readlines() # Le arquivo txt
                                    lastLine = text[-1] # Ultima linha do txt
                                    values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                                    ag.append(float(values[-1])) # Adiciona FO ao vetor
                if ((populations[p] == 50 and filhosGerados[fi] == 1) or (populations[p] == 25 and filhosGerados[fi] == 100)):
                    # Imprime apenas as relações de tamPop 50 e filhos 1 e tamPop 25 e filhos 100, para apm e padrao
                    printLatex(de,es_0_1,es_0_2,es_1_1,es_1_2,ag)


'''
for func in range (totalFunctions):
    de = [] # DE
    es_0_1 = [] # ES 0 sigma global (+ 1)
    es_0_2 = [] # ES 0 sigma indivudal (+ 2)
    es_1_1 = [] # ES 1 sigma global (, 1)
    es_1_2 = [] # ES 1 sigma individual (, 2)
    ag = [] # AG
    for m in range(totalMethods):
        for es in range(totalTipoES):
            for sig in range(totalSigma):
                for s in range (totalSeeds):
                    if methods[m] == 'DE' and es == 0 and sig == 0:
                        file = open("/home/pedrohen/Documentos/PedroIC/Resultados/DE/DE_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe]),"r")
                        text = file.readlines() # Le arquivo txt
                        lastLine = text[-1] # Ultima linha do txt
                        values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                        de.append(float(values[-1])) # Adiciona FO ao vetor
                    elif methods[m] == 'ES': #ES_FUNC _ SEED _ POP _ X _ FILHOS _ FE _ TIPOES _ SIGMAGLOBAL
                        if tipoES[es] == 0: # Es0 +
                            if sigmas[sig] == 1: # Sigma Global
                                file = open("/home/pedrohen/Documentos/PedroIC/Resultados/ES/ES0/sigmaGlobal/ES_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe],tipoES[es],sigmas[sig]),"r")
                                text = file.readlines() # Le arquivo txt
                                lastLine = text[-1] # Ultima linha do txt
                                values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                                es_0_1.append(float(values[-1])) # Adiciona FO ao vetor
                            elif sigmas[sig] == 2: # Sigma Individual
                                file = open("/home/pedrohen/Documentos/PedroIC/Resultados/ES/ES0/sigmaIndividual/ES_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe],tipoES[es],sigmas[sig]),"r")
                                text = file.readlines() # Le arquivo txt
                                lastLine = text[-1] # Ultima linha do txt
                                values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                                es_0_2.append(float(values[-1])) # Adiciona FO ao vetor
                        elif tipoES[es] == 1: # Es1 ,
                            if sigmas[sig] == 1: # Sigma Global
                                file = open("/home/pedrohen/Documentos/PedroIC/Resultados/ES/ES1/sigmaGlobal/ES_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe],tipoES[es],sigmas[sig]),"r")
                                text = file.readlines() # Le arquivo txt
                                lastLine = text[-1] # Ultima linha do txt
                                values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                                es_1_1.append(float(values[-1])) # Adiciona FO ao vetor
                            elif sigmas[sig] == 2: # Sigma Individual
                                file = open("/home/pedrohen/Documentos/PedroIC/Resultados/ES/ES1/sigmaIndividual/ES_{}_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe],tipoES[es],sigmas[sig]),"r")
                                text = file.readlines() # Le arquivo txt
                                lastLine = text[-1] # Ultima linha do txt
                                values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                                es_1_2.append(float(values[-1])) # Adiciona FO ao vetor
                    elif methods[m] == 'AG' and es == 0 and sig == 0: #AG_FUNC _ SEED _ POP _ X _ FILHOS _ FE _ PROBCROSSOVER
                        file = open("/home/pedrohen/Documentos/PedroIC/Resultados/AG/AG_{}_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe],probCrossovers[c]),"r")
                        text = file.readlines() # Le arquivo txt
                        lastLine = text[-1] # Ultima linha do txt
                        values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                        ag.append(float(values[-1])) # Adiciona FO ao vetor
    printLatex(de,es_0_1,es_0_2,es_1_1,es_1_2,ag)
'''



'''
texto = file.readlines() # Le arquivo txt
line = texto[-1] # Pega ultima linha
#file = open("/home/pedro/Área de Trabalho/Diversos/IC/MinimizarFuncao/Resultado-{:02d}-{}-{:02d}.txt".format(j,tam,k),'r')
print(texto)
print(line)
values = line.split("\t")
print(values)
fo1 = values[-1]
fo2 = values[1]
print(fo2)
print(fo1)

res = []
res.append(float(fo1))
res.append(float(fo2))
print(res)
print('%e' %(np.mean(res)))'''