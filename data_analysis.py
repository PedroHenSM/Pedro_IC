# -*- coding: utf-8 -*-
import sys
#print (sys.version_info) # Versão python
import  numpy as np

def cabecalho():
    print("C01\tMean\t\tStd\t\tBest\t\tWorst")

def imprime(de,es01,es02,es11,es12,ag):
    print("AG\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(ag),np.std(ag),np.amin(ag),np.amax(ag)))
    print("DE\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(de),np.std(de),np.amin(de),np.amax(de)))
    print("ES + G\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es01),np.std(es01),np.amin(es01),np.amax(es01)))
    print("ES + I\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es02),np.std(es02),np.amin(es02),np.amax(es02)))
    print("ES , G\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es11),np.std(es11),np.amin(es11),np.amax(es11)))
    print("ES , I\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es12),np.std(es12),np.amin(es12),np.amax(es12)))

totalMethods = 3 
methods = ['DE','ES','AG']
m = 0

totalFunctions = 1 # Funcoes
functions = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
func = 0


totalSeeds = 5 # Seeds
seeds = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30] 
s = 0

totalPopulation=1 # Tamanho Populacao
#populations=(50 25)
populations = [50,25]
p = 0

totalDimensions = 1 # TamX
#dimensions=(10 30)
dimensions=[10,30]
d = 0

totalFilhosGerados = 1
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


#DE FUNC _ SEED _ POP _ X _ FILHOS _ FE  ** IMPRESSAO TXT **
deFo= []
deV = []

cabecalho()
'''
for func in range (totalFunctions):
    for s in range (totalSeeds):
        for p in range (totalPopulation):
            for d in range(totalDimensions):
                for fi in range (totalFilhosGerados):
                    for fe in range (totalFE):
                        file = open("/home/pedrohen/Documentos/PedroIC/Resultados/DE/DE_{}_{}_{}_{}_{}_{}.txt".format(functions[func],seeds[s],populations[p],dimensions[d],filhosGerados[fi],functionEvaluations[fe]),"r")
                        text = file.readlines() # Le arquivo txt
                        lastLine = text[-1] # Ultima linha do txt
                        values = lastLine.split("\t") # Separa ultima linha em valores e salva em array
                        deFo.append(float(values[-1])) # Adiciona FO ao vetor
                        deV.append(float(values[1])) # Adiciona violacao ao vetor
                        print("DE\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(deFo),np.std(deFo),np.amin(deFo),np.amax(deFo)))
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
                        
                #print("DE\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(deFo),np.std(deFo),np.amin(deFo),np.amax(deFo)))
    imprime(de,es_0_1,es_0_2,es_1_1,es_1_2,ag)
    #print("DE\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(de),np.std(de),np.amin(de),np.amax(de)))
    #print("AG\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(ag),np.std(ag),np.amin(ag),np.amax(ag)))
    #print("ES01\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es01),np.std(es01),np.amin(es01),np.amax(es01)))
    #print("ES02\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es02),np.std(es02),np.amin(es02),np.amax(es02)))
    #print("ES11\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es11),np.std(es11),np.amin(es11),np.amax(es11)))
    #print("ES12\t{:e}\t{:e}\t{:e}\t{:e}".format(np.mean(es12),np.std(es12),np.amin(es12),np.amax(es12)))

print("DE\t\\textbf{{{:e}}}\t{:e}\t{:e}\t{:e}".format(np.mean(de),np.std(de),np.amin(de),np.amax(de)))
#\textbf
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