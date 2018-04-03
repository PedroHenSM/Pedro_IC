#AG N 'X' POP 'Y' PROB_CROSSOVER 'Z'
#Ex: ; 2 100
import numpy as np
import pandas as pd
#X = pd.read_csv("teste.txt", sep='\t')
X = pd.read_table("saida.txt",header=-1,sep='\t')
X = np.array(X)
#print(X)
#print(X[0,0])
valores = []
contSeed=0 # numero de seeds
print ("TipoAG\tMediana\tMédia\tDesvioPadrão\tMínimo\tMáximo")
for i in range (len(X)):
	#print("ENtrou for")
	#print(i)
	#print("Linha:",i,"Valor:",X[i,0])
	if X[i,1] == 100000: # Num max de contaObj
		contSeed=contSeed+1	
		valores.append(X[i,-1])
	if contSeed == 30:
		tipoAG= X[i,0].split()
		print("AG N",tipoAG[0] + " POP",tipoAG[1] + " PROB",tipoAG[2],end = '\t')
		print('{:f}'.format(np.median(valores)),'{:f}'.format(np.mean(valores)),'{:f}'.format(np.std(valores)),'{:f}'.format(np.amin(valores)),'{:f}'.format(np.amax(valores)),sep = '\t')
		contSeed=0
		valores=[]
#\textbf{numero} Negrito em latex
'''
print(X[0,0][-1])
splitado = X[0,0].split()
print (splitado)
print("ES",splitado[0] + " N",splitado[1]+ " POP",splitado[2])
'''
