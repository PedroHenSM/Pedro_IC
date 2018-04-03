#TipoES ';/+' N 'x' TAM' 'y'
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
print ("TipoES\tMediana\tMédia\tDesvioPadrão\tMínimo\tMáximo")
for i in range (len(X)):
	#print("ENtrou for")
	#print(i)
	#print("Linha:",i,"Valor:",X[i,0])
	if X[i,1] == 100100: # Num max de contaObj
		contSeed=contSeed+1	
		valores.append(X[i,-1])
	if contSeed == 30:
		tipoES= X[i,0].split()
		print("ES 1",tipoES[0] + " N",tipoES[1]+ " POP",tipoES[2],end = '\t')
		print('{:f}'.format(np.median(valores)),'{:f}'.format(np.mean(valores)),'{:f}'.format(np.std(valores)),'{:f}'.format(np.amin(valores)),'{:f}'.format(np.amax(valores)),sep = '\t')
		#print("%f\t%f\t%f\t%f\t%f", %(np.median(valores),np.mean(valores),np.std(valores),np.amin(valores),np.amax(valores)))
		#print("%s\t%s",%(a,b))
#{:.2f}.format( valor ) 
		contSeed=0
		valores=[]
#\textbf{numero} Negrito em latex
'''
print(X[0,0][-1])
splitado = X[0,0].split()
print (splitado)
print("ES",splitado[0] + " N",splitado[1]+ " POP",splitado[2])
'''
