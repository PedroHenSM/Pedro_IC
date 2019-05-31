#!/bin/bash

<<COMMENT
dont need it

totalPenaltyMethods=1 # PenaltyMethods | 1 Deb Penalty and 2 for APM
penaltyMethods=(1 2)

totalOffsprings=1 # offspringsSize 4
offsprings=(50 100 200 400)

totalFE=1 # Max funtions evaluations
functionEvaluations=(20000)

totalProbCrossover=2 # Probabilidade Crossover 2
probCrossovers=(80 100)

totalTipoES=2 # Tipos de ES (+ ,) 2
tipoES=(0 1) # 0 + | 1 ,

totalSigma=2 # Sigma global 2
sigmas=(1 2) # 1 sigmaGlobal | (!=1)sigmaIndividual

totalPopulation=1 # parentsSize 2
populations=(50)

COMMENT

# totalAlgorithms=3 # 3
# algorithms=(DE ESCMA)

totalFunctions=2 # Funcoes 
functions=(81 83 84 85 86 87 88 89 810 811 812 813 814 815 816 817 818 819 820 821 822 823 824 825 826 827 828 829 830)

totalSeeds=3 # Seeds 30
seeds=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30) 

totalDimensions=2 
dimensions=(2 10 50 100)

n=0
while(($n<$totalDimensions))
do
	s=0
	while(($s<$totalSeeds))
	do
		f=0
		while(($f<$totalFunctions))
		do
			echo "Executing function ${functions[f]} with seed ${seeds[s]} and dimension ${dimensions[n]}"
			python3 /home/pedrohen/Documents/Pedro_IC/individual.py -f ${functions[f]} -s ${seeds[s]} -n ${dimensions[n]} >> data_results/f${functions[f]}n${dimensions[n]}
			f=$((f+1))
		done
		s=$((s+1))
	done
	n=$((n+1))
done

