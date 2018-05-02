#!/bin/bash

#totalFunctions=18 # Funcoes
#functios=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18)

#./Otimizacao -f 1 -s 1 -p 50 -x 10 -g 1 -m 20000 -c 0 -e 0 -t 0 >> ../../../Resultados/testin.txt 

totalMethods=3 
methods=(DE ES AG)

totalFunctions=2 # Funcoes
functions=(1 2)

totalSeeds=5 # Seeds
seeds=(1 2 3 4 5) 

totalPopulation=2 # Tamanho Populacao
#populations=(50 25)
populations=(50)

totalDimensions=2 # TamX
#dimensions=(10 30)
dimensions=(10)

totalFilhosGerados=2
#filhosGerados=(1 100)
filhosGerados=(1)

totalFE=1 # Max funtions evaluations (inicalmente 1 - apenas 20000)
functionEvalutions=(20000)

totalProbCrossover=2 # Probabilidade Crossover
probCrossovers=(100 80)

totalTipoES=2 # Tipos de ES (+ ,)
tipoES=(1 2) # 1 + | 2 ,

totalSigma=2 # Sigma global
sigmas=(1 2) # 1 sigmaGlobal | (!=1)sigmaIndividual


m=0
while(($m<$totalMethods))
do	
	func=0
	while(($func<$totalFunctions))
	do
		s=0
		while(($s<$totalSeeds))
		do
			p=0
			while(($p<$totalPopulation))
			do
				d=0
				while(($d<$totalDimensions))
				do
					fi=0
					while(($fi<$totalFilhosGerados))
					do
						fe=0
						while(($fe<$totalFE))
						do
							c=2 # Apenas para AG ** Inicialize com 2 para rodar testes DE e ES **
							while(($c<$totalProbCrossover)) 
							do
								es=2 # Apenas para ES ** Inicialize com 2 para rodar testes AG e DE **
								while(($es<$totalTipoES))
								do
									sig=2 # Apenas para ES ** Inicialize com 2 para rodar testes AG e DE **
									while(($sig<$totalSigma))
									do
										# FUNCAO SEED POP N FILHOSGERADOS FE PROBCROSSOVER TIPOES SIGMAGLOBAL
										if [ ${methods[m]} = DE ];
											then
												/home/pedrohen/Documentos/PedroIC/Otimizacao/bin/Debug/Otimizacao -n ${methods[m]} -f ${functions[func]} -s ${seeds[s]} -p ${populations[p]} -x ${dimensions[d]} -g ${filhosGerados[fi]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -t ${sigmas[sig]} >> /home/pedrohen/Documentos/PedroIC/Resultados/DE_${functions[func]}_${seeds[s]}_${populations[p]}_${dimensions[d]}_${filhosGerados[fi]}_${functionEvaluations[fe]}_${probCrossovers[c]}_${tipoES[es]}_${sigmas[sig]}.txt
											elif [ ${methods[m]} = ES ];
											then
												/home/pedrohen/Documentos/PedroIC/Otimizacao/bin/Debug/Otimizacao -n ${methods[m]} -f ${functions[func]} -s ${seeds[s]} -p ${populations[p]} -x ${dimensions[d]} -g ${filhosGerados[fi]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -t ${sigmas[sig]} >> /home/pedrohen/Documentos/PedroIC/Resultados/ES_${functions[func]}_${seeds[s]}_${populations[p]}_${dimensions[d]}_${filhosGerados[fi]}_${functionEvaluations[fe]}_${probCrossovers[c]}_${tipoES[es]}_${sigmas[sig]}.txt
											elif [ ${methods[m]} = AG ];
											then
												/home/pedrohen/Documentos/PedroIC/Otimizacao/bin/Debug/Otimizacao -n ${methods[m]}-f ${functions[func]} -s ${seeds[s]} -p ${populations[p]} -x ${dimensions[d]} -g ${filhosGerados[fi]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -t ${sigmas[sig]} >> /home/pedrohen/Documentos/PedroIC/Resultados/AG_${functions[func]}_${seeds[s]}_${populations[p]}_${dimensions[d]}_${filhosGerados[fi]}_${functionEvaluations[fe]}_${probCrossovers[c]}_${tipoES[es]}_${sigmas[sig]}.txt
											else
												echo "Metodo nao encontrado"
												exit 1
										fi
										sig=$((sig+1))
									done
									es=$((es+1))
								done
								c=$((c+1))
							done
							fe=$((fe+1))
						done
						fi=$((fi+1))
					done
					d=$((d+1))
				done
				p=$((p+1))
			done
			s=$((s+1))
		done
		func=$((func+1))
	done
	m=$((m+1))
done


