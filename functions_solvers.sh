#!/bin/bash

#./Otimizacao -f 1 -s 1 -p 50 -x 10 -g 1 -m 20000 -c 0 -e 0 -t 0 >> ../../../Resultados/testin.txt 

totalMethods=3 
methods=(DE ES AG)

totalFunctions=2 # Funcoes
functions=(1 2)

totalSeeds=5 # Seeds
seeds=(1 2 3 4 5) 

totalPopulation=1 # Tamanho Populacao
#populations=(50 25)
populations=(50)

totalDimensions=1 # TamX
#dimensions=(10 30)
dimensions=(10)

totalFilhosGerados=1
#filhosGerados=(1 100)
filhosGerados=(1)

totalFE=1 # Max funtions evaluations (inicalmente 1 - apenas 20000)
functionEvaluations=(20000)

totalProbCrossover=1 # Probabilidade Crossover
#probCrossovers=(100 80)
probCrossovers=(100)

totalTipoES=2 # Tipos de ES (+ ,)
tipoES=(0 1) # 0 + | 1 ,

totalSigma=2 # Sigma global
sigmas=(1 2) # 1 sigmaGlobal | (!=1)sigmaIndividual

echo "Rodou antes dos whiles"
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
							c=0 # Apenas para AG ** Inicialize com 2 para rodar testes DE e ES **
							while(($c<$totalProbCrossover)) 
							do
								es=0 # Apenas para ES ** Inicialize com 2 para rodar testes AG e DE **
								while(($es<$totalTipoES))
								do
									sig=0 # Apenas para ES ** Inicialize com 2 para rodar testes AG e DE **
									while(($sig<$totalSigma))
									do
										# FUNCAO SEED POP N FILHOSGERADOS FE PROBCROSSOVER TIPOES SIGMAGLOBAL
										if [ ${methods[m]} = DE -a $c -eq 0 -a $es -eq 0 -a $sig -eq 0 ];
										then
											# DE FUNC _ SEED _ POP _ X _ FILHOS _ FE  ** IMPRESSAO TXT **
											echo "Executando DE"
											/home/pedrohen/Documentos/PedroIC/Otimizacao/bin/Debug/Otimizacao -n ${methods[m]} -f ${functions[func]} -s ${seeds[s]} -p ${populations[p]} -x ${dimensions[d]} -g ${filhosGerados[fi]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -t ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/DE/DE_${functions[func]}_${seeds[s]}_${populations[p]}_${dimensions[d]}_${filhosGerados[fi]}_${functionEvaluations[fe]}.txt
											#/home/pedrohen/Documentos/PedroIC/Otimizacao/bin/Debug/Otimizacao -n ${methods[m]} -f ${functions[func]} -s ${seeds[s]} -p ${populations[p]} -x ${dimensions[d]} -g ${filhosGerados[fi]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -t ${sigmas[sig]} >> /home/pedrohen/Documentos/PedroIC/Resultados/DE/DE_${functions[func]}_${seeds[s]}_${populations[p]}_${dimensions[d]}_${filhosGerados[fi]}_${functionEvaluations[fe]}_${probCrossovers[c]}_${tipoES[es]}_${sigmas[sig]}.txt
										elif [ ${methods[m]} = ES -a $c -eq 0 ];
										then
											# es FUNC _ SEED _ POP _ X _ FILHOS _ FE _ TIPOES _ SIGMAGLOBAL  ** IMPRESSAO TXT **
											echo "Executando ES"
											if [ $es -eq 0 -a $sig -eq 0 ]; # ES + e Sigma Global
											then
												/home/pedrohen/Documentos/PedroIC/Otimizacao/bin/Debug/Otimizacao -n ${methods[m]} -f ${functions[func]} -s ${seeds[s]} -p ${populations[p]} -x ${dimensions[d]} -g ${filhosGerados[fi]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -t ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/ES/ES0/sigmaGlobal/ES_${functions[func]}_${seeds[s]}_${populations[p]}_${dimensions[d]}_${filhosGerados[fi]}_${functionEvaluations[fe]}_${tipoES[es]}_${sigmas[sig]}.txt
											elif [ $es -eq 0 -a $sig -eq 1 ];	# ES + e Sigma Individual
											then
												/home/pedrohen/Documentos/PedroIC/Otimizacao/bin/Debug/Otimizacao -n ${methods[m]} -f ${functions[func]} -s ${seeds[s]} -p ${populations[p]} -x ${dimensions[d]} -g ${filhosGerados[fi]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -t ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/ES/ES0/sigmaIndividual/ES_${functions[func]}_${seeds[s]}_${populations[p]}_${dimensions[d]}_${filhosGerados[fi]}_${functionEvaluations[fe]}_${tipoES[es]}_${sigmas[sig]}.txt									
											elif [ $es -eq 1 -a $sig -eq 0 ];
											then
												/home/pedrohen/Documentos/PedroIC/Otimizacao/bin/Debug/Otimizacao -n ${methods[m]} -f ${functions[func]} -s ${seeds[s]} -p ${populations[p]} -x ${dimensions[d]} -g ${filhosGerados[fi]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -t ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/ES/ES1/sigmaGlobal/ES_${functions[func]}_${seeds[s]}_${populations[p]}_${dimensions[d]}_${filhosGerados[fi]}_${functionEvaluations[fe]}_${tipoES[es]}_${sigmas[sig]}.txt
											elif [ $es -eq 1 -a $sig -eq 1 ];
											then
												/home/pedrohen/Documentos/PedroIC/Otimizacao/bin/Debug/Otimizacao -n ${methods[m]} -f ${functions[func]} -s ${seeds[s]} -p ${populations[p]} -x ${dimensions[d]} -g ${filhosGerados[fi]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -t ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/ES/ES1/sigmaIndividual/ES_${functions[func]}_${seeds[s]}_${populations[p]}_${dimensions[d]}_${filhosGerados[fi]}_${functionEvaluations[fe]}_${tipoES[es]}_${sigmas[sig]}.txt
											fi
										elif [ ${methods[m]} = AG -a $es -eq 0 -a $sig -eq 0 ];
										then
											# AG FUNC _ SEED _ POP _ X _ FILHOS _ FE _ PROBCROSSOVER  ** IMPRESSAO TXT **
											echo "Executando AG"
											/home/pedrohen/Documentos/PedroIC/Otimizacao/bin/Debug/Otimizacao -n ${methods[m]} -f ${functions[func]} -s ${seeds[s]} -p ${populations[p]} -x ${dimensions[d]} -g ${filhosGerados[fi]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -t ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/AG/AG_${functions[func]}_${seeds[s]}_${populations[p]}_${dimensions[d]}_${filhosGerados[fi]}_${functionEvaluations[fe]}_${probCrossovers[c]}.txt
										else
											echo "Metodo nao encontrado"
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

echo "Finalizou Execucoes"

