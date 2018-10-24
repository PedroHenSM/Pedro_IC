#!/bin/bash

totalAlgorithms=3 # 3
algorithms=(DE ES GA)

totalFunctions=1 # Funcoes 5
functions=(210 225 260 272 2942)

totalSeeds=3 # Seeds 30
seeds=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30) 

totalPenaltyMethods=1 # PenaltyMethods | 1 Deb Penalty and 2 for APM
penaltyMethods=(1 2)

totalPopulation=1 # parentsSize 2
populations=(50 100)

<<COMMENT
dont need it
totalDimensions=1 # will be subscribed on execution TODO: Dont really need it?!
dimensions=(10)
COMMENT

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


a=0
while(($a<$totalAlgorithms))
do	
	func=0
	while(($func<$totalFunctions))
	do
		s=0
		while(($s<$totalSeeds))
		do
			p=0
			while(($p<$totalPenaltyMethods))
			do
				u=0
				while(($u<$totalPopulation))
				do
					l=0
					while(($l<$totalOffsprings))
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
										echo "Executando SIG:$sig ES:$es Cross:$c FE:$fe Filhos:$fg N:$d POP:$u SEED:$s FUNC:$func Methods:$a"
										# FUNCAO SEED POP N FILHOSGERADOS FE PROBCROSSOVER TIPOES SIGMAGLOBAL
										if [ ${algorithms[a]} = DE -a $c -eq 0 -a $es -eq 0 -a $sig -eq 0 ];
										then
											# DE FUNC _ SEED _ PenMethod _ POP _ FILHOS _ FE  ** IMPRESSAO TXT **
											echo "Executando DE"
											python /home/pedrohen/Documentos/PedroIC/individual.py -a ${algorithms[a]} -f ${functions[func]} -s ${seeds[s]} -p ${penaltyMethods[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -g ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/trusses/DE/DE_truss_${functions[func]}_${seeds[s]}_${penaltyMethods[p]}_${populations[u]}_${offsprings[l]}_${functionEvaluations[fe]}.txt
										elif [ ${algorithms[a]} = ES -a $c -eq 0 ];
										then
											# es FUNC _ SEED _ PenMethod _ POP _ FILHOS _ FE _ TIPOES _ SIGMAGLOBAL  ** IMPRESSAO TXT **
											echo "Executando ES"
											if [ $es -eq 0 -a $sig -eq 0 ]; # ES + e Sigma Global
											then
												python /home/pedrohen/Documentos/PedroIC/individual.py -a ${algorithms[a]} -f ${functions[func]} -s ${seeds[s]} -p ${penaltyMethods[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -g ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/trusses/ES/ES0/sigmaGlobal/ES_truss_${functions[func]}_${seeds[s]}_${penaltyMethods[p]}_${populations[u]}_${offsprings[l]}_${functionEvaluations[fe]}_${tipoES[es]}_${sigmas[sig]}.txt
											elif [ $es -eq 0 -a $sig -eq 1 ];	# ES + e Sigma Individual
											then
												python /home/pedrohen/Documentos/PedroIC/individual.py -a ${algorithms[a]} -f ${functions[func]} -s ${seeds[s]} -p ${penaltyMethods[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -g ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/trusses/ES/ES0/sigmaIndividual/ES_truss_${functions[func]}_${seeds[s]}_${penaltyMethods[p]}_${populations[u]}_${offsprings[l]}_${functionEvaluations[fe]}_${tipoES[es]}_${sigmas[sig]}.txt
											elif [ $es -eq 1 -a $sig -eq 0 ];
											then
												python /home/pedrohen/Documentos/PedroIC/individual.py -a ${algorithms[a]} -f ${functions[func]} -s ${seeds[s]} -p ${penaltyMethods[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -g ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/trusses/ES/ES1/sigmaGlobal/ES_truss_${functions[func]}_${seeds[s]}_${penaltyMethods[p]}_${populations[u]}_${offsprings[l]}_${functionEvaluations[fe]}_${tipoES[es]}_${sigmas[sig]}.txt
											elif [ $es -eq 1 -a $sig -eq 1 ];
											then
												python /home/pedrohen/Documentos/PedroIC/individual.py -a ${algorithms[a]} -f ${functions[func]} -s ${seeds[s]} -p ${penaltyMethods[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -g ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/trusses/ES/ES1/sigmaIndividual/ES_truss_${functions[func]}_${seeds[s]}_${penaltyMethods[p]}_${populations[u]}_${offsprings[l]}_${functionEvaluations[fe]}_${tipoES[es]}_${sigmas[sig]}.txt
											fi
										elif [ ${algorithms[a]} = GA -a $es -eq 0 -a $sig -eq 0 ];
										then
											# AG FUNC _ SEED _ PenMethod _ POP _ FILHOS _ FE _ PROBCROSSOVER  ** IMPRESSAO TXT **
											echo "Executando AG"
											python /home/pedrohen/Documentos/PedroIC/individual.py -a ${algorithms[a]} -f ${functions[func]} -s ${seeds[s]} -p ${penaltyMethods[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -g ${sigmas[sig]} > /home/pedrohen/Documentos/PedroIC/Resultados/trusses/GA/GA_truss_${functions[func]}_${seeds[s]}_${penaltyMethods[p]}_${populations[u]}_${offsprings[l]}_${functionEvaluations[fe]}_${probCrossovers[c]}.txt
										else
											echo "Metodo nao encontrado ou já executado. Algorithm: ${algorithms[a]}"
											#echo "-a ${algorithms[a]} -f ${functions[func]} -s ${seeds[s]} -p ${penaltyMethods[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -c ${probCrossovers[c]} -e ${tipoES[es]} -g ${sigmas[sig]}"
										fi
										sig=$((sig+1))
									done
									es=$((es+1))
								done
								c=$((c+1))
							done
							fe=$((fe+1))
						done
						l=$((l+1))
					done
					u=$((u+1))
				done
				p=$((p+1))			
			done
			s=$((s+1))
		done
		func=$((func+1))
	done
	a=$((a+1))
done


