import numpy as np
def rosenbrockFunction(n):
    sum = 0
    for i in range(n):
        sum = sum + (100*(x[i+1] - x**2) + (1-x[i])**2)
    return sum

def ESCMA(function, seed, penaltyMethod, parentsSize, nSize, offspringsSize, maxFE, crossoverProb, esType, globalSigma):  # Evolution Strategy
    strFunction = str(function)
    seed = 1
    nSize = 10
    crossoverProb = -1
    np.random.seed(seed)
    functionEvaluations = 0
    # User defined parameters
    sigma = 0.5
    # xmean = np.random.randn(nSize)  # np.random.randn(nSize, 1)
    xmean = [0.134364, 0.847434, 0.763775, 0.255069, 0.495435, 0.449491, 0.651593, 0.788723, 0.093860, 0.028347] ## just for tests
    """
    λ ≥ 2, population size, sample size, number of offspring, see (5).
    µ ≤ λ parent number, number of (positively) selected search points in the population, number
    of strictly positive recombination weights, see (6).
    """
    # Strategy parameters setting: Selection
    parentsSize = 4 + np.floor(3 * np.log(nSize))  # parentsSize is biased on nSize
    mu = parentsSize / 2  # mu is NOT offspringsSize.
    muList = [i + 1 for i in range(int(mu))]
    weights = np.log(mu+1/2)-np.log(muList)  # muXone recombination weights
    mu = np.floor(mu)
    weights = weights/np.sum(weights)
    mueff = np.power(sum(weights),2) / np.sum(np.power(weights,2))

    # Strategy parameter setting: Adaptation
    cc = (4+mueff / nSize) / (nSize+4 + 2*mueff/nSize) # time constant for cumulation for C
    cs = (mueff+2) / (nSize+mueff+5)  # t-const for cumulation for sigma control
    c1 = 2 / (np.power(nSize + 1.3, 2) + mueff)  # learning rate for rank one update of C
    cmu = np.minimum(1 - c1, 2 * (mueff - 2 + 1/mueff)/(np.power(nSize+2, 2) + 2*mueff/2))  # and for rank-mu update
    damps = 1 + 2*np.maximum(0, np.sqrt((mueff -1) / (nSize + 1)) -1 ) + cs  # damping for sigma

    # Initiliaze dynamic (internal) strategy parameters  and constants
    pc = [0] * nSize  # evolutions paths for C
    ps = [0] * nSize  # evolutions paths for sigma
    B = np.eye(nSize) # B defines de coordinate system
    D = np.eye(nSize) # diagonal matrix D defines the scaling
    AUX = (B * D)  # auxiliar tranpose matrix
    AUX = AUX.transpose()
    C = B * D * AUX  # covariance matrix
    eigenval = 0  # B and D update at counteval == 0
    chinN = nSize**(0.5) * (1-1/(4*nSize)+1 / (21*np.power(nSize, 2)))  # expectation of ||N(0,I)|| == norm(randn(N,1))

    # CODIGO ANTIGO
    """
    generatedOffspring = int(offspringsSize / parentsSize)
    lowerBound = upperBound = truss = 0
    if strFunction[0] == "2":  # solving trusses
        truss, lowerBound, upperBound = initializeTruss(function)
        nSize = truss.getDimension()
        gSize, hSize, constraintsSize = initializeConstraintsTrusses(truss)
        penaltyCoefficients = [-1 for i in range(constraintsSize)]
        avgObjFunc = -1  # will be subscribed on 'calculatePenaltyCoefficients'
        parents = Population(parentsSize, nSize, function, lowerBound, upperBound)
        offsprings = Population(offspringsSize, nSize, function, lowerBound, upperBound)
        functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize, functionEvaluations, truss)
    else:
        gSize, hSize, constraintsSize = initializeConstraints(function)  # Initialize constraints
        penaltyCoefficients = [-1 for i in range(constraintsSize)]
        avgObjFunc = 0  # will be subscribed on 'calculatePenaltyCoefficients'
        parents = Population(parentsSize, nSize, function)  # Initialize parents population
        offsprings = Population(offspringsSize, nSize, function)  # Initialize offsprings population
        functionEvaluations = parents.evaluate(parentsSize, function, nSize, gSize, hSize,
                                               functionEvaluations)  # Evaluate parents

    if penaltyMethod == 1:  # Padrao?  (not apm)
        parents.sumViolations(parentsSize, gSize, hSize)
    elif penaltyMethod == 2:  # // Adaptive Penalty Method ( APM )
        parents.uniteConstraints(parentsSize, gSize, hSize)
        avgObjFunc = parents.calculatePenaltyCoefficients(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
        parents.calculateAllFitness(parentsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
    else:
        print("Penalthy method not encountered")
        sys.exit("Penalty method not encountered")
    parents.initializeEvolutionStrategy(offsprings, nSize, parentsSize, offspringsSize, globalSigma)
    """
    # FIM CODIGO ANTIGO
    while functionEvaluations < maxFE:
        # Generate and evaluate lambda offspring
        arzAuxList = []
        arxAuxList = []
        for k in range(parentsSize):
            arz = np.random.randn(nSize)  # standard normally distributed vector  (38)
            arzAuxList.append(arz)
            arz = [-2.66652, -0.73817, 1.50790,  0.60194, -0.45066, -0.70544, -0.42442, 0.54571, 1.69134, 0.37733]
            arz1 = [0.018461, -1.164203, 1.184779, -1.219738, 0.107310, -0.580791, 0.829065, -0.083036, 1.213385, -0.17495 ]
            # np.dot(a,b) If a is an N-D array and b is a 1-D array, it is a sum product over the last axis of a and b.
            arx = xmean + sigma * (np.dot(B*D, arz))  # add mutation N(m,sigma²C) (40)
            arxAuxList.append(arx)
            # arx é  individuals[i].n
            """
            for i in range(nSize):  # Copies arx to individual[].n TODO This can be done as offsprings.individuals[k].n = [i for i in arx] ?!
                offsprings.individuals[k].n[i] = arx[i]
            """
            # TODO evaluates individuals (objective function call)
            # avalia funcao aqui com o individuo k
        # JUST FOR TEST, VALUES FROM FOR ABOVE
        arzAuxList = [[ -2.666522,  0.018461,  0.262710,  0.114931,  0.103913, -0.704036, -0.943168,  1.154129, -0.251690, -1.034194],
        [-0.738172, -1.164203,  0.334877,  0.925357,  1.364233, -1.927049, -1.136317, -0.736228, -0.053444, -0.757933],
        [1.507904,  1.184779, -2.409754, -0.331684, -2.577779, -0.435071,  0.377811, -0.108602, -1.272071, -1.859290],
        [0.601943, -1.219738, -0.422966,  0.313458,  1.460150, -1.389020, -0.836552, -0.986827,  0.671340,  0.511160],
        [-0.450661,  0.107310, -0.408185, -0.274270, -0.845093, -0.832287, -1.713540,  1.162427, -0.412784, -0.117363],
        [-0.705443, -0.580791, -0.600879, -0.732635, -0.012908,  1.837349,  0.093211,  0.532659, -0.102758,  0.297341],
        [-0.424425,  0.829065, -0.045379, -0.939494,  0.730737, -1.374715, -0.970766, -1.981590,  0.561319, -0.954346],
        [0.545705, -0.083036,  0.856319,  0.010801, -1.597652,  1.105937,  1.214798, -0.532279, -1.006279,  0.842970],
        [1.691342,  1.213385,  0.205616, -0.658486, -0.986162, -1.665159,  1.199242, -0.170409, -1.582087,  0.432230],
        [0.377328, -0.174952, -0.260771,  1.454996,  0.462593, -0.563949,  2.110802,  0.948742, -0.838144,  2.346825 ]]

        arxAuxList = [[ -1.1988966,  0.1435945,  0.2657193,  0.1918297,  0.1863209, -0.2176538, -0.3372198,  0.7114285,  0.0085193, -0.3827326],
        [0.4783477,  0.2653321,  1.0148725,  1.3101121,  1.5295500, -0.1160909,  0.2792752,  0.4793197,  0.8207119,  0.4684674],
        [1.5177266,  1.3561643, -0.4411025,  0.5979328, -0.5251148,  0.5462393,  0.9526800,  0.7094737,  0.1277391, -0.1658703],
        [0.5560404, -0.3547997,  0.0435862,  0.4117982,  0.9851442, -0.4394412, -0.1632072, -0.2383446,  0.5907389,  0.5106490],
        [0.2701045,  0.5490903,  0.2913427,  0.3583002,  0.0728884,  0.0792915, -0.3613348,  1.0766484,  0.2890430,  0.4367535],
        [0.0967695,  0.1590957,  0.1490518,  0.0831734,  0.4430372,  1.3681657,  0.4960967,  0.7158207,  0.3981118,  0.5981618],
        [0.4393806,  1.0661255,  0.6289034,  0.1818460,  1.0169612, -0.0357643,  0.1662101, -0.3392018,  0.9322524,  0.1744202],
        [1.0615760,  0.7472055,  1.2168827,  0.7941236, -0.0101029,  1.3416918,  1.3961222,  0.5225840,  0.2855838,  1.2102082],
        [0.9395308,  0.7005519,  0.1966676, -0.2353835, -0.3992216, -0.7387198,  0.6934805,  0.0086550, -0.6971837,  0.3099745],
        [0.2170117, -0.0591285, -0.1020378,  0.7558454,  0.2596442, -0.2536270,  1.0837483,  0.5027183, -0.3907246,  1.2017602 ]]


        bestIndividuals =  [[  0.2657193,   0.1863209,   0.1435945,   0.7114285,  -0.2176538],
                            [ 1.0148725,   1.5295500,   0.2653321,   0.4793197,  -0.1160909],
                            [ -0.4411025,  -0.5251148,   1.3561643,   0.7094737,   0.5462393],
                            [ 0.0435862,   0.9851442,  -0.3547997,  -0.2383446,  -0.4394412],
                            [ 0.2913427,   0.0728884,   0.5490903,   1.0766484,   0.0792915],
                            [ 0.1490518,   0.4430372,   0.1590957,   0.7158207,   1.3681657],
                            [ 0.6289034,   1.0169612,   1.0661255,  -0.3392018,  -0.0357643],
                            [ 1.2168827,  -0.0101029,   0.7472055,   0.5225840,   1.3416918],
                            [ 0.1966676,  -0.3992216,   0.7005519,   0.0086550,  -0.7387198],
                            [ -0.1020378,   0.2596442,  -0.0591285,   0.5027183,  -0.2536270]]
        # TODO Avaliar tudo ao final
        arzMatrix = np.vstack(arzAuxList)  # matrix nd.array with all values calculated above
        arxMatrix = np.vstack(arxAuxList)  # matrix nd.array with all values calculated above

        # TODO sort by fitness (TODO)
        ###  Até aqui está funcionando corretamente.
        # linha abaixo arz ordenado
        arz = [[0.262710,0.334877,-2.409754,-0.422966,-0.408185,-0.600879,-0.045379,0.856319,0.205616,-0.260771],
               [0.103913,1.364233,-2.577779,1.460150,-0.845093,-0.012908,0.730737,-1.597652,-0.986162,0.462593],
               [0.018461,-1.164203,1.184779,-1.219738,0.107310,-0.580791,0.829065,-0.083036,1.213385,-0.174952],
               [1.154129,-0.736228,-0.108602,-0.986827,1.162427,0.532659,-1.981590,-0.532279,-0.170409,0.948742],
               [-0.704036,-1.927049,-0.435071,-1.389020,-0.832287,1.837349,-1.374715,1.105937,-1.665159,-0.563949],
               [-0.251690,-0.053444,-1.272071,0.671340,-0.412784,-0.102758,0.561319,-1.006279,-1.582087,-0.838144],
               [-2.666522,-0.738172,1.507904,0.601943,-0.450661,-0.705443,-0.424425,0.545705,1.691342,0.377328],
               [0.114931,0.925357,-0.331684,0.313458,-0.274270,-0.732635,-0.939494,0.010801,-0.658486,1.454996],
               [-0.943168,-1.136317,0.377811,-0.836552,-1.713540,0.093211,-0.970766,1.214798,1.199242,2.110802],
               [-1.034194,-0.757933,-1.859290,0.511160,-0.117363,0.297341,-0.954346,0.842970,0.432230,2.346825]]

        zmean = [0.231408, 0.221386, -1.625593, -0.115071, -0.319800, -0.279611, 0.107678, -0.072482, -0.033343, 0.044362]

        y = np.subtract(xmean, xold)
        ps = (1-cs) * ps + (np.sqrt(cs*(2-cs)*mueff)) * np.dot(B, zmean)
        pc = np.array([0.292072, 0.279423, -2.051747, -0.145237, -0.403636, -0.352912, 0.135906, -0.091483, -0.042084, 0.055992])
        hsig = True



"""
CfromOCtave
9.6910e-01, 1.9975e-03, -1.6429e-02, -2.3695e-03, -6.5479e-04, -2.6784e-03, -2.4696e-03, -7.0619e-04, 8.7800e-05, 1.9568e-03
1.9975e-03, 9.8421e-01, -3.9315e-02, 1.6214e-02, -1.0326e-02, -3.7369e-03, 6.5928e-03, -9.7537e-03, -9.6414e-03, 2.9046e-03
-1.6429e-02, -3.9315e-02, 1.1233e+00, -1.0842e-02, 3.3975e-02, 2.1803e-02, -9.6463e-03, 5.8966e-03, 1.5740e-02, -3.2129e-03
-2.3695e-03, 1.6214e-02, -1.0842e-02, 9.8569e-01, -6.0540e-03, 3.1189e-03, 6.7317e-03, -1.5414e-02, -1.1925e-02, 4.0680e-03
-6.5479e-04, -1.0326e-02, 3.3975e-02, -6.0540e-03, 9.7520e-01, 4.5658e-03, -7.1157e-03, 3.1525e-03, 4.8335e-03, 5.7456e-04
-2.6784e-03, -3.7369e-03, 2.1803e-02, 3.1189e-03, 4.5658e-03, 9.7311e-01, -5.2199e-03, -3.4103e-03, -4.8728e-03, 1.7739e-03
-2.4696e-03, 6.5928e-03, -9.6463e-03, 6.7317e-03, -7.1157e-03, -5.2199e-03, 9.7774e-01, -6.1129e-03, 9.4061e-04, -1.2355e-03
-7.0619e-04, -9.7537e-03, 5.8966e-03, -1.5414e-02, 3.1525e-03, -3.4103e-03, -6.1129e-03, 9.8650e-01, 9.1551e-03, -7.3053e-03
8.7800e-05, -9.6414e-03, 1.5740e-02, -1.1925e-02, 4.8335e-03, -4.8728e-03, 9.4061e-04, 9.1551e-03, 9.7657e-01, -3.5075e-03
1.9568e-03, 2.9046e-03, -3.2129e-03, 4.0680e-03, 5.7456e-04, 1.7739e-03, -1.2355e-03, -7.3053e-03, -3.5075e-03, 9.6821e-01

    # First iteration
    B (autvetores Octave)
    -0.72718   0.17242  -0.31273  -0.27214   0.10685   0.26866  -0.35019   0.21399  -0.10109   0.08593
    0.17281   0.62453   0.08878   0.04641   0.19643   0.55325   0.23564  -0.09095   0.30960   0.25641
    -0.18405   0.13093   0.12107   0.03900   0.13285   0.09143   0.05074  -0.12633   0.23132  -0.91365
    -0.22408  -0.51201  -0.05023  -0.15262  -0.37105   0.28749   0.29538  -0.07776   0.58021   0.10746
    0.51624  -0.28968  -0.53142  -0.06404   0.18136   0.43713  -0.18610   0.24291  -0.05423  -0.20620
    0.11429   0.43677  -0.43646  -0.13089  -0.48047  -0.35394   0.19342   0.37530   0.19949  -0.11979
    -0.03829   0.06324  -0.57153   0.17280   0.10644  -0.21469  -0.11950  -0.73127   0.15994   0.07738
    -0.21446  -0.04404  -0.24470   0.28927  -0.07578   0.20765   0.66967  -0.04251  -0.55179  -0.07006
    0.13630   0.14251   0.13308  -0.17458  -0.67471   0.33211  -0.30546  -0.36873  -0.32556  -0.11013
    -0.10795  -0.00000   0.00000   0.85658  -0.24211   0.12191  -0.32236   0.22704   0.15810   0.02604
    
    B (autovetores Python)
    -0.085933	0.101089	0.213985	0.350195	0.268657	-0.861275	0.020060	0.017073	-0.128712	0.120797
    -0.256409	-0.309597	-0.090954	-0.235644	0.553249	0.043409	-0.441009	0.528551	0.213959	0.111192
    0.913647	-0.231316	-0.126331	-0.050740	0.091432	-0.141806	0.050845	0.154445	0.078247	0.177395
    -0.107464	-0.580207	-0.077765	-0.295382	0.287492	-0.107124	0.462104	-0.386375	-0.235140	-0.285180
    0.206201	0.054230	0.242906	0.186102	0.437131	0.258165	-0.428189	-0.561856	-0.036835	-0.085654
    0.119787	-0.199493	0.375300	-0.193422	-0.353937	-0.131171	-0.478148	0.277394	-0.371524	-0.466915
    -0.077383	-0.159943	-0.731273	0.119501	-0.214692	-0.211118	-0.405529	-0.180229	-0.355367	0.185504
    0.070061	0.551793	-0.042510	-0.669673	0.207647	-0.160305	-0.023279	-0.067838	-0.398077	0.195518
    0.110132	0.325564	-0.368731	0.305459	0.332115	0.163408	0.085772	0.294033	-0.158742	-0.624473
    -0.026038	-0.158103	0.227039	0.322362	0.121914	0.209550	0.068264	0.176165	-0.653380	0.415086
    
    D (autovalores Octave)
    0.98212         0         0         0         0         0         0         0         0         0
    0   0.98212         0         0         0         0         0         0         0         0
    0         0   0.98212         0         0         0         0         0         0         0
    0         0         0   0.98212         0         0         0         0         0         0
    0         0         0         0   0.98212         0         0         0         0         0
    0         0         0         0         0   0.98487         0         0         0         0
    0         0         0         0         0         0   0.98775         0         0         0
    0         0         0         0         0         0         0   0.99191         0         0
    0         0         0         0         0         0         0         0   1.00816         0
    0         0         0         0         0         0         0         0         0   1.07280
    
    D (autovalores Python)
    1.072803	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000
    0.000000	1.008161	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000
    0.000000	0.000000	0.991906	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000
    0.000000	0.000000	0.000000	0.987753	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000
    0.000000	0.000000	0.000000	0.000000	0.984874	0.000000	0.000000	0.000000	0.000000	0.000000
    0.000000	0.000000	0.000000	0.000000	0.000000	0.982121	0.000000	0.000000	0.000000	0.000000
    0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.982121	0.000000	0.000000	0.000000
    0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.982121	0.000000	0.000000
    0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.982121	0.000000
    0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.982121
    
    # Second iteration
    B (autovetores Octave)
       0.8125186  -0.1691631   0.1675565   0.3036075  -0.1988456   0.1318318  -0.1036725  -0.2877685   0.1459485   0.1384410
      -0.0739690  -0.3296476  -0.0329443   0.6577079  -0.0396555  -0.1959856  -0.1202157   0.5352523  -0.3141237   0.1122162
       0.0860094  -0.0239913   0.1101355   0.1209071  -0.1053128   0.0914022   0.0361878   0.1891868   0.1766665  -0.9367609
       0.1710587   0.2926066   0.2159082   0.0375735   0.5360150  -0.0165138  -0.0518435   0.4741602   0.5406532   0.1722959
      -0.2360260   0.0704926   0.1786761   0.3642173  -0.0273500  -0.6633435   0.2029216  -0.4181460   0.3421923  -0.0291834
       0.3219345  -0.1709058  -0.6090203  -0.2932113   0.1209225  -0.5504422  -0.2586601   0.0761662   0.0657677  -0.1250207
       0.2978594   0.2397491   0.1315881  -0.2290375  -0.3545198  -0.2913323   0.6339612   0.3568469  -0.1951680   0.0782984
       0.2180601   0.3777284   0.1375072   0.1244026   0.5351562  -0.1354133  -0.0336620  -0.2216937  -0.6187924  -0.1935786
       0.0427945  -0.1413315  -0.5286527   0.2671639   0.3311497   0.2788639   0.6465257  -0.1108081   0.0977846  -0.0091029
      -0.0152181   0.7226418  -0.4399604   0.3182367  -0.3512522   0.1076815  -0.2041696   0.0257117   0.0673136   0.0294391
    
    B (autovetores python)
        -0.1384259   -0.1459358   -0.2877570   -0.1036655   -0.1318742   0.1988636   -0.3035690   0.8125861   -0.1674470   -0.1690090
        -0.1122053   0.3140650   0.5352581   -0.1201945   0.1960611   0.0397105   -0.6576593   -0.0738250   0.0332257   -0.3297553
        0.9367781   -0.1766552   0.1891549   0.0361811   -0.0913873   0.1052953   -0.1209273   0.0860089   -0.1100735   -0.0239841
        -0.1722673   -0.5406829   0.4741142   -0.0518599   0.0165789   -0.5360392   -0.0377038   0.1709855   -0.2159188   0.2926103
        0.0291819   -0.3421753   -0.4182125   0.2029998   0.6633460   0.0273967   -0.3642544   -0.2359662   -0.1785075   0.0703491
        0.1250075   -0.0657565   0.0761354   -0.2586457   0.5504205   -0.1208912   0.2935396   0.3219849   0.6089145   -0.1707660
        -0.0782899   0.1951337   0.3568705   0.6340030   0.2912502   0.3544898   0.2290191   0.2977878   -0.1316752   0.2398369
        0.1935509   0.6188264   -0.2216695   -0.0336583   0.1354607   -0.5351428   -0.1245001   0.2179817   -0.1375005   0.3777189
        0.0091046   -0.0977679   -0.1108045   0.6464700   -0.2789075   -0.3311802   -0.2669541   0.0428545   0.5287858   -0.1413232
        -0.0294343   -0.0673134   0.0257076   -0.2041704   -0.1076399   0.3512333   -0.3181570   -0.0153531   0.4400222   0.7226519
    
    D (autovalores Octave)
       0.96460         0         0         0         0         0         0         0         0         0
         0   0.96515         0         0         0         0         0         0         0         0
         0         0   0.96589         0         0         0         0         0         0         0
         0         0         0   0.96719         0         0         0         0         0         0
         0         0         0         0   0.97089         0         0         0         0         0
         0         0         0         0         0   0.97380         0         0         0         0
         0         0         0         0         0         0   0.97977         0         0         0
         0         0         0         0         0         0         0   1.00104         0         0
         0         0         0         0         0         0         0         0   1.00995         0
         0         0         0         0         0         0         0         0         0   1.09349
    D (autovalores Python)
    1.09350   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
    0.00000   1.00995   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
    0.00000   0.00000   1.00104   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
    0.00000   0.00000   0.00000   0.97977   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
    0.00000   0.00000   0.00000   0.00000   0.97380   0.00000   0.00000   0.00000   0.00000   0.00000
    0.00000   0.00000   0.00000   0.00000   0.00000   0.97088   0.00000   0.00000   0.00000   0.00000
    0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.96719   0.00000   0.00000   0.00000
    0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.96460   0.00000   0.00000
    0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.96589   0.00000
    0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.96515
"""



        ### Individuo em COLUNAS, CONTA FUNCIONANDO
        Caux1 = (1-c1-cmu) * C
        Caux2 = np.outer(pc, pc) + ((1-hsig) * cc * (2-cc) * C)
        BDARGZ = np.matmul(np.matmul(B,D), muBestZ)
        Caux3 = c1*Caux2
        Caux4 = (cmu * BDARGZ)
        Caux5 = np.matmul(np.diag(weights), BDARGZ.transpose())
        Caux6 = np.matmul(Caux4, Caux5)
        C = Caux1 + Caux3 + Caux6


        ###########################################
        # Individuo em LINHAS, CONTA FUNCIONANDO
        Caux1 = (1-c1-cmu) * C
        Caux2 = np.outer(pc, pc) + ((1-hsig) * cc * (2-cc) * C)
        BDARGZ = np.matmul(np.matmul(B,D), muBestZ.T)
        Caux3 = c1*Caux2
        Caux4 = (cmu * BDARGZ)
        Caux5 = np.matmul(np.diag(weights), BDARGZ.transpose())
        Caux6 = np.matmul(Caux4, Caux5)
        CLinha = Caux1 + Caux3 + Caux6


        if strFunction[0] == "2":
            offsprings.bounding(nSize, function, offspringsSize, lowerBound, upperBound)
            functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize,
                                                      functionEvaluations, truss)
        else:
            offsprings.bounding(nSize, function, offspringsSize)
            functionEvaluations = offsprings.evaluate(offspringsSize, function, nSize, gSize, hSize,
                                                      functionEvaluations)
        if penaltyMethod == 1:
            offsprings.sumViolations(offspringsSize, gSize, hSize)
        else:
            offsprings.uniteConstraints(offspringsSize, gSize, hSize)
            avgObjFunc = offsprings.calculatePenaltyCoefficients(offspringsSize, constraintsSize, penaltyCoefficients,
                                                                 avgObjFunc)
            offsprings.calculateAllFitness(offspringsSize, constraintsSize, penaltyCoefficients, avgObjFunc)
        parents.elitismES(offsprings, parentsSize, offspringsSize, nSize, gSize, hSize, constraintsSize, globalSigma,
                          esType, generatedOffspring, penaltyMethod, strFunction, truss, lowerBound, upperBound)
        parents.printBest(nSize, parentsSize, penaltyMethod)

if __name__ == '__main__':
    ESCMA(210, 1, 1, 10, 10, 10, 2000, 0, 0, 0)