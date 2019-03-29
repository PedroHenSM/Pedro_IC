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

        xmean = np.dot(arx[0:mu], weights[:mu], transpose=True)  # eq 42 TODO VERIFICAR TRANPOSE
        y = np.subtract(xmean, xold)
        ps = (1-cs) * ps + (np.sqrt(cs*(2-cs)*mueff)) * (B *zmean)

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