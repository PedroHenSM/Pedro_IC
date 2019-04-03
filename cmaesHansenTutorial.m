function xmin=purecmaes
    % CMA-ES: Evolution Strategy with Covariance Matrix Adaptation for
    % nonlinear function minimization.
    %
    % This code is an excerpt from cmaes.m and implements the key parts
    % of the algorithm. It is intendend to be used for READING and
    % UNDERSTANDING the basic flow and all details of the CMA *algorithm*.
    % Computational efficiency is sometimes disregarded.
    
    % -------------------- Initialization --------------------------------
    % User defined input parameters (need to be edited)
    % Next 2 lines fix the seed, just for personal tests
    rand('state', 1); % just for personal tests
    randn('state',1);  
    strfitnessfct = 'felli'; % name of objective/fitness function
    N = 10; % number of objective variables/problem dimension
    xmean = rand(N,1); % objective variables initial point
    sigma = 0.5; % coordinate wise standard deviation (step-size)
    stopfitness = 1e-10; % stop if fitness < stopfitness (minimization)
    stopeval = 1e3*N^2; % stop after stopeval number of function evaluations
    
    % Strategy parameter setting: Selection
    lambda = 4+floor(3*log(N)); % population size, offspring number
    mu = lambda/2; % lambda=12; mu=3; weights = ones(mu,1); would be (3_I,12)-ES
    % duvida1 = (1:mu);
    % duvida2 = (1:mu)';  % Python code: muList = [i+1 for i in range(mu)]
    % peso1 = log(mu+1/2);
    % peso2 = log(1:mu)';
    % pesos = peso1 - peso2;
    % PythonCode: weights = np.log(mu+1/2)-np.log(muList)
    weights = log(mu+1/2)-log(1:mu)'; % muXone recombination weights
    mu = floor(mu); % number of parents/points for recombination
    weights = weights/sum(weights); % normalize recombination weights array
    mueff=sum(weights)^2/sum(weights.^2); % variance-effective size of mu
    
    % Strategy parameter setting: Adaptation
    cc = (4+mueff/N) / (N+4 + 2*mueff/N); % time constant for cumulation for C
    cs = (mueff+2)/(N+mueff+5); % t-const for cumulation for sigma control
    c1 = 2 / ((N+1.3)^2+mueff); % learning rate for rank-one update of C
    cmu = 2 * (mueff-2+1/mueff) / ((N+2)^2+2*mueff/2); % and for rank-mu update
    damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma
    
    % Initialize dynamic (internal) strategy parameters and constants
    pc = zeros(N,1); ps = zeros(N,1); % evolution paths for C and sigma
    B = eye(N); % B defines the coordinate system
    D = eye(N); % diagonal matrix D defines the scaling
    C = B*D*(B*D)'; % covariance matrix
    eigeneval = 0; % B and D updated at counteval == 0
    chiN=N^0.5*(1-1/(4*N)+1/(21*N^2)); % expectation of
    % ||N(0,I)|| == norm(randn(N,1))
    
    % -------------------- Generation Loop --------------------------------
    
    counteval = 0; % the next 40 lines contain the 20 lines of interesting code
    while counteval < stopeval
    % randn: Normally distributed random numbers
        % Generate and evaluate lambda offspring
        for k=1:lambda,
          % pythonCode: provavelmente vou ter que criar um arz dentro de individual
          arz(:,k) = randn(N,1); % standard normally distributed vector
          arx(:,k) = xmean + sigma * (B*D * arz(:,k)); % add mutation % Eq. 40
          arfitness(k) = feval(strfitnessfct, arx(:,k)); % objective function call
          counteval = counteval+1;
        end

        % Sort by fitness and compute weighted mean into xmean
        [arfitness, arindex] = sort(arfitness); % minimization
        % arindex(1:mu) pega os melhores mu indices (após ordenados)
        % arx(:,arindex(1:mu)) seleciona os mu melhores individuos (ordenados)
        % arz(:,arindex(1:mu)) seleciona os mu melhores distribuicoes (dos individuos) NAO TEM NO MEU CODIGO
        % xmean se torna um vetor 10x1
        % pythonCode: xmean = np.dot(muBestIndividuals, weights)
        % pythonCode: zmean = np.dot(muBestIndividuals, weights)
        xmean = arx(:,arindex(1:mu))*weights; % recombination % Eq. 42
        zmean = arz(:,arindex(1:mu))*weights; % == Dˆ-1*B’*(xmean-xold)/sigma

        % Cumulation: Update evolution paths
        ps = (1-cs)*ps + (sqrt(cs*(2-cs)*mueff)) * (B * zmean); % Eq. 43
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4+2/(N+1);
        pc = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (B*D*zmean); % Eq. 45
        
        % Adapt covariance matrix C
        C = (1-c1-cmu) * C ... % regard old matrix % Eq. 47
        + c1 * (pc*pc' ... % plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ... % minor correction
        + cmu ... % plus rank mu update
        * (B*D*arz(:,arindex(1:mu))) ...
        * diag(weights) * (B*D*arz(:,arindex(1:mu)))';
        
        % Adapt step-size sigma
        sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); % Eq. 44

        % Update B and D from C
        if counteval - eigeneval > lambda/(c1+cmu)/N/10 % to achieve O(Nˆ2)
            eigeneval = counteval;
            C=triu(C)+triu(C,1)'; % enforce symmetry
            [B,D] = eig(C); % eigen decomposition, B==normalized eigenvectors
            D = diag(sqrt(diag(D))); % D contains standard deviations now
        end

        % Break, if fitness is good enough
        if arfitness(1) <= stopfitness
            break;
        end

        % Escape flat fitness, or better terminate?
        if arfitness(1) == arfitness(ceil(0.7*lambda))
            sigma = sigma * exp(0.2+cs/damps);
            disp('warning: flat fitness, consider reformulating the objective');
        end

        disp([num2str(counteval) ': ' num2str(arfitness(1))]);
        
    end % while, end generation loop

    % -------------------- Final Message ---------------------------------
    
    disp([num2str(counteval) ': '  num2str(arfitness(1))]);
    xmin = arx(:, arindex(1)); % Return best point of last generation.
    % Notice that xmean is expected to be even
    % better.

% ---------------------------------------------------------------
function f=felli(x)
    N = size(x,1); if N < 2 error('dimension must be greater one'); end
    f=1e6.^((0:N-1)/(N-1)) * x.^2; % condition number 1e6