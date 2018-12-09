#-------------------------------#
#--- Ajuste do Metamodelo GP ---#
#-------------------------------#
import numpy as np
import pandas as pd
import sys

from sklearn import model_selection
from sklearn.gaussian_process  import  GaussianProcessRegressor 
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C


import random

thetaGlobal = []

def ranked (pred, y):
    n = len(pred)
    aux = pd.DataFrame(np.arange(0, n), columns = ['posto'])
    
    pred = pd.concat([aux, pred], axis = 1)
    pred = pred.sort_values(pred.columns.values[1]).reset_index(drop=True)
    
    y = pd.concat([aux, y], axis = 1)
    y = y.sort_values(y.columns.values[1]).reset_index(drop=True)
    
    count = 0
    
    for i in range(n):
        if pred[pred.columns.values[1]][i] != y[y.columns.values[1]][i]:
            if y.posto[i] == pred.posto[i]:
                count = count + 1
        
        else:
            count = count + 1
            
    return(count/n)

def hyperPar (x,y):
    aux = pd.DataFrame(np.zeros(shape=(20, 3)), columns = ['TAX', 'sig', 'l'])
    
    for i in range(20):
        x_train, x_test, y_train, y_test = model_selection.train_test_split(x, y)
        aux.iloc[i,1] = np.random.uniform(1.0, 10.0, 1)
        aux.iloc[i,2] = np.random.uniform(2.0, 10.0, 1)
                
        kernel=(C(aux.iloc[i,1], (1e-3, 1e3)))*RBF(aux.iloc[i,2],(1e-3, 1e3))
        gpr  =  GaussianProcessRegressor( kernel = kernel, alpha=1e-12,  normalize_y = False).fit(x_train,y_train)
        GaussianProcessRegressor()

        pred = gpr.predict(x_test)
        pred = pd.DataFrame(pred, columns = ['pred'])
        
        pred[pred < 0.0] = 0.0
        # print("y test --> {}".format(y_test))
        # print("pred --> {}".format(pred))
        y_test = pd.DataFrame(np.array(y_test).reshape(len(y_test),1)).reset_index(drop=True)
        # print("y test --> {}".format(y_test))
        #sys.exit("uai")
        # print("pred --> {}".format(pred))
        aux.iloc[i,0] = ranked(pred, y_test)

    # print(aux)
    aux1 = list(aux.TAX)
    aux2 = aux1.index(max(aux.TAX))
    return([aux.sig[aux2], aux.l[aux2]])

def pesos (x, y, p0, p1):
    for i in range(x.shape[0]):
        if y[i] == 0:
            x.iloc[i,:] = p0*x.iloc[i,:]
        else:
            x.iloc[i,:] = p1*x.iloc[i,:]
    
    return(x)

def surGPR_training (data, n, p, functionEvaluations):
    #argumentos:
    #data = dados em formato de lista
    #n = número de indivíduos
    #p = dimensão da estrutura (número de X) + 1 (csum)
    df = pd.DataFrame(np.array(data).reshape(n,p))
    x = np.log(df.iloc[:,0:(p-1)])
    csum = df.iloc[:,(p-1)]
    if len(csum[csum == 0]) == 0 or len(csum[csum == 0]) == len(csum):
        pass
    else:
        p0 = 1/ (len(csum[csum == 0])/n)
        p1 = 1/ (1-(len(csum[csum == 0])/n))
        x = pesos(x, csum, p0, p1)

    """
    if functionEvaluations == 200 or functionEvaluations == 7500:  # only calculates hyperparameters if on first iteration or in the middle of the process
        # print("CALCULOU HYPERPARAMETROS")
        theta = hyperPar(x, csum)
        thetaGlobal.clear()
        thetaGlobal.append(theta[0])
        thetaGlobal.append(theta[1])
        # list.append(theta[0])
        # list.append(theta[1])
        # theta = theta1
        # print(type(theta))
    """

    thetaGlobal = hyperPar(x,csum)  # TODO Descomentar para rodar atualizando parâmetro toda hora, e comentar if acima.
    # print("{}\t{}".format(thetaGlobal[0],thetaGlobal[1]))
    # print("thetaGlobal: {}".format(thetaGlobal))

    kernel = (C(thetaGlobal[0], (1e-3, 1e3))) * RBF(thetaGlobal[1], (1e-3, 1e3))
    # kernel = (C(thetaGlobal[0], (1e-4, 1e4))) * RBF(thetaGlobal[1], (1e-4, 1e4))
    gpr  =  GaussianProcessRegressor( kernel=kernel, alpha=1e-12, normalize_y = False).fit(x,csum)
    # aux = gpr.kernel_
    # print(aux)
    # thetaGlobal[0] = aux.get_params()['k1__constant_value']
    # thetaGlobal[1] = aux.get_params()['k2__length_scale']
    # print("{}\t{}".format(thetaGlobal[0], thetaGlobal[1]))
    crossV = model_selection.cross_val_score(gpr, x, csum, cv=5, scoring='explained_variance')
    # print ( "% 0.2f, % 0.2f "  %  ( crossV.mean(),  crossV.std()))
    #escrever o valor do CV e desvio sempre que rodar a função    
    return gpr, crossV.mean(), crossV.std()

#Como proceder:
    #(A) surGPR_training é a função de ajuste/treino do modelo. Ela retorna
    #   a estrutura do modelo ajustado. Ex.: model_GPR = surGPR_training (lista, 200, 11)
    
    #(B) após o treino/ajuste do modelo, o mesmo pode ser usado para as predições
    #   de csum com a extenção .predict() da função. Ex.: model_GPR.preditc(np.log(novos_individuos))
    
    #(C) no item (A) é necessário passa o csum, em (B) não, pois estmos predizendo-os
    
    #(D) a etapa (A) deve ser usada SEMPRE que precisar atualizar o modelo

#ATENÇÃO:
    #(1) Em um primeiro momento passar os dados transformados em escala logaritmica.
    #   Logo, ao usar o predict aplicar o np.log() nos dados antes, como segue no
    #   exemplo (B).
     