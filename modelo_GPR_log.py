#-------------------------------#
#--- Ajuste do Metamodelo GP ---#
#-------------------------------#
import numpy as np
import pandas as pd

from sklearn import model_selection
from sklearn.gaussian_process  import  GaussianProcessRegressor 
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C


import random

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
                
        kernel=(C(aux.iloc[i,1], (1e-4, 1e3)))*RBF(aux.iloc[i,2],(1e-4, 1e3))
        gpr  =  GaussianProcessRegressor( kernel = kernel, n_restarts_optimizer=9, normalize_y = False).fit(x_train,y_train) 
        pred = gpr.predict(x_test)
        pred = pd.DataFrame(pred, columns = ['pred'])
        
        pred[pred < 0.0] = 0.0
        y_test = pd.DataFrame(y_test, columns = ['csum']).reset_index(drop=True)
        
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

def surGPR_training (data, n, p):
    #argumentos:
    #data = dados em formato de lista
    #n = número de indivíduos
    #p = dimensão da estrutura (número de X) + 1 (csum)
    df = pd.DataFrame(np.array(data).reshape(n,p))
    # df = data
    # print(type(df))
    x = np.log(df.iloc[:,0:(p-1)])
    csum = df.iloc[:,(p-1)]
    
    p0 = 1/ (len(csum[csum == 0])/n)
    p1 = 1/ (1-(len(csum[csum == 0])/n))
    x = pesos(x, csum, p0, p1)
        
    theta = hyperPar(x,csum)  # TODO: comentado, voltar?!
    kernel = (C(theta[0], (1e-4, 1e3))) * RBF(theta[1], (1e-4, 1e3))
    gpr  =  GaussianProcessRegressor( kernel=kernel, n_restarts_optimizer=9, normalize_y = False).fit(x,csum)
    
    crossV = model_selection.cross_val_score(gpr, x, csum, cv=5, scoring='explained_variance')
    # print ( "% 0.2f, % 0.2f "  %  ( crossV.mean(),  crossV.std()))
    #escrever o valor do CV e desvio sempre que rodar a função    
    return gpr, crossV.mean(), crossV.std()

print("Rodou")
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
     