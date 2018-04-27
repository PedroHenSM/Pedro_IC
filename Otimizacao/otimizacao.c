/*
Titulo: Computação Evolucionista Aplicada à Engenharia
Autor: Pedro Henrique Santos
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "functions.h"
#define TAM_POPULACAO 50
#define NUM_FILHOS_GERADOS 1
#define TAM_POPULACAO_FILHOS TAM_POPULACAO*NUM_FILHOS_GERADOS
#define MAX_CALC_OBJ 20000// 100000 original
#define TAM_X 10
#define TAM_X_MAXIMO 30
#define MEDIA 0
#define DESVIO_PADRAO 1
#define SEED 1   //BUG SEED 1,2,3 com N EM BETA 2 PROBLEMA PRA HOLDER SBX (gerabeta) seed 7 bug pasa rosenbrock,com ES
#define ELITISMO 5
#define GERA_ARQUIVO 1
#define PROB_CROSSOVER 100// 100 = 100% , 80 = 80%
#define TIPO_CROSSOVER 1 // 0 - Simples , 1 - SBX
#define TIPO_PROBLEMA 0 // 0 - Rosenbrock , 1 - Holder, 2 - Rastrigin
#define N_EM_BETA 2 // Normalmente entre 2 e 5. Quanto maior , maior a chance de gerar filhos mais 'pr�ximo' do pai
#define N_EM_MUTACAOREAL 4 // Para mutação Real (slide SBX)
#define DELTA_MAX 10 // Para mutação Real (slide SBX)
#define TIPO_FUNCAO 1 // Funcoes com restricao 14 e 13 = 3g,0h | 01 = 2g,0h
#define TIPO_ES 0 // Se 0: es+ Se 1:es, Se 2: es, modificado
#define SIGMA_GLOBAL 2 // Se sigma global == 1, sigma global, se nao, vetor de sigmas[]
#define NUM_RESTRICOES_G 2
#define NUM_RESTRICOES_H 0
#define EPSILON 0.0001 // Modulo da restricao de igualdade tem que ser menor que epsilon
#define CR 0.9 // Crossover Probability , DE [0,1] 0.9
#define F 0.5 // Differential Weight , DE [0,2] 0.7

typedef struct Individuo {
    float x[TAM_X_MAXIMO];
    float sigma[TAM_X_MAXIMO]; // Strategy Parameters ( Parametros de Controle )
    float funcaoObjetivo[TAM_X_MAXIMO];
    float g[NUM_RESTRICOES_G];
    float h[NUM_RESTRICOES_H];
    float v[NUM_RESTRICOES_G+NUM_RESTRICOES_H]; // Vetor de violacoes
    float violacao; // Soma das violacoes // Soma das violacoes
}Individuo;

void imprimeContaObj(int contaObj){
    printf("%i\t",contaObj);
}

void imprimeIndividuo(Individuo ind,int tamanhoX){ // Imprimindo em arquivo(1), ou cmd(0)
    int i;
    for(i=0;i<tamanhoX;i++){
        //printf("x: %f\n",ind.x[i]);
    }
    printf("Violacao: %f\tFO: %f\n",ind.violacao,ind.funcaoObjetivo[0]);
}

void imprimePopulacao(Individuo pop[],int tamanhoX,int tamanhoPop){
    int i=0;
    for(i=0;i<tamanhoPop;i++){
        imprimeIndividuo(pop[i],tamanhoX);

    }
    printf("\n");
}

void avaliaFuncaoRosenbrock(Individuo pop[],int *contaObj,int tamanhoX,int tamanhoPop){
    int i,j;
    for(i=0;i<tamanhoPop;i++){
        (*contaObj)++;
        pop[i].funcaoObjetivo[0]=0;
        for(j=0;j<tamanhoX-1;j++){
            float aux =pop[i].x[j] * pop[i].x[j]; // xi^2
            //float aux = pow(pop[i].x[j],2);
            //pop[i].funcaoObjetivo[0]+=100*(pow(pop[i].x[j+1] - aux,2)) + pow(pop[i].x[j]-1,2);
            pop[i].funcaoObjetivo[0]+=100*((pop[i].x[j+1]-aux) * (pop[i].x[j+1]-aux)) + ((1-pop[i].x[j])*(1-pop[i].x[j]));
            //printf("%i  %i\n",i,j);
        }
    }
}

void avaliaFuncaoRastrigin(Individuo pop[],int *contaObj,int tamanhoX,int tamanhoPop){
    int i,j;
    for(i=0;i<tamanhoPop;i++){
        (*contaObj)++;
        pop[i].funcaoObjetivo[0]=10*tamanhoX;
        for(j=0;j<tamanhoX;j++){
            pop[i].funcaoObjetivo[0]+= ((pop[i].x[j] * pop[i].x[j])) - (10*cos(2*M_PI*pop[i].x[j]));
        }
    }
}

void avaliaFuncaoHolder(Individuo pop[],int *contaObj){ // f(x,y) = TAM_X = 2
    int i,j=0;
    for(i=0;i<TAM_POPULACAO;i++){
        (*contaObj)++;
        float x,y,res;
        x = pop[i].x[j];
        y = pop[i].x[j+1];
        res = -fabs(sin(x)*cos(y)*exp(fabs(1-(sqrt(x*x+y*y)/M_PI))));
        pop[i].funcaoObjetivo[0]=res;
    }
}

void avaliaFuncao (Individuo populacao[],int *contaObj,int tamanhoX,int tamanhoPop){
    if(TIPO_PROBLEMA == 0){ // Rosenbrock
        avaliaFuncaoRosenbrock(populacao,contaObj,tamanhoX,tamanhoPop);
    }
    else if (TIPO_PROBLEMA == 1){ // Holder
        avaliaFuncaoHolder(populacao,contaObj);
    }
    else if (TIPO_PROBLEMA == 2){ // Rastrigin
        avaliaFuncaoRastrigin(populacao,contaObj,tamanhoX,tamanhoPop);
    }
}

int comparaFuncao(const void *a, const void *b){
    return (*(float*)a - *(float*)b);
}

int comparaFuncaoObjetivo0(const void *a, const void *b){ // Ordena crescente ROSENBROCK
    if((*(Individuo*)a).funcaoObjetivo[0] == (*(Individuo*)b).funcaoObjetivo[0])
        return 0;
    else{
        if((*(Individuo*)a).funcaoObjetivo[0] < (*(Individuo*)b).funcaoObjetivo[0]){
            return -1; // vem antes
        }
        else
            return 1; // vem depois
    }
}

int comparaFuncaoObjetivo1(const void *a, const void *b){ // Ordena os mais pr�ximos de -19 HOLDER TABLE FUNCTION
    if((*(Individuo*)a).funcaoObjetivo[0] == (*(Individuo*)b).funcaoObjetivo[0])
        return 0;
    else{
        if(fabs((*(Individuo*)a).funcaoObjetivo[0]+19) < fabs((*(Individuo*)b).funcaoObjetivo[0]+19)){
            return -1; // vem antes
        }
        else
            return 1; // vem depois
    }
}

int comparaFuncaoObjetivo2(const void *a, const void *b){ // Ordena os mais proximos de 0 RASTRIGIN FUNCTION
    if((*(Individuo*)a).funcaoObjetivo[0] == (*(Individuo*)b).funcaoObjetivo[0])
        return 0;
    else{
        if(fabs((*(Individuo*)a).funcaoObjetivo[0]) < fabs((*(Individuo*)b).funcaoObjetivo[0])){
            return -1; // vem antes
        }
        else
            return 1; // vem depois
    }
}

float randNormal (float media, float desvioPadrao){ // M�todo de Marsaglia (Polar Method)

  float U1, U2, W, mult;
  static float X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (media + desvioPadrao * (float) X2);
    }

  do
    {
      U1 = -1 + ((float) rand () / RAND_MAX) * 2;
      U2 = -1 + ((float) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (media + desvioPadrao * (float) X1);
}

void inicializaPopulacao(Individuo populacao[],int tamanhoX){
    int i,j;
    for(i=0;i<TAM_POPULACAO;i++){ // Inicializa Popula��o
        for(j=0;j<tamanhoX;j++){
            float eq = (rand()/(float)(RAND_MAX)); // gera valores entre [0,1]
            //float val =eq * (rand () % 21 - 10); // Valores decimais entre [-10,10]
            float val =eq * (rand () % 11); // gera valores deciamis entre [0,10]
            populacao[i].x[j] = val;
        }
    }
}

void selecao(Individuo populacao[],Individuo filhos[],int tamanhoX){ /// Torneio NOTE: NECESSARIO GERAR MAIS FILHOS AQUI?!
    int i,j; // Passa os ''melhores'' para filhos
    for(i=0;i<TAM_POPULACAO_FILHOS;i++){ // Sele��o
        int indice1,indice2;
        indice1= rand () % TAM_POPULACAO;
        indice2= rand () % TAM_POPULACAO;
        //printf("indice1: %i //// indice2: %i  \n ",indice1,indice2);
        if(populacao[indice1].funcaoObjetivo[0] <= populacao[indice2].funcaoObjetivo[0]){
                for(j=0;j<tamanhoX;j++){
                    filhos[i].x[j] = populacao[indice1].x[j];
                }
                filhos[i].funcaoObjetivo[0]=populacao[indice1].funcaoObjetivo[0];
        }
        else{
            for(j=0;j<tamanhoX;j++){
                filhos[i].x[j] = populacao[indice2].x[j];
            }
            filhos[i].funcaoObjetivo[0]=populacao[indice2].funcaoObjetivo[0];
        }
    }
}

void selecaoRestricao(Individuo populacao[],Individuo filhos[],int tamanhoX,int tamanhoPop,int tamanhoPopFilhos){ /// Torneio NOTE: NECESSARIO GERAR MAIS FILHOS AQUI?!
    int i,j; // Passa os ''melhores'' para filhos
    for(i=0;i<tamanhoPopFilhos;i++){ // Selecao
        int indice1,indice2;
        indice1 = rand () % tamanhoPop;
        indice2 = rand () % tamanhoPop;
        //printf("indice1: %i //// indice2: %i  \n ",indice1,indice2);

        if(populacao[indice1].violacao < populacao[indice2].violacao){
            for(j=0;j<tamanhoX;j++){
                filhos[i].x[j] = populacao[indice1].x[j];
            }
            filhos[i].funcaoObjetivo[0]=populacao[indice1].funcaoObjetivo[0];
        }
        else if(populacao[indice1].violacao > populacao[indice2].violacao){
            for(j=0;j<tamanhoX;j++){
                filhos[i].x[j] = populacao[indice2].x[j];
            }
            filhos[i].funcaoObjetivo[0]=populacao[indice2].funcaoObjetivo[0];
        }
        else if(populacao[indice1].funcaoObjetivo[0] < populacao[indice2].funcaoObjetivo[0]){
            for(j=0;j<tamanhoX;j++){
                filhos[i].x[j] = populacao[indice1].x[j];
            }
            filhos[i].funcaoObjetivo[0]=populacao[indice1].funcaoObjetivo[0];
        }
        else{
            for(j=0;j<tamanhoX;j++){
                filhos[i].x[j] = populacao[indice2].x[j];
            }
            filhos[i].funcaoObjetivo[0]=populacao[indice2].funcaoObjetivo[0];
        }
        /// O bloco abaixo faz a mesma coisa que os ifs de cima, porem eh mais lento.
        /*
        if(populacao[indice1].violacao < populacao[indice2].violacao){
            filhos[i] = populacao[indice1];
        }
        else if(populacao[indice1].violacao > populacao[indice2].violacao){
            filhos[i] = populacao[indice2];
        }
        else if(populacao[indice1].funcaoObjetivo[0] < populacao[indice2].funcaoObjetivo[0]){
            filhos[i] = populacao[indice1];
        }
        else{
            filhos[i] = populacao[indice2];
        }*/
    }
}

void crossover (Individuo filhos[],int tamanhoX,int tamanhoPopFilhos){ /// NOTE: NECESSARIO GERAR MAIS FILHOS AQUI?!
    int i,l;
    for(i=0;i<tamanhoPopFilhos;i+=2){ // Crossover
        float x1,x2;
        float v1,v2;
        float eq;
        for(l=0;l<tamanhoX;l++){
            //printf("Equacao: %f \n",eq);
            x1=filhos[i].x[l];
            x2=filhos[i+1].x[l];
            //printf("x1: %f \n x2: %f \n",x1,x2);
            if(x1 <=  x2){ // Ex 15-30
                eq =(rand()/(float)(RAND_MAX));
                v1 = (eq * (x2-x1))+ x1;
                eq =(rand()/(float)(RAND_MAX));
                v2 = (eq * (x2-x1))+ x1;
                //printf("V1: %f \n V2: %f \n",v1,v2);
                filhos[i].x[l]=v1; // Preenche o vetor completo em X
                filhos[i+1].x[l]=v2;
            }
            else{ // 30-15
                eq =(rand()/(float)(RAND_MAX));
                v1 = (eq * (x1-x2))+ x2;
                eq =(rand()/(float)(RAND_MAX));
                v2 = (eq * (x1-x2))+ x2;
                //printf("V1: %f \n V2: %f \n",v1,v2);
                filhos[i].x[l]=v1; // Preenche o vetor completo em X
                filhos[i+1].x[l]=v2;
            }
        }
    }
}

float geraBeta(){
    float u,beta;
    u = (rand()/(RAND_MAX+1.0));
    //printf("u = %f\n",u); // gera valores entre [0,1]
	if(u < 0.5){
		beta = pow(2.0*u,(1.0/(N_EM_BETA+1.0)));
	}
	else{
		beta = pow((0.5/(1.0-u)),(1.0/(N_EM_BETA+1.0)));
	}
	//printf("beta = %f\n",beta);
	return beta;
}

void crossoverSBX(Individuo populacao[],Individuo filhos[],int tamanhoX,int tamPop,int numFilhosGerados){ /// NOTE: NECESSARIO GERAR MAIS FILHOS AQUI?!
    int i,j,k;
    for(i=0;i<tamPop;i+=2){
        for(k=0;k<numFilhosGerados;k++){
            for(j=0;j<tamanhoX;j++){
                float c1,c2,p1,p2,media,beta;
                if(populacao[i].x[j] < populacao[i+1].x[j]){ // p2  > p1 (sempre)
                    p1 = populacao[i].x[j];
                    p2 = populacao[i+1].x[j];
                }
                else{
                    p1 = populacao[i+1].x[j];
                    p2 = populacao[i].x[j];
                }
                /*c1 = filhos[i].x[j];
                c2 = filhos[i+1].x[j];
                beta = abs((c1-c2)/(p1-p2)); // conferir
                */
                beta = geraBeta(); // Outra forma de gerar o beta
                media = (p1+p2)/2;
                // c1 e c2 'definitivos'
                c1 = media - 0.5*beta*(p2-p1);
                c2 = media + 0.5*beta*(p2-p1);
                filhos[i].x[j] = c1;
                filhos[i+1].x[j] = c2;
            }
        }
    }
}

int probabilidadeCrossover(int prob){
    if(rand() % 100 < prob ){ // 'Prob' chance de ocorrer
        return 1;
    }
    else{
        return 0;
    }
}

void mutacaoReal(Individuo populacao[],Individuo filhos[],int tamanhoX){ // Slide do SBX
    int i,j;
    for(i=0;i<TAM_POPULACAO;i++){
        for(j=0;j<tamanhoX;j++){
            if(rand() % 10 < 3){ // 30% chance
                float val,eq;
                val = (rand()/(float)(RAND_MAX)) * (rand() % 3 - 1); // Gera valor entre [-1,1]
                eq = 0.5*(N_EM_MUTACAOREAL+1)*pow(1-abs(val),N_EM_MUTACAOREAL);
                if(rand() % 2 == 1){ // 50% chance de ser positivo ou negativo
                    filhos[i].x[j]= populacao[i].x[j] + eq * DELTA_MAX;
                }
                else{
                    filhos[i].x[j]= populacao[i].x[j] - eq * DELTA_MAX;
                }
            }
        }
    }
}

void mutacao(Individuo filhos[],int tamanhoX){
    int i,j;
    for(i=0;i<TAM_POPULACAO_FILHOS;i++){// Muta��o
        for(j=0;j<tamanhoX;j++){
            filhos[i].x[j] += randNormal(MEDIA,DESVIO_PADRAO);//Preenche o vetor com os valores de x com  a distruibuicao normal
        }
    }
}

int individuoMelhorou(Individuo populacao,Individuo filho,int tipoProblema){ // Retorna 1 se filho melhor que pai
    if(TIPO_PROBLEMA == 0 || TIPO_PROBLEMA == 2) { // Rosenbrock ou Rastrigin
    //if(tipoProblema == 0){
        if(filho.funcaoObjetivo[0] < populacao.funcaoObjetivo[0]){ // Filho melhor que pai
            return 1;
        }
        else{
            return 0;
        }
    }
    else{
        if(fabs(filho.funcaoObjetivo[0]+19) < fabs(populacao.funcaoObjetivo[0]+19)){
            return 1;
        }
        else{
            return 0;
        }
    }
}

int geracaoMelhorou(Individuo populacao[],Individuo filhos[]){ // Testa se os os filhos são melhores que os pais
    int i,contaMelhores=0;
    for(i=0;i<TAM_POPULACAO;i++){
        if(filhos[i].funcaoObjetivo[0] < populacao[i].funcaoObjetivo[0]){ // Geração Melhorou
            contaMelhores++;
        }
    }
    //printf("\nMelhorou?: %f \n",contaMelhores/(float)TAM_POPULACAO);
    if(contaMelhores/(float)TAM_POPULACAO > 0.5){ // Mais da metade da populacao melhorou
        return 1;
    }
    else{
        return 0;
    }
}

void ordenaMelhores(Individuo populacao[],Individuo filhos[],int tipoProblema){
    // qsort(vetor,tamanhoVet,tamanho do tipo de dado,funcao comparacao)
    if(TIPO_PROBLEMA == 0){ // Rosenbrock
        qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaFuncaoObjetivo0); //Ordena populacao
        qsort(filhos,TAM_POPULACAO_FILHOS,sizeof(Individuo),comparaFuncaoObjetivo0);
    }
    else if(TIPO_PROBLEMA == 1){ // Holder
        qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaFuncaoObjetivo1);
        qsort(filhos,TAM_POPULACAO_FILHOS,sizeof(Individuo),comparaFuncaoObjetivo1);
    }
    else if(TIPO_PROBLEMA == 2){ // Rastrigin
        qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaFuncaoObjetivo2);
        qsort(filhos,TAM_POPULACAO_FILHOS,sizeof(Individuo),comparaFuncaoObjetivo2);
    }
}

void imprimeMelhores(Individuo populacao[],Individuo filhos[],int tamanhoX){ // Ser� "a mesma" se popula��o e filhos estiverem ordenados
    if(populacao[0].funcaoObjetivo[0] < filhos[0].funcaoObjetivo[0]){
        imprimeIndividuo(populacao[0],tamanhoX);
    }
    else{
        imprimeIndividuo(filhos[0],tamanhoX);
    }
}

void copiaIndividuos(Individuo filhos[],Individuo populacao[],int inicioDaCopia,int tamanhoX,int tamanhoPop){ // Copia de filhos para populacao
    // Inicio da copia: A partir de onde ser� copiado
    int i,j,k;
    for(i=0,k=inicioDaCopia;k<tamanhoPop;i++,k++){
        populacao[k].funcaoObjetivo[0]=filhos[i].funcaoObjetivo[0];
        for(j=0;j<tamanhoX;j++){
            populacao[k].x[j]=filhos[i].x[j];
        }
    }
}

void copiaIndividuosRestricao(Individuo filhos[],Individuo populacao[],int inicioDaCopia,int tamanhoX,int tamanhoPop){ // Copia de filhos para populacao
    int i,k;
    for(i=0,k=inicioDaCopia;k<tamanhoPop;i++,k++){
        populacao[k]=filhos[i];
    }
}

void elitismo(Individuo filhos[],Individuo populacao[],int tamanhoX,int tamanhoPop){
    copiaIndividuos(filhos,populacao,ELITISMO,tamanhoX,tamanhoPop);
}

void elitismoRestricao(Individuo filhos[],Individuo populacao[],int tamanhoX,int tamanhoPop){
    copiaIndividuosRestricao(filhos,populacao,ELITISMO,tamanhoX,tamanhoPop);
}

void geraDados(){
    int s,n;
    for(n=2;n<11;n++){ // TamanhoProblema
        for(s=1;s<31;s++){// seed
            //algoritmoGeneticoGeraDados(n,s);
        }
    }
}

void algoritmoGenetico(){
    int contaObj=0;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(SEED);
    inicializaPopulacao(populacao,TAM_X);
    avaliaFuncao(populacao,&contaObj,TAM_X,TAM_POPULACAO);
    while(contaObj < MAX_CALC_OBJ){
        selecao(populacao,filhos,TAM_X);
        if(probabilidadeCrossover(PROB_CROSSOVER) == 1){
            if(TIPO_CROSSOVER == 0){
                crossover(filhos,TAM_X,TAM_POPULACAO_FILHOS);
            }
            else if(TIPO_CROSSOVER == 1){
                crossoverSBX(populacao,filhos,TAM_X,TAM_POPULACAO,NUM_FILHOS_GERADOS);
            }
        }
        mutacao(filhos,TAM_X);
        avaliaFuncao(filhos,&contaObj,TAM_X,TAM_POPULACAO_FILHOS);
        ordenaMelhores(populacao,filhos,TIPO_PROBLEMA);
        elitismo(filhos,populacao,TAM_X,TAM_POPULACAO);
        // Pegar melhor individuo
        imprimeContaObj(contaObj);
        imprimeMelhores(populacao,filhos,TAM_X);
    }
}

int comparaViolacaoRestricaoMinimizacao(const void *a, const void *b){
    if((*(Individuo*)a).violacao < (*(Individuo*)b).violacao){
        return -1; // vem antes
    }
    else if((*(Individuo*)a).violacao > (*(Individuo*)b).violacao){
        return 1; // vem depois
    }
    else if((*(Individuo*)a).funcaoObjetivo[0] < (*(Individuo*)b).funcaoObjetivo[0]){
        return -1;
    }
    else if ((*(Individuo*)a).funcaoObjetivo[0] > (*(Individuo*)b).funcaoObjetivo[0]){
        return 1;
    }
    else{ // funcao objetivo igual
        return 0;
    }
}

void inicializaPopulacaoRestricao(Individuo populacao[],int tamanhoX,int tamanhoPop,int tipoFuncao){ ///TODO TROCAR IF POR SWITCH
    int i,j;
    for(i=0;i<tamanhoPop;i++){
        for(j=0;j<tamanhoX;j++){
            float eq = (rand()/(float)(RAND_MAX)); // gera valores entre [0,1]
            float val;
            if (tipoFuncao == 1){
                val = eq * (rand() % 11); // Valores decimais entre [0,10]
            }
            else if(tipoFuncao == 2){
                val = eq * (rand() % 11 - 5); // Gera valores decimais entre [-5,5]
            }
            else if(tipoFuncao == 3 || tipoFuncao == 12 || tipoFuncao == 14 || tipoFuncao == 15){
                val = eq * (rand() % 2001 - 1000); // Valores decimais entre [-1000,1000]
            }
            else if(tipoFuncao == 4 || tipoFuncao == 18){
                val = eq * (rand() % 101 - 50); // Valores decimais entre [-50,50]
            }
            else if(tipoFuncao == 5 || tipoFuncao == 6){
                val = eq * (rand() % 1201 - 600); // Valores decimais entre [-600,600]
            }
            else if(tipoFuncao == 7 || tipoFuncao == 8){
                val = eq * (rand() % 281 - 140); // Valores decimais entre [-140,140]
            }
            else if(tipoFuncao == 9 || tipoFuncao == 10 || tipoFuncao == 13){
                val = eq * (rand() % 1001 - 500); // valores decimais entre [-500,500]
            }
            else if(tipoFuncao == 11){
                val = eq * (rand() % 201 - 100); // valores decimais entre [-100,100]
            }
            else if(tipoFuncao == 16 || tipoFuncao == 17){
                val = eq * (rand() % 21 - 10); // valores decimais entre [-10,10]
            }
            else{
                printf("Funcao nao encontrada\n");
                exit(0);
            }
            populacao[i].x[j] = val;
        }
    }
}

void ordenaMelhoresRestricao(Individuo populacao[],Individuo filhos[],int tamanhoPop,int tamanhoPopFilhos){ // Ordena com base nas violacoes
    qsort(populacao,tamanhoPop,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
    qsort(filhos,tamanhoPopFilhos,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
}

void C01 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C02 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C03 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C03 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C04 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C05 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C06 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C07 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C08 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C09 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C10 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C11 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C12 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C13 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C14 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C15 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C16 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C17 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void C18 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh);

void inicializaRestricoes(int tipoFuncao,int *numRestricoesG,int *numRestricoesH){ ///Inicializa tipos de restricoes
    switch(tipoFuncao){
        case 1:
            printf("entrou caso 1 inicializa restricao\n");
            (*numRestricoesG) = 2;
            (*numRestricoesH) = 0;
            break;
        case 2:
            (*numRestricoesG) = 2;
            (*numRestricoesH) = 1;
            break;
        case 3:
            (*numRestricoesG) = 0;
            (*numRestricoesH) = 1;
            break;
        case 4:
            (*numRestricoesG) = 0;
            (*numRestricoesH) = 4;
            break;
        case 5:
            (*numRestricoesG) = 0;
            (*numRestricoesH) = 2;
            break;
        case 6:
            (*numRestricoesG) = 0;
            (*numRestricoesH) = 2;
            break;
        case 7:
            (*numRestricoesG) = 1;
            (*numRestricoesH) = 0;
            break;
        case 8:
            (*numRestricoesG) = 1;
            (*numRestricoesH) = 0;
            break;
        case 9:
            (*numRestricoesG) = 0;
            (*numRestricoesH) = 1;
            break;
        case 10:
            (*numRestricoesG) = 0;
            (*numRestricoesH) = 1;
            break;
        case 11:
            (*numRestricoesG) = 0;
            (*numRestricoesH) = 1;
            break;
        case 12:
            (*numRestricoesG) = 1;
            (*numRestricoesH) = 1;
            break;
        case 13:
            (*numRestricoesG) = 3;
            (*numRestricoesH) = 0;
            break;
        case 14:
            (*numRestricoesG) = 3;
            (*numRestricoesH) = 0;
            break;
        case 15:
            (*numRestricoesG) = 3;
            (*numRestricoesH) = 0;
            break;
        case 16:
            (*numRestricoesG) = 2;
            (*numRestricoesH) = 2;
            break;
        case 17:
            (*numRestricoesG) = 2;
            (*numRestricoesH) = 1;
            break;
        case 18:
            (*numRestricoesG) = 1;
            (*numRestricoesH) = 1;
            break;
        default:
            printf("Funcao nao encontrada\n");
            exit(1);
    }
}

void avaliaFuncaoRestricao(Individuo populacao[],int tipoFuncao,int tamanhoPop,int *numG,int *numH,int *fe){
    int i;
    printf("tipoFuncao em avaliaFuncao: %i\n",tipoFuncao);
    printf("tamPOp em avaliaFuncao: %i\n",tamanhoPop);
    printf("numG: %i\n",(*numG));
    printf("numH: %i\n",(*numH));
    printf("wtt");
    for(i=0;i<tamanhoPop;i++){ // parametros tam_x,tamFuncObj,tam_g,tam_h
        (*fe)++;
        //printf("entrou for\n");
        switch(tipoFuncao){
        case 1:
            //printf("entrou caso 1\n");
            C01(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 2:
            C02(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 3:
            C03(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 4:
            C04(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 5:
            C05(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 6:
            C06(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 7:
            C07(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 8:
            C08(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 9:
            C09(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 10:
            C10(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 11:
            C11(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 12:
            C12(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 13:
            C13(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 14:
            C14(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 15:
            C15(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 16:
            C16(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 17:
            C17(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        case 18:
            C18(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,(*numG),(*numH));
            break;
        default:
            printf("Funcao nao encontrada\n");
            exit(1);
        }
    }
        /*if(tipoFuncao == 1){
            C01(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,2,0);
        }
        else if(tipoFuncao == 2){
            C02(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,2,1);
        }
        else if(tipoFuncao == 3){
            C03(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,0,1);
        }
        else if(tipoFuncao == 4){
            C04(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,0,4);
        }
        else if(tipoFuncao == 10){
            C10(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,0,1);
        }
        else if(tipoFuncao == 13){
            C13(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,3,0);
        }
        else if(tipoFuncao == 14){
            C14(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,3,0);
        }
    }*/
}

float maiorValorArray(float g[],int tam){ // Pega maior valor do array, se for maior que 0
    float maior = g[0];
    int i;
    for(i=1;i<tam;i++){
        if(g[i] > maior){
            maior = g[i];
        }
    }
    if(maior < 0){ // Não tem violação, retorna 0
        maior = 0;
    }
    return maior;
}

float somaValoresArray(float v[],int tam){
    int i;
    float soma=0;
    for(i=0;i<tam;i++){
        if(v[i] > 0){ // Soma apenas os valores maiores que 0, violacao negativa nao eh somada
            soma+= v[i];
        }
    }
    return soma;
}

void corrigeLimitesX(Individuo populacao[],int tamanhoX,int tipoFuncao,int tamanhoPop){ // Arruma limites para os valores x
    float xmin,xmax;
    int i,j;
    if(tipoFuncao == 1){
        xmin = 0;
        xmax = 10;
    }
    else if(tipoFuncao == 2){
        xmin = -5.12;
        xmax = 5.12;
    }
    else if(tipoFuncao == 4){
        xmin = -50;
        xmax = 50;
    }
    else if(tipoFuncao == 13 || tipoFuncao == 10){
        xmin = -500;
        xmax = 500;
    }
    else if(tipoFuncao == 3 || tipoFuncao == 14){
        xmin = -1000;
        xmax = 1000;
    }
    for(i=0;i<tamanhoPop;i++){
        for(j=0;j<tamanhoX;j++){
            if(populacao[i].x[j] > xmax){
                populacao[i].x[j] = xmax;
            }
            else if(populacao[i].x[j] < xmin){
                populacao[i].x[j] = xmin;
            }
        }
    }
}

void consertaIgualdades(Individuo populacao[],int tamanhoPop){ // Faz o calculo das restricoes com igualdades e transforma em desigualdade
    if(NUM_RESTRICOES_H != 0){ // Possui igualdade.
        int i,indH;
        for(i=0;i<tamanhoPop;i++){
            for(indH=0;indH<NUM_RESTRICOES_H;indH++){
                populacao[i].h[indH] = fabs(populacao[i].h[indH]) - 0.0001;
            }
        }
    }
}

void somaViolacoes(Individuo populacao[],int tamanhoPop,int *numG,int *numH){ // Coloca as violacoes g e h no vetor 'v' de violacoes
    int i,j;
    int indG,indH;
    for(i=0;i<tamanhoPop;i++){
        indG=0,indH=0; // Indice do vetor g e h
        for(j=0;j<(*numG)+(*numH);j++){ // Percorrer vetor 'v' de violacao
            if(j < numG){
                populacao[i].v[j] = populacao[i].g[indG];
                indG++;
            }
            else{
                populacao[i].v[j] = fabs(populacao[i].h[indH]) - EPSILON; // Transforma h(x)=0 em g(x) <= 0
                indH++;
            }
        }
        populacao[i].violacao = somaValoresArray(populacao[i].v,(*numG)+(*numH));
    }
}

void imprimeInformacoesIndividuo(Individuo ind){
    int j;
    for(j=0;j<TAM_X;j++){
        printf("Os valores de x[] : %f\n",ind.x[j]);
    }
    printf("Os valores de funcaoObj: %e\n",ind.funcaoObjetivo[0]);
    for(j=0;j<NUM_RESTRICOES_G;j++){
        printf("Os valores de g[]: %f\n",ind.g[j]);
    }
    for(j=0;j<NUM_RESTRICOES_H;j++){
        printf("Os valores de h[]: %f\n",ind.h[j]);
    }
    for(j=0;j<NUM_RESTRICOES_G+NUM_RESTRICOES_H;j++){
        printf("O valor de v[]: %f\n",ind.v[j]);
    }
    printf("O valor de violacao: %f\n",ind.violacao);
    printf("--------------------ACABOU A IMPRESSAO---------------------\n");
}

void imprimeVioleFO(Individuo ind){
    printf("V\t%f\tFO\t%f\n",ind.violacao,ind.funcaoObjetivo[0]);
}

void imprimeInformacoesPopulacao(Individuo populacao[],int tamanhoPop){
    int i=0;
    for(i=0;i<tamanhoPop;i++){
        imprimeInformacoesIndividuo(populacao[i]);
    }
}

Individuo melhorIndividuoRestricao(Individuo populacao[],int tamanhoPop){
    int i;
    Individuo melhor = populacao[0];
    for(i=1;i<tamanhoPop;i++){
        if(populacao[i].violacao < melhor.violacao){
            melhor = populacao[i];
        }
        else if(populacao[i].violacao == melhor.violacao){
            if(populacao[i].funcaoObjetivo[0] < melhor.funcaoObjetivo[0]){
                melhor = populacao[i];
            }
        }
    }
    return melhor;
}

int selecionaPopulacaoDE(int solucao, int ch[], int tamanhoPop){
    int indice;
    int igual;
    int a;
    while(1){
        igual = 0;
        indice = rand() % tamanhoPop;
        if(indice != solucao){
            for(a=0;a<3;a++){
                if(indice == ch[a]){
                    igual = 1;
                    break;
                }
            }
            if(igual == 0){
                break;
            }
        }
    }
    return indice;
}

void selecaoDE(Individuo populacao[],Individuo filhos[],int tamanhoPop){
    int i; // No DE implementado, cada pai gera apenas um filho.
    for(i=0;i<tamanhoPop;i++){
        if(filhos[i].violacao < populacao[i].violacao){ // Filho melhor que pai
            populacao[i] = filhos[i];
        }
        else if(filhos[i].violacao == populacao[i].violacao){
            if(filhos[i].funcaoObjetivo[0] < populacao[i].funcaoObjetivo[0]){ // Filho melhor que pai
                populacao[i] = filhos[i];
            }
        }
    }
}

void inicalizaEstrategiaEvolutiva(Individuo populacao[], Individuo filhos[],int tamanhoX,int sigmaGlobal){
    int i,j;
    for(i=0;i<TAM_POPULACAO;i++){// Para Estratégia Evolutiva (ES)
        if(sigmaGlobal == 1){
            populacao[i].sigma[0] = DESVIO_PADRAO;
        }
        else{
            for(j=0;j<tamanhoX;j++){
                populacao[i].sigma[j] = DESVIO_PADRAO;
            }
        }
    }
    for(i=0;i<TAM_POPULACAO_FILHOS;i++){
        if(sigmaGlobal == 1){
            filhos[i].sigma[0] = DESVIO_PADRAO;
        }
        else{
            for(j=0;j<tamanhoX;j++){
                filhos[i].sigma[j]= DESVIO_PADRAO;
            }
        }
    }
}

void autoAdaptacaoSigma(Individuo populacao[],Individuo filhos[],int tamanhoX,int sigmaGlobal){
    /*int i,j;
    float epsilon,tau; // global step-size
    for(i=0;i<TAM_POPULACAO;i++){
        tau=1/sqrt(tamanhoX);
        epsilon = tau*randNormal(0,1);
        //printf("EPSILON: %f \n",epsilon);
        for(j=0;j<tamanhoX;j++){
            filhos[i].sigma[j] = populacao[i].sigma[j] * exp(epsilon) * exp(epsilon); //NOTE: Adicionar ao quadrado?
            filhos[i].x[j] = populacao[i].x[j] + filhos[i].sigma[j] * randNormal(0,1);
        }
    }*/
    int i,j,k,l=0;
    float epsilon,tau;
    for(i=0;i<TAM_POPULACAO;i++){
        for(j=0;j<NUM_FILHOS_GERADOS;j++){ // Cada pai deve gerar X filhos
            tau = 1/sqrt(tamanhoX);
            epsilon = tau*randNormal(0,1);
            if(sigmaGlobal == 1){
                filhos[l].sigma[0] = populacao[i].sigma[0] * exp(epsilon);
            }
            for(k=0;k<tamanhoX;k++){
                if(sigmaGlobal == 1){
                    filhos[l].x[k] = populacao[i].x[k] + filhos[l].sigma[0] * randNormal(0,1);
                }
                else{
                    filhos[l].sigma[k] = populacao[i].sigma[k] * exp(epsilon) * exp(epsilon); // NOTE: O quadrado eh necessario?
                    filhos[l].x[k] = populacao[i].x[k] + filhos[l].sigma[k] * randNormal(0,1);

                }
            }
            l++;
        }
    }
}

void selecionaMelhoresRestricao(Individuo populacao[],Individuo filhos[],int tipoES){
    if(tipoES == 0){ // Es + Junta pais e filhos, e pega os melhores.
        int i,k=0;
        Individuo aux[TAM_POPULACAO+TAM_POPULACAO_FILHOS];
        for(i=0;i<TAM_POPULACAO+TAM_POPULACAO_FILHOS;i++){
            if(i<TAM_POPULACAO){ // Copia populacao para aux
                aux[i]=populacao[i];
            }
            else{ // Copia filhos para aux
                aux[i]=filhos[k];
                k++;
            }
        }
        qsort(aux,(TAM_POPULACAO+TAM_POPULACAO_FILHOS),sizeof(Individuo),comparaViolacaoRestricaoMinimizacao); // Ordena auxiliar
        for(i=0;i<TAM_POPULACAO;i++){
            populacao[i]=aux[i]; // Copia os melhores de aux para populacao
        }
    }
    else if (tipoES == 1){ // Es , O período de vida de cada indivíduo está restrito a uma geração. Pega o melhor filho de cada pai
        int i,j=0;
        for(i=0;i<TAM_POPULACAO;i++){
            Individuo melhor = filhos[j];
            while(j<NUM_FILHOS_GERADOS*(i+1)){ // Percorre os X filhos de cada pai
                if(filhos[j].violacao < melhor.violacao){ // Seleciona o melhor filho
                    melhor = filhos[j];
                }
                else if(filhos[j].violacao == melhor.violacao){
                    if(filhos[j].funcaoObjetivo[0] < melhor.funcaoObjetivo[0]){ // Violacao eh igual, confere funcObj
                        melhor = filhos[j];
                    }
                    // Se nao, melhor ja eh melhor
                }
                j++;
            }
            populacao[i] = melhor;
        }
        qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
    }
    else if(tipoES == 2){ // Es , Adaptado (Seleciona melhor filho de cada pai, e se for melhor que o pai, substitui)
        int i,j=0;
        for(i=0;i<TAM_POPULACAO;i++){
            Individuo melhor = filhos[j];
            while(j<NUM_FILHOS_GERADOS*(i+1)){ // Pega o melhor filho de cada pai
                if(filhos[j].funcaoObjetivo[0] < melhor.funcaoObjetivo[0]){
                    melhor = filhos[j];
                }
                j++;
            }
            if(melhor.violacao < populacao[i].violacao){ // Filho melhor que pai
                populacao[i] = melhor; // Filho eh colocado no lugar do pai
            }
            if(melhor.violacao == populacao[i].violacao){
                if(melhor.funcaoObjetivo[0] < populacao[i].funcaoObjetivo[0]){ // Filho melhor que pai
                    populacao[i] = melhor; // Filho eh colocado no lugar do pai
                }
            }
        }
        qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
    }
    else{
        printf("Es nao reconhecido\n");
        exit(0);
    }
}

void AG(int tipoFuncao,int seed,int tamPopulacao,int tamX,int numFilhosGerados,int maxFE,int tipoES,int sigmaGlobal,int probCrossover){
    int fe=0; // Function Evaluation
    int tamPopulacaoFilhos = tamPopulacao*numFilhosGerados;
    int numG,numH;
    Individuo populacao[tamPopulacao];
    Individuo filhos[tamPopulacaoFilhos];
    srand(seed);
    inicializaRestricoes(tipoFuncao,&numG,&numH);
    inicializaPopulacaoRestricao(populacao,tamX,tamPopulacao,tipoFuncao);
    avaliaFuncaoRestricao(populacao,tipoFuncao,tamPopulacao,&fe,numG,numH);
    somaViolacoes(populacao,tamPopulacao,numG,numH);
    while(fe < maxFE){
        selecaoRestricao(populacao,filhos,tamX,tamPopulacao,tamPopulacaoFilhos);
        if(probabilidadeCrossover(PROB_CROSSOVER) == 1){
            if(TIPO_CROSSOVER == 0){
                crossover(filhos,tamX,tamPopulacaoFilhos);
            }
            else if(TIPO_CROSSOVER == 1){
                crossoverSBX(populacao,filhos,tamX,tamPopulacao,numFilhosGerados);
            }
        }
        mutacao(filhos,tamX);
        corrigeLimitesX(filhos,tamX,tipoFuncao,tamPopulacaoFilhos);
        avaliaFuncaoRestricao(populacao,tipoFuncao,tamPopulacao,&fe,numG,numH);
        somaViolacoes(filhos,tamPopulacaoFilhos,numG,numH);
        ordenaMelhoresRestricao(populacao,filhos,tamPopulacao,tamPopulacaoFilhos);
        elitismoRestricao(filhos,populacao,tamX,tamPopulacao);
        // Pegar melhor individuo
        ordenaMelhoresRestricao(populacao,filhos,tamPopulacao,tamPopulacaoFilhos);
        //printf("FE: %i\n",fe);
        //imprimeInformacoesIndividuo(populacao[0]);
    }
    //(*melhorAG) = melhorIndividuoRestricao(populacao,TAM_POPULACAO);
}

void DE(int tipoFuncao,int seed,int tamPopulacao,int tamX,int numFilhosGerados,int maxFE,int tipoES,int sigmaGlobal,int probCrossover){
    int fe=0;
    int a,i,j;
    int tamPopulacaoFilhos = tamPopulacao*numFilhosGerados;
    int numG;
    int numH;
    printf("pop filhos: %i\n",tamPopulacaoFilhos);
    printf("numG: %i\nnumH:%i\n",numG,numH);
    Individuo populacao[tamPopulacao];
    Individuo filhos[tamPopulacaoFilhos];
    srand(seed);
    inicializaRestricoes(tipoFuncao,&numG,&numH);
    printf("numG: %i\nnumH:%i\n",numG,numH);
    inicializaPopulacaoRestricao(populacao,tamX,tamPopulacao,tipoFuncao);
        printf("opa\n");

    avaliaFuncaoRestricao(populacao,tipoFuncao,tamPopulacao,&numG,&numH,&fe);
        printf("opa\n");

    somaViolacoes(populacao,tamPopulacao,&numG,&numH);
        printf("opa\n");
    exit(1);
    while(fe < maxFE){
    //while(iteracoes< MAX_CALC_OBJ){
        selecaoRestricao(populacao,filhos,tamX,tamPopulacao,tamPopulacaoFilhos);
        int ch[3] = {-1,-1,-1}; // Vetor de indices, Poderia ser de tamanho 3?
        for(i=0;i<tamPopulacao;i++){
            for(a=0;a<3;++a){ // Preenche vetor de indices
                ch[a] = selecionaPopulacaoDE(i,ch,tamPopulacao);
            }
            int R = rand() % tamX; // Indice aleatorio, baseado na dimensao do problema
            for(j=0;j<tamX;j++){ // Altera valor de x
                float Ri = (rand() / (float)(RAND_MAX));
                if(Ri < CR || j == R){
                    filhos[i].x[j] = populacao[ch[0]].x[j] + F * (populacao[ch[1]].x[j] - populacao[ch[2]].x[j]);
                }
                else{
                    filhos[i].x[j] = populacao[i].x[j];
                }
            }
        }
        corrigeLimitesX(filhos,tamX,tipoFuncao,tamPopulacaoFilhos);
        avaliaFuncaoRestricao(populacao,tipoFuncao,tamPopulacao,&numG,&numH,&fe);
        somaViolacoes(filhos,tamPopulacaoFilhos,&numG,&numH);
        selecaoDE(populacao,filhos,tamPopulacaoFilhos);
        //imprimeIndividuo(melhorIndividuoDE(populacao,TAM_POPULACAO),TAM_X);
        //printf("FE: %i\n",fe);
        //printf("Iteracao: %i\n",iteracoes);
        imprimeInformacoesIndividuo(melhorIndividuoRestricao(populacao,tamPopulacao));
    }
    //(*melhorDE) = melhorIndividuoRestricao(populacao,tamPopulacao);
}

void ES(int tipoFuncao,int seed,int tamPopulacao,int tamX,int numFilhosGerados,int maxFE,int tipoES,int sigmaGlobal,int probCrossover){
    int i,j=0;
    int fe=0;
    int tamPopulacaoFilhos = tamPopulacao*numFilhosGerados;
    int numG,numH;
    Individuo populacao[tamPopulacao];
    Individuo filhos[tamPopulacaoFilhos];
    srand(seed);
    inicializaRestricoes(tipoFuncao,&numG,&numH);
    inicializaPopulacaoRestricao(populacao,tamX,tamPopulacao,tipoFuncao);
    corrigeLimitesX(populacao,tamX,tipoFuncao,tamPopulacao); // No caso de C01 x[-10,10]
    // Avalia funcao
    avaliaFuncaoRestricao(populacao,tipoFuncao,tamPopulacao,&fe,&numG,&numH);
    somaViolacoes(populacao,tamPopulacao,numG,numH); // Preenche vetor 'v' de violacoes e seta variavel violacao, que eh a soma das violacoes
    inicalizaEstrategiaEvolutiva(populacao,filhos,tamX,sigmaGlobal);
    /// TODO: MELHORAR FOR ABAIXO
    for(i=0;i<tamPopulacaoFilhos;i++,j++){ // Copia populacao para filhos
        filhos[i]=populacao[j];
        if(j >=tamPopulacao){ // Acessaria lixo
            j=0; // reinicializa j, e continua copiando.
        }
    }
    //imprimeInformacoesPopulacao(filhos,TAM_POPULACAO_FILHOS);
    while(fe < MAX_CALC_OBJ){
    //while(numIteracoes <MAX_CALC_OBJ){
        autoAdaptacaoSigma(populacao,filhos,tamX,sigmaGlobal);
        // Avalia funcao
        corrigeLimitesX(filhos,tamX,tipoFuncao,tamPopulacaoFilhos);
        avaliaFuncaoRestricao(populacao,tipoFuncao,tamPopulacao,&fe,numG,numH);
        somaViolacoes(filhos,tamPopulacaoFilhos,numG,numH);
        ordenaMelhoresRestricao(populacao,filhos,tamPopulacao,tamPopulacaoFilhos);
        //ordenaMelhores(populacao,filhos);// NOTE e necessario?
        selecionaMelhoresRestricao(populacao,filhos,tipoES);
        //selecionaMelhores(populacao,filhos,TIPO_ES);
        // Pegar melhor individuo
        //imprimeContaObj(contaObj);
        //imprimeIndividuo(populacao[0],TAM_X);
        //printf("FE: %i\n",fe);
        //printf("iteracao: %i\n",numIteracoes);
        //printf("O melhor individuo: %f\n",populacao[0].funcaoObjetivo[0]);
        imprimeInformacoesIndividuo(populacao[0]); //Melhor Individuo
        //imprimeMelhores(populacao,filhos,TAM_X);
    }
    //(*melhorES) = populacao[0];
    //qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaFuncaoObjetivo0);
    //imprimeInformacoesIndividuo(populacao[0]);
}

int main(){
    DE(1,1,50,10,1,20000,0,0,0);
    return 0;
}


/*teste para ver se copia direta funciona. Funciona!
    Individuo a,b;
    a.funcaoObjetivo[0]=99;
    a.x[0]=1;
    a.x[1]=2;
    a.v[0]=5;
    a.v[1]=5;
    a.violacao=10;
    b.funcaoObjetivo[0]=33;
    b.violacao=5;
    printf("a\n");
    imprimeIndividuo(a,TAM_X);
    printf("b\n");
    imprimeIndividuo(b,TAM_X);
    b=a;
    printf("\n");
    imprimeIndividuo(a,TAM_X);
    printf("b\n");
    imprimeIndividuo(b,TAM_X);
*/
