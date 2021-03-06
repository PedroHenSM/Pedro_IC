/*
Titulo: Computação Evolucionista Aplicada à Engenharia
Autor: Pedro Henrique Santos
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include "functions.h"
#include "apm.c"
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
#define TIPO_ES 0 // Se 0: es+ Se 1:es, Se 2: es, modificado
#define SIGMA_GLOBAL 2 // Se sigma global == 1, sigma global, se nao, vetor de sigmas[]
#define NUM_RESTRICOES 10
#define EPSILON 0.0001 // Modulo da restricao de igualdade tem que ser menor que epsilon
#define CR 0.9 // Crossover Probability , DE [0,1] 0.9
#define F 0.5 // Differential Weight , DE [0,2] 0.7


typedef struct Individuo {
    float x[TAM_X_MAXIMO];
    float sigma[TAM_X_MAXIMO]; // Strategy Parameters ( Parametros de Controle )
    float funcaoObjetivo[TAM_X_MAXIMO];
    float g[NUM_RESTRICOES];
    float h[NUM_RESTRICOES];
    //float penaltyCoefficients[NUM_RESTRICOES]; // For APM
    float fitness; // Fitness for APM
    //float avgObjFunc; // For APM
    float v[2*NUM_RESTRICOES]; // Vetor de violacoes
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

void selecaoRestricaoAPM(Individuo populacao[],Individuo filhos[],int tamanhoX,int tamanhoPop,int tamanhoPopFilhos){
    int i,j;
    for(i=0;i<tamanhoPopFilhos;i++){ // Selecao
        int indice1,indice2;
        indice1 = rand () % tamanhoPop;
        indice2 = rand () % tamanhoPop;
        //printf("indice1: %i //// indice2: %i  \n ",indice1,indice2);

        if(populacao[indice1].fitness < populacao[indice2].fitness){
            for(j=0;j<tamanhoX;j++){
                filhos[i].x[j] = populacao[indice1].x[j];
            }
            filhos[i].funcaoObjetivo[0]=populacao[indice1].funcaoObjetivo[0];
        }
        else if(populacao[indice1].fitness >= populacao[indice2].fitness){
            for(j=0;j<tamanhoX;j++){
                filhos[i].x[j] = populacao[indice2].x[j];
            }
            filhos[i].funcaoObjetivo[0]=populacao[indice2].funcaoObjetivo[0];
        }
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

void mutacao(Individuo filhos[],int tamanhoX,int tamanhoPopFilhos){
    int i,j;
    for(i=0;i<tamanhoPopFilhos;i++){// Muta��o
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

int comparaFitnesseObjFunc(const void *a, const void *b){
    if((*(Individuo*)a).fitness < (*(Individuo*)b).fitness){
        return -1; // vem antes
    }
    else if((*(Individuo*)a).fitness > (*(Individuo*)b).fitness){
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
            switch(tipoFuncao){
                case 1:
                    val = eq * (rand() % 11); // Valores decimais entre [0,10]
                    break;
                case 2:
                    val = eq * (rand() % 11 - 5); // Gera valores decimais entre [-5,5]
                    break;
                case 3:
                case 12:
                case 14:
                case 15:
                    val = eq * (rand() % 2001 - 1000); // Valores decimais entre [-1000,1000]
                    break;
                case 4:
                case 18:
                    val = eq * (rand() % 101 - 50); // Valores decimais entre [-50,50]
                    break;
                case 5:
                case 6:
                    val = eq * (rand() % 1201 - 600); // Valores decimais entre [-600,600]
                    break;
                case 7:
                case 8:
                    val = eq * (rand() % 281 - 140); // Valores decimais entre [-140,140]
                    break;
                case 9:
                case 10:
                case 13:
                    val = eq * (rand() % 1001 - 500); // valores decimais entre [-500,500]
                    break;
                case 11:
                    val = eq * (rand() % 201 - 100); // valores decimais entre [-100,100]
                    break;
                case 16:
                case 17:
                    val = eq * (rand() % 21 - 10); // valores decimais entre [-10,10]
                    break;
                default:
                    printf("Tipo funcao nao encontrado\n");
                    exit(1);
                    break;
            }
            populacao[i].x[j] = val;
        }
    }
}

void ordenaMelhoresRestricao(Individuo populacao[],Individuo filhos[],int tamanhoPop,int tamanhoPopFilhos,int penaltyMethod){ // Ordena com base nas violacoes
    if(penaltyMethod == 1){ // Not APM
        qsort(populacao,tamanhoPop,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
        qsort(filhos,tamanhoPopFilhos,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
    }
    else if(penaltyMethod == 2){ // APM
        qsort(populacao,tamanhoPop,sizeof(Individuo),comparaFitnesseObjFunc);
        qsort(filhos,tamanhoPopFilhos,sizeof(Individuo),comparaFitnesseObjFunc);
    }
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

void avaliaFuncaoRestricao(Individuo populacao[],int tipoFuncao,int tamanhoPop,int numG,int numH,int *fe){
    int i;
    for(i=0;i<tamanhoPop;i++){ // parametros tam_x,tamFuncObj,tam_g,tam_h
        (*fe)++;
        switch(tipoFuncao){
        case 1:
            C01(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 2:
            C02(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 3:
            C03(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 4:
            C04(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 5:
            C05(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 6:
            C06(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 7:
            C07(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 8:
            C08(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 9:
            C09(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 10:
            C10(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 11:
            C11(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 12:
            C12(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 13:
            C13(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 14:
            C14(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 15:
            C15(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 16:
            C16(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 17:
            C17(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        case 18:
            C18(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,numG,numH);
            break;
        default:
            printf("Funcao nao encontrada\n");
            exit(1);
        }
    }
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
    switch(tipoFuncao){
        case 1:
            xmin = 0;
            xmax = 10;
            break;
        case 2:
            xmin = -5.12;
            xmax = 5.12;
            break;
        case 3:
        case 12:
        case 14:
        case 15:
            xmin = -1000;
            xmax = 1000;
            break;
        case 4:
        case 18:
            xmin = -50;
            xmax = 50;
            break;
        case 5:
        case 6:
            xmin = -600;
            xmax = 600;
            break;
        case 7:
        case 8:
            xmin = -140;
            xmax = 140;
            break;
        case 9:
        case 10:
        case 13:
            xmin = -500;
            xmax = 500;
            break;
        case 11:
            xmin = -100;
            xmax = 100;
            break;
        case 16:
        case 17:
            xmin = -10;
            xmax = 10;
            break;
        default:
            printf("Tipo funcao nao encontrado\n");
            exit(1);
            break;
    }
    /*
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
    }*/
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

void somaViolacoes(Individuo populacao[],int tamanhoPop,int numG,int numH){ // Coloca as violacoes g e h no vetor 'v' de violacoes e soma
    int i,j;
    int indG,indH;
    for(i=0;i<tamanhoPop;i++){
        indG=0,indH=0; // Indice do vetor g e h
        for(j=0;j<numG+numH;j++){ // Percorrer vetor 'v' de violacao
            if(j < numG){
                populacao[i].v[j] = populacao[i].g[indG];
                indG++;
            }
            else{
                populacao[i].v[j] = fabs(populacao[i].h[indH]) - EPSILON; // Transforma h(x)=0 em g(x) <= 0
                indH++;
            }
        }
        populacao[i].violacao = somaValoresArray(populacao[i].v,numG+numH);
    }
}

void uniteConstraints(Individuo populacao[],int tamanhoPop,int numG,int numH){ // Coloca as violacoes g e h no vetor 'v' de violacoes
    int i,j;
    int indG,indH;
    for(i=0;i<tamanhoPop;i++){
        indG=0,indH=0; // Indice do vetor g e h
        for(j=0;j<numG+numH;j++){ // Percorrer vetor 'v' de violacao
            if(j < numG){
                populacao[i].v[j] = populacao[i].g[indG];
                indG++;
            }
            else{
                populacao[i].v[j] = fabs(populacao[i].h[indH]) - EPSILON; // Transforma h(x)=0 em g(x) <= 0
                indH++;
            }
        }
    }
}

void imprimeInformacoesIndividuo(Individuo ind,int tamX,int numG,int numH){
    int j;
    for(j=0;j<tamX;j++){
        printf("Os valores de x[] : %f\n",ind.x[j]);
    }
    printf("Os valores de funcaoObj: %e\n",ind.funcaoObjetivo[0]);
    for(j=0;j<numG;j++){
        printf("Os valores de g[]: %f\n",ind.g[j]);
    }
    for(j=0;j<numH;j++){
        printf("Os valores de h[]: %f\n",ind.h[j]);
    }
    for(j=0;j<numG+numH;j++){
        printf("O valor de v[]: %f\n",ind.v[j]);
    }
    printf("O valor de violacao: %f\n",ind.violacao);
    printf("--------------------ACABOU A IMPRESSAO---------------------\n");
}

void imprimeVioleFO(Individuo ind){
    printf("\nV\t%e\tFO\t%e",ind.violacao,ind.funcaoObjetivo[0]);
}

void imprimeFiteFO(Individuo ind){
    if(ind.fitness == ind.funcaoObjetivo[0]){
        printf("\nF\t1\tFO\t%e",ind.funcaoObjetivo[0]);
    }
    else{
        printf("\nF\t1\tFO\t%e",ind.funcaoObjetivo[0]);
    }
    //printf("\nF\t%e\tFO\t%e",ind.fitness,ind.funcaoObjetivo[0]);
}

void imprimeInformacoesPopulacao(Individuo populacao[],int tamanhoPop,int tamX,int numG,int numH){
    int i=0;
    for(i=0;i<tamanhoPop;i++){
        imprimeInformacoesIndividuo(populacao[i],tamX,numG,numH);
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

Individuo melhorIndividuoRestricaoAPM(Individuo populacao[],int tamanhoPop){
    int i=0;
    int bestIdx=0;
    for(i=1;i<tamanhoPop;i++){
        if(populacao[i].fitness < populacao[bestIdx].fitness){
            bestIdx = i;
        }
        else if(populacao[i].fitness == populacao[bestIdx].fitness){ /// TODO: Verificar esse if
            if(populacao[i].funcaoObjetivo[0] < populacao[bestIdx].funcaoObjetivo[0]){
                bestIdx = i;
            }
        }
    }
    return populacao[bestIdx];
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

void selecaoDEAPM(Individuo populacao[],Individuo filhos[],int tamanhoPop){
    int i;
    for(i=0;i<tamanhoPop;i++){
        if(filhos[i].fitness < populacao[i].fitness){ // Filho melhor que pai
            populacao[i] = filhos[i];
        }
        else if(filhos[i].fitness == populacao[i].fitness){ /// TODO: Verificar esse if
            if(filhos[i].funcaoObjetivo[0] < populacao[i].funcaoObjetivo[0]){
                populacao[i] = filhos[i]; // Filho melhor que pai
            }
        }
    }
}

void inicalizaEstrategiaEvolutiva(Individuo populacao[], Individuo filhos[],int tamanhoX,int tamanhoPop,int tamanhoPopFilhos,int sigmaGlobal){
    int i,j;
    for(i=0;i<tamanhoPop;i++){// Para Estratégia Evolutiva (ES)
        if(sigmaGlobal == 1){
            populacao[i].sigma[0] = DESVIO_PADRAO;
        }
        else{
            for(j=0;j<tamanhoX;j++){
                populacao[i].sigma[j] = DESVIO_PADRAO;
            }
        }
    }
    for(i=0;i<tamanhoPopFilhos;i++){
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

void autoAdaptacaoSigma(Individuo populacao[],Individuo filhos[],int tamanhoX,int tamanhoPop,int numFilhosGerados,int sigmaGlobal){
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
    for(i=0;i<tamanhoPop;i++){
        for(j=0;j<numFilhosGerados;j++){ // Cada pai deve gerar X filhos
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

void selecionaMelhoresRestricao(Individuo populacao[],Individuo filhos[],int tipoES,int tamanhoPop,int tamanhoPopFilhos,int numFilhosGerados,int penaltyMethod){
    if(tipoES == 0){ // Es + Junta pais e filhos, e pega os melhores.
        int i,k=0;
        Individuo aux[tamanhoPop+tamanhoPopFilhos];
        for(i=0;i<tamanhoPop+tamanhoPopFilhos;i++){
            if(i<tamanhoPop){ // Copia populacao para aux
                aux[i]=populacao[i];
            }
            else{ // Copia filhos para aux
                aux[i]=filhos[k];
                k++;
            }
        }
        if (penaltyMethod == 1){
            qsort(aux,(tamanhoPop+tamanhoPopFilhos),sizeof(Individuo),comparaViolacaoRestricaoMinimizacao); // Ordena auxiliar
        }
        else if (penaltyMethod == 2){ // APM
            qsort(aux,(tamanhoPop+tamanhoPopFilhos),sizeof(Individuo),comparaFitnesseObjFunc); // Ordena auxiliar
        }
        for(i=0;i<tamanhoPop;i++){
            populacao[i]=aux[i]; // Copia os melhores de aux para populacao
        }
    }
    else if (tipoES == 1){ // Es , O período de vida de cada indivíduo está restrito a uma geração. Pega o melhor filho de cada pai
        int i,j=0;
        for(i=0;i<tamanhoPop;i++){
            Individuo melhor = filhos[j];
            while(j<numFilhosGerados*(i+1)){ // Percorre os X filhos de cada pai
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
        if (penaltyMethod == 1){
            qsort(populacao,tamanhoPop,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
        }
        else if(penaltyMethod == 2){ // APM
            qsort(populacao,tamanhoPop,sizeof(Individuo),comparaFitnesseObjFunc);
        }
    }
    else if(tipoES == 2){ // Es , Adaptado (Seleciona melhor filho de cada pai, e se for melhor que o pai, substitui)
        int i,j=0;
        for(i=0;i<tamanhoPop;i++){
            Individuo melhor = filhos[j];
            while(j<numFilhosGerados*(i+1)){ // Pega o melhor filho de cada pai
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
        if (penaltyMethod == 1){
            qsort(populacao,tamanhoPop,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
        }
        else if (penaltyMethod == 2){ // APM
            qsort(populacao,tamanhoPop,sizeof(Individuo),comparaFitnesseObjFunc);
        }
    }
    else{
        printf("Es nao reconhecido\n");
        exit(1);
    }
}

void calculatePenaltyCoefficients(Individuo populacao[],int populationSize,/*float* objectiveFunctionValue/,float* constraintViolationValues*/int numberOfConstraints,float* penaltyCoefficients,float* averageObjectiveFunctionValues){

	int i;
	int j;
	int l;
	float sumObjectiveFunction = 0;
	//foreach candidate solution
	for (i = 0; i < populationSize; i++) {

		///sumObjectiveFunction += objectiveFunctionValues[ i ];

		sumObjectiveFunction += populacao[i].funcaoObjetivo[0];

	}
	//the absolute of the sumObjectiveFunction
	if (sumObjectiveFunction < 0) {
		sumObjectiveFunction = -sumObjectiveFunction;
	}

	//the average of the objective function values
	*averageObjectiveFunctionValues = sumObjectiveFunction / populationSize;



	//the denominator of the equation of the penalty coefficients
	float denominator = 0;
	//the sum of the constraint violation values
	//these values are recorded to be used in the next situation
	float* sumViolation = (float*) malloc(numberOfConstraints * sizeof ( float));
	for (l = 0; l < numberOfConstraints; l++) {

		sumViolation[ l ] = 0;
		for (i = 0; i < populationSize; i++) {

			///sumViolation[ l ] += constraintViolationValues[ i ][ l ] > 0? constraintViolationValues[ i ][ l ]: 0;

			sumViolation[ l ] += populacao[i].v[l] > 0? populacao[i].v[l]: 0;

		}

		denominator += sumViolation[ l ] * sumViolation[ l ];
	}

	//the penalty coefficients are calculated
	for (j = 0; j < numberOfConstraints; j++) {

		penaltyCoefficients[ j ] = denominator == 0? 0: (sumObjectiveFunction / denominator) * sumViolation[ j ];



	}

	//remove auxiliary variables
	free(sumViolation);

}

void calculateAllFitness(Individuo populacao[],/*float* fitnessValues,*/int populationSize,/*float* objectiveFunctionValues,float** constraintViolationValues,*/int numberOfConstraints,float* penaltyCoefficients,float averageObjectiveFunctionValues){

	//indicates if the candidate solution is infeasible
	_Bool infeasible;
	int i;
	int j;
	//the penalty value
	float penalty;
	for (i = 0; i < populationSize; i++) {

		//the candidate solutions are assumed feasibles
		infeasible = 0;
		penalty = 0;

		for (j = 0; j < numberOfConstraints; j++) {

			///if ( constraintViolationValues[ i ][ j ] > 0 ) {
			if(populacao[i].v[j] > 0){
				//the candidate solution is infeasible if some constraint is violated
				infeasible = 1;
				//the penalty value is updated
				///penalty += penaltyCoefficients[ j ] * constraintViolationValues[ i ][ j ];
                penalty += penaltyCoefficients [j] * populacao[i].v[j];
			}

		}

		//the fitness is the sum of the objective function and penalty values
		//if the candidate solution is infeasible and just the objective function value,
		//otherwise
		/**fitnessValues[ i ] = infeasible ?
			(objectiveFunctionValues[ i ] > averageObjectiveFunctionValues? objectiveFunctionValues[ i ] + penalty: averageObjectiveFunctionValues + penalty) :
			objectiveFunctionValues[ i ];
        **/
        populacao[i].fitness = infeasible ? (populacao[i].funcaoObjetivo[0] > averageObjectiveFunctionValues? populacao[i].funcaoObjetivo[0] + penalty: averageObjectiveFunctionValues + penalty) :populacao[i].funcaoObjetivo[0];

	}

}

void AG(int tipoFuncao,int seed,int penaltyMethod,int tamPopulacao,int tamX,int numFilhosGerados,int maxFE,int probCrossover,int tipoES,int sigmaGlobal){
    int fe=0; // Function Evaluation
    int tamPopulacaoFilhos = tamPopulacao*numFilhosGerados;
    int numG,numH,numConstraints;
    Individuo populacao[tamPopulacao];
    Individuo filhos[tamPopulacaoFilhos];
    srand(seed);
    inicializaRestricoes(tipoFuncao,&numG,&numH);
    numConstraints = numG+numH;
    float avgObjFunc;
    float penaltyCoefficients[numConstraints];
    inicializaPopulacaoRestricao(populacao,tamX,tamPopulacao,tipoFuncao);
    avaliaFuncaoRestricao(populacao,tipoFuncao,tamPopulacao,numG,numH,&fe);
    if(penaltyMethod == 1){ // Padrao?  (not apm)
        somaViolacoes(populacao,tamPopulacao,numG,numH);
    }
    else if(penaltyMethod == 2){ // Adaptive Penalty Method ( APM )
        uniteConstraints(populacao,tamPopulacao,numG,numH);
        calculatePenaltyCoefficients(populacao,tamPopulacao,numConstraints,penaltyCoefficients,&avgObjFunc);
        calculateAllFitness(populacao,tamPopulacao,numConstraints,penaltyCoefficients,avgObjFunc);
    }
    else{
        printf("Penalthy mehthod doesn't exist");
        exit(2);
    }
    while(fe < maxFE){
        if(penaltyMethod == 1){ // Padrao?  (not apm)
            selecaoRestricao(populacao,filhos,tamX,tamPopulacao,tamPopulacaoFilhos);
        }
        else if(penaltyMethod == 2){ // Adaptive Penalty Method ( APM )
            selecaoRestricaoAPM(populacao,filhos,tamX,tamPopulacao,tamPopulacaoFilhos);
        }
        if(probabilidadeCrossover(probCrossover) == 1){
            if(TIPO_CROSSOVER == 0){
                crossover(filhos,tamX,tamPopulacaoFilhos);
            }
            else if(TIPO_CROSSOVER == 1){
                crossoverSBX(populacao,filhos,tamX,tamPopulacao,numFilhosGerados);
            }
        }
        mutacao(filhos,tamX,tamPopulacaoFilhos);
        corrigeLimitesX(filhos,tamX,tipoFuncao,tamPopulacaoFilhos);
        avaliaFuncaoRestricao(filhos,tipoFuncao,tamPopulacaoFilhos,numG,numH,&fe);

        if(penaltyMethod == 1){ // Padrao?  (not apm)
            somaViolacoes(filhos,tamPopulacaoFilhos,numG,numH);
        }
        else if(penaltyMethod == 2){ // Adaptive Penalty Method ( APM )
            uniteConstraints(filhos,tamPopulacaoFilhos,numG,numH);
            calculatePenaltyCoefficients(filhos,tamPopulacaoFilhos,numConstraints,penaltyCoefficients,&avgObjFunc);
            calculateAllFitness(filhos,tamPopulacaoFilhos,numConstraints,penaltyCoefficients,avgObjFunc);
        }
        ordenaMelhoresRestricao(populacao,filhos,tamPopulacao,tamPopulacaoFilhos,penaltyMethod);
        elitismoRestricao(filhos,populacao,tamX,tamPopulacao);
        // Pegar melhor individuo
        ordenaMelhoresRestricao(populacao,filhos,tamPopulacao,tamPopulacaoFilhos,penaltyMethod);
        //printf("FE: %i\n",fe);
        if(penaltyMethod == 1){ // Padrao?  (not apm)
            imprimeVioleFO(populacao[0]);
        }
        else if(penaltyMethod == 2){ // Adaptive Penalty Method ( APM )
            imprimeFiteFO(populacao[0]);
        }
    }
}

void DE(int tipoFuncao,int seed,int penaltyMethod,int tamPopulacao,int tamX,int numFilhosGerados,int maxFE,int probCrossover,int tipoES,int sigmaGlobal){
    int fe=0;
    int a,i,j;
    numFilhosGerados=1; // Sempre gera apenas um filho
    int tamPopulacaoFilhos = tamPopulacao*numFilhosGerados;
    int numG,numH,numConstraints;
    Individuo populacao[tamPopulacao];
    Individuo filhos[tamPopulacaoFilhos];
    srand(seed);
    inicializaRestricoes(tipoFuncao,&numG,&numH);
    numConstraints = numG+numH;
    float avgObjFunc;
    float penaltyCoefficients[numConstraints];
    inicializaPopulacaoRestricao(populacao,tamX,tamPopulacao,tipoFuncao);
    avaliaFuncaoRestricao(populacao,tipoFuncao,tamPopulacao,numG,numH,&fe);

    if(penaltyMethod == 1){ // Padrao?  (not apm)
        somaViolacoes(populacao,tamPopulacao,numG,numH);
    }
    else if(penaltyMethod == 2){ // Adaptive Penalty Method ( APM )
        uniteConstraints(populacao,tamPopulacao,numG,numH);
        calculatePenaltyCoefficients(populacao,tamPopulacao,numConstraints,penaltyCoefficients,&avgObjFunc);
        calculateAllFitness(populacao,tamPopulacao,numConstraints,penaltyCoefficients,avgObjFunc);
    }
    else{
        printf("Penalthy mehthod doesn't exist");
        exit(2);
    }

    while(fe < maxFE){
        if(penaltyMethod == 1){ // Padrao?  (not apm)
            selecaoRestricao(populacao,filhos,tamX,tamPopulacao,tamPopulacaoFilhos);
        }
        else if(penaltyMethod == 2){ // Adaptive Penalty Method ( APM )
            selecaoRestricaoAPM(populacao,filhos,tamX,tamPopulacao,tamPopulacaoFilhos);
        }

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
        avaliaFuncaoRestricao(filhos,tipoFuncao,tamPopulacaoFilhos,numG,numH,&fe);
        if(penaltyMethod == 1){ // Padrao?  (not apm)
            somaViolacoes(filhos,tamPopulacaoFilhos,numG,numH);
            selecaoDE(populacao,filhos,tamPopulacaoFilhos);
            imprimeVioleFO(melhorIndividuoRestricao(populacao,tamPopulacao));
        }
        else if(penaltyMethod == 2){ // Adaptive Penalty Method ( APM )
            uniteConstraints(filhos,tamPopulacaoFilhos,numG,numH);
            calculatePenaltyCoefficients(filhos,tamPopulacaoFilhos,numConstraints,penaltyCoefficients,&avgObjFunc);
            calculateAllFitness(filhos,tamPopulacaoFilhos,numConstraints,penaltyCoefficients,avgObjFunc);
            selecaoDEAPM(populacao,filhos,tamPopulacaoFilhos);
            imprimeFiteFO(melhorIndividuoRestricaoAPM(populacao,tamPopulacao));
        }
        //imprimeInformacoesIndividuo(melhorIndividuoRestricao(populacao,tamPopulacao),tamX,numG,numH);
    }
}

void ES(int tipoFuncao,int seed,int penaltyMethod,int tamPopulacao,int tamX,int numFilhosGerados,int maxFE,int probCrossover,int tipoES,int sigmaGlobal){
    int i,j=0;
    int fe=0;
    int tamPopulacaoFilhos = tamPopulacao*numFilhosGerados;
    int numG,numH,numConstraints;
    Individuo populacao[tamPopulacao];
    Individuo filhos[tamPopulacaoFilhos];
    srand(seed);
    inicializaRestricoes(tipoFuncao,&numG,&numH);
    numConstraints = numG+numH;
    float avgObjFunc;
    float penaltyCoefficients[numConstraints];
    inicializaPopulacaoRestricao(populacao,tamX,tamPopulacao,tipoFuncao);
    corrigeLimitesX(populacao,tamX,tipoFuncao,tamPopulacao); // No caso de C01 x[-10,10]
    // Avalia funcao
    avaliaFuncaoRestricao(populacao,tipoFuncao,tamPopulacao,numG,numH,&fe);
    if(penaltyMethod == 1){ // Padrao?  (not apm)
        somaViolacoes(populacao,tamPopulacao,numG,numH);
    }
    else if(penaltyMethod == 2){ // Adaptive Penalty Method ( APM )
        uniteConstraints(populacao,tamPopulacao,numG,numH);
        calculatePenaltyCoefficients(populacao,tamPopulacao,numConstraints,penaltyCoefficients,&avgObjFunc);
        calculateAllFitness(populacao,tamPopulacao,numConstraints,penaltyCoefficients,avgObjFunc);
    }
    else{
        printf("Penalthy mehthod doesn't exist");
        exit(2);
    }
    inicalizaEstrategiaEvolutiva(populacao,filhos,tamX,tamPopulacao,tamPopulacaoFilhos,sigmaGlobal);
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
        autoAdaptacaoSigma(populacao,filhos,tamX,tamPopulacao,numFilhosGerados,sigmaGlobal);
        // Avalia funcao
        corrigeLimitesX(filhos,tamX,tipoFuncao,tamPopulacaoFilhos);
        avaliaFuncaoRestricao(filhos,tipoFuncao,tamPopulacaoFilhos,numG,numH,&fe);

        if(penaltyMethod == 1){ // Padrao?  (not apm)
            somaViolacoes(filhos,tamPopulacaoFilhos,numG,numH);
        }
        else if(penaltyMethod == 2){ // Adaptive Penalty Method ( APM )
            uniteConstraints(filhos,tamPopulacao,numG,numH);
            calculatePenaltyCoefficients(filhos,tamPopulacao,numConstraints,penaltyCoefficients,&avgObjFunc);
            calculateAllFitness(filhos,tamPopulacao,numConstraints,penaltyCoefficients,avgObjFunc);
        }
        ordenaMelhoresRestricao(populacao,filhos,tamPopulacao,tamPopulacaoFilhos,penaltyMethod);
        //ordenaMelhores(populacao,filhos);// NOTE e necessario?
        selecionaMelhoresRestricao(populacao,filhos,tipoES,tamPopulacao,tamPopulacaoFilhos,numFilhosGerados,penaltyMethod);
        if(penaltyMethod == 1){ // Padrao?  (not apm)
            imprimeVioleFO(populacao[0]);
        }
        else if(penaltyMethod == 2){ // Adaptive Penalty Method ( APM )
            imprimeFiteFO(populacao[0]);
        }
    }
    //(*melhorES) = populacao[0];
    //qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaFuncaoObjetivo0);
    //imprimeInformacoesIndividuo(populacao[0]);
}

int main(int argc,char* argv[]){
    /// Se 0: es+ Se 1:es, Se 2: es, modificado
    /// Se 1: Penalty method padrao , Se 2: APM
    /// Se sigma global == 1, sigma global, se nao, vetor de sigmas[]
    /// Parametros tipoFuncao|seed|tamPopulacao|tamX|numFilhosGerados|maxFE|probCrossover|tipoES|sigmaGlobal
    char* method = "null";
    int tipoFuncao = 1;
    int seed = 1;
    int penaltyMethod = 1;
    int tamPopulacao = 50;
    int tamX = 10;
    int numFilhosGerados = 10;
    int maxFE = 200000;
    int probCrossover = 0;
    int tipoES = 0;
    int sigmaGlobal = 1;
    //DE(1,1,50,10,1,20000,0,0,0);
    int opt;
    while ((opt = getopt(argc, argv, "n:f:s:k:p:x:g:m:c:e:t:")) != -1) {
    	switch (opt) {
    	    case 'n': method = optarg; break;
    	    case 'f': tipoFuncao = atoi(optarg); break;
			case 's': seed = atoi(optarg); break;
            case 'k': penaltyMethod = atoi(optarg); break;
			case 'p': tamPopulacao = atoi(optarg); break;
			case 'x': tamX = atoi(optarg); break;
			case 'g': numFilhosGerados = atoi(optarg); break;
			case 'm': maxFE = atoi(optarg); break;
            case 'c': probCrossover = atoi(optarg); break;
			case 'e': tipoES = atoi(optarg); break;
			case 't': sigmaGlobal = atoi(optarg); break;
			default: abort();
    	}
    }
    //method = "DE";
    if(!strcmp(method,"DE")){ // Retorna 0 se str1 == str2
        DE(tipoFuncao,seed,penaltyMethod,tamPopulacao,tamX,numFilhosGerados,maxFE,probCrossover,tipoES,sigmaGlobal);
    }
    else if(!strcmp(method,"ES")){
        ES(tipoFuncao,seed,penaltyMethod,tamPopulacao,tamX,numFilhosGerados,maxFE,probCrossover,tipoES,sigmaGlobal);
    }
    else if(!strcmp(method,"AG")){
        AG(tipoFuncao,seed,penaltyMethod,tamPopulacao,tamX,numFilhosGerados,maxFE,probCrossover,tipoES,sigmaGlobal);
    }
    else{
        printf("Metodo nao encontrado, impresso no codigo em C\n");
        exit(1);
    }
    //AG(1,1,50,10,1,20000,100,0,0);
    //ES(1,1,25,10,100,20000,0,1,0);
    //ES(1,1,3,10,9,20000,0,1,0)
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
