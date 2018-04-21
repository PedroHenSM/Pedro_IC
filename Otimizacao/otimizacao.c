/*
Titulo: Computação Evolucionista Aplicada à Engenharia
Autor: Pedro Henrique Santos
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define TAM_POPULACAO 10
#define NUM_FILHOS_GERADOS 50
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
#define TIPO_FUNCAO 2 // Funcoes com restricao 14 e 13 = 3g,0h | 01 = 2g,0h
#define TIPO_ES 0 // Se 0: es+ Se 1:es, Se 2: es, modificado
#define SIGMA_GLOBAL 1 // Se sigma global == 1, sigma global, se nao, vetor de sigmas[]
#define NUM_RESTRICOES_G 2
#define NUM_RESTRICOES_H 1
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

void selecaoRestricao(Individuo populacao[],Individuo filhos[],int tamanhoX){ /// Torneio NOTE: NECESSARIO GERAR MAIS FILHOS AQUI?!
    int i,j; // Passa os ''melhores'' para filhos
    for(i=0;i<TAM_POPULACAO_FILHOS;i++){ // Sele��o
        int indice1,indice2;
        indice1 = rand () % TAM_POPULACAO;
        indice2 = rand () % TAM_POPULACAO;
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

void crossover (Individuo filhos[],int tamanhoX){ /// NOTE: NECESSARIO GERAR MAIS FILHOS AQUI?!
    int i,l;
    for(i=0;i<TAM_POPULACAO_FILHOS;i+=2){ // Crossover
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

void crossoverSBX(Individuo populacao[],Individuo filhos[],int tamanhoX){ /// NOTE: NECESSARIO GERAR MAIS FILHOS AQUI?!
    int i,j,k;
    for(i=0;i<TAM_POPULACAO;i+=2){
        for(k=0;k<NUM_FILHOS_GERADOS;k++){
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

void algoritmoGeneticoGeraDados(int tamanhoX,int seed){
    int contaObj=0,contador;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO];
    srand(seed);
    inicializaPopulacao(populacao,tamanhoX);
    avaliaFuncao(populacao,&contaObj,tamanhoX,TAM_POPULACAO);
    while(contaObj < MAX_CALC_OBJ){
        contador=10-tamanhoX;
        selecao(populacao,filhos,tamanhoX);
        if(probabilidadeCrossover(PROB_CROSSOVER) == 1){
            if(TIPO_CROSSOVER == 0){
                crossover(filhos,tamanhoX);
            }
            else if(TIPO_CROSSOVER == 1){
                crossoverSBX(populacao,filhos,tamanhoX);
            }
        }
        mutacao(filhos,tamanhoX);
        avaliaFuncao(filhos,&contaObj,tamanhoX,TAM_POPULACAO_FILHOS);
        ordenaMelhores(populacao,filhos,0);
        elitismo(filhos,populacao,tamanhoX,TAM_POPULACAO);
        // Pegar melhor individuo
        printf("%i %i %i\t",tamanhoX,TAM_POPULACAO,PROB_CROSSOVER);
        imprimeContaObj(contaObj);
        while(contador>0){
            //printf("1.234567\t");
            printf("---------\t");
            contador--;
        }
        imprimeMelhores(populacao,filhos,tamanhoX);
    }
}

void geraDados(){
    int s,n;
    for(n=2;n<11;n++){ // TamanhoProblema
        for(s=1;s<31;s++){// seed
            algoritmoGeneticoGeraDados(n,s);
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
                crossover(filhos,TAM_X);
            }
            else if(TIPO_CROSSOVER == 1){
                crossoverSBX(populacao,filhos,TAM_X);
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

void inicializaPopulacaoRestricao(Individuo populacao[],int tamanhoX,int tamanhoPop,int tipoFuncao){
    int i,j;
    for(i=0;i<tamanhoPop;i++){
        for(j=0;j<tamanhoX;j++){
            float eq = (rand()/(float)(RAND_MAX)); // gera valores entre [0,1]
            float val;
            if (tipoFuncao == 1){
                val = eq * (rand() % 11); // Valores decimais entre [0,10] C01
            }
            else if(tipoFuncao == 2){
                val = eq * (rand() % 11 - 5); // Gera valores decimais entre [-5,5] C02
            }
            else if(tipoFuncao == 3 || tipoFuncao == 14){
                val = eq * (rand() % 2001 - 1000); // Valores decimais entre [-1000,1000]
            }
            else if(tipoFuncao == 4){
                val = eq * (rand() % 101 - 50); // Valores decimais entre [-50,50}
            }
            else if(tipoFuncao == 13 || tipoFuncao == 10){
                val = eq * (rand() % 1001 - 500); // valores decimais entre [-500,500]
            }
            else{
                printf("Funcao nao encontrada\n");
                exit(0);
            }
            populacao[i].x[j] = val;
        }
    }
}

void ordenaMelhoresRestricao(Individuo populacao[],Individuo filhos[]){ // Ordena com base nas violacoes
    qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
    qsort(filhos,TAM_POPULACAO_FILHOS,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
}

void C01 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh){ // x[0,10]
  int j;
  float f1, f2, f3, g1, g2, e[nx];
  float o[30] = {0.030858718087483,	-0.078632292353156,	0.048651146638038,	-0.069089831066354,	-0.087918542941928,	0.088982639811141,	0.074143235639847,	-0.086527593580149,	-0.020616531903907,	0.055586106499231,	0.059285954883598,	-0.040671485554685,	-0.087399911887693,	-0.01842585125741,	-0.005184912793062,	-0.039892037937026,	0.036509229387458,	0.026046414854433,	-0.067133862936029,	0.082780189144943,	-0.049336722577062,	0.018503188080959,	0.051610619131255,	0.018613117768432,	0.093448598181657,	-0.071208840780873,	-0.036535677894572,	-0.03126128526933,	0.099243805247963,	0.053872445945574};
  for (j = 0; j < nx; j++){
	  e[j]= x[j] - o[j];
  }
  /* objective function */
  f1 = 0.;
  f2 = 1.;
  f3 = 0.;
  g1 = 1.;
  g2 = 0.;
  for (j = 0; j < nx; j++)
    {
      f1 = f1 + pow (cos (e[j]), 4);
      f2 = f2 * cos (e[j]) * cos (e[j]);
      f3 = f3 + ((float) (j + 1)) * e[j] * e[j];
      g1 = g1 * e[j];
      g2 = g2 + e[j];
    }
  f[0] = fabs ((f1 - 2 * f2) / sqrt (f3));
  f[0] = -f[0];
  g[0] = 0.75 - g1;
  g[1] = g2 - 7.5 * ((float) nx);
}

void C02 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh){ // x[-5.12,5.12]
  int j;
  float f1, g1, g2, h1, e[nx];
  float PI = 4.0 * atan (1.0);
  float o[30] = {-0.066939099286697,	0.470966419894494,	-0.490528349401176,	-0.312203454689423,	-0.124759576300523,	-0.247823908806285,	-0.448077079941866,	0.326494954650117,	0.493435908752668,	0.061699778818925,	-0.30251101183711,	-0.274045146932175,	-0.432969960330318,	0.062239193145781,	-0.188163731545079,	-0.100709842052095,	-0.333528971180922,	-0.496627672944882,	-0.288650116941944,	0.435648113198148,	-0.348261107144255,	0.456550427329479,	-0.286843419772511,	0.145639015401174,	-0.038656025783381,	0.333291935226012,	-0.293687524888766,	-0.347859473554797,	-0.089300971656411,	0.142027393193559};
  for (j = 0; j < nx; j++)
  {
	  e[j]= x[j] - o[j];
  }
  /* objective function */
  f1 = e[0];
  g1 = 0.;
  g2 = 0.;
  h1 = 0.;
  for (j = 0; j < nx; j++)
    {
      if (e[j] > f1)
	  {
		  f1 = e[j];
	  }

	  g1 = g1 + (e[j] * e[j] - 10 * cos(2 * PI * e[j]) + 10);
	  g2 = g2 + (e[j] * e[j] - 10 * cos(2 * PI * e[j]) + 10);
	  h1 = h1 + ((e[j] - 0.5) * (e[j] - 0.5) - 10 * cos(2 * PI * (e[j] - 0.5)) + 10);
    }


  f[0] = f1;
  g[0] = 10.0 - g1/((float) nx);
  g[1] = -15.0 + g2/((float) nx);
  h[0] = -20.0 + h1/((float) nx);
}

void C03 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh){ // x[-1000,1000]
    int j;
    float f1, h1, e[nx];
    float o[30] = {111.17633500088529,	92.07880492633424,	417.9818592609036,	253.16188128024302,	363.5279986597767,	314.334093889305,	187.32739056163342,	240.4363027535162,	422.60090880560665,	327.63042902581515,	62.04762897064405,	25.435663968682125,	360.56773191905114,	154.9226721156832,	33.161292034425806,	177.8091733067186,	262.58198940407755,	436.9800562237075,	476.6400624069227,	331.2167787340325,	75.205948242522,	484.33624811710115,	258.4696246506982,	419.8919566566751,	357.51468895930395,	166.3771729386268,	47.59455935830133,	188.20606700809785,	184.7964918401363,	267.9201349178807};
    for (j = 0; j < nx; j++){
        e[j]= x[j] - o[j];
    }

  /* objective function */
    f1 = 0.;
    h1 = 0.;
    for (j = 0; j < (nx - 1); j++){
        f1 = f1 + (100 * pow((e[j] * e[j] - e[j+1]),2) + pow((e[j] - 1),2));
        h1 = h1 +  pow((e[j] - e[j+1]),2);
    }
    f[0] = f1;
    h[0] = h1;
}

void C04 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh){
  int j;
  float f1, h1, h2, h3, h4, e[nx];
  float o[30] = {0.820202353727904, 5.260154140335203,	-1.694610371739177,	-5.589298730330406,	-0.141736605495543,	9.454675508078164,	8.795744608532939,	9.687346331423548,	-3.246522827444976,	6.647399971577617,	1.434490229836026,	-0.506531215086801,	0.558594225280784,	7.919942423520642,	1.383716002673571,	-1.520153615528276,	-2.266737465474915,	6.48052999726508,	-8.893207968949003,	-3.528743044935322,	6.063486037065154,	-4.51585211274229,	7.320477892009357,	-8.990263774675665,	9.446412007392851,	-6.41068985463494,	-9.135251626491991,	2.07763837492787,	8.051026378030816,	-1.002691032064544};
  for (j = 0; j < nx; j++)
  {
	  e[j]= x[j] - o[j];
  }

  /* objective function */
  f1 = e[0];
  h1 = 0.;
  h2 = 0.;
  h3 = 0.;
  h4 = 0.;
  for (j = 0; j < nx; j++)
    {
      if (e[j] > f1)
	  {
		  f1 = e[j];
	  }

	  h1 = h1 + e[j] * cos(sqrt(fabs(e[j])));
      h4 = h4 + e[j];
    }

  for (j = 0; j < ( nx/2 - 1); j++)
    {

	  h2 = h2 +  pow((e[j] - e[j+1]),2);

    }
  for (j = nx/2  ; j < (nx - 1); j++)
    {

	  h3 = h3 +  pow((e[j] * e[j] - e[j+1]),2);

    }
  f[0] = f1;
  h[0] = h1/((float) nx);
  h[1] = h2;
  h[2] = h3;
  h[3] = h4;
}

void C10 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh){
  int j, i, k, l;
  float f1, h1, y[nx], M[nx][nx], e[nx];
  float o[30] = {-41.03250252873486,	-35.70280591875908,	-48.66938576680659,	94.51946988004894,	31.68700466174738,	99.69508270219342,	30.778279925351967,	-31.041222172110807,	-46.21010370947247,	27.26190010072706,	-2.093622677920422,	22.246274570582585,	-42.887366421312436,	89.88377145577851,	-6.731523713182725,	97.86439204258224,	49.49993772881544,	23.210695390854696,	-81.36716857155828,	-20.15688556597543,	36.692155371634726,	44.37408948075327,	-15.984549833405907,	-49.68391424581281,	98.3715576810595,	0.127593155843627,	61.709914317965655,	-84.0189999580673,	-35.39565398431638,	-5.143979333218638};

  float M1[10][10] = {-0.040348736997622,	-0.57547478301341,	-0.389765935131311,	 0.140885590781398,	-0.144670835908614,	 0.115450535800564,	-0.518304793509814,	 0.149470047989466,	 0.386138410795004,	 0.145230906380518,
                      -0.065866355431932,	-0.257122057117828,	 0.298076832606074,	-0.339838329245853,	 0.36335902728584,	 0.69140513648482,	-0.028135015592,	 0.214592766458436,	-0.21337789362989,	 0.150888625576118,
                      -0.154282201350543,	 0.233990300992969,	 0.196071677766854,	-0.516011967110319,	-0.265721705059327,	-0.060905234855331,	 0.017107816544837,	 0.459602625258773,	 0.548125243387881,	-0.174490925899168,
                       0.363857337521185,	 0.622123871321307,	 0.051207493885242,	 0.281214293456459,	 0.077610534696611,	 0.38475834924149,	-0.410410905293451,	-0.035659107285799,	 0.245258564359423,	 0.122242429111897,
                       0.050671935660746,	 0.081648277010661,	-0.034403618638332,	-0.604379533928053,	 0.029715665982435,	-0.255779368776897,	-0.321832198773717,	-0.454789288464763,	-0.04205042244723,	 0.49580649289475,
                      -0.580165394752893,	 0.155095038022222,	 0.210065625307939,	 0.265731671126386,	-0.481942328396859,	 0.107406464626288,	-0.031135086609527,	 0.076014134025855,	-0.182793563791557,	 0.490558019238156,
                       0.525119522958187,	-0.007317154673771,	-0.078860784151735,	-0.118095618729855,	-0.388828033497333,	-0.11024307003898,	-0.140281776355321,	 0.496093855143328,	-0.517431971433034,	 0.084650093567346,
                      -0.221934029584865,	 0.08937666134435,	-0.161579047653623,	-0.202383925865339,	-0.358394121525824,	 0.305489267359872,	-0.319803594563591,	-0.329350210528224,	-0.279367850055518,	-0.604257992994783,
                      -0.14024729892381,	-0.085950562359788,	 0.562192927749687,	 0.15621897836716,	 0.259175107289986,	-0.401255905586421,	-0.564318107548709,	 0.125407586464797,	-0.135691838706447,	-0.227388734136487,
                      -0.392433919094433,	 0.338788001793312,	-0.569620515414919,	-0.064134166835177,	 0.436408519589215,	-0.115926326668514,	-0.12994296127687,	 0.36903024504463,	-0.213067396847211,	 0.016735182121047};

  float M2[30][30] = {0.210498733952046,	0.044384555963486,	0.039557467453157,	0.322626150860574,	0.090728331050138,	-0.046559513284792,	0.279058079389439,	0.216086043131122,	0.001910966340658,	0.005653042988124,	0.055129408353881,	-0.157826900224939,	0.126895165086409,	0.21964276539513,	-0.030777527838514,	-0.174795699173724,	-0.071684319926197,	-0.435002022420347,	0.300879515782067,	0.196965690878292,	0.015408513092551,	-0.088767209425526,	-0.056721005316589,	0.004563253603378,	-0.411103445073825,	-0.098798851715865,	-0.017178459384374,	0.223973235768929,	-0.038630814336146,	0.153324298640425,
                     -0.123853898371964,	0.006208844977804,	-0.196968610000223,	-0.297886862654356,	0.337304789619692,	0.151336922415474,	-0.190737243250796,	0.371272076371965,	0.05154807935738,	-0.089924110414285,	-0.147609151282806,	-0.083802701820601,	0.017052091139248,	0.207360614681,	  -0.014730657062048,	0.030357000089238,	-0.38493302429957,	-0.098471468276102,	-0.24989513579197,	0.265831026621616,	0.039327358400047,	0.102952558500907,	0.112835401766515,	0.075442015984702,	-0.050054414064672,	0.152147800707443,	0.275356821244483,	-0.075104525195784,	0.130811841464527,	0.128154556261642,
                     -0.247552776888796,	-0.170900469984963,	-0.025609465501807,	-0.074857181493926,	0.06402729466598,	-0.249589575845155,	0.291141404096232,	-0.151955001483896,	-0.182101792706067,	0.177528281451576,	-0.11080027107009,	-0.101265961425373,	0.05885945181901,	0.282705041568106,	-0.075385821998282,	-0.183086148754718,	-0.03660323576297,	0.095170386364688,	0.004340336079074,	0.000556182984035,	0.321477491995385,	0.031769482356446,	-0.111377507159358,	0.532122152715201,	0.088882029439331,	-0.098799193147706,	0.101408217806575,	-0.110636423867126,	-0.211274352044764,	-0.155079279400244,
					  0.025663646265114,	0.102509721975959,	0.309345681084264,	0.070334750076729,	0.097241604458889,	-0.273266320148514,	-0.239907149554807,	0.150740659614263,	0.069994908090163,	-0.039122490587809,	0.125893036315579,	-0.340093749710236,	0.154051275851787,	0.186934934550147,	-0.25121208048304,	-0.196077026994064,	0.01583365232178,	0.55909380072088,	0.112181997526451,	-0.008728203934299,	-0.217367189352171,	-0.106529307174376,	0.050536269466879,	-0.006551233713265,	-0.045903154193401,	0.167087870245007,	0.002122466601521,	-0.012682827041684,	-0.002595974426726,	0.093386409335028,
                     -0.096757718998012,	-0.039650722973667,	-0.056221933017481,	0.247299467670382,	0.12500615760951,	-0.129535632819352,	0.037407381782829,	-0.049321379462225,	0.022699926035277,	0.192965409321344,	0.104996464498579,	0.085577579883565,	0.078833431392308,	0.19343561654434,	0.253727682682793,	0.025915952751218,	-0.084868665088399,	0.037854288550344,	-0.393242682132416,	-0.312008732104883,	0.022684882865586,	0.127011588278839,	-0.048595771065928,	-0.309159369729801,	-0.363341449374355,	0.16665805071768,	-0.056515345289132,	-0.298609174881293,	-0.294119559817135,	0.099285471881439,
					  0.043186472508942,	0.001771226886988,	0.14263426111345,	-0.026446427851299,	-0.042202844595988,	0.229858934929016,	0.192708259881167,	0.454739275512023,	0.078993562199931,	0.012272017371985,	0.208615477072962,	0.0286479351207,	-0.078639302217307,	-0.271449720411795,	0.183880260881944,	-0.221305977752936,	0.090741804187672,	0.221668479563272,	-0.18125553211922,	0.202526121479001,	0.143690682699279,	-0.009428483121244,	-0.138122988901157,	-0.0087984112929,	0.180843519453226,	-0.167186223826231,	0.059308954954177,	0.044175554473605,	-0.451367400847036,	0.155461979668354,
					  0.208459638473921,	-0.192313493572805,	0.042799453883608,	0.056319917449526,	-0.064050781971354,	0.336443527835829,	-0.024913686145588,	-0.135392681180607,	0.084208982811827,	0.19368687774262,	0.273506736028856,	0.161025819114185,	-0.145878857750618,	0.175226077989248,	0.197067036614125,	0.01612829628023,	-0.266892552385116,	0.294523656999125,	0.111195525470887,	-0.201292974303023,	0.007989116024197,	-0.133381335708404,	-0.05978698001368,	0.163532731543214,	-0.166898927792151,	-0.069407809244429,	0.410416043355545,	0.195220789227849,	0.206696988670138,	-0.030109638400114,
					 -0.161888866215218,	-0.477833679054076,	0.090201553328432,	0.197276529552719,	0.257525091427261,	0.234933865087851,	0.111554961171882,	0.039671046122422,	0.1024028037035,	-0.275023229663029,	-0.179789894147703,	-0.149367349549056,	0.097199867696211,	0.038487537865377,	-0.102826054337925,	0.041005419958829,	0.071428957044728,	0.133047965356434,	-0.271733080060097,	-0.07745323185521,	-0.074759189714073,	0.047497229049009,	-0.165726777593199,	-0.084379708142798,	-0.046302009416911,	-0.250215681879312,	-0.226883004983199,	0.246816548042778,	0.174568318097348,	-0.215049387958273,
					  0.163861486371016,	-0.148103013020227,	-0.198506174180913,	-0.04219384123679,	0.118302213986748,	0.022695570191254,	0.206707310041736,	0.156232002478571,	0.093575716314624,	0.487161837511012,	-0.277679678378239,	-0.05898824271195,	-0.101340463707747,	-0.197187024207523,	-0.356221119308529,	0.175907621024382,	0.025724776145916,	0.177555844621935,	0.186546393359421,	-0.011515740909936,	-0.063280757601228,	0.108572545320813,	0.250175342159778,	-0.217188740737622,	-0.118288673330692,	-0.018858445496624,	0.125813534590454,	-0.005766830114285,	-0.173770830530608,	-0.204076434641454,
					  0.076726374752626,	-0.059157318075244,	-0.109120969315065,	0.281735964889155,	0.113701355600812,	0.057934126040283,	0.056930714392195,	-0.260311629272444,	0.066131852169431,	0.254150636558942,	-0.123595487021829,	-0.068770016151727,	-0.084370591918571,	-0.210531372886378,	-0.058314753673756,	-0.373308427361442,	-0.015030712250411,	0.013086120163669,	-0.210199938641698,	0.044179063442584,	-0.090836039745107,	0.110863603159718,	0.197449444025479,	0.163833987692559,	0.172463370094917,	-0.02624766621562,	-0.075091307788291,	-0.035678286150929,	0.258868590341026,	0.540009857641777,
					 -0.11448537731886,	     0.158421318146668,	-0.119175728122649,	0.201903441501447,	-0.009316901256506,	0.173542979966456,	-0.234634060687318,	-0.135893056643097,	0.271225318340463,	0.108809056996355,	-0.016309462418078,	0.118496233811428,	0.304052823739158,	-0.007722437001817,	0.094448421400381,	0.034726127685042,	-0.116005525261297,	0.078981290570637,	0.205206605235178,	0.30369188618657,	-0.312746271441503,	0.393986794526704,	-0.195899452039796,	0.20494103596369,	-0.018949028418705,	-0.184962367629022,	-0.079402591539207,	-0.121729951631301,	-0.181982754615953,	-0.156530779671702,
					  0.162253760728648,	-0.153789438331871,	0.00868657051539,	-0.131911170907241,	0.420552377770407,	-0.216321029385655,	-0.133791045486745,	-0.270383692536091,	-0.113959519847022,	-0.182238016801448,	0.191497999470082,	-0.054919230163731,	0.162994884573224,	-0.222520719869188,	0.021240982443138,	0.090745017688942,	-0.291288801107928,	-0.047608464418265,	0.152396018978015,	-0.070526756136424,	0.048450789083476,	-0.045567922857345,	0.102798772566309,	-0.159656473982117,	0.118251926660962,	-0.430675707603667,	0.10458086074914,	-0.029209095590176,	-0.242559117346188,	0.146082148041515,
					 -0.054769889978099,	-0.102998778221723,	0.064729032089379,	0.181703388132909,	-0.238051881438163,	0.125211106412204,	-0.301322664395015,	-0.055687348687259,	-0.039615420818266,	-0.10893447819021,	-0.225903659932295,	-0.168752128472236,	0.050214382872752,	-0.209623736459376,	-0.021152472752855,	-0.530138693435876,	-0.134544018857329,	-0.145183011013022,	0.052406331352567,	-0.066953538369538,	0.210839727968681,	-0.005790710509523,	0.040900628077853,	-0.219106039109813,	-0.058779545307731,	0.037152749763887,	0.306738266478922,	-0.155410315881546,	-0.059669333754521,	-0.307937539980894,
                      0.188677164396131,	0.102763453891941,	-0.392205809826425,	0.023274604305339,	-0.241062713234987,	-0.134335764348244,	-0.156920883882407,	0.118394352690385,	-0.242783527480707,	0.148199620429226,	0.038598382023358,	-0.373777385965017,	-0.079990650844341,	0.109379832985269,	0.276235130978874,	0.000424181313047,	-0.1581803206613,	0.194537934564581,	-0.112037113974561,	0.092669956020601,	0.151181460492394,	0.064779413736991,	0.03298441610138,	-0.105713078293658,	0.026477488720875,	-0.301872238837486,	-0.304572846125097,	0.129303116825904,	0.150275305655397,	-0.143819273917584,
					 -0.121136819961029,	-0.280376857590378,	-0.142884824324944,	-0.234336987516145,	-0.27318548804741,	-0.062906552289143,	0.117125592225031,	-0.111153019559446,	0.170636020457398,	-0.049474387767955,	0.250337970843838,	-0.015737111315197,	0.303098653510777,	-0.031695841766391,	-0.018294649882712,	-0.004853390742363,	0.288923483741941,	0.102932886997556,	-0.030117891708216,	0.348873109957197,	0.107864455557603,	-0.053217559061395,	0.155554070063422,	-0.131389140271268,	-0.306220365198726,	-0.153022880982604,	0.128394915701578,	-0.254920058505309,	0.228391386551728,	0.124793864045053,
					 -0.085458356559741,	-0.06749024827502,	-0.243684986433717,	0.195223151576616,	-0.056144203577952,	-0.107022717570903,	-0.119766313023139,	0.21180214818079,	0.286641525002163,	-0.256687343688466,	-0.016698212375744,	0.168676200561914,	-0.300510419700627,	0.3859169863873,	-0.025097350370833,	-0.064878242312045,	0.223840463303938,	-0.013692765327541,	0.191347293229108,	-0.221476004310888,	0.014539395850687,	0.108928515315896,	0.2763380292334,	-0.022071088238979,	0.190680276913542,	-0.305290419919566,	0.107640319223234,	-0.104538877590241,	-0.084705533603217,	0.108822927690725,
					  0.359256748293047,	0.154976198577344,	-0.116153265766918,	0.13931677738968,	0.374496005165856,	0.178417709752365,	-0.020215045948477,	0.0511359434473,	-0.30173528700312,	-0.216545463965152,	-0.113137661147475,	0.186762444883836,	0.040082745946815,	-0.053364283388025,	0.065488033743007,	-0.033154988554239,	0.317998692979354,	0.247815755943833,	0.126243283202638,	0.038451263306162,	0.256900869660439,	0.043331197504268,	-0.051702015143296,	0.115290293077177,	-0.162878255904405,	0.042941144091639,	-0.0273285865418,	-0.345298257763101,	0.15936874911257,	-0.076890125653739,
					 -0.28099424578992,	    0.107396814515739,	-0.084563105000033,	0.251259020314969,	-0.046244227618743,	0.194849917362354,	0.105731360526473,	-0.116729435444227,	-0.125871612130733,	-0.004143188817115,	0.031303723499165,	0.193842070959737,	-0.071047584389393,	0.163798237107731,	-0.33103883422486,	0.083042404997121,	-0.314536638066646,	0.1850656428806,	0.139857282060989,	0.258138857521682,	0.241402804816099,	-0.267773012822818,	-0.136618598358345,	-0.352889529882469,	0.178761358680642,	0.004675588386463,	-0.130523685747766,	-0.11010510100568,	0.012337514111924,	0.133039188998738,
					 -0.136214748186763,	0.061446441408168,	0.124443001628158,	-0.384194680239079,	0.161386673633437,	0.05151848059208,	0.050237897866986,	0.092332554624287,	0.025301015752734,	0.144267170372807,	0.177218756700414,	0.330110661903483,	-0.032723550578364,	0.032862846981569,	-0.039329546595698,	-0.479455648611786,	-0.142384919826993,	-0.00804128544337,	0.153527223636046,	-0.147592550588846,	0.01671488256034,	0.153461668336777,	0.143391962131602,	-0.1032423296903,	-0.127806420761444,	-0.112212283814506,	-0.419393741016172,	0.011554752050117,	0.174489436742579,	-0.152583142805368,
					 -0.183525414815911,	0.084264423476034,	-0.513582573163132,	-0.060641305770964,	0.096244159556264,	0.049715888666474,	0.311001693872411,	0.001259440460654,	0.036651495187505,	-0.187364828692437,	-0.017792472149735,	-0.108962196548542,	0.106923649625955,	-0.137418200183617,	0.231554866808663,	-0.211864650236558,	-0.056173351885998,	0.108205208097382,	0.152251721125532,	-0.172564127758839,	-0.379707703841598,	-0.399495209658996,	-0.017301373440465,	0.016739053799758,	0.035848677241745,	0.150058535685014,	0.019465311724744,	-0.030054112339037,	-0.049980094722222,	-0.078081830537862,
					  0.061058425410664,	0.042041509977416,	-0.211003741431923,	-0.126388800908604,	0.039219152461988,	0.225495863827954,	-0.366358041819373,	-0.190018182245224,	-0.135422744295414,	0.106462362604256,	-0.028597656153562,	0.151442873263496,	0.108762194743749,	0.202895035408, 	-0.220198514387979,	-0.145116392843135,	0.287960518621411,	-0.046429895952516,	-0.236435128849732,	0.064337841563859,	-0.044861474741858,	-0.323541240367847,	0.116777979642275,	0.058194390691809,	-0.143859430803049,	-0.058357266848021,	-0.06663691339304,	0.323975481839132,	-0.363821240073463,	0.025044920981486,
					 -0.189525638545991,	-0.006677446894394,	0.040135587178237,	0.030562197515064,	-0.081557175019313,	-0.35988601916995,	-0.046748366964196,	0.052796424347975,	-0.293965333147152,	-0.107256786396328,	-0.266428993150317,	0.331112621223704,	-0.3240691448606,	-0.183930905848232,	0.033931266862913,	-0.056733070481899,	-0.071906584046074,	0.165910393063134,	-0.066596131882292,	0.187535095996836,	-0.282573695097042,	-0.012822431812993,	-0.130628684496259,	0.019603565202997,	-0.36530186945985,	-0.218948877174082,	0.15928161680357,	0.117469132881907,	0.02739594774293,	0.086724839893498,
					 -0.046415954902964,	-0.160644004446126,	-0.033524440660352,	0.07500744681351,	0.016322145454599,	-0.107804845719038,	-0.239974742220443,	0.045234085507995,	0.410025327436819,	-0.041785955831231,	-0.088824232512646,	0.019324971910103,	-0.164179563969319,	-0.274883654548913,	0.08344816321592,	0.134361131744622,	-0.168293304036166,	0.050119827090831,	0.107998968576299,	0.024222465405286,	0.333568438641683,	-0.314273943234102,	0.088627173413873,	0.307987145180008,	-0.25957029953184,	0.139693331842536,	-0.366318566467329,	-0.024315207215651,	-0.097866140810133,	0.031984120313315,
					  0.414611908814294,	-0.37112136873643,	-0.073866553818035,	-0.29087757946088,	-0.217072217367287,	-0.029123302474621,	0.010661588246943,	-0.056648186984039,	0.052448395670158,	-0.125771363235671,	-0.323377762528554,	0.103065249994115,	0.03990538947911,	0.20892679696363,	0.025498721902253,	-0.133437428932213,	-0.116769147892009,	0.075000348696883,	0.147043486128104,	0.012254394305512,	-0.070162283643009,	0.104480173528188,	-0.343231962242353,	-0.124041674266431,	0.079443082346753,	0.192704701695008,	-0.169380089268582,	-0.054285719705745,	-0.145500023209769,	0.255395103787066,
					  0.207119315322088,	-0.320073115712564,	0.0559314734117,	0.172649138634129,	-0.011542461945528,	-0.130620764160256,	-0.10888053008216,	0.267564406009138,	-0.120890663229249,	0.202202465644977,	0.149425179881808,	0.216934380895458,	0.065653150449446,	0.041794472366999,	-0.048128520065833,	-0.017549459873998,	-0.077842831218385,	-0.198170596730465,	-0.149865167620775,	0.12634718117617,	-0.30097962077431,	-0.356324559219511,	0.015984934796442,	0.112407601721805,	0.237966910014296,	-0.056653284700555,	-0.084496688943403,	-0.372456205768623,	0.097560810635958,	-0.241013145399317,
					  0.042993730292494,	0.004956761653254,	0.053895052194156,	-0.088854185716889,	0.326910901526066,	-0.087794341394476,	-0.020802090219748,	-0.298708928515842,	0.25242600994189,	0.101197681470221,	0.07152866676548,	-0.181921564652412,	-0.505149592427631,	0.143193361309159,	0.186981040042569,	-0.133652828493287,	0.200265248087893,	-0.073752726703129,	-0.000238503991611,	0.374550360408052,	-0.050281630947108,	-0.092458775511721,	-0.205541771113341,	-0.211350594454595,	0.048295879169752,	0.061913849844292,	0.046521349189447,	-0.065574597414797,	0.000727202776555,	-0.214369198817042,
					 -0.372708113027754,	-0.273111422565572,	-0.028044496936687,	-0.031315545207917,	0.128692442270306,	0.103912301156068,	-0.317409754881132,	0.193362522820359,	-0.263946940397186,	0.328539161755314,	0.014885444611456,	-0.135231025676715,	0.012770730354147,	-0.056470973842145,	0.187010128145666,	0.06551126700651,	0.264488080744978,	-0.092362300447396,	0.374463443107909,	-0.114265316590068,	0.014321931394336,	-0.014287975486322,	-0.251977653040875,	-0.068814134318263,	0.037272928895606,	0.036847910078949,	-0.004344868007102,	0.000682914575515,	0.075677153681488,	0.26759596689543,
					  0.064818689916455,	0.284053411520265,	-0.063930667567762,	-0.055888531296039,	0.084971324574045,	-0.248553864049376,	-0.053882957221272,	0.109974358415086,	0.361665795220151,	0.212622412400492,	-0.222105881590724,	0.159407186211676,	0.277766857748719,	-0.040024212790115,	-0.045419517646103,	-0.049565963430445,	0.075309844691921,	-0.01462002387698,	-0.125901695229838,	-0.179360862877057,	0.201763873122918,	-0.220450927984847,	-0.414575162457279,	-0.112162765618009,	0.098487134259606,	-0.239965053948616,	0.170255100801187,	0.078980100492884,	0.243403068580785,	0.021368523436668,
					 -0.054279013736138,	0.018659420049998,	0.330359950213943,	0.023123143363581,	0.035090878248845,	0.025527653285977,	0.082763577820282,	-0.044575644917522,	-0.004626086242742,	0.103520181686413,	-0.432497376883502,	0.119032674935153,	0.236464940348734,	0.185294902781396,	0.510006376075268,	0.054475637131684,	0.021019631752974,	0.107026280159219,	0.091553865611036,	0.169591926725838,	0.02088417330554,	-0.183147177835376,	0.427381645911547,	-0.13187711537175,	0.138348964310852,	-0.037276342733556,	-0.010728013977158,	0.085070861966194,	-0.008940851260341,	0.001071606497694,
					 -0.001220559638941,	-0.207365392796594,	-0.23185510790523,	0.208482953784606,	0.122000548923194,	-0.322922674036945,	-0.041150767726355,	0.035241090632951,	-0.012254747396228,	-0.030675366498019,	0.216247246460075,	0.297863251529541,	0.168847156819068,	-0.082183953429559,	0.065957904193529,	-0.09572234076031,	0.040748898076806,	0.06592096841622,	0.028275231975731,	0.180261847132564,	0.167295789767305,	0.208499007714482,	0.06073181729704,	-0.09136260314961,	0.195634463938498,	0.405036010088252,	0.090356343227339,	0.443504939597125,	0.041746547462442,	-0.143522086623263};

 for (j = 0; j < nx; j++)
  {
	  e[j]= x[j] - o[j];
	  y[j] = 0.;
	  if (nx == 10)
	  {
	   for ( l = 0; l < nx; l++)
		    M[l][j] = M1[l][j];
	  }
	  if (nx == 30)
	  {
	   for ( l = 0; l < nx; l++)
		    M[l][j] = M2[l][j];
	  }
  }

  for (i = 0; i < nx; i++)
  {
	  for (k = 0; k < nx; k++)
	  {
		  y[i] = y[i] + e[k] * M[k][i];
	  }
  }

  /* objective function */
  f1 = 0.;
  h1 = 0.;

  for (j = 0; j < (nx - 1); j++)
    {

	  f1 = f1 + (100 * pow(((e[j]+1) * (e[j]+1) - (e[j+1]+1)),2) + pow(((e[j]+1) - 1),2));

    }
  for (j = 0; j < nx; j++)
    {
      h1 = h1 + y[j] * sin(sqrt(fabs(y[j])));

    }
  f[0] = f1;
  h[0] = h1;

}

void C13 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh){ // x[-500,500]
    int j;
    float f1, g1, g2, g3, e[nx];
    float PI = 4.0 * atan (1.0);
    float o[30] = {69.69311714880897,	1.509803311435702,	67.6746198312362,	80.43173609273597,	80.47622449424348,	51.21092936019716,	52.7723719926014,	17.248465789326257,	52.40150903116374,	39.64846247456716,	89.86375903333635,	32.079301315169474,	43.192499277837946,	70.79294586561508,	1.48440984483988,	19.8566700417119,	29.502667246412756,	34.256788127976684,	12.643016541338264,	78.57234385195876,	26.51647349482587,	97.06430708087798,	10.180504722002471,	82.90799886855778,	63.540231382573154,	74.78243308676124,	87.20817289266436,	50.779655804893764,	43.05412185616204,	33.862234518700916};
    for (j = 0; j < nx; j++)
    {
        e[j]= x[j] - o[j];
    }
  /* objective function */
    f1 = 0.;
    g1 = 0.;
    g2 = 0.;
    g3 = 1.;

    for (j = 0; j < nx; j++){
        f1 = f1 - e[j] * sin(sqrt(fabs(e[j])));
        g1 = g1 + e[j] * e[j];
        g2 = g2 + sin((1/50.0) * PI * e[j]);
        g3 = g3 * cos(e[j]/sqrt(((float) (j+1))));
    }

    f[0] = f1/((float) nx);
    g[0] = -50.0 + (1/(100.0 * ((float) nx))) * g1;
    g[1] = (50.0/((float) nx)) * g2;
    g[2] = 75.0 - 50.0 * (g1/4000.0 - g3 + 1.0);

}

void C14 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh){ // x[-1000,1000]
    int j;
    float f1, g1, g2, g3, e[nx];
    float o[30] = {-31.718907007204272,	-39.536680684207184, -46.033718058035944,	-42.2004014684422,	-28.331307546159135,	-38.64403177375364,	-11.313371899853626,	-11.717383190039943,	-43.345049558717875, -31.46016185891229,	-35.57742732758397,	-45.49638850141341,	-4.177473725277878,	-26.974808661067083,	-46.30991533784743,	-45.997883193212814,	-29.479673271045964,	-4.336542960830036,	-43.66244285780764,	-22.43896852522004,	-25.89273808052249,	-24.221450510218993,	-30.3952886350567,	-31.170730638052895,	-9.859463575974534,	-16.727846507426452,	-44.35226340706524,	-33.10843069426064,	-7.175153678947718,	-4.601421202670486};
    for (j = 0; j < nx; j++){
        e[j]= x[j] - o[j];
    }

  /* objective function */
    f1 = 0.;
    g1 = 0.;
    g2 = 0.;
    g3 = 0.;

    for (j = 0; j < (nx - 1); j++){
	  f1 = f1 + (100 * pow(((e[j]+1) * (e[j]+1) - (e[j+1]+1)),2) + pow(((e[j]+1) - 1),2));
    }
    for (j = 0; j < nx; j++){
      g1 = g1 - e[j] * cos(sqrt(fabs(e[j])));
	  g2 = g2 + e[j] * cos(sqrt(fabs(e[j])));
	  g3 = g3 + e[j] * sin(sqrt(fabs(e[j])));
    }
    f[0] = f1;
    g[0] = g1 -((float) nx);
    g[1] = g2 -((float) nx);
    g[2] = g3 - 10 *((float) nx);
}

void avaliaFuncaoRestricao(Individuo populacao[],int tipoFuncao,int tamanhoPop,int *fe){
    int i;
    for(i=0;i<tamanhoPop;i++){ // parametros tam_x,tamFuncObj,tam_g,tam_h
        (*fe)++;
        if(tipoFuncao == 1){
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

void somaViolacoes(Individuo populacao[],int tamanhoPop){ // Coloca as violacoes g e h no vetor 'v' de violacoes
    int i,j;
    int indG,indH;
    for(i=0;i<tamanhoPop;i++){
        indG=0,indH=0; // Indice do vetor g e h
        for(j=0;j<NUM_RESTRICOES_G+NUM_RESTRICOES_H;j++){ // Percorrer vetor 'v' de violacao
            if(j < NUM_RESTRICOES_G){
                populacao[i].v[j] = populacao[i].g[indG];
                indG++;
            }
            else{
                populacao[i].v[j] = fabs(populacao[i].h[indH]) - EPSILON; // Transforma h(x)=0 em g(x) <= 0
                indH++;
            }
        }
        populacao[i].violacao = somaValoresArray(populacao[i].v,NUM_RESTRICOES_G+NUM_RESTRICOES_H);
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

void algoritmoGeneticoRestricaoDEBUG(Individuo *melhorAG){
    int fe=0;
    int numIteracoes=0;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(SEED);
    //inicializaPopulacao(populacao,TAM_X);
    //avaliaFuncao(populacao,&contaObj,TAM_X,TAM_POPULACAO);
    inicializaPopulacaoRestricao(populacao,TAM_X,TAM_POPULACAO,TIPO_FUNCAO);
    avaliaFuncaoRestricao(populacao,TIPO_FUNCAO,TAM_POPULACAO,&fe);
    somaViolacoes(populacao,TAM_POPULACAO);
    printf("POP\n");
    imprimePopulacao(populacao,TAM_X,TAM_POPULACAO);
    //printf("POP\n");
    while(numIteracoes < MAX_CALC_OBJ){
        printf("POP\n");
        imprimePopulacao(populacao,TAM_X,TAM_POPULACAO);
        //selecao(populacao,filhos,TAM_X);
        selecaoRestricao(populacao,filhos,TAM_X);
        if(probabilidadeCrossover(PROB_CROSSOVER) == 1){
            if(TIPO_CROSSOVER == 0){
                crossover(filhos,TAM_X);
            }
            else if(TIPO_CROSSOVER == 1){
                crossoverSBX(populacao,filhos,TAM_X);
            }
        }
        mutacao(filhos,TAM_X);
        //avaliaFuncao(filhos,&contaObj,TAM_X,TAM_POPULACAO_FILHOS);
        corrigeLimitesX(filhos,TAM_X,TIPO_FUNCAO,TAM_POPULACAO_FILHOS);
        avaliaFuncaoRestricao(filhos,TIPO_FUNCAO,TAM_POPULACAO_FILHOS,&fe);
        somaViolacoes(filhos,TAM_POPULACAO_FILHOS);
        printf("FILHOS\n");
        imprimePopulacao(filhos,TAM_X,TAM_POPULACAO_FILHOS);
        ordenaMelhoresRestricao(populacao,filhos);
        printf("POP APOS ORDENAR\n");
        imprimePopulacao(populacao,TAM_X,TAM_POPULACAO);
        printf("FILHOS APOS ORDENAR\n");
        imprimePopulacao(filhos,TAM_X,TAM_POPULACAO_FILHOS);
        //exit(1);
        //printf("pos elitismo\n");
        //elitismo(filhos,populacao,TAM_X);
        elitismoRestricao(filhos,populacao,TAM_X,TAM_POPULACAO);
        printf("POP FINAL, FILHOS E PAIS:\n");
        imprimePopulacao(populacao,TAM_X,TAM_POPULACAO);
        //exit(1);
        // Pegar melhor individuo
        ordenaMelhoresRestricao(populacao,filhos);
        printf("POP FINAL, FILHOS E PAIS ORDENADOS\n");
        imprimePopulacao(populacao,TAM_X,TAM_POPULACAO);
        printf("Iteracao AG: %i\n",numIteracoes);
        imprimeInformacoesIndividuo(populacao[0]);
        numIteracoes++;
    }
    (*melhorAG) = melhorIndividuoRestricao(populacao,TAM_POPULACAO);
}

void AGRestricao(Individuo *melhorAG,int seed){
    int fe=0; // Function Evaluation
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(seed);
    //inicializaPopulacao(populacao,TAM_X);
    //avaliaFuncao(populacao,&contaObj,TAM_X,TAM_POPULACAO);
    inicializaPopulacaoRestricao(populacao,TAM_X,TAM_POPULACAO,TIPO_FUNCAO);
    avaliaFuncaoRestricao(populacao,TIPO_FUNCAO,TAM_POPULACAO,&fe);
    somaViolacoes(populacao,TAM_POPULACAO);
    while(fe < MAX_CALC_OBJ){
        selecaoRestricao(populacao,filhos,TAM_X);
        if(probabilidadeCrossover(PROB_CROSSOVER) == 1){
            if(TIPO_CROSSOVER == 0){
                crossover(filhos,TAM_X);
            }
            else if(TIPO_CROSSOVER == 1){
                crossoverSBX(populacao,filhos,TAM_X);
            }
        }
        mutacao(filhos,TAM_X);
        corrigeLimitesX(filhos,TAM_X,TIPO_FUNCAO,TAM_POPULACAO_FILHOS);
        avaliaFuncaoRestricao(filhos,TIPO_FUNCAO,TAM_POPULACAO_FILHOS,&fe);
        somaViolacoes(filhos,TAM_POPULACAO_FILHOS);
        ordenaMelhoresRestricao(populacao,filhos);
        elitismoRestricao(filhos,populacao,TAM_X,TAM_POPULACAO);
        // Pegar melhor individuo
        ordenaMelhoresRestricao(populacao,filhos);
        //printf("FE: %i\n",fe);
        //imprimeInformacoesIndividuo(populacao[0]);
    }
    (*melhorAG) = melhorIndividuoRestricao(populacao,TAM_POPULACAO);
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

void DERestricao(Individuo *melhorDE,int seed){
    int fe=0;
    int a,i,j;
    int iteracoes=0;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(seed);
    inicializaPopulacaoRestricao(populacao,TAM_X,TAM_POPULACAO,TIPO_FUNCAO);
    avaliaFuncaoRestricao(populacao,TIPO_FUNCAO,TAM_POPULACAO,&fe);
    somaViolacoes(populacao,TAM_POPULACAO);
    while(fe < MAX_CALC_OBJ){
    //while(iteracoes< MAX_CALC_OBJ){
        selecaoRestricao(populacao,filhos,TAM_X);
        int ch[3] = {-1,-1,-1}; // Vetor de indices, Poderia ser de tamanho 3?
        for(i=0;i<TAM_POPULACAO;i++){
            for(a=0;a<3;++a){ // Preenche vetor de indices
                ch[a] = selecionaPopulacaoDE(i,ch,TAM_POPULACAO);
            }
            int R = rand() % TAM_X; // Indice aleatorio, baseado na dimensao do problema
            for(j=0;j<TAM_X;j++){ // Altera valor de x
                float Ri = (rand() / (float)(RAND_MAX));
                if(Ri < CR || j == R){
                    filhos[i].x[j] = populacao[ch[0]].x[j] + F * (populacao[ch[1]].x[j] - populacao[ch[2]].x[j]);
                }
                else{
                    filhos[i].x[j] = populacao[i].x[j];
                }
            }
        }
        corrigeLimitesX(filhos,TAM_X,TIPO_FUNCAO,TAM_POPULACAO_FILHOS);
        avaliaFuncaoRestricao(filhos,TIPO_FUNCAO,TAM_POPULACAO_FILHOS,&fe);
        somaViolacoes(filhos,TAM_POPULACAO_FILHOS);
        selecaoDE(populacao,filhos,TAM_POPULACAO_FILHOS);
        //imprimeIndividuo(melhorIndividuoDE(populacao,TAM_POPULACAO),TAM_X);
        //printf("FE: %i\n",fe);
        iteracoes++;
        //printf("Iteracao: %i\n",iteracoes);
        //imprimeInformacoesIndividuo(melhorIndividuoRestricao(populacao,TAM_POPULACAO));
    }
    (*melhorDE) = melhorIndividuoRestricao(populacao,TAM_POPULACAO);
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

void ESRestricao(Individuo *melhorES,int tipoES,int sigmaGlobal,int seed){
    int i=0,j=0;
    int fe=0;
    //int numIteracoes=0;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(seed);
    inicializaPopulacaoRestricao(populacao,TAM_X,TAM_POPULACAO,TIPO_FUNCAO);
    corrigeLimitesX(populacao,TAM_X,TIPO_FUNCAO,TAM_POPULACAO); // No caso de C01 x[-10,10]
    // Avalia funcao
    avaliaFuncaoRestricao(populacao,TIPO_FUNCAO,TAM_POPULACAO,&fe);
    somaViolacoes(populacao,TAM_POPULACAO); // Preenche vetor 'v' de violacoes e seta variavel violacao, que eh a soma das violacoes
    inicalizaEstrategiaEvolutiva(populacao,filhos,TAM_X,sigmaGlobal);
    for(i=0;i<TAM_POPULACAO_FILHOS;i++,j++){ // Copia populacao para filhos
        filhos[i]=populacao[j];
        if(j >=TAM_POPULACAO){ // Acessaria lixo
            j=0; // reinicializa j, e continua copiando.
        }
    }
    //imprimeInformacoesPopulacao(filhos,TAM_POPULACAO_FILHOS);
    while(fe < MAX_CALC_OBJ){
    //while(numIteracoes <MAX_CALC_OBJ){
        autoAdaptacaoSigma(populacao,filhos,TAM_X,sigmaGlobal);
        // Avalia funcao
        corrigeLimitesX(filhos,TAM_X,TIPO_FUNCAO,TAM_POPULACAO_FILHOS);
        avaliaFuncaoRestricao(filhos,TIPO_FUNCAO,TAM_POPULACAO_FILHOS,&fe);
        somaViolacoes(filhos,TAM_POPULACAO_FILHOS);
        ordenaMelhoresRestricao(populacao,filhos);
        //ordenaMelhores(populacao,filhos);// NOTE e necessario?
        selecionaMelhoresRestricao(populacao,filhos,tipoES);
        //selecionaMelhores(populacao,filhos,TIPO_ES);
        // Pegar melhor individuo
        //imprimeContaObj(contaObj);
        //imprimeIndividuo(populacao[0],TAM_X);
        //printf("FE: %i\n",fe);
        //printf("iteracao: %i\n",numIteracoes);
        //printf("O melhor individuo: %f\n",populacao[0].funcaoObjetivo[0]);
        //imprimeInformacoesIndividuo(populacao[0]);
        //imprimeMelhores(populacao,filhos,TAM_X);
    }
    (*melhorES) = populacao[0];
    //qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaFuncaoObjetivo0);
    //imprimeInformacoesIndividuo(populacao[0]);
}

int main()
{
        //AG N 'X' POP 'Y' PROB CROSSOVER 'Z'
    //geraDados();
    //algoritmoGenetico(); // 0.001934sbx 0.000280normal
    int i,j,k;
    /*
    /// AG Restricao
    for(i=1;i<6;i++){
        Individuo melhorAG;
        AGRestricao(&melhorAG,i);
        printf("Melhor individuo do AG | seed: %i \n",i);
        //imprimeInformacoesIndividuo(melhorAG);
        imprimeVioleFO(melhorAG);
    }
    /// DE Restricao
    for(i=1;i<6;i++){
        Individuo melhorDE;
        DERestricao(&melhorDE,i);
        printf("Melhor individuo do DE | seed: %i \n",i);
        //imprimeInformacoesIndividuo(melhorDE);
        imprimeVioleFO(melhorDE);
    }

    /// ES Restricao

    for(i=0;i<1;i++){ // tipoES
        for(j=0;j<2;j++){ // SIGMA
            for(k=1;k<6;k++){ // SEED
                Individuo melhorES;
                ESRestricao(&melhorES,i,j,k);
                if(i == 0){ // ES +
                    if(j == 1){
                        printf("Melhor individuo do ES + Sigma Global | seed: %i \n",k);
                    }
                    else{
                        printf("Melhor individuo do ES + Sigma Individual | seed: %i \n",k);
                    }
                }
                if(i == 1){ // ES ,
                    if(j == 1){
                        printf("Melhor individuo do ES , Sigma Global | seed: %i \n",k);
                    }
                    else{
                        printf("Melhor individuo do ES , Sigma Individual | seed: %i \n",k);
                    }
                }
                //imprimeInformacoesIndividuo(melhorDE);
                imprimeVioleFO(melhorES);
            }
        }
    }
    */
    /// PARA SIGMA , GERAR MAIS FILHOS
    for(i=1;i<2;i++){ // tipoES
        for(j=0;j<2;j++){ // SIGMA
            for(k=1;k<6;k++){ // SEED
                Individuo melhorES;
                ESRestricao(&melhorES,i,j,k);
                if(i == 0){ // ES +
                    if(j == 1){
                        printf("Melhor individuo do ES + Sigma Global | seed: %i \n",k);
                    }
                    else{
                        printf("Melhor individuo do ES + Sigma Individual | seed: %i \n",k);
                    }
                }
                if(i == 1){ // ES ,
                    if(j == 1){
                        printf("Melhor individuo do ES , Sigma Global | seed: %i \n",k);
                    }
                    else{
                        printf("Melhor individuo do ES , Sigma Individual | seed: %i \n",k);
                    }
                }
                //imprimeInformacoesIndividuo(melhorDE);
                imprimeVioleFO(melhorES);
            }
        }
    }
    /*
    Individuo melhorES;
    printf("teste\n");
    ESRestricao(&melhorES,TIPO_ES,SIGMA_GLOBAL,SEED);
    imprimeVioleFO(melhorES);
    */
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
