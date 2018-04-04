/*
Titulo: Computação Evolucionista Aplicada à Engenharia
Autor: Pedro Henrique Santos
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define TAM_POPULACAO 50
#define NUM_FILHOS_GERADOS 1
#define TAM_POPULACAO_FILHOS TAM_POPULACAO*NUM_FILHOS_GERADOS
#define MAX_CALC_OBJ 15000// 100000 original
#define TAM_X 10
#define TAM_X_MAXIMO 30
#define MEDIA 0
#define DESVIO_PADRAO 1
#define SEED 1 //BUG SEED 1,2,3 com N EM BETA 2 PROBLEMA PRA HOLDER SBX (gerabeta) seed 7 bug pasa rosenbrock,com ES
#define ELITISMO 5
#define GERA_ARQUIVO 1
#define PROB_CROSSOVER 100// 100 = 100% , 80 = 80%
#define TIPO_CROSSOVER 1 // 0 - Simples , 1 - SBX
#define TIPO_PROBLEMA 0 // 0 - Rosenbrock , 1 - Holder, 2 - Rastrigin
#define N_EM_BETA 2 // Normalmente entre 2 e 5. Quanto maior , maior a chance de gerar filhos mais 'pr�ximo' do pai
#define N_EM_MUTACAOREAL 4 // Para mutação Real (slide SBX)
#define DELTA_MAX 10 // Para mutação Real (slide SBX)
#define TIPO_FUNCAO 1 // Funcoes com restricao 14 e 13 = 3g,0h | 01 = 2g,0h
#define NUM_RESTRICOES_G 2
#define NUM_RESTRICOES_H 0
#define CR 0.8 // Crossover Probability , DE [0,1]
#define F 0.7 // Differential Weight , DE [0,2]

typedef struct {
    float x[TAM_X_MAXIMO];
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
        printf("%f\t",ind.x[i]);
    }
    printf("FO:%f\n",ind.funcaoObjetivo[0]);
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
            else if(tipoFuncao == 13){
                val = eq * (rand() % 1001 - 500); // valores decimais entre [-500,500]
            }
            else{
                printf("Funcao nao encontrada\n");
                exit(0);
            }
            //float val = eq * (rand() % 11 - 5); // Gera valores decimais entre [-5,5]
            //float val =eq * (rand () % 21 - 10); // Valores decimais entre [-10,10]
            //float val = eq * (rand() % 2001 - 1000); // Valores decimais entre [-1000,1000] C03
            populacao[i].x[j] = val;
        }

        for(j=0;j<NUM_RESTRICOES_G;j++){
            populacao[i].g[j]=0;
        }
        for(j=0;j<NUM_RESTRICOES_H;j++){
            populacao[i].h[j]=0;
        }
        for(j=0;j<NUM_RESTRICOES_G+NUM_RESTRICOES_H;j++){
            populacao[i].v[j]=0;
        }
        populacao[i].violacao=0;
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

void avaliaFuncaoRestricao(Individuo populacao[],int tipoFuncao,int tamanhoPop){
    int i;
    for(i=0;i<tamanhoPop;i++){ // parametros tam_x,tamFuncObj,tam_g,tam_h
        if(tipoFuncao == 1){
            C01(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,2,0);
        }
        else if(tipoFuncao == 2){
            C02(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,2,1);
        }
        else if(tipoFuncao == 3){
            C03(populacao[i].x,populacao[i].funcaoObjetivo,populacao[i].g,populacao[i].h,TAM_X,1,0,1);
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
    else if(tipoFuncao == 13){
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

void somaViolacoes(Individuo populacao[],int tamanhoPop){ // Coloca as violacoes g e h no vetor 'v' de violacoes
    int i,j;
    int indG,indH;
    for(i=0;i<tamanhoPop;i++){
        indG=0,indH=0; // Indice do vetor g e h
        for(j=0;j<NUM_RESTRICOES_G+NUM_RESTRICOES_H;j++){ // Percorrer vetor 'v' de violacao
            if(j < NUM_RESTRICOES_G){
                //printf("O valor de g(x): %f",populacao[i].g[indG]);
                populacao[i].v[j] = populacao[i].g[indG];
                indG++;
                //printf("indG: %i",indG);
            }
            else{
                //printf("entrou else\n");
                populacao[i].v[j] = populacao[i].h[indH];
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
    printf("Os valores de funcaoObj: %f\n",ind.funcaoObjetivo[0]);
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

void algoritmoGeneticoRestricao(Individuo *melhorAG){
    //int contaObj=0;
    int numIteracoes=0;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(SEED);
    //inicializaPopulacao(populacao,TAM_X);
    //avaliaFuncao(populacao,&contaObj,TAM_X,TAM_POPULACAO);
    inicializaPopulacaoRestricao(populacao,TAM_X,TAM_POPULACAO,TIPO_FUNCAO);
    avaliaFuncaoRestricao(populacao,TIPO_FUNCAO,TAM_POPULACAO);
    somaViolacoes(populacao,TAM_POPULACAO);
    printf("POP\n");
    imprimePopulacao(populacao,TAM_X,TAM_POPULACAO);
    //printf("POP\n");
    //exit(1);
    while(numIteracoes < MAX_CALC_OBJ){
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
        avaliaFuncaoRestricao(filhos,TIPO_FUNCAO,TAM_POPULACAO_FILHOS);
        somaViolacoes(filhos,TAM_POPULACAO_FILHOS);
        //printf("FILHOS\n");
        //imprimePopulacao(filhos,TAM_X,TAM_POPULACAO_FILHOS);
        //ordenaMelhores(populacao,filhos,TIPO_PROBLEMA);
        ordenaMelhoresRestricao(populacao,filhos);
        /*printf("POP\n");
        imprimePopulacao(populacao,TAM_X,TAM_POPULACAO);
        printf("FILHOS\n");
        imprimePopulacao(filhos,TAM_X,TAM_POPULACAO_FILHOS);*/
        //exit(1);
        //printf("pos elitismo\n");
        //elitismo(filhos,populacao,TAM_X);
        elitismoRestricao(filhos,populacao,TAM_X,TAM_POPULACAO);
        //imprimePopulacao(populacao,TAM_X,TAM_POPULACAO);
        //exit(1);
        // Pegar melhor individuo
        ordenaMelhoresRestricao(populacao,filhos); // eh necessario ordenar de novo? *podem ser copiados filhos melhores que pais*
        printf("Iteracao AG: %i\n",numIteracoes);
        imprimeInformacoesIndividuo(populacao[0]);
        numIteracoes++;
    }
    (*melhorAG) = melhorIndividuoRestricao(populacao,TAM_POPULACAO);
}

int selecionaPopulacao(int solucao,int ch [],int tamanhoPop){
    int indice;
    int achou = 1;
    while(1){
        indice = rand() % tamanhoPop;
        if(indice != solucao){
            int a = 0;
            while(ch[a] != -1){
                if(indice == ch[a]){ // Já está na lista
                    achou = 0;
                    break;
                }
                ++a;
            }
            if(achou == 1){
                break;
            }
            else{
                achou = 1;
            }
        }
    }
    return indice;
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

void algoritmoDE(Individuo *melhorDE){
    int a,i,j,numIteracoes = 0;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(SEED);
    inicializaPopulacaoRestricao(populacao,TAM_X,TAM_POPULACAO,TIPO_FUNCAO);
    avaliaFuncaoRestricao(populacao,TIPO_FUNCAO,TAM_POPULACAO);
    somaViolacoes(populacao,TAM_POPULACAO);
    while(numIteracoes < MAX_CALC_OBJ){
        selecaoRestricao(populacao,filhos,TAM_X);
        int ch[3] = {-1,-1,-1}; // Vetor de indices, Poderia ser de tamanho 3?
        for(i=0;i<TAM_POPULACAO;i++){
            for(a=0;a<3;++a){ // Preenche vetor de indices
                //ch[a] = selecionaPopulacao(i,ch,TAM_POPULACAO);
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
        avaliaFuncaoRestricao(filhos,TIPO_FUNCAO,TAM_POPULACAO_FILHOS);
        somaViolacoes(filhos,TAM_POPULACAO_FILHOS);
        selecaoDE(populacao,filhos,TAM_POPULACAO_FILHOS);
        //imprimeIndividuo(melhorIndividuoDE(populacao,TAM_POPULACAO),TAM_X);
        imprimeInformacoesIndividuo(melhorIndividuoRestricao(populacao,TAM_POPULACAO));
        printf("Iteracao DE: %i\n",numIteracoes);
        numIteracoes++;
    }
    (*melhorDE) = melhorIndividuoRestricao(populacao,TAM_POPULACAO);
}

int main()
{
        //AG N 'X' POP 'Y' PROB CROSSOVER 'Z'
    //geraDados();
    //algoritmoGenetico(); // 0.001934sbx 0.000280normal

    /// AG Restricao
    /*
    Individuo melhorAG;
    algoritmoGeneticoRestricao(&melhorAG);
    printf("Melhor individuo do AG: \n");
    imprimeInformacoesIndividuo(melhorAG);
    */
    /// DE Restricao

    Individuo melhorDE;
    algoritmoDE(&melhorDE);
    printf("Melhor individuo do DE: \n");
    imprimeInformacoesIndividuo(melhorDE);


    /*algoritmoDE(&melhorDE);
    printf("Melhor individuo do AG: \n")
    imprimeInformacoesIndividuo(melhorAG);
    printf("UAIPORRA");
    printf("Melhor Individuo do DE: \n");
    imprimeInformacoesIndividuo(melhorDE);
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
