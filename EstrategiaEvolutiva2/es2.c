/*
Titulo: Computação Evolucionista Aplicada à Engenharia
Autor: Pedro Henrique Santos
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define TAM_POPULACAO 10
#define NUM_FILHOS_GERADOS 30
#define TAM_POPULACAO_FILHOS TAM_POPULACAO*NUM_FILHOS_GERADOS
#define MAX_CALC_OBJ 100000// 100000
#define TAM_X 10
#define TAM_X_MAXIMO 30
#define MEDIA 0
#define DESVIO_PADRAO 1
#define SEED 1 // BUG 21 E 30
#define ELITISMO 5
#define GERA_ARQUIVO 1
#define PROB_CROSSOVER 0// 100 = 100% , 80 = 80%
#define TIPO_CROSSOVER 1 // 0 - Simples , 1 - SBX
#define TIPO_PROBLEMA 2 // 0 - Rosenbrock , 1 - Holder, 2 - Rastrigin
//#define TAU 1/sqrt(TAM_X)
#define TIPO_ES 1 // Se 0: es+ Se 1:es, Se 2: es, adaptado
#define TIPO_FUNCAO 13 // Funcoes com restricao
#define NUM_RESTRICOES_G 3
#define NUM_RESTRICOES_H 0

typedef struct {
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
    for(i=0;i<tamanhoX;i++){ // Imrime no cmd
        printf("%f\t",ind.x[i]);
    }
    printf("%f\n",ind.funcaoObjetivo[0]);
}

void imprimePopulacao(Individuo pop[],int tamanhoX){
    int i=0;
    for(i=0;i<TAM_POPULACAO;i++){
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

int comparaViolacaoRestricaoMinimizacao(const void *a, const void *b){
    if((*(Individuo*)a).violacao < (*(Individuo*)b).violacao){
        return -1; // vem antes
    }
    else if((*(Individuo*)a).violacao >(*(Individuo*)b).violacao){
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
            float val = eq * (rand() % 11 - 5); // Gera valores decimais entre [-5,5]
            //float val =eq * (rand () % 21 - 10); // Valores decimais entre [-10,10]
            populacao[i].x[j] = val;
        }
    }
}

void inicalizaEstrategiaEvolutiva(Individuo populacao[], Individuo filhos[],int tamanhoX){
    int i,j;
    for(i=0;i<TAM_POPULACAO;i++){// Para Estratégia Evolutiva (ES)
        for(j=0;j<tamanhoX;j++){
            populacao[i].sigma[j] = DESVIO_PADRAO;
        }
    }
    for(i=0;i<TAM_POPULACAO_FILHOS;i++){
        for(j=0;j<tamanhoX;j++){
            filhos[i].sigma[j]= DESVIO_PADRAO;
        }
    }
}

int individuoMelhorou(Individuo populacao,Individuo filho,int tipoProblema){ // Retorna 1 se filho melhor que pai
    //if(TIPO_PROBLEMA == 0 || TIPO_PROBLEMA == 2) { // Rosenbrock ou Rastrigin
    if(tipoProblema == 0){
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

void autoAdaptacaoSigma(Individuo populacao[],Individuo filhos[],int tamanhoX){
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
            for(k=0;k<tamanhoX;k++){
                filhos[l].sigma[k] = populacao[i].sigma[k] * exp(epsilon) * exp(epsilon); // NOTE: O quadrado eh necessario?
                filhos[l].x[k] = populacao[i].x[k] + filhos[l].sigma[k] * randNormal(0,1);
            }
            l++;
        }
    }
}

void ordenaMelhores(Individuo populacao[],Individuo filhos[]){
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

void selecionaMelhores(Individuo populacao[],Individuo filhos[],int tipoES){ // "Elitismo"
    if(tipoES == 0){ // Es +
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
        qsort(aux,TAM_POPULACAO+TAM_POPULACAO_FILHOS,sizeof(Individuo),comparaFuncaoObjetivo0); // Ordena auxiliar
        for(i=0;i<TAM_POPULACAO;i++){
            populacao[i]=aux[i]; // Copia os melhores de aux para populacao
        }
    }
    else if(tipoES == 1){ // Es ,
        int i,j=0;
        for(i=0;i<TAM_POPULACAO;i++){
            Individuo melhor = filhos[j];
            while(j<NUM_FILHOS_GERADOS*(i+1)){ // Percorre os X filhos de cada pai
                if(filhos[j].funcaoObjetivo[0] < melhor.funcaoObjetivo[0]){ // Seleciona o melhor filho
                    melhor = filhos[j];
                }
                j++;
            }
            populacao[i] = melhor;
        }
        qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaFuncaoObjetivo0);
    }
    else{ // Es , Adaptado
        int i,j=0;
        for(i=0;i<TAM_POPULACAO;i++){
            Individuo melhor = filhos[j];
            while(j<NUM_FILHOS_GERADOS*(i+1)){ // Pega o melhor filho de cada pai
                if(filhos[j].funcaoObjetivo[0] < melhor.funcaoObjetivo[0]){
                    melhor = filhos[j];
                }
                j++;
            }
            if(melhor.funcaoObjetivo[0] < populacao[i].funcaoObjetivo[0]){ // Filho melhor que pai
                populacao[i]=melhor; // Filho eh colocado no lugar do pai
            }
        }
        qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaFuncaoObjetivo0);
    }
}

void estrategiaEvolutivaGera(int seed,int tamanhoX,int tipoEs){
    int contaObj=0,i,j=0;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(seed);
    inicializaPopulacao(populacao,tamanhoX);
    avaliaFuncao(populacao,&contaObj,tamanhoX,TAM_POPULACAO);
    inicalizaEstrategiaEvolutiva(populacao,filhos,tamanhoX);
    for(i=0;i<TAM_POPULACAO_FILHOS;i++,j++){ // Copia populacao para filhos
        filhos[i]=populacao[j];
        if(j >= TAM_POPULACAO){
            j=0;
        }
    }
    while(contaObj < MAX_CALC_OBJ){
        autoAdaptacaoSigma(populacao,filhos,tamanhoX);
        avaliaFuncao(filhos,&contaObj,tamanhoX,TAM_POPULACAO_FILHOS);
        ordenaMelhores(populacao,filhos);
        selecionaMelhores(populacao,filhos,tipoEs);
        imprimeContaObj(contaObj);
        imprimeIndividuo(populacao[0],tamanhoX);
    }
}

void estrategiaEvolutiva(){
    int contaObj=0,i,j=0;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(SEED);
    inicializaPopulacao(populacao,TAM_X);
    //inicializaPopulacao(filhos,TAM_X);
    avaliaFuncao(populacao,&contaObj,TAM_X,TAM_POPULACAO);
    inicalizaEstrategiaEvolutiva(populacao,filhos,TAM_X);
    for(i=0;i<TAM_POPULACAO_FILHOS;i++,j++){ // Copia populacao para filhos
        filhos[i]=populacao[j];
        if(j >= TAM_POPULACAO){
            j=0;
        }
    }
    while(contaObj < MAX_CALC_OBJ){
        autoAdaptacaoSigma(populacao,filhos,TAM_X);
        avaliaFuncao(filhos,&contaObj,TAM_X,TAM_POPULACAO_FILHOS);
        ordenaMelhores(populacao,filhos);// NOTE e necessario?
        selecionaMelhores(populacao,filhos,TIPO_ES);
        // Pegar melhor individuo
        imprimeContaObj(contaObj);
        imprimeIndividuo(populacao[0],TAM_X);
    }
}

void estrategiaEvolutivaGeraDados(int seed,int tamanhoX,int tipoES){
    int contaObj=0,i,contador,j=0;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(seed);
    inicializaPopulacao(populacao,tamanhoX);
    avaliaFuncao(populacao,&contaObj,tamanhoX,TAM_POPULACAO);
    inicalizaEstrategiaEvolutiva(populacao,filhos,tamanhoX);
    for(i=0;i<TAM_POPULACAO_FILHOS;i++,j++){ // Copia populacao para filhos
        filhos[i]=populacao[j];
        if(j >= TAM_POPULACAO){
            j=0;
        }
    }
    while(contaObj < MAX_CALC_OBJ){
        contador = TAM_X_MAXIMO - tamanhoX;
        autoAdaptacaoSigma(populacao,filhos,tamanhoX);
        avaliaFuncao(filhos,&contaObj,tamanhoX,TAM_POPULACAO_FILHOS);
        ordenaMelhores(populacao,filhos);// NOTE e necessario?
        selecionaMelhores(populacao,filhos,tipoES);
        if(tipoES ==0){ // ; 2 100
            printf("+ %i %i\t",tamanhoX,TAM_POPULACAO);
        }
        else{
            printf(", %i %i\t",tamanhoX,TAM_POPULACAO);
        }
        imprimeContaObj(contaObj);
        while(contador>0){
            //printf("1.234567\t");
            printf("---------\t");
            contador--;
        }
        // Pegar melhor individuo
        imprimeIndividuo(populacao[0],tamanhoX);
    }
}

void geraDados(){
    int s,n;
    for(n=2;n<=11;n+=3){ // TamanhoProblema
        for(s=1;s<31;s++){// ES +
            estrategiaEvolutivaGeraDados(s,n,0);
        }
        for(s=1;s<31;s++){ // ES ,
            estrategiaEvolutivaGeraDados(s,n,1);
        }
    }
}

int comparaViolacao(const void *a, const void *b){ // Ordena em ordem crescentes vilacao
    if((*(Individuo*)a).violacao == (*(Individuo*)b).violacao)
        return 0;
    else{
        if((*(Individuo*)a).violacao < (*(Individuo*)b).violacao){
            return -1; // vem antes
        }
        else
            return 1; // vem depois
    }
}

void inicializaPopoulacaoRestricao(Individuo populacao[],int tamanhoX,int tamanhoPop,int tipoFuncao){
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
    qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaViolacao);
    qsort(filhos,TAM_POPULACAO_FILHOS,sizeof(Individuo),comparaViolacao);
}

void selecionaMelhoresRestricao(Individuo populacao[],Individuo filhos[],int tipoES){ // "Elitismo"
    if(tipoES == 0){ // Es +
        int i,k=0;
        int numViolacao=0;
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
                }
                // Se nao, melhor ja eh melhor
                j++;
            }
            populacao[i] = melhor;
        }
        qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaViolacaoRestricaoMinimizacao);
    }
    ///Possivel if para substituir IF porco do Es ,
        /*
                if(filhos[j].violacao < melhor.violacao){ // Seleciona o melhor filho
                    melhor = filhos[j];
                }
                else if(filhos[j].violacao = melhor.violacao){
                    if(filhos[j].funcaoObjetivo[0] < melhor.funcaoObjetivo[0]){ // Violacao eh igual, confere funcObj
                        melhor = filhos[j];
                    }
                    // Se nao, melhor ja eh melhor
                }*/
    else if(tipoES == 1){ // Es ,
        int i,j=0;
        for(i=0;i<TAM_POPULACAO;i++){
            Individuo melhor = filhos[j];
            while(j<NUM_FILHOS_GERADOS*(i+1)){ // Percorre os X filhos de cada pai
                if(filhos[j].funcaoObjetivo[0] < melhor.funcaoObjetivo[0]){ // Seleciona o melhor filho
                    melhor = filhos[j];
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

void C02 (float *x, float *f, float *g, float *h, int nx, int nf, int ng, int nh){
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
	  g3 = g3 + e[j] * sin(sqrt(fabs(e[j])));;
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

void estrategiaEvolutivaRestricao(){
    int contaObj=0,i=0,j=0;
    int numIteracoes=0;
    Individuo populacao[TAM_POPULACAO];
    Individuo filhos[TAM_POPULACAO_FILHOS];
    srand(SEED);
    inicializaPopoulacaoRestricao(populacao,TAM_X,TAM_POPULACAO,TIPO_FUNCAO);
    corrigeLimitesX(populacao,TAM_X,TIPO_FUNCAO,TAM_POPULACAO); // No caso de C01 x[-10,10]
    // Avalia funcao
    avaliaFuncaoRestricao(populacao,TIPO_FUNCAO,TAM_POPULACAO);
    somaViolacoes(populacao,TAM_POPULACAO); // Preenche vetor 'v' de violacoes e seta variavel violacao, que eh a soma das violacoes
    inicalizaEstrategiaEvolutiva(populacao,filhos,TAM_X);
    for(i=0;i<TAM_POPULACAO_FILHOS;i++,j++){ // Copia populacao para filhos
        filhos[i]=populacao[j];
        if(j >=TAM_POPULACAO){ // Acessaria lixo
            j=0; // reinicializa j, e continua copiando.
        }
    }
    while(numIteracoes < MAX_CALC_OBJ){
        autoAdaptacaoSigma(populacao,filhos,TAM_X);
        // Avalia funcao
        corrigeLimitesX(filhos,TAM_X,TIPO_FUNCAO,TAM_POPULACAO_FILHOS);
        avaliaFuncaoRestricao(filhos,TIPO_FUNCAO,TAM_POPULACAO_FILHOS);
        somaViolacoes(filhos,TAM_POPULACAO_FILHOS);
        ordenaMelhoresRestricao(populacao,filhos);
        //ordenaMelhores(populacao,filhos);// NOTE e necessario?
        selecionaMelhoresRestricao(populacao,filhos,TIPO_ES);
        //selecionaMelhores(populacao,filhos,TIPO_ES);
        // Pegar melhor individuo
        //imprimeContaObj(contaObj);
        //imprimeIndividuo(populacao[0],TAM_X);
        printf("Iteracao %i\n",numIteracoes);
        //printf("O melhor individuo: %f\n",populacao[0].funcaoObjetivo[0]);
        imprimeInformacoesIndividuo(populacao[0]); /// NOTE: VERIFICAR IMPRESSAO DO MELHOR
        //imprimeMelhores(populacao,filhos,TAM_X);
        numIteracoes++;
    }
    //qsort(populacao,TAM_POPULACAO,sizeof(Individuo),comparaFuncaoObjetivo0);
    //imprimeInformacoesIndividuo(populacao[0]);
}

int main(){
    //geraDados();
    //estrategiaEvolutiva();
    estrategiaEvolutivaRestricao();
    //estrategiaEvolutivaGeraDados(1,2,1);
    return 0;
}
