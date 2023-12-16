# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>


typedef struct matriceS{
    double **matrix;
    int32_t size;
}matrice_t;
matrice_t* prodmat(matrice_t* A, matrice_t* B);
matrice_t* transpose(matrice_t* A);
matrice_t* initmat(int32_t size);
void freemat(matrice_t* A);
matrice_t* inverser1(matrice_t* A);
matrice_t* inverser2(matrice_t* AtA, matrice_t* At);
matrice_t* partitionmat(matrice_t *AtA, uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2);
matrice_t* lirefichiermatrice(FILE* fichier);
matrice_t* inverser(matrice_t* mat);
matrice_t* submat(matrice_t* A,matrice_t*B);
matrice_t* addmat(matrice_t *A,matrice_t* B);
matrice_t* recoller(matrice_t* A, matrice_t* B, matrice_t* C, matrice_t* D);
matrice_t* negativ(matrice_t* A);
void printmat(matrice_t *A);
matrice_t* strassen(matrice_t* A, matrice_t* B);

int main( int argc , const char * argv []) {
    FILE* f=fopen("matrice.txt","r");
    if(f==NULL){
        printf("mauvais fichier");
        return 1;
    }
    
    matrice_t* matrice = lirefichiermatrice(f);
    fclose(f);
    if(matrice==NULL){
        printf("erreur lecture matrice\n");
        return EXIT_FAILURE;
    }
    matrice_t* res=initmat(matrice->size);
    res=inverser(matrice);
    printmat(res);

    freemat(matrice);
    freemat(res);

    
    return EXIT_SUCCESS;
}

matrice_t *prodmat(matrice_t *A, matrice_t *B)
{
    uint32_t i=0, j=0,k=0;
    matrice_t* res=initmat(A->size);
    for(i=0; i<A->size;i++){
        for(j=0;j<A->size;j++){
            res->matrix[i][j]=0;
            for(k=0;k<A->size;k++){
                res->matrix[i][j]+= A->matrix[i][k]*B->matrix[k][j];
            }
        }
    }
    return res;
}

matrice_t *transpose(matrice_t *A) {
    matrice_t* At = initmat(A->size); 
    uint32_t i, j;

    for (i = 0; i < A->size; i++) {
        for (j = 0; j < A->size; j++) {
            At->matrix[i][j] = A->matrix[j][i]; 
        }
    }

    return At;
}


matrice_t *initmat(int32_t size) {
    matrice_t* A = malloc(sizeof(matrice_t)); 
    if (A == NULL) {
        printf("alloc flop\n");
        return NULL;
    }

    
    A->matrix =(double**)malloc(size * sizeof(double *)); 
    if (A->matrix == NULL) {
        printf("ALLOC FLOP\n");
        free(A); 
        return NULL;
    }

    for (int i = 0; i < size; i++) {
        A->matrix[i] = (double*)malloc(size * sizeof(double)); 
    }
    A->size = size;
    return A;
}


void freemat(matrice_t* A)
{ 
    uint32_t k=0;
    for(k=0;k<A->size;k++){
        free(A->matrix[k]); 
    }
    free(A->matrix);
    free(A);
    return;
}


matrice_t* inverser1(matrice_t *A) {
    matrice_t *At = transpose(A); 
    matrice_t *AtA = prodmat(At, A); 

    
    freemat(At);


    return AtA; 
}

matrice_t* inverser2(matrice_t *AtA, matrice_t* At)
{   
    if(AtA->size==1){
        AtA->matrix[0][0]=(1/(AtA->matrix[0][0]))*At->matrix[0][0];
        return AtA;
    }
    uint32_t taille = AtA->size;
    uint32_t demi_taille = taille / 2;

    matrice_t* B= partitionmat(AtA,0,0,demi_taille-1,demi_taille-1);
    matrice_t* Ct= partitionmat(AtA,0,demi_taille,demi_taille-1,taille-1);
    matrice_t* C=transpose(Ct);
    matrice_t* D= partitionmat(AtA,demi_taille,demi_taille,taille-1,taille-1);
    matrice_t* Bi=inverser(B);
    matrice_t* CBi= prodmat(C,Bi);  
    matrice_t* BiCt=prodmat(Bi,Ct);
    matrice_t* CBiCt=prodmat(CBi,Ct);
    matrice_t* S=submat(D,CBiCt);
    matrice_t* Si=inverser(S);
    matrice_t* SiCBi=prodmat(Si,CBi); 
    matrice_t* BiCtSiCBi=prodmat(BiCt,SiCBi);
    matrice_t* BiPBiCtSiCBi=addmat(Bi,BiCtSiCBi);
    matrice_t* BiCtSi=prodmat(BiCt,Si);
    matrice_t* mBiCtSi=negativ(BiCtSi);
    matrice_t* mSiCBi=negativ(SiCBi);
    matrice_t* res= recoller(BiPBiCtSiCBi,mBiCtSi,mSiCBi,Si);
    res= prodmat(res,At);
    freemat(B);
    freemat(C);
    freemat(Ct);
    freemat(D);
    freemat(Bi);
    freemat(CBi);
    freemat(BiCt);
    freemat(CBiCt);
    freemat(S);
    freemat(Si);
    freemat(SiCBi);
    freemat(BiCtSiCBi);
    freemat(BiPBiCtSiCBi);
    freemat(BiCtSi);
    freemat(mBiCtSi);
    freemat(mSiCBi);

    return res;
    

}



matrice_t *partitionmat(matrice_t *AtA, uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2) {
    matrice_t *res = initmat(x2 - x1 + 1);
    uint32_t i, j, a = 0, b = 0;

    for (i = x1; i <= x2; i++) {
        b = 0;
        for (j = y1; j <= y2; j++) {
            res->matrix[a][b] = AtA->matrix[i][j];
            b++;
        }
        a++;
    }

    return res;
}


matrice_t *inverser(matrice_t *mat)
{
    matrice_t* At=transpose(mat);
    
    mat=inverser1(mat);
    mat=inverser2(mat, At);
    freemat(At);
    return mat;
}

matrice_t *submat(matrice_t *A,matrice_t* B)
{ 
    matrice_t* res=initmat(A->size);
    uint32_t i,j;
    for(i=0;i<A->size;i++){
        for(j=0;j<A->size;j++){
            res->matrix[i][j]=A->matrix[i][j]-B->matrix[i][j];
        }        
    }
    return res;
}
matrice_t *addmat(matrice_t *A,matrice_t* B)
{ 
    matrice_t* res=initmat(A->size);
    uint32_t i,j;
    for(i=0;i<A->size;i++){
        for(j=0;j<A->size;j++){
            res->matrix[i][j]=A->matrix[i][j]+B->matrix[i][j];
        }        
    }
    return res;
}


matrice_t *recoller(matrice_t *A, matrice_t *B, matrice_t *C, matrice_t *D) {
    matrice_t *res = initmat((A->size) * 2);
    uint32_t i, j;

    for (i = 0; i < A->size; i++) {
        for (j = 0; j < A->size; j++) {
            res->matrix[i][j] = A->matrix[i][j];
            res->matrix[i][j + A->size] = B->matrix[i][j];
            res->matrix[i + A->size][j] = C->matrix[i][j];
            res->matrix[i + A->size][j + A->size] = D->matrix[i][j];
        }
    }


    return res;
}


matrice_t *negativ(matrice_t *A) {
    matrice_t *res = initmat(A->size);
    uint32_t i, j;

    for (i = 0; i < A->size; i++) {
        for (j = 0; j < A->size; j++) {
            res->matrix[i][j] = -A->matrix[i][j];
        }
    }

    return res;
}

matrice_t *lirefichiermatrice(FILE *fichier)
{
    int32_t size;
    fscanf(fichier,"%d",&size);
    matrice_t* res= initmat(size);
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            fscanf(fichier, "%lf", &res->matrix[i][j]);
        }
    }
    return res;
}

void printmat(matrice_t *A)
{
    uint32_t i=0, j=0;
    for(i=0; i<A->size;i++){
        for(j=0;j<A->size;j++){
            printf("%lf  ",A->matrix[i][j]);
        }
        printf("\n");
    }
}

matrice_t *strassen(matrice_t *A, matrice_t *B)
{
    int32_t n= A->size;
    if(n==1){
        matrice_t* res=initmat(1);
        res->matrix[0][0]=A->matrix[0][0]*B->matrix[0][0];
        return res;
    }

    matrice_t* A11=partitionmat(A, 0, 0, (n/2)-1,(n/2)-1);
    matrice_t* A12=partitionmat(A, 0, (n/2)-1, (n/2)-1,(n-1));


}
