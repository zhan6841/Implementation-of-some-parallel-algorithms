#include <malloc.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#define  MAX(m,n)    (m>n?m:n)

typedef struct{
  	int pedlen;
  	int psuffixlen;
  	int pednum;
}pntype;


/*对模式串进行周期分析，并计算相应的new和newval值*/
void Next(char *W,int patlen,int *nextval,pntype *pped)
{
	int i,j,plen;
   	int *next;

   	if((next=(int *)malloc((patlen+1)*sizeof(int)))==NULL){
     	printf("no enough memory\n");
       	exit(1);
   	}

  	/*计算next和nextval*/    
    // set the next[0] as -1
    next[0]=nextval[0]=-1;
  	j=1;
   	while(j<=patlen){
   		i=next[j-1];
      	while(i!=(-1) && W[i]!=W[j-1]) i=next[i];
      	next[j]=i+1;

      	if(j!=patlen){
          	if( W[j]!=W[i+1])
              	nextval[j]=i+1;
           	else
    			nextval[j]=nextval[i+1]; 
      	}
        j++;
  	} 
    // get the nextval[patlen], when finish one match
    // it can regard as the next character of pattern W missmatch to continue next match
    pped->pedlen=patlen-next[patlen]; 
    // the times that can rematch
 	pped->pednum=(int)(patlen/pped->pedlen); 
    // length of the last several characters in the pattern W
   	pped->psuffixlen=patlen%pped->pedlen; 

  	free(next);
}

/*改进的KMP算法*/
// prefix_flag means if this time has to prefix for the next prcocessor
void kmp(char *T,char*W,int textlen,int patlen,int *nextval,pntype *pped,int prefix_flag,int matched_num,int *match,int *prefixlen)
{
	int i,j;

  	i=matched_num;              
  	j=matched_num;           

    while(i<textlen)
    {
      // if don't have to prefix and the left length of text is not enough
    	if((prefix_flag==1)&&((patlen-j)>(textlen-i))) {
          	break;
		}
        // when mismatch, shift back the pattern W
        while(j!=(-1) && W[j]!=T[i])  j=nextval[j];
        // if match successful, set match[i-(patlen-1)] as 1
        if(j==(patlen-1)) {    
			match[i-(patlen-1)]=1;
          // if the pattern can only rematch once
        	if(pped->pednum+pped->psuffixlen==1)
            j = -1 ; 
      // if there are other previous suffix can match, continue rematch
      // from the longest presuffix to the shortest one
   		else                                   
          	j=patlen-1-pped->pedlen; 
        }
   		j++;
      	i++; 
   	}
    // record the length of the prefix pattern
    // the index of the character in the pattern that need to match next for the next processor
   	(*prefixlen)=j;
}

/*重构模式串以及next函数*/
void Rebuild_info(int patlen,pntype *pped,int *nextval,char *W) 
{ 
	int i; 
    // if the pattern can only rematch once
   	if (pped->pednum == 1) 
   		memcpy(W+pped->pedlen,W,pped->psuffixlen); 
	else {  
       	memcpy(W+pped->pedlen,W,pped->pedlen);
       	for (i=3; i<=pped->pednum; i++){ 
        	memcpy(W+(i-1)*pped->pedlen,W,pped->pedlen);
           	memcpy(nextval+(i-1)*pped->pedlen,nextval+pped->pedlen,pped->pedlen*sizeof(int));
      	} 

       	if(pped->psuffixlen!=0){
       		memcpy(W+(i-1)*pped->pedlen,W,pped->psuffixlen);
           	memcpy(nextval+(i-1)*pped->pedlen,nextval+pped->pedlen,pped->psuffixlen*sizeof(int));            
       	}
 	} 
} 
 
/*生成文本串*/
void gen_string(int strlen,int pedlen,char *string,int seed)
{
	int suffixlen,num,i,j;

   	srand(seed);
   	for(i=0;i<pedlen;i++){
    	num=rand()%26;
        string[i]='a'+num;
   	}
   	for(j=1;j<(int)(strlen/pedlen);j++)
    	strncpy(string+j*pedlen,string,pedlen);
   	if((suffixlen=strlen%pedlen)!=0)
    	strncpy(string+j*pedlen,string,suffixlen);
}  

/*从文件读入模式串信息*/ 
void GetFile(char *filename,char **place, int *number) 
{ 
	FILE *fp;
    struct stat statbuf;  

    if ((fp=fopen(filename,"rb"))==NULL) {
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 

    fstat(fileno(fp),&statbuf);  
    if(((*place)=(char *)malloc(sizeof(char)*statbuf.st_size)) == NULL){
		printf ("Error alloc memory\n");
        exit(1);
 	}
     
   	if (fread((*place),1,statbuf.st_size,fp)!=statbuf.st_size){
		printf ("Error in reading num\n"); 
        exit(0); 
	} 
    fclose (fp); 
    (*number)=statbuf.st_size; 
} 

/*打印运行参数信息*/
// Modify by Zhanganlan
void PrintFile_info(char *filename,char *T,int id,int strlen)
{ 
	FILE *fp; 
	int i;

	if ((fp=fopen(filename,"a"))==NULL){
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 

 	// fprintf (fp,"The Text on node %d is %s .\n",id,T); 
  fprintf (fp,"The Text on node %d is ",id); 
  for(i = 0; i < strlen; i++){
    if(T[i] != '\0'){
      fprintf(fp, "%c", T[i]);
    }
  }
  fprintf(fp, " .\n");
   	
	fclose (fp); 
} 


/*打印匹配结果*/
void PrintFile_res(char *filename,int *t,int len,int init,int id)
{ 
	FILE *fp; 
	int i;

	if ((fp=fopen(filename,"a"))==NULL){
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 

 	fprintf (fp,"This is the match result on node %d \n",id); 
   	for (i=0; i<=len-1; i++) 
    	if(t[i]==1)
 			fprintf (fp,"(%d)  +\n",i+init); 
    	else
  			fprintf (fp,"(%d)  -\n",i+init);
	fclose (fp); 
} 

void main(int argc,char *argv[]) 
{ 
	char *T,*W; 
	int	*nextval,*match; 
  	int	textlen,patlen,pedlen,nextlen_send; 
   	pntype pped; 
 	int	i,myid,numprocs,prefixlen,ready; 
  	MPI_Status  status; 

	MPI_Init(&argc,&argv); 
   	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
   	MPI_Comm_rank(MPI_COMM_WORLD,&myid); 

   	nextlen_send=0;
   	ready=1;
   	/*读如文本串长度*/
   	textlen=atoi(argv[1]);
    // modify by Zhanganlan
    // consider the situation that texlen % numprocs != 0
    int totallen = textlen;
    if(textlen%numprocs){
      textlen = textlen/numprocs + 1;
    }
    else{
      textlen=textlen/numprocs;
    }

  	pedlen=atoi(argv[2]);
   	if((T=(char *)malloc(textlen*sizeof(char)))==NULL){
     	printf("no enough memory\n");
       	exit(1);
   	}
   	gen_string(textlen,pedlen,T,myid);
    // add by Zhanganlan
    if(myid == numprocs-1){
      if(totallen%numprocs){
        for(i = totallen - textlen*(numprocs-1); i < textlen; i++){
          T[i] = '\0';
        }
      }
    }

  // every processor send a ready signal to its next processor
	if(myid==0){
		PrintFile_info("match_result",T,myid,textlen);
		if(numprocs>1)
			MPI_Send(&ready,1,MPI_INT,1,0,MPI_COMM_WORLD);
	}
   	else{
  		MPI_Recv(&ready,1,MPI_INT,myid-1,myid-1,MPI_COMM_WORLD,&status);
		PrintFile_info("match_result",T,myid,textlen);
		if(myid!=numprocs-1)
			MPI_Send(&ready,1,MPI_INT,myid+1,myid,MPI_COMM_WORLD);
	}

	printf("\n");
   	
	if((match=(int *)malloc(textlen*sizeof(int)))==NULL){
		printf("no enough memory\n");
		exit(1);
	}
  
  	/*处理器0读入模式串，并记录运行参数*/
   	if(myid==0){ 
  		printf("processor num = %d \n",numprocs);
      // Modify by Zhanganlan
    	// printf("textlen = %d\n",textlen*numprocs); 
      printf("textlen = %d\n",totallen);  

        GetFile("pattern.dat",&W,&patlen); 
    	printf("patlen= %d\n",patlen); 

    	if((nextval=(int *)malloc(patlen*sizeof(int)))==NULL){
        	printf("no enough memory\n");
           	exit(1);
        }
		/*对模式串进行分析，对应于算法14.6步骤（1）*/

      	Next(W,patlen,nextval,&pped);
        // 
        if(numprocs>1){
        	if (pped.pednum==1) 
           		nextlen_send = patlen;
            else 
        		nextlen_send = pped.pedlen*2;
        }
    }

	/*向各个处理器播送模式串的信息，对应于算法14.6步骤（2）*/
  	if(numprocs>1){
     	MPI_Bcast(&patlen, 1, MPI_INT, 0, MPI_COMM_WORLD);  
  		if(myid!=0)
    		if(((nextval=(int *)malloc(patlen*sizeof(int)))==NULL)
				||((W=(char *)malloc(patlen*sizeof(char)))==NULL)){
           		printf("no enough memory\n");
            	exit(1);
            }

 	 	MPI_Barrier(MPI_COMM_WORLD);
    // this is to save the cost of communication?
    	MPI_Bcast(&pped,3,MPI_INT,0,MPI_COMM_WORLD);  
    	MPI_Bcast(&nextlen_send,1,MPI_INT,0,MPI_COMM_WORLD);
    	MPI_Bcast(nextval,nextlen_send,MPI_INT,0,MPI_COMM_WORLD); 
    	MPI_Bcast(W,pped.pedlen,MPI_CHAR,0,MPI_COMM_WORLD);
   	}

    MPI_Barrier(MPI_COMM_WORLD);

   	/*调用修改过的KMP算法进行局部串匹配，对应于算法14.6步骤（3）*/
  	if(numprocs==1) {
  		kmp(T,W,textlen,patlen,nextval,&pped,1,0,match+patlen-1,&prefixlen);
   	}
   	else { 
    	if(myid!=0)
    		/*各个处理器分别根据部分串数据以及周期信息重构模式串*/
        	Rebuild_info(patlen,&pped,nextval,W); 
      // the first numprocs-1 processors need to prefix the last patlen-1 text for the next processor
    	if(myid!=numprocs-1)
  			kmp(T,W,textlen,patlen,nextval,&pped,0,0,match+patlen-1,&prefixlen);
		else
  			kmp(T,W,textlen,patlen,nextval,&pped,1,0,match+patlen-1,&prefixlen);

   		MPI_Barrier(MPI_COMM_WORLD);

		/*各个处理器进行段间匹配，对应于算法14.6步骤（4）*/
    	if(myid<numprocs-1) 
        	MPI_Send(&prefixlen,1,MPI_INT,myid+1,99,MPI_COMM_WORLD); 

    	if(myid>0) 
    		MPI_Recv(&prefixlen,1,MPI_INT,myid-1,99,MPI_COMM_WORLD,&status); 

    	MPI_Barrier(MPI_COMM_WORLD);

      // the last numprocs-1 processors need to try to match the pattern between the adjcent two parts
    	if((myid>0)&&(prefixlen!=0))  
   			kmp(T-prefixlen,W,prefixlen+patlen-1,patlen,nextval,&pped,1,prefixlen,match+patlen-1-prefixlen,&prefixlen); 

   		MPI_Barrier(MPI_COMM_WORLD);
   	}

	/*输出匹配结果*/
  // MPI_Send & Recv make sure the output is in order
	if(myid==0){
		PrintFile_res("match_result",match+patlen-1,textlen-patlen+1,0,myid);
		if(numprocs>1)
			MPI_Send(&ready,1,MPI_INT,1,0,MPI_COMM_WORLD);
	}
   	else{
  		MPI_Recv(&ready,1,MPI_INT,myid-1,myid-1,MPI_COMM_WORLD,&status);
		PrintFile_res("match_result",match,textlen,myid*textlen-patlen+1,myid);
		if(myid!=numprocs-1)
			MPI_Send(&ready,1,MPI_INT,myid+1,myid,MPI_COMM_WORLD);
	}

	free(T);
    free(W);
    free(nextval);
    MPI_Finalize(); 
 } 
