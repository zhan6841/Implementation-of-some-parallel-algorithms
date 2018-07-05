#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define INIT_TYPE 10
#define ALLTOONE_TYPE 100
#define ONETOALL_TYPE 200
#define MULTI_TYPE 300
#define RESULT_TYPE 400
#define RESULT_LEN 500
#define MULTI_LEN 600

int Spt;
long DataSize;
int *arr,*arr1;
int mylength;
int *index;
int *temp1;

main(int argc,char* argv[])
{
    long BaseNum = 1;
    int PlusNum;
    int MyID;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);

    // PlusNum=60;
    // for test convenient
    PlusNum=60;
    DataSize = BaseNum*PlusNum;

    // print the total length of the dataset
    if (MyID==0)
        printf("The DataSize is : %lu\n",PlusNum);
    // call psra algorithm
    Psrs_Main();

    if (MyID==0)
        printf("\n");

    MPI_Finalize();
}


Psrs_Main( )
{
    int i,j;
    int MyID,SumID;
    int n,c1,c2,c3,c4,k,l;
    FILE * fp;
    int ready;
    MPI_Status status[32*32*2];
    MPI_Request request[32*32*2];

    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    MPI_Comm_size(MPI_COMM_WORLD,&SumID);

    Spt = SumID-1;

	/*初始化参数*/
    // why double data size?
    arr = (int *)malloc(2*DataSize*sizeof(int));
    if (arr==0) merror("malloc memory for arr error!");
    arr1 = &arr[DataSize];
    // if there are more than one process
    if (SumID>1)
    {
        temp1 = (int *)malloc(sizeof(int)*SumID*Spt);
        if (temp1==0) merror("malloc memory for temp1 error!");
        index = (int *)malloc(sizeof(int)*2*SumID);
        if (index==0) merror("malloc memory for index error!");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // modify by Zhanganlan
    if(MyID == SumID - 1){
        // the last process may not have the same length of data
        mylength = DataSize - (DataSize / SumID) * (SumID - 1);
    }
    else{
        mylength = DataSize / SumID;
    }
    
    srand(MyID);
    // if(MyID==2){
    printf("This is node %d \n",MyID);
    printf("On node %d the input data is:\n",MyID);
    // }
    for (i=0;i<mylength;i++)
    {
        arr[i] = (int)rand() % 1000;
        // if(MyID==2)
        printf("%d : ",arr[i]);
    }
    // if(MyID==2)
    printf("\n");

	/*每个处理器将自己的n/P个数据用串行快速排序(Quicksort)，得到一个排好序的序列，对应于算法13.5步骤（1）*/
    MPI_Barrier(MPI_COMM_WORLD);
    // call quicksort to sort local data in every process
    quicksort(arr,0,mylength - 1);
    /*
    if(MyID==0){
    printf("This is node %d \n after quicksort",MyID);
    for (i=0;i<mylength;i++)
    {
        printf("%d : ",arr[i]);
    }
    printf("\n");
    }
    */
    MPI_Barrier(MPI_COMM_WORLD);

	/*每个处理器从排好序的序列中选取第w，2w，3w，…，(P-1)w个共P-1个数据作为代表元素，其中w=n/P*P，对应于算法13.5步骤（2）*/
    if (SumID>1)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        n = (int)(mylength/(Spt+1));
        for (i=0;i<Spt;i++)
            temp1[i] = arr[(i+1)*n-1];

        MPI_Barrier(MPI_COMM_WORLD);

        if (MyID==0)
        {
			/*每个处理器将选好的代表元素送到处理器P0中，对应于算法13.5步骤（3） */
            j = 0;
            for (i=1;i<SumID;i++)
                /*
                int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
                OUT request: request handle for future inquiry
                 */
                // why sizeof(int)*Spt? why MPI_CHAR?
                MPI_Irecv(&temp1[i*Spt], sizeof(int)*Spt, MPI_CHAR, i, ALLTOONE_TYPE+i, MPI_COMM_WORLD, &request[j++]);
            /*
            int MPI_Waitall(int Count, MPI_Request *array_of_requests, MPI_Status *array_of_status)
            INOUT array_of_requests: the array of all request handles
            INOUT array_of_status: the array of all the status of all requests
            wait until after the communication of all the requests are done
             */
            MPI_Waitall(SumID-1, request, status);
            /*
            printf("This is node %d \n after choose from every process\n",MyID);
            for (i=0;i<SumID*Spt;i++)
            {
                printf("%d : ",temp1[i]);
            }
            printf("\n");
            */

			/* 处理器P0将上一步送来的P段有序的数据序列做P路归并，再选择排序后的第P-1，2(P-1)，…，(P-1)(P-1)个共P-1个主元，，对应于算法13.5步骤（3）*/
            MPI_Barrier(MPI_COMM_WORLD);
            quicksort(temp1,0,SumID*Spt-1);
            /*
            printf("This is node %d \n after quicksort\n",MyID);
            for (i=0;i<SumID*Spt;i++)
            {
                printf("%d : ",temp1[i]);
            }
            printf("\n");
            */
            MPI_Barrier(MPI_COMM_WORLD);

            for (i=1;i<Spt+1;i++)
                temp1[i] = temp1[i*Spt-1];
			/*处理器P0将这P-1个主元播送到所有处理器中，对应于算法13.5步骤（4）*/
            // although the broadcast number is P, P-1 will be used, temp1[0] won't be use
            MPI_Bcast(temp1, sizeof(int)*(1+Spt), MPI_CHAR, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        else
        {
            // other process excpet 0, first send the representative numbers to process 0
            MPI_Send(temp1,sizeof(int)*Spt,MPI_CHAR,0,ALLTOONE_TYPE+MyID, MPI_COMM_WORLD);
            MPI_Barrier( MPI_COMM_WORLD);
            MPI_Barrier( MPI_COMM_WORLD);
            // then receive the broadcast from process 0
            MPI_Bcast(temp1, sizeof(int)*(1+Spt), MPI_CHAR, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }

		/*每个处理器根据上步送来的P-1个主元把自己的n/P个数据分成P段，记为处理器Pi的第j+1段，其中i=0,…,P-1，j=0,…,P-1，对应于算法13.5步骤（5）*/
        n = mylength;
        index[0] = 0;
        i = 1;
        // part i begin from index[2*i], its length = index[2*i+1] - index[2*i], i < SumID
        // find the fisrt number from process 0 that is larger than the first number in the process
        while ((arr[0]>=temp1[i])&&(i<SumID))
        {
            // the first several parts are all from index 0 to index 0
            index[2*i-1] = 0;
            index[2*i] = 0;
            i++;
        }
        // if in the process all the numbers are larger than the number get from process 0
        // the final part is from index 0 to index n (mylength)
        if (i==SumID) index[2*i-1] = n;
        c1 = 0;
        while (i<SumID)
        {
            c4 = temp1[i];
            c3 = n;
            c2 = (int)((c1+c3)/2);
            // find the index of the first number in processor k that larger than temp1[i]
            while ((arr[c2]!=c4)&&(c1<c3))
            {
                if (arr[c2]>c4)
                {
                    c3 = c2-1;
                    c2 = (int)((c1+c3)/2);
                }
                else
                {
                    c1 = c2+1;
                    c2 = (int)((c1+c3)/2);
                }
            }
            while ((arr[c2]<=c4)&&(c2<n)) c2++;
            // if reach the end of all the number in processor k
            if (c2==n)
            {
                index[2*i-1] = n;
                for (k=i;k<SumID;k++)
                {
                    index[2*k] = 0;
                    index[2*k+1] = 0;
                }
                i = SumID;
            }
            else
            {
                index[2*i] = c2;
                index[2*i-1] = c2;
            }
            c1 = c2;
            c2 = (int)((c1+c3)/2);
            i++;
        }
        if (i==SumID) index[2*i-1] = n;

        MPI_Barrier( MPI_COMM_WORLD);

		/*每个处理器送它的第i段给处理器Pi，从而使得第i个处理器含有所有处理器的第i段数据(i=0,…,P-1)，，对应于算法13.5步骤（6）*/

        j = 0;
        for (i=0;i<SumID;i++)
        {
            // every processor first choose the i th part of itself, than send other parts to other processors
            if (i==MyID)
            {
                // the length
                temp1[i] = index[2*i+1]-index[2*i];
                for (n=0;n<SumID;n++)
                    if (n!=MyID)
                {
                    // the length
                    k = index[2*n+1]-index[2*n];
                    MPI_Send(&k, sizeof(int), MPI_CHAR, n, MULTI_LEN+MyID,MPI_COMM_WORLD);
                }

            }
            // every processor receive data from other processors
            else
            {
                MPI_Recv(&temp1[i], sizeof(int), MPI_CHAR, i,MULTI_LEN+i, MPI_COMM_WORLD, &status[j++]);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        j = 0;
        k = 0;
        l = 0;

        // use arr1 to store all the recieved data
        for (i=0;i<SumID;i++)
        {
            MPI_Barrier(MPI_COMM_WORLD);

            if (i==MyID)
            {
                for (n=index[2*i];n<index[2*i+1];n++)
                    arr1[k++] = arr[n];
            }

            MPI_Barrier(MPI_COMM_WORLD);

            if (i==MyID)
            {
                for (n=0;n<SumID;n++)
                    if (n!=MyID)
                {
                    MPI_Send(&arr[index[2*n]], sizeof(int)*(index[2*n+1]-index[2*n]),MPI_CHAR, n, MULTI_TYPE+MyID, MPI_COMM_WORLD);
                }

            }
            else
            {
                l = temp1[i];
                MPI_Recv(&arr1[k], l*sizeof(int), MPI_CHAR, i ,MULTI_TYPE+i, MPI_COMM_WORLD, &status[j++]);
                k=k+l;
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }
        // renew mylength
        mylength = k;
        /*
        if(MyID==2){
            for(i = 0; i < mylength; i++){
                printf("%d : ", arr1[i]);
            }
            printf("\n");
            for(i = 0; i < SumID; i++){
                printf("%d : ", temp1[i]);
            }
            printf("\n");
        }
        */
        MPI_Barrier(MPI_COMM_WORLD);

		/*每个处理器再通过P路归并排序将上一步的到的数据排序；从而这n个数据便是有序的，，对应于算法13.5步骤（7） */
        k = 0;
        multimerge(arr1,temp1,arr,&k,SumID);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    // when k (iter) is odd, the result is stored in arr, else it is stored in arr1
    // if(MyID==2){
    // Modify by Zhanganlan
    printf("On node %d the sorted data is : \n",MyID);
    if(k%2){
        for (i=0;i<mylength;i++)
            printf("%d : ",arr[i]);
        printf("\n");
    }
    else{
        for (i=0;i<mylength;i++)
            printf("%d : ",arr1[i]);
        printf("\n");
    }
    // }
}


/*输出错误信息*/
merror(char* ch)
{
    printf("%s\n",ch);
    exit(1);
}


/*串行快速排序算法*/
quicksort(int *datas,int bb,int ee)
{
    int tt,i,j;
    tt = datas[bb];
    i = bb;
    j = ee;

    if (i<j)
    {
        while(i<j)
        {
            while ((i<j)&&(tt<=datas[j])) j--;
            if (i<j)
            {
                datas[i] = datas[j];
                i++;
                while ((i<j)&&(tt>datas[i])) i++;
                if (i<j)
                {
                    datas[j] = datas[i];
                    j--;
                    if (i==j) datas[i] = tt;
                }
                else datas[j] = tt;
            } else datas[i] = tt;
        }

        quicksort(datas,bb,i-1);
        quicksort(datas,i+1,ee);
    }
}


/*串行多路归并算法*/
multimerge(int *data1,int *ind,int *data,int *iter,int SumID)
{
    int i,j,n;

    j = 0;
    for (i=0;i<SumID;i++)
        if (ind[i]>0)
    {
        ind[j++] = ind[i];
        if (j<i+1) ind[i] = 0;
    }
    if ( j>1 )
    {
        n = 0;
        for (i=0;i<j,i+1<j;i=i+2)
        {
            merge(&(data1[n]),ind[i],ind[i+1],&(data[n]));
            ind[i] += ind[i+1];
            ind[i+1] = 0;
            n += ind[i];
        }
        // Modify by Zhanganlan, can not use data[n++]=data1[n] !
        if (j%2==1 )
            for (i=0;i<ind[j-1];i++){
                 data[n]=data1[n];
                 n++;
            }
        (*iter)++;
        multimerge(data,ind,data1,iter,SumID);
    }
}


merge(int *data1,int s1,int s2,int *data2)
{
    int i,l,m;

    l = 0;
    m = s1;
    for (i=0;i<s1+s2;i++)
    {
        if (l==s1)
            data2[i]=data1[m++];
        else
        if (m==s2+s1)
            data2[i]=data1[l++];
        else
        if (data1[l]>data1[m])
            data2[i]=data1[m++];
        else
            data2[i]=data1[l++];
    }
}
