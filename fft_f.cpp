/********************************************************************
	    Time:	2013/07/28 12:32
	filename: 	fft_f.cpp
	  author:   huiming
	
	 purpose:	FFT�任 Ƶ���2��ȡ
*********************************************************************/
#include<iostream>
#include<math.h>
#include<string.h>
#define PI 3.14159
using namespace std;
double x1[100000],x2[100000],w1[100000],w2[100000];//�Ǳ�1����ʵ������,2��������
int visited[100000];
int s[1000];    
void wnp(int M,int L,int N) //��(WN)^p
{
	int i,k;
	int t1=(int)(pow(2.0,M-L))-1;
	double t2=-2*PI/N;
	double a,b;
	memset(w1,0,sizeof(w1));
	memset(w2,0,sizeof(w2));
	k=(int)(pow(2.0,L-1));
	a=cos(t2*k);
	b=sin(t2*k);
	w1[0]=1; //��0�����
	w2[0]=0;
	for(i=1;i<=t1;i++)
	{
		w1[i]=a*w1[i-1]-b*w2[i-1];
		w2[i]=a*w2[i-1]+b*w1[i-1];
	}
}
void FFT(int N,int M)
{
	int i,j,n=N;
	double a,b;
	for(i=1;i<=M;i++)  //�ӵ�1����M��
	{
		n/=2;
		memset(visited,0,sizeof(visited));//�����������
		wnp(M,i,N);  //����wnp����    
		for(j=0;j<N;j++)   //�ֱ���ÿһ���ĵ�ǰʵ�����ֺ���������
		{
			if(!visited[j])
			{
				visited[j]=1;
				visited[j+n]=1;
				a=x1[j]-x1[j+n]; //�˴������� Ӧ�ȼ������ֵ ��Ϊ���� x1[j]��x2[j]��ֵ��ı�
				b=x2[j]-x2[j+n];
				x1[j]=x1[j]+x1[j+n];
				x2[j]=x2[j]+x2[j+n];                

				int t=j%(N/(int)(pow(2.0,i*1.0)));
				x1[j+n]=a*w1[t]-b*w2[t];
				x2[j+n]=a*w2[t]+b*w1[t];
			}
		}    
	}
}
void solve(double *x,int N,int M) //��λ����
{
	int a,i,j,k;
	double t;
	for(k=0;k<N/2;k++)
	{
		i=k;
		a=0;
		memset(s,0,sizeof(s));    
		for(j=0;j<M;j++)
		{
			s[j]=i%2;
			i/=2;
		}
		for(j=0;j<M;j++)
		{
			a=a+s[j]*(int)(pow(2.0,M-1-j));
		}
		t=x[a];
		x[a]=x[k];
		x[k]=t;
	}    
}
int main()
{
	//freopen("d:\\1.txt","r",stdin);
	int N,i,M;
	cout<<"���������γ���N(N����2�������η�): ";
	cin>>N;
	M=floor(log10(N*1.0)/log10(2.0)+0.5);
	cout<<"��ֱ�����N������ֵ���и�����ʵ�����鲿: "<<endl;
	for(i=0;i<N;i++)
	{
		printf("ʵ��x1[%d]=",i);
		cin>>x1[i];
		printf("�鲿x2[%d]=",i);
		cin>>x2[i];
	}
	FFT(N,M);
	solve(x1,N,M);
	solve(x2,N,M);
	printf("�õ���Ƶ��ֵΪ:\n");
	for(i=0;i<N;i++)
		printf("X[%d]=(%.2lf)+(%.2lf)j\n",i,x1[i],x2[i]);
	return 0;
}