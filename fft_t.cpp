/********************************************************************
	    Time:	2013/07/25 7:59
	filename: 	fft.cpp
	  author:   huiming
	
	 purpose:	���ٸ���Ҷ�任 (tʱ���2��ȡ)
*********************************************************************/
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;
#define  PI acos(-1)

//������
class complex
{
public:
	complex(double x=0.0,double y=0.0):real(x),img(y){}
	complex operator + (const complex Rhs)
	{
		return complex(this->real+Rhs.real,this->img+Rhs.img);
	}
	complex operator -(const complex Rhs)
	{
		return complex(this->real-Rhs.real,this->img-Rhs.img);
	}
	complex operator *(const complex Rhs)
	{
		//�����ĳ˷������ƶ���ʽ�˷�
		return complex(this->real*Rhs.real-this->img*Rhs.img,this->real*Rhs.img+this->img*Rhs.real);
	}
	double getReal() const
	{
		return real;
	}
	double getImg() const
	{
		return img;
	}
private:
	double real;
	double img;
};

int log2(int value)
{
	int i=0;
	while (value>1)
	{
		value=value>>1;
		i++;
	}
	return i;
}

//λ��ת�û�
void BRC(vector<complex> &v)
{
	int k,i,j,t;
	for (i=0;i<v.size()-1;i++)
	{
		k=i;
		j=0;
		t=log2(v.size());
		while((t--)>0)
		{
			j=j<<1;
			j|=(k&1);
			k=k>>1;
		}
		if (j>i)
		{
			swap(v[i],v[j]);
		}
	}
}

//
void FFT(vector<complex> &v,int type)
{
	int i,j,k,l;
	BRC(v);
	//һ����������
	for (i=0;i<log2(v.size());i++)
	{

		l=1<<i;
		//һ���������
		for (j=0;j<v.size();j+=l*2)
		{
			//һ����������
			for (k=0;k<l;k++)
			{
				complex w(cos(type*2*PI*k/(2*l)),-1*sin(type*2*PI*k/(2*l)));
				complex v1=v[j+k]+w*v[j+k+l];
				complex v2=v[j+k]-w*v[j+k+l];
				v[j+k]=v1;
				v[j+k+l]=v2;
			}

		}
	}
	//fft��任
	if (type==-1)
	{
		for (i=0;i<v.size();i++)
		{
			v[i]=complex(v[i].getReal()/v.size(),v[i].getImg()/v.size());
		}
	}
}

int main()
{
	int lenght;
	vector<complex> com_v;
	cout<<"please input the input vector size: \n";
	cin>>lenght;
	if (lenght<0||lenght%2)
	{
		return;
	}
	cout<<"log2(length)="<<log2(lenght)<<"PI="<<PI<<"real:img"<<endl;
	for (int i=0;i<lenght;i++)
	{
// 		double a,b;
// 		cin>>a>>b;
		complex temp(i,0);
		com_v.push_back(temp);
	}
	FFT(com_v,-1);
	for ( i=0;i<lenght;i++)
	{
		cout<<com_v[i].getReal()<<":"<<com_v[i].getImg()<<endl;
	}
}







