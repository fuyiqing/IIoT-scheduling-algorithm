#ifndef _GENETIC_H_
#define _GENETIC_H_
#include"libxl.h"
using namespace std;

#define  PI    3.14159265358979323846


//machine�б���λ��MFLops��
//�Ͷ� machine
constexpr double machineLow[20] = { 80,64,15,30,50,80,72,25,48,64,63,24,16,80,32,48,32,48,70,24 };
//�ж�machine
constexpr double machineMid[20] = { 550,280,120,180,150,200,300,252,500,600,800,100,960,266,240,200,100,300,264,480 };
//�߶�machine
constexpr double machineHigh[20] = { 2500,1000,2200,1800,1400,1500,1200,1100,1176,1900,2300,2500,2600,1600,1700,1800,1300,2200,3200,3000 };

//�Ŵ��㷨��������Ⱥ��ģ��0~100������ֳ��������������������������ʡ��������
# define GROUP_SCALE    10   
# define MAX_GENS       500
# define N_VARS         510
# define P_MATING       0.8
# define P_MUTATION     0.25

//�Ͷ˻�����50-100M
//�ж˻�����
//�߶˻�����
#define TASK_MIN 50
#define TASK_MAX 100
//machine��λ�� TFLOPS ��T����������/�룩 0.015-0.025
//#define MACHINE_MAX 10
//#define MACHINE_MIN 4
//Ŀ�ģ���������Ӧ��ֵ������0-1֮��

//��Ӧ�����ֵ
#define FITNESS_MAX 10000000

struct Individual
{
	int Xn[N_VARS];      //��ű���ֵ(����������)
	double Fitness;         //��Ӧֵ
	double ReFitness;       //��Ӧֵ�����ܶ�
	double SumFitness;      //�ۼӷֲ���Ϊ����ת
};
struct X_Range
{
	double Upper;           //�������Ͻ�ȡֵ
	double Lower;           //�������½�ȡֵ
};

template<typename T>
T randT(T Lower, T Upper); //���������������������

void crossover(int &seed);
void elitist();        //������
void evaluate();

void initGroup(int &seed);

void selectBest();
void mutate(int &seed);

double r8_uniform_ab(double a, double b, int &seed);
int i4_uniform_ab(int a, int b, int &seed);

void report(int Xnration, libxl::Sheet* sheet);
void selector(int &seed);
void showTime();
void Xover(int one, int two, int &seed);

//*@ param task      ����ʼ����task����
//*@ param length  ����ĳ���
//*@ note                ��task��������������ֵ
void taskInit(double *task, int length);

//*@ param machine      ����ʼ����machine����
//*@ param length         ����ĳ���
//*@ note                       ��machine��������������ֵ
void machineInit(double *machine, int length);
#endif // !_GENETIC_H_
