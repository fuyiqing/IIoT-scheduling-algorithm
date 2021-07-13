#ifndef _GENETIC_H_
#define _GENETIC_H_
#include"libxl.h"
using namespace std;

#define  PI    3.14159265358979323846

//�Ŵ��㷨��������Ⱥ��ģ��0~100������ֳ��������������������������ʡ��������
# define GROUP_SCALE    10   
# define MAX_GENS       350
# define N_VARS         510
# define P_MATING       0.8
# define P_MUTATION     0.25

//task ���������ֵ����Сֵ�������ٶ����ֵ��Сֵ
#define TASK_MIN 10
#define TASK_MAX 20
#define MACHINE_MAX 5
#define MACHINE_MIN 3

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

void report(int Xnration, libxl::Sheet* sheet, int nodenn);
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
