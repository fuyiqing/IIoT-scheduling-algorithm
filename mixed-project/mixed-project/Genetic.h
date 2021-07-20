#ifndef _GENETIC_H_
#define _GENETIC_H_
#include"libxl.h"
using namespace std;

#define  PI    3.14159265358979323846


//machine列表（单位：MFLops）
//低端 machine
constexpr double machineLow[20] = { 80,64,15,30,50,80,72,25,48,64,63,24,16,80,32,48,32,48,70,24 };
//中端machine
constexpr double machineMid[20] = { 550,280,120,180,150,200,300,252,500,600,800,100,960,266,240,200,100,300,264,480 };
//高端machine
constexpr double machineHigh[20] = { 2500,1000,2200,1800,1400,1500,1200,1100,1176,1900,2300,2500,2600,1600,1700,1800,1300,2200,3200,3000 };

//遗传算法参数，种群规模（0~100）、繁殖代数、函数变量个数、交叉概率、编译概率
# define GROUP_SCALE    10   
# define MAX_GENS       500
# define N_VARS         510
# define P_MATING       0.8
# define P_MUTATION     0.25

//低端机器：50-100M
//中端机器：
//高端机器：
#define TASK_MIN 50
#define TASK_MAX 100
//machine单位： TFLOPS （T个浮点运算/秒） 0.015-0.025
//#define MACHINE_MAX 10
//#define MACHINE_MIN 4
//目的：将最优适应度值控制在0-1之间

//适应度最大值
#define FITNESS_MAX 10000000

struct Individual
{
	int Xn[N_VARS];      //存放变量值(即机器索引)
	double Fitness;         //适应值
	double ReFitness;       //适应值概率密度
	double SumFitness;      //累加分布，为轮盘转
};
struct X_Range
{
	double Upper;           //变量的上界取值
	double Lower;           //变量的下界取值
};

template<typename T>
T randT(T Lower, T Upper); //产生任意类型随机数函数

void crossover(int &seed);
void elitist();        //基因保留
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

//*@ param task      待初始化的task数组
//*@ param length  数组的长度
//*@ note                给task数组进行随机化赋值
void taskInit(double *task, int length);

//*@ param machine      待初始化的machine数组
//*@ param length         数组的长度
//*@ note                       给machine数组进行随机化赋值
void machineInit(double *machine, int length);
#endif // !_GENETIC_H_
