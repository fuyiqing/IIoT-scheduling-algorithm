#ifndef _GENETIC_H_
#define _GENETIC_H_
#include"libxl.h"
using namespace std;

#define  PI    3.14159265358979323846

//遗传算法参数，种群规模（0~100）、繁殖代数、函数变量个数、交叉概率、编译概率
# define GROUP_SCALE    10   
# define MAX_GENS       350
# define N_VARS         510
# define P_MATING       0.8
# define P_MUTATION     0.25

//task 任务量最大值、最小值、机器速度最大值最小值
#define TASK_MIN 10
#define TASK_MAX 20
#define MACHINE_MAX 5
#define MACHINE_MIN 3

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

void report(int Xnration, libxl::Sheet* sheet, int nodenn);
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
