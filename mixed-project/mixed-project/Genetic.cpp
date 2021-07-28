# include <cstdlib>
#include<cstdio>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include"MyGraph.h"
#include "Genetic.h"
#pragma warning(disable : 4996)
#include "libxl.h"

using namespace libxl;

extern int seed;
//extern  constexpr double machineLow[20];
//extern  constexpr double machineMid[20];
//extern double machineHigh[20];

//定义一个graph
MyGraph graph;

//申请节点对应的机器的内存（种群内存）
struct Individual nodeMachine[GROUP_SCALE + 1] = {};

//申请种群内存，其中多加1个是放置上一代中最优秀个体
struct Individual Population[GROUP_SCALE + 1];

//申请task machine内存
double task[NODE_MAX + 2] = {};
double machine[NODE_MAX + 2] = {};

//确定 task machine 数量
int taskNum = 10;
int machineNum = 10;
X_Range  XnRange[N_VARS] = { { -3.0,12.1}, {4.1,5.8} };

//节点个数 即 每个种群的个体个数
int nodeNum = 0;

//每个图中的nodes数组
double nodes[NODE_MAX] = {};

//申请写入xlsx文件
Book* book = xlCreateBookW(); // xlCreateXMLBook() for xlsx
Sheet* sheet = book->addSheet(L"Sheet1");

//有交配权的所有父代进行交叉
void crossover(int &seed)
{
	const double a = 0.0;
	const double b = 1.0;
	int mem;
	int one;
	int first = 0;
	double x;

	for (mem = 0; mem < GROUP_SCALE; ++mem)
	{
		x = r8_uniform_ab(0.0, 1.0, seed);
		//x = r8_uniform_ab(a, b, seed);//产生交配概率

		if (x < P_MATING)
		{
			++first;

			if (first % 2 == 0)//交配
			{
				Xover(one, mem, seed);
			}
			else
			{
				one = mem;
			}

		}
	}
	return;
}

//对最差的一代和最优的一代的处理，起到优化的目的
void elitist()
{
	int i;
	double best;
	int best_mem;
	double worst;
	int worst_mem;

	best = nodeMachine[0].Fitness;
	worst = nodeMachine[0].Fitness;

	for (i = 0; i < GROUP_SCALE - 1; ++i)
	{
		//如果后一代比现在的差
		if (nodeMachine[i + 1].Fitness > nodeMachine[i].Fitness)
		{
			//如果最好的比现在的差
			if (best >= nodeMachine[i].Fitness)
			{
				best = nodeMachine[i].Fitness;
				best_mem = i;
			}
			//如果最坏的比后一代好
			if (nodeMachine[i + 1].Fitness >= worst)
			{
				worst = nodeMachine[i + 1].Fitness;
				worst_mem = i + 1;
			}

		}
		else
		{
			//如果最差的比现在的好
			if (nodeMachine[i].Fitness >= worst)
			{
				worst = nodeMachine[i].Fitness;
				worst_mem = i;
			}
			//如果最好的比后一代差
			if (best >= nodeMachine[i + 1].Fitness)
			{
				best = nodeMachine[i + 1].Fitness;
				best_mem = i + 1;
			}

		}

	}

	//对于当前代的最优值的处理，如果当前的最优值小于上一代则将上一代的值最优个体取代当前的最弱个体
	//基因保留
	if (nodeMachine[GROUP_SCALE].Fitness >= best)
	{
		for (i = 0; i < nodeNum; i++)
		{
			nodeMachine[GROUP_SCALE].Xn[i] = nodeMachine[best_mem].Xn[i];
		}
		nodeMachine[GROUP_SCALE].Fitness = nodeMachine[best_mem].Fitness;
	}
	else
	{
		for (i = 0; i < nodeNum; i++)
		{
			nodeMachine[worst_mem].Xn[i] = nodeMachine[GROUP_SCALE].Xn[i];
		}
		nodeMachine[worst_mem].Fitness = nodeMachine[GROUP_SCALE].Fitness;
	}
	return;
}
//计算适应度值
void evaluate()
{
	int member;
	int i;
	double x[N_VARS + 1];

	//计算机器累加量，如果两个节点机器相同那么时间累加（解决后面三位迭代图机器数量增加而时间不变的问题）
//	double machineaccelerate[100] = {};
	



	for (member = 0; member < GROUP_SCALE; member++)
	{
		int times[100] = {};

		for (int i = 0; i < 100; i++)
		{
			times[i] = 1;

		}
	
			for (int j = 0; j < nodeNum; j++)
			{
				times[nodeMachine[member].Xn[j]]++;
			}
		


		for (i = 0; i < nodeNum; i++)
		{
			int index;
			double t, m;
			try
			{
				index = nodeMachine[member].Xn[i];
				t = task[i], m = machine[index];
			
				if (fabs(m - 0) <= 0.0001)
				{
					throw 0;
				}
			}
			catch (...)
			{
				cerr << "Genetic 174行 发现异常 machine 为0" << endl;
			}
			//计算机器累加量，如果两个节点机器相同那么时间累加（解决后面三位迭代图机器数量增加而时间不变的问题）
			//就是累加时间（这个算法之后再聊）
			
			nodes[i] = t / m * times[nodeMachine[member].Xn[i]];
		}

		criticalPath(graph, nodes, nodeNum);
		nodeMachine[member].Fitness = graph.InstNum;
	}
	return;
}


//产生整形的随机数
int i4_uniform_ab(int a, int b, int &seed)
{
	int c;
	const int i4_huge = 2147483647;
	int k;
	float r;
	int value;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "I4_UNIFORM_AB - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}
	//保证a小于b
	if (b < a)
	{
		c = a;
		a = b;
		b = c;
	}

	k = seed / 127773;
	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed = seed + i4_huge;
	}

	r = (float)(seed)* 4.656612875E-10;
	//
	//  Scale R to lie between A-0.5 and B+0.5.
	//
	r = (1.0 - r) * ((float)a - 0.5)
		+ r * ((float)b + 0.5);
	//
	//  Use rounding to convert R to an integer between A and B.
	//
	value = round(r);//四舍五入
	//保证取值不越界
	if (value < a)
	{
		value = a;
	}
	if (b < value)
	{
		value = b;
	}

	return value;
}

//初始化种群个体
void initGroup(int &seed)
{
	int i;
	int j;
	double lbound;
	double ubound;
	// 
	//  initGroup variables within the bounds 
	//
//	for (i = 0; i < N_VARS; i++)
	//{
		//input >> lbound >> ubound;

	for (j = 0; j < GROUP_SCALE; j++)
	{
		nodeMachine[j].Fitness = FITNESS_MAX;
		nodeMachine[j].ReFitness = 0;
		nodeMachine[j].SumFitness = 0;
		for (i = 0; i < nodeNum; i++)
		{
			nodeMachine[j].Xn[i] = i4_uniform_ab(0, machineNum - 1, seed);
		}
		//Population[j].Xn[i] = r8_uniform_ab(XnRange[i].Lower, XnRange[i].Upper, seed);
	}
	nodeMachine[GROUP_SCALE].Fitness = FITNESS_MAX;
	nodeMachine[GROUP_SCALE].ReFitness = 0;
	nodeMachine[GROUP_SCALE].SumFitness = 0;
	//	}
	return;
}


//挑选出最小值，保存在种群数组的最后一个位置
void selectBest()
{
	int cur_best;
	int mem;
	int i;

	cur_best = 0;

	for (mem = 0; mem < GROUP_SCALE; mem++)
	{
		if (nodeMachine[GROUP_SCALE].Fitness > nodeMachine[mem].Fitness)
		{
			cur_best = mem;
			nodeMachine[GROUP_SCALE].Fitness = nodeMachine[mem].Fitness;
		}
	}

	for (i = 0; i < nodeNum; i++)
	{
		nodeMachine[GROUP_SCALE].Xn[i] = nodeMachine[cur_best].Xn[i];
	}

	return;
}

//个体变异
void mutate(int &seed)
{
	const double a = 0.0;
	const double b = 1.0;
	int i;
	int j;
	int lbound;
	int ubound;
	double x;

	for (i = 0; i < GROUP_SCALE; i++)
	{
		for (j = 0; j < nodeNum; j++)
		{
			//x = r8_uniform_ab(a, b, seed);
			x = r8_uniform_ab(a, b, seed);//突变概率
			if (x < P_MUTATION)
			{
				lbound = 0;
				ubound = machineNum - 1;
				nodeMachine[i].Xn[j] = i4_uniform_ab(lbound, ubound, seed);
				//Population[i].Xn[j] = r8_uniform_ab(lbound, ubound, seed);
			}
		}
	}

	return;
}

//模板函数，用于生成各种区间上的数据类型
template<typename T>
T randT(T Lower, T Upper)
{
	return rand() / (double)RAND_MAX *(Upper - Lower) + Lower;
}

//产生小数随机数
double r8_uniform_ab(double a, double b, int &seed)

{
	int i4_huge = 2147483647;
	int k;
	double value;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "R8_UNIFORM_AB - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}

	k = seed / 127773;
	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed = seed + i4_huge;
	}

	value = (double)(seed)* 4.656612875E-10;

	value = a + (b - a) * value;

	return value;
}

//输出每一代进化的结果
void report(int Xnration, libxl::Sheet* sheet)
{
	double avg;
	double best_val;
	int i;
	double square_sum;
	double stddev;
	double sum;
	double sum_square;



	//	if (Xnration == 0)
	//	{
	//		cout << "\n";
	//		cout << "  Xnration       Best            Average       Standard \n";
	//		cout << "  number           value           Fitness       deviation \n";
	//		cout << "\n";
	//	}
	sum = 0.0;
	sum_square = 0.0;

	for (i = 0; i < GROUP_SCALE; i++)
	{
		sum = sum + nodeMachine[i].Fitness;
		sum_square = sum_square + nodeMachine[i].Fitness * nodeMachine[i].Fitness;
	}

	avg = sum / (double)GROUP_SCALE;
	square_sum = avg * avg * GROUP_SCALE;
	stddev = sqrt((sum_square - square_sum) / (GROUP_SCALE - 1));
	best_val = nodeMachine[GROUP_SCALE].Fitness;

	//	cout << "  " << setw(8) << Xnration
//	cout << "  " << setw(14) << best_val << endl;
	sheet->writeNum(machineNum-3, nodeNum-4, best_val);

	return;
}

//选择有交配权的父代
void selector(int &seed)
{
	struct Individual NewPopulation[GROUP_SCALE + 1];//临时存放挑选的后代个体
	const double a = 0.0;
	const double b = 1.0;
	int i;
	int j;
	int mem;
	double p;
	double sum;

	sum = 0.0;
	for (mem = 0; mem < GROUP_SCALE; mem++)
	{
		sum = sum + nodeMachine[mem].Fitness;
	}
	//计算概率密度
	for (mem = 0; mem < GROUP_SCALE; mem++)
	{
		nodeMachine[mem].ReFitness = nodeMachine[mem].Fitness / sum;
	}
	// 计算累加分布，思想是轮盘法
	nodeMachine[0].SumFitness = nodeMachine[0].ReFitness;
	for (mem = 1; mem < GROUP_SCALE; mem++)
	{
		nodeMachine[mem].SumFitness = nodeMachine[mem - 1].SumFitness +
			nodeMachine[mem].ReFitness;
	}
	// 选择个体为下一代繁殖，选择优秀的可能性大，这是轮盘法的奥秘之处
	for (i = 0; i < GROUP_SCALE; i++)
	{
		p = r8_uniform_ab(a, b, seed);
		if (p < nodeMachine[0].SumFitness)
		{
			NewPopulation[i] = nodeMachine[0];
		}
		else
		{
			for (j = 0; j < GROUP_SCALE; j++)
			{
				if (nodeMachine[j].SumFitness <= p && p < nodeMachine[j + 1].SumFitness)
				{
					NewPopulation[i] = nodeMachine[j + 1];
					break;
				}
			}
		}
	}
	//更新后代个体 
	for (i = 0; i < GROUP_SCALE; i++)
	{
		nodeMachine[i] = NewPopulation[i];
	}
	return;
}

//显示系统时间
void showTime()
{
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	size_t len;
	time_t now;

	now = time(NULL);
	tm = localtime(&now);

	len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

	cout << time_buffer << "\n";

	return;
# undef TIME_SIZE
}

//交叉产生子代
void Xover(int one, int two, int &seed)
{
	int i;
	int point;
	int t;
	//随机选择交叉点，这里的点是以变量的整个长度为单位
	point = i4_uniform_ab(0, nodeNum - 1, seed);
	//point = i4_uniform_ab(0, N_VARS - 1, seed);
	//交叉
	for (i = 0; i < point; i++)
	{
		t = nodeMachine[one].Xn[i];
		nodeMachine[one].Xn[i] = nodeMachine[two].Xn[i];
		nodeMachine[two].Xn[i] = t;
	}
	return;
}

void taskInit(double *task, int length)
{
	for (int i = 0; i < length; i++)
	{
		if (i == 0 || i == length - 1)
		{
			task[i] = 0;
			continue;
		}
		//首先写低端的
		task[i] = r8_uniform_ab(TASK_MIN,TASK_MAX,seed);
	}
}

void machineInit(double *machine, int length)
{
	for (int i = 0; i < length; i++)
	{
		while (machine[i] == 0)
		{
			int m = i4_uniform_ab(0,19,seed);
			//printf("fffffff%dfffffff\n",m);
			machine[i] = machineLow[m];
			//printf("fffffff%lffffffff\n", machine[i]);
		}
	}
}