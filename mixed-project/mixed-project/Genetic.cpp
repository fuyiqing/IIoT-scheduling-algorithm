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

//����һ��graph
MyGraph graph;

//����ڵ��Ӧ�Ļ������ڴ棨��Ⱥ�ڴ棩
struct Individual nodeMachine[GROUP_SCALE + 1] = {};

//������Ⱥ�ڴ棬���ж��1���Ƿ�����һ�������������
struct Individual Population[GROUP_SCALE + 1];

//����task machine�ڴ�
double task[NODE_MAX + 2] = {};
double machine[NODE_MAX + 2] = {};

//ȷ�� task machine ����
int taskNum = 10;
int machineNum = 10;
X_Range  XnRange[N_VARS] = { { -3.0,12.1}, {4.1,5.8} };

//�ڵ���� �� ÿ����Ⱥ�ĸ������
int nodeNum = 0;

//ÿ��ͼ�е�nodes����
double nodes[NODE_MAX] = {};

//����д��xlsx�ļ�
Book* book = xlCreateBookW(); // xlCreateXMLBook() for xlsx
Sheet* sheet = book->addSheet(L"Sheet1");

//�н���Ȩ�����и������н���
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
		//x = r8_uniform_ab(a, b, seed);//�����������

		if (x < P_MATING)
		{
			++first;

			if (first % 2 == 0)//����
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

//������һ�������ŵ�һ���Ĵ������Ż���Ŀ��
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
		//�����һ�������ڵĲ�
		if (nodeMachine[i + 1].Fitness > nodeMachine[i].Fitness)
		{
			//�����õı����ڵĲ�
			if (best >= nodeMachine[i].Fitness)
			{
				best = nodeMachine[i].Fitness;
				best_mem = i;
			}
			//�����ıȺ�һ����
			if (nodeMachine[i + 1].Fitness >= worst)
			{
				worst = nodeMachine[i + 1].Fitness;
				worst_mem = i + 1;
			}

		}
		else
		{
			//������ı����ڵĺ�
			if (nodeMachine[i].Fitness >= worst)
			{
				worst = nodeMachine[i].Fitness;
				worst_mem = i;
			}
			//�����õıȺ�һ����
			if (best >= nodeMachine[i + 1].Fitness)
			{
				best = nodeMachine[i + 1].Fitness;
				best_mem = i + 1;
			}

		}

	}

	//���ڵ�ǰ��������ֵ�Ĵ��������ǰ������ֵС����һ������һ����ֵ���Ÿ���ȡ����ǰ����������
	//������
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
//������Ӧ��ֵ
void evaluate()
{
	int member;
	int i;
	double x[N_VARS + 1];

	//��������ۼ�������������ڵ������ͬ��ôʱ���ۼӣ����������λ����ͼ�����������Ӷ�ʱ�䲻������⣩
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
				cerr << "Genetic 174�� �����쳣 machine Ϊ0" << endl;
			}
			//��������ۼ�������������ڵ������ͬ��ôʱ���ۼӣ����������λ����ͼ�����������Ӷ�ʱ�䲻������⣩
			//�����ۼ�ʱ�䣨����㷨֮�����ģ�
			
			nodes[i] = t / m * times[nodeMachine[member].Xn[i]];
		}

		criticalPath(graph, nodes, nodeNum);
		nodeMachine[member].Fitness = graph.InstNum;
	}
	return;
}


//�������ε������
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
	//��֤aС��b
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
	value = round(r);//��������
	//��֤ȡֵ��Խ��
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

//��ʼ����Ⱥ����
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


//��ѡ����Сֵ����������Ⱥ��������һ��λ��
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

//�������
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
			x = r8_uniform_ab(a, b, seed);//ͻ�����
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

//ģ�庯�����������ɸ��������ϵ���������
template<typename T>
T randT(T Lower, T Upper)
{
	return rand() / (double)RAND_MAX *(Upper - Lower) + Lower;
}

//����С�������
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

//���ÿһ�������Ľ��
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

//ѡ���н���Ȩ�ĸ���
void selector(int &seed)
{
	struct Individual NewPopulation[GROUP_SCALE + 1];//��ʱ�����ѡ�ĺ������
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
	//��������ܶ�
	for (mem = 0; mem < GROUP_SCALE; mem++)
	{
		nodeMachine[mem].ReFitness = nodeMachine[mem].Fitness / sum;
	}
	// �����ۼӷֲ���˼�������̷�
	nodeMachine[0].SumFitness = nodeMachine[0].ReFitness;
	for (mem = 1; mem < GROUP_SCALE; mem++)
	{
		nodeMachine[mem].SumFitness = nodeMachine[mem - 1].SumFitness +
			nodeMachine[mem].ReFitness;
	}
	// ѡ�����Ϊ��һ����ֳ��ѡ������Ŀ����Դ��������̷��İ���֮��
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
	//���º������ 
	for (i = 0; i < GROUP_SCALE; i++)
	{
		nodeMachine[i] = NewPopulation[i];
	}
	return;
}

//��ʾϵͳʱ��
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

//��������Ӵ�
void Xover(int one, int two, int &seed)
{
	int i;
	int point;
	int t;
	//���ѡ�񽻲�㣬����ĵ����Ա�������������Ϊ��λ
	point = i4_uniform_ab(0, nodeNum - 1, seed);
	//point = i4_uniform_ab(0, N_VARS - 1, seed);
	//����
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
		//����д�Ͷ˵�
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