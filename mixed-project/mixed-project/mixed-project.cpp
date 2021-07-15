// mixed-project.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#pragma warning(disable : 4996)
#include <iostream>
#include "Genetic.h"
#include "MyGraph.h"
#include "libxl.h"
#include <time.h>
#include<cmath>
#include<string>
#include<string.h>
#include<fstream>
#include<vector>
#include<random>

using namespace libxl;

//遗传算法的参数----------------------------------------------------------//
int seed = 123456789;
extern Individual Population[GROUP_SCALE + 1];

extern MyGraph graph;

extern Individual nodeMachine[GROUP_SCALE + 1];

extern double task[NODE_MAX + 2];
extern double machine[NODE_MAX + 2];

extern int taskNum;
extern int machineNum;
extern int nodeNum;
extern double nodes[NODE_MAX];
//---------------------------------------------------------------------------//

//布谷鸟算法的参数---------------------------------------------------------//
//每次迭代的鸟巢数
constexpr int NestNum = 10;
//最大迭代次数
constexpr int MaxIterationTimes = 500;
//被宿主发现的概率
constexpr double Pa = 0.55;
//最大种群规模
# define N_VARS     510
//自变量结构体
struct Nest {
	int Xn[N_VARS];      //存放变量值(即机器索引)

	double fitness;
};
#define pi acos(-1)
void fitFunc(Nest& nest);
int findBetterNest(std::vector<Nest>&);
std::vector<Nest> levy(std::vector<Nest> OldNestPop, std::size_t bestNum);
std::vector<Nest> RandomAbandonPaNestPop(std::vector<Nest> OldNestPop);
//随机数引擎
static std::default_random_engine e(time(0));
static std::uniform_real_distribution<double> u(0, 1);

//产生整形的随机数
int my_i4_uniform_ab(int a, int b, int &seed);
//--------------------------------------------------------------------------//


Book* bookGenetic = xlCreateBookW();
Sheet* sheetGen = bookGenetic->addSheet(L"Sheet1");

Book* bookCuckoo = xlCreateBookW();
Sheet* sheetCuc = bookCuckoo->addSheet(L"Sheet1");

//*@ param diedai 节点个数
//*@ note    对全局变量graph进行操作：初始化，随机化 
//*@ note    对全局变量task和machine进行随机化
void taskMachineGraphFunction(int diedai);

//*@ param diedai 节点个数
//*@ param 没有其他输入参数，输入参数都是全局变量
//*@ note    genetic-algorithm算法  对每一个图以及对应的task和machine生成对应的最优值 
void genetic(int diedai);

//*@ param diedai 节点个数
//*@ note    cuckoo-algorithm算法 对每一个图以及对应的task和machine生成对应的最优值 
void cuckoo(int diedai);

//*@ param 没有参数
//*@ note    输出当前时间（年月日 时分秒）
void myShowtime();

int main()
{
	
	//输入genetic对应的适应度
	sheetGen->writeStr(1, 0, L"节点个数");
	for (int nodeG = 5; nodeG <= 50; nodeG++)
		sheetGen->writeNum(1, nodeG - 4, nodeG);
	for (int iiG = 1; iiG <= 500; iiG++)
		sheetGen->writeNum(iiG + 1, 0, iiG);

	//输入cuckoo对应的适应度
	sheetCuc->writeStr(1,0,L"节点个数");
	for (int i = 1; i <= 500; i++)
		sheetCuc->writeNum(i + 1, 0, i);
	for (int i = 5; i <= 50; i++)
		sheetCuc->writeNum(1, i - 4, i);

	srand(time(nullptr));
	

	//---------------------------------------主要代码-----------------------------------------------
	
	for (int diedai = 5; diedai <= 50; diedai++) 
	{
		//先对机器、任务、图进行初始化并存在txt文件中
		taskMachineGraphFunction(diedai);

		//首先开始进行遗传算法
		genetic(diedai);

		//其次开始布谷鸟算法
		cuckoo(diedai);

	}
	//-------------------------------------主要代码-----------------------------------------------------

	//保存并输出BookGen
	bookGenetic->save(L"genetic_5-50.xlsx");
	bookGenetic->release();

	//保存并输出BookCuc
	bookCuckoo->save(L"cuckoo_5-50.xlsx");
	bookCuckoo->release();

	//输出时间：年月日 时分秒
	myShowtime();
}


void taskMachineGraphFunction(int diedai)
{
	//先试试解决八个
	nodeNum = diedai;

	//先初始化task和machine
	taskNum = diedai;
	machineNum = ((int)(sqrt(taskNum))) * 2;
	taskInit(task, taskNum);
	machineInit(machine, machineNum);

	//再初始化随机化graph
	graphInit(graph, nodeNum);
	graphRandom(graph);

	//将graph输出到graph_diedai.txt
	ofstream output;
	string fileName = "graph" + to_string(diedai)+".txt";
	const char*  p = fileName.data();

	output.open(p);
	output << "digraph D" <<diedai<<"{"<< endl;
	for (int i = 0; i < nodeNum; i++)
	{
		Arc* arc = graph.nodes[i].edge;
		while (arc)
		{
			string x = to_string(i) + "->" + to_string(arc->go) + ";";
			output << x << endl;
			arc = arc->next;
		}
	}
	output << "}" << endl;

	//输出task数组
	output << endl << endl;
	output << "task:" <<taskNum<< "个" << endl;
	for (int i = 0; i < taskNum; i++)
		output << task[i] << endl;

	//输出machine数组
	output << endl << endl;
	output << "machine:" << machineNum << "个" << endl;
	for (int i = 0; i < machineNum; i++)
		output << machine[i] << endl;

	output.close();
}

void genetic(int diedai)
{
	clock_t start, finish;
	double  duration;
	start = clock();
	//先试试解决八个
	nodeNum = diedai;

	//先初始化task和machine
	taskNum = diedai;
	machineNum = ((int)(sqrt(taskNum))) * 2;
	taskInit(task, taskNum);
	machineInit(machine, machineNum);

	//再初始化随机化graph
	graphInit(graph, nodeNum);
	graphRandom(graph);

	//初始化种群
	initGroup(seed);
	//计算适应度
	evaluate();
	//选择最优的个体放入nodeMachine[GROUP_SCALE]中 
	selectBest();

	//这上面的代码都测试成功了
	//记得把Population数组替换成nodeMachine



	//开始遗传迭代了
	int Xnration = 0; int i = 0;

//	showTime();
	initGroup(seed);
	evaluate();
	selectBest();
	for (Xnration = 0; Xnration < MAX_GENS; Xnration++)
	{
		selector(seed);
		crossover(seed);
		mutate(seed);
		report(Xnration, sheetGen, diedai);
		evaluate();
		elitist();
	}


	finish = clock();
	duration = 1000 * (finish - start) / CLOCKS_PER_SEC;
	sheetGen->writeStr(508, diedai - 4, L"总共用时(毫秒)：");
	sheetGen->writeNum(509, diedai - 4, duration);
	//printf("\n总共用时%lf毫秒\n", duration);
}


void cuckoo(int index)
{
	//记录运行时间
	clock_t start, finish;
	double duration;
	start = clock();
		//现在的鸟巢群
	std::vector<Nest> current_nestPop;
	//迭代次数
	int num = 0;
	//最优鸟巢
	int BestNestCurrent;


	//初始化
	for (int i = 0; i < NestNum; ++i)
	{
		Nest nestinitial;
		for (int j = 0; j < nodeNum; j++)
		{
			nestinitial.Xn[j] = i4_uniform_ab(0, machineNum - 1, seed);
		}
		fitFunc(nestinitial);
		//		printf("%lf\n",nestinitial.fitness);
		current_nestPop.push_back(nestinitial);
	}
	BestNestCurrent = findBetterNest(current_nestPop);


	//开始循环
	double xOld[NODE_MAX] = {};
	while (num < MaxIterationTimes)
	{
		//储存上次的最优解的X,Y
		for (int i = 0; i < nodeNum; i++)
			xOld[i] = current_nestPop[BestNestCurrent].Xn[i];

		//levy飞行--位置更新
		std::vector<Nest> NewNestPop = levy(current_nestPop, BestNestCurrent);

		//用适应值较好的鸟窝位置替换适应值较差的鸟窝位置
		for (decltype(NewNestPop.size()) i = 0; i < NewNestPop.size(); ++i)
		{
			if (i != BestNestCurrent && NewNestPop[i].fitness < current_nestPop[i].fitness)
			{
				current_nestPop[i] = NewNestPop[i];
			}
		}//此时得到更优的鸟窝位置


		//存安去险 保留鸟窝中被发现概率较小的鸟窝位置，并随机改变发现概率较大的鸟窝位置
		NewNestPop = RandomAbandonPaNestPop(current_nestPop);

		for (decltype(NewNestPop.size()) i = 0; i < NewNestPop.size(); ++i)
		{
			if (i != BestNestCurrent && NewNestPop[i].fitness < current_nestPop[i].fitness)
			{
				current_nestPop[i] = NewNestPop[i];
			}
		}//此时得到更优的鸟窝位置

		BestNestCurrent = findBetterNest(current_nestPop);//现在的最优鸟巢位置
		int nodeNew[NODE_MAX] = {};
		bool flag = false;
		for (int j = 0; j < nodeNum; j++) {
			nodeNew[j] = current_nestPop[BestNestCurrent].Xn[j];
		}

		sheetCuc->writeNum(num + 2, index - 4, current_nestPop[BestNestCurrent].fitness);
		//			printf("%lf\n", current_nestPop[BestNestCurrent].fitness);
		num++;
	}

	finish = clock();
	double time = 1000 * (finish - start) / CLOCKS_PER_SEC;
	sheetCuc->writeStr(509, index - 4, L"运行时间(毫秒)：");
	sheetCuc->writeNum(510, index - 4, time);

//	std::cout << "运行时间" << time << "毫秒" << std::endl;
}



/*
*@ param: 鸟巢
*@ note:    适应度函数，计算每个鸟巢的适应度
*/
void fitFunc(Nest& nest)
{
	int member;
	int i;

	for (i = 0; i < nodeNum; i++)
	{
		int index;
		double t, m;
		try
		{
			index = nest.Xn[i];
			t = task[i], m = machine[index];
			if (fabs(m - 0) <= 0.0001)
			{
				throw 0;
			}
		}
		catch (...)
		{
			std::cerr << "Genetic 174行 发现异常 machine 为0" << std::endl;
		}
		nodes[i] = t / m;
	}

	criticalPath(graph, nodes, nodeNum);
	nest.fitness = graph.InstNum;

	//nest.fitness = -(nest.x - 1) * (nest.x - 1) + 1;
}

int findBetterNest(std::vector<Nest>& nestPop)
{
	int BestNum = 0;
	for (decltype(nestPop.size()) i = 0; i < nestPop.size(); ++i)
	{
		if (nestPop[i].fitness < nestPop[BestNum].fitness)
		{
			BestNum = i;
		}
	}

	return BestNum;
}

std::vector<Nest> levy(std::vector<Nest> OldNestPop, std::size_t bestNum)
{
	double beta = 1.5;
	//	double alpha = 0.01 * R(e);//有的论文写
	double alpha = 0.4;
	double sigma_u = pow((tgamma(1 + beta) * sin(pi * beta / 2)) / (beta * tgamma((1 + beta) / 2) * pow(2, (beta - 1) / 2)), 1 / beta);
	double sigma_v = 1;

	static std::normal_distribution<double> R(0, sigma_u);
	static std::normal_distribution<double> R1(0, sigma_v);

	for (auto& i : OldNestPop)
	{
		double step[NODE_MAX] = {};
		//前面的系数是保证最优鸟巢不会进行levy飞行
		for (int j = 0; j < nodeNum; j++)
		{
			step[j] = (i.Xn[j] - OldNestPop[bestNum].Xn[j])*R(e) * R(e) / (pow(abs(R1(e)), 1 / beta));
		}
		//按范围更新NODE
		for (int j = 0; j < nodeNum; j++)
		{
			if (i.Xn[j] + alpha * step[j] > machineNum - 1)
			{
				i.Xn[j] = machineNum - 1;
			}
			else if (i.Xn[j] + alpha * step[j] < 0) {
				i.Xn[j] = 0;
			}
			else
				i.Xn[j] = (int)(i.Xn[j] + alpha * step[j]);
		}

		fitFunc(i);
	}

	return OldNestPop;
}

std::vector<Nest> RandomAbandonPaNestPop(std::vector<Nest> OldNestPop)
{
	double step_sizeX = 0;
	double step_sizeY = 0;
	static std::uniform_int_distribution<int> randomInt(0, OldNestPop.size() - 1);
	for (decltype(OldNestPop.size()) i = 0; i < OldNestPop.size(); ++i)
	{
		if (u(e) < Pa)//被宿主发现了，要重新寻找新巢
		{
			double step[NODE_MAX] = {};
			for (int j = 0; j < nodeNum; j++)
			{
				step[j] = u(e) * (OldNestPop[randomInt(e)].Xn[j] - OldNestPop[randomInt(e)].Xn[j]);
			}
			for (int j = 0; j < nodeNum; j++)
			{
				if (OldNestPop[i].Xn[j] + step[j] > machineNum - 1)
				{
					OldNestPop[i].Xn[j] = machineNum - 1;
				}
				else if (OldNestPop[i].Xn[j] + step[j] < 0)
				{
					OldNestPop[i].Xn[j] = 0;
				}
				else
				{
					OldNestPop[i].Xn[j] += (int)step[j];
				}
			}
			fitFunc(OldNestPop[i]);
		}
	}

	return OldNestPop;
}



//产生整形的随机数
int my_i4_uniform_ab(int a, int b, int &seed)
{
	int c;
	const int i4_huge = 2147483647;
	int k;
	float r;
	int value;

	if (seed == 0)
	{
		std::cerr << "\n";
		std::cerr << "I4_UNIFORM_AB - Fatal error!\n";
		std::cerr << "  Input value of SEED = 0.\n";
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
	r = (1.0 - r) * ((float)a - 0.5)
		+ r * ((float)b + 0.5);
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


void myShowtime() {
	time_t now = time(nullptr);
	char buffer[100];
	struct tm* tm = localtime(&now);
	size_t len = strftime(buffer, 100, "%d %B %Y %I:%M:%S %p", tm);
//	std::cout << buffer << std::endl;
}