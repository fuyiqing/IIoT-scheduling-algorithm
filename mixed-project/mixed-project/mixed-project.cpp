// mixed-project.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include <iostream>
#include "Genetic.h"
#include "MyGraph.h"
#include "libxl.h"
#include <time.h>
#include<cmath>
#include<string>
#include<string.h>
#include<fstream>

using namespace libxl;

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

int main()
{
	/*
	//输入genetic对应的适应度
	sheetGen->writeStr(1, 0, L"节点个数");
	for (int nodeG = 5; nodeG <= 50; nodeG++)
		sheetGen->writeNum(1, nodeG - 4, nodeG);
	for (int iiG = 1; iiG <= 350; iiG++)
		sheetGen->writeNum(iiG + 1, 0, iiG);

	//输入cuckoo对应的适应度
	sheetCuc->writeStr(1,0,L"节点个数");
	for (int i = 1; i <= 350; i++)
		sheetCuc->writeNum(i + 1, 0, i);
	for (int i = 5; i <= 50; i++)
		sheetCuc->writeNum(1, i - 4, i);

	*/

	//---------------------------------------主要代码-----------------------------------------------
	
	for (int diedai=5;diedai<=50;diedai++)
	{
		taskMachineGraphFunction(diedai);
	}


















	//-------------------------------------主要代码-----------------------------------------------------



	/*

	//保存并输出BookGen
	bookGenetic->save(L"genetic_5-50");
	bookGenetic->release();

	//保存并输出BookCuc
	bookCuckoo->save(L"cuckoo_5-50");
	bookCuckoo->release();*/
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
	/*ofstream output;
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

	output.close();*/
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

	showTime();
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
	sheetGen->writeStr(358, diedai - 4, L"总共用时(毫秒)：");
	sheetGen->writeNum(359, diedai - 4, duration);
	//printf("\n总共用时%lf毫秒\n", duration);
}
