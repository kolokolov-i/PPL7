#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cfloat>
#include <streambuf>
#include <string>
#include <iomanip>
#include <chrono>
#include <vector>
#include <thread>

using namespace std;

typedef struct
{
	float x;
	float y;
} Point;

void readGraph();
void writeMST();
void release();
void buildMST();
void runBuildKey(bool* key, bool* mstSet);
void runMinKey(bool* key, bool* mstSet);
void buildKey(bool* key, bool* mstSet, int pfrom, int count);
void minKey(int& u, float& m, int& pj, bool* key, bool* mstSet, int pfrom, int count);
void writeLog();

const int REPEAT_COUNT = 10;
const int THREAD_COUNT = 3;
int* timeMeter;
thread threads[THREAD_COUNT];

int n; // vertex count
int pathCount;
float* graph;
int* result;
int* minKeys;
int* jj;
float* mm;
Point* points;

int main(int argc, char* argv[])
{
	int ierr;
	ierr = MPI_Init(&argc, &argv);
	if (ierr != MPI_SUCCESS)
		return ierr;

	readGraph();
	timeMeter = new int[REPEAT_COUNT];
	minKeys = new int[THREAD_COUNT];
	mm = new float[THREAD_COUNT];
	jj = new int[THREAD_COUNT];
	result = new int[n];
	chrono::time_point<chrono::system_clock> start, end;
	for (int w = 0; w < REPEAT_COUNT; w++)
	{
		start = chrono::system_clock::now();
		buildMST();
		end = chrono::system_clock::now();
		timeMeter[w] = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		cout << "Round " << w + 1 << " finished" << endl;
	}
	writeMST();
	writeLog();
	release();

	MPI_Finalize();
}


void readGraph()
{
	cout << "Reading graph ..." << endl;
	ifstream inf("graph.bin", ios::in | ios::binary);
	if (!inf)
	{
		cerr << "Graph reading error!" << endl;
		exit(1);
	}
	inf.read((char*)&n, sizeof(n));
	pathCount = (n * (n - 1)) / 2;
	float* paths = new float[pathCount];
	inf.read((char*)paths, pathCount * sizeof(float));
	points = new Point[n];
	inf.read((char*)points, n * sizeof(Point));
	inf.close();
	graph = new float[n * n];
	int p = 0;
	for (int i = 0; i < n; i++)
	{
		int j = 0;
		for (; j < i; j++)
		{
			graph[i * n + j] = graph[j * n + i];
		}
		j = i;
		graph[i * n + j] = 0;
		for (j++; j < n; j++)
		{
			graph[i * n + j] = paths[p++];
		}
	}
	delete[] paths;
}

void buildMST()
{
	bool* key = new bool[n];
	bool* mstSet = new bool[n];
	for (int i = 0; i < n; i++)
	{
		key[i] = false;
		mstSet[i] = false;
		result[i] = -1;
	}
	mstSet[0] = true;
	result[0] = -1;
	for (int i = 1; i < n; i++)
	{
		runBuildKey(key, mstSet);
		runMinKey(key, mstSet);
		int u = 0;
		int ji = 0;
		float mmKey = FLT_MAX;
		for (int j = 0; j < THREAD_COUNT; j++) {
			if (mm[j] < mmKey) {
				u = minKeys[j];
				mmKey = mm[j];
				ji = jj[j];
			}
		}
		result[u] = ji;
		mstSet[u] = true;
	}
	delete[] key;
	delete[] mstSet;
}

void runBuildKey(bool* key, bool* mstSet)
{
	int nominalCount = n / THREAD_COUNT;
	for (int startIndex = 0, i = 0; startIndex < n; startIndex += nominalCount, i++)
	{
		int count = nominalCount;
		if (startIndex + nominalCount < n && startIndex + nominalCount * 2 >= n)
		{
			count = n - startIndex;
		}
		threads[i] = thread(buildKey, key, mstSet, startIndex, count);
		if (count != nominalCount)
		{
			break;
		}
	}
	for (int i = 0; i < THREAD_COUNT; i++)
	{
		threads[i].join();
	}
}

void runMinKey(bool* key, bool* mstSet)
{
	int nominalCount = n / THREAD_COUNT;
	for (int startIndex = 0, i = 0; startIndex < n; startIndex += nominalCount, i++)
	{
		int count = nominalCount;
		if (startIndex + nominalCount < n && startIndex + nominalCount * 2 >= n)
		{
			count = n - startIndex;
		}
		threads[i] = thread(minKey, ref(minKeys[i]), ref(mm[i]), ref(jj[i]), key, mstSet, startIndex, count);
		if (count != nominalCount)
		{
			break;
		}
	}
	for (int i = 0; i < THREAD_COUNT; i++)
	{
		threads[i].join();
	}
}

void buildKey(bool* key, bool* mstSet, int pfrom, int count)
{
	int k = 0;
	for (int i = pfrom; k < count; i++, k++)
	{
		if (mstSet[i])
		{
			key[i] = false;
		}
		else
		{
			bool f = false;
			for (int j = 0; j < n; j++)
			{
				if (!mstSet[j])
				{
					continue;
				}
				if (graph[i * n + j] > 0)
				{
					f = true;
					break;
				}
			}
			key[i] = f;
		}
	}
}

void minKey(int& v, float& m, int& pj, bool* key, bool* mstSet, int pfrom, int count)
{
	float min = FLT_MAX;
	int minIndex = 0;
	int jetIndex = 0;
	int k = 0;
	for (int i = pfrom; k < count; i++, k++)
	{
		if (key[i])
		{
			for (int j = 0; j < n; j++)
			{
				if (mstSet[j])
				{
					float t = graph[i * n + j];
					if (t < min)
					{
						min = t;
						minIndex = i;
						jetIndex = j;
					}
				}
			}
		}
	}
	v = minIndex;
	m = min;
	pj = jetIndex;
}

void writeMST()
{
	cout << "Writing mst ..." << endl;
	ofstream outf("mst.bin", ios::out | ios::binary);
	if (!outf)
	{
		cerr << "MST writing error!" << endl;
		exit(1);
	}
	outf.write((char*)(&n), sizeof(n));
	outf.write((char*)result, n * sizeof(int));
	outf.close();
}

void writeLog()
{
	cout << "Writing time log ..." << endl;
	ofstream outf("log/time_par.txt", ios::out);
	if (!outf)
	{
		cerr << "Time log writing error!" << endl;
		exit(1);
	}
	int timeMin, timeMax, timeSum;
	timeMin = timeMeter[0];
	timeMax = timeMeter[0];
	timeSum = 0 + timeMeter[0];
	float timeAverage;
	for (int i = 1; i < REPEAT_COUNT; i++)
	{
		int t = timeMeter[i];
		if (t < timeMin)
		{
			timeMin = t;
		}
		if (t > timeMax)
		{
			timeMax = t;
		}
		timeSum += t;
	}
	timeAverage = ((float)timeSum) / REPEAT_COUNT;
	outf << "Time test \"Parallel\"" << endl;
	outf << "Threads: " << THREAD_COUNT << endl;
	outf << "Repeat: " << REPEAT_COUNT << endl;
	outf << "Time Average (ms): " << timeAverage << endl;
	outf << "Time Min (ms): " << timeMin << endl;
	outf << "Time Max (ms): " << timeMax << endl;
	for (int i = 0; i < REPEAT_COUNT; i++)
	{
		outf << "[" << i << "] " << timeMeter[i] << endl;
	}
	outf.close();
}

void release()
{
	delete[] graph;
	delete[] result;
	delete[] points;
	delete[] timeMeter;
	delete[] minKeys;
	delete[] mm;
	delete[] jj;
	cout << "Completed" << endl;
}