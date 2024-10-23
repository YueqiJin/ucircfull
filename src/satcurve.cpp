#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <thread>
#include <mutex>
#include <string>
#include <cstring>
#include <map>
#include <thread>
#include <random>
#include <queue>
#include <algorithm>
#include <filesystem>
#include <boost/timer/progress_display.hpp>
#include "circfullUtils.hpp"
#include "version.hpp"
#include "argparse/argparse.hpp"

int cntCirc = 0, cntRead = 0;
std::string *circID;
std::map<std::string, std::string> readIDcirc;
std::map<std::string, std::string> readIDumi;
std::set<std::string> circIDs;
std::string *readID;
std::map<std::string, int> readLen;
int *seqLen;
std::queue<std::pair<std::pair<int, unsigned long long int>, int>> res;
std::mutex mt;
int interval, rep, nthread;
bool quiet, help, circfullOutput, haveUmi = false;
std::string oPathString, fastqString, circfullFileNameString, readsInfoFileString, umiFileString, seqkit;

void readinCircfull(std::filesystem::path tableFile, std::filesystem::path seqInfoFile)
{
	cntCirc = circfull::getFileNrow(tableFile);
	circID = new std::string[cntCirc];
	std::ifstream tableIn(tableFile, std::ios::in);
	int tmpCnt;
	std::string tmpReadID;
	for (int i = 0; i < cntCirc; i++)
	{
		tableIn >> circID[i];
		tableIn >> tmpCnt;
		for (int j = 0; j < tmpCnt; j++)
		{
			tableIn >> tmpReadID;
			readIDcirc[tmpReadID] = circID[i];
		}
	}
	tableIn.close();
	cntRead = circfull::getFileNrow(seqInfoFile);
	readID = new std::string[cntRead];
	seqLen = new int[cntRead];
	std::ifstream seqIn(seqInfoFile, std::ios::in);
	for (int i = 0; i < cntRead; i++)
	{
		seqIn >> readID[i];
		seqIn >> seqLen[i];
	}
	seqIn.close();
}

int countCircCnt(std::filesystem::path tableFile)
{
	FILE *stream;
	std::string cmdLine = "cut -f 1 " + tableFile.string() + " | sort | uniq | wc -l";
	std::string result;
	char buf[128];
	int n = -1;
	if ((stream = popen(cmdLine.c_str(), "r")) != NULL)
	{
		while (fgets(buf, 128, stream) != NULL)
		{
			result += buf;
		}
		sscanf(result.c_str(), "%d", &n);
		pclose(stream);
		return n;
	}
	else
	{
		return -1;
	}
}

void readUmi(std::filesystem::path umiFile)
{
	std::ifstream umiIn(umiFile, std::ios::in);
	std::string tmpReadID, tmpUmi;
	while (!umiIn.eof())
	{
		umiIn >> tmpUmi >> tmpReadID;
		readIDumi[tmpReadID] = tmpUmi;
	}
	umiIn.close();
}

void readin(std::filesystem::path tableFile, std::filesystem::path seqInfoFile)
{
	cntCirc = countCircCnt(tableFile);
	circID = new std::string[cntCirc];
	std::ifstream tableIn(tableFile, std::ios::in);
	std::string tmpCircID, tmpReadID;
	int circCnt = 0;
	while (!tableIn.eof())
	{
		tableIn >> tmpReadID >> tmpCircID;
		if (circIDs.count(tmpCircID) == 0)
		{
			circID[circCnt++] = tmpCircID;
			circIDs.insert(tmpCircID);
		}
		readIDcirc[tmpReadID] = tmpCircID;
	}
	tableIn.close();
	cntRead = circfull::getFileNrow(seqInfoFile);
	readID = new std::string[cntRead];
	seqLen = new int[cntRead];
	std::ifstream seqIn(seqInfoFile, std::ios::in);
	for (int i = 0; i < cntRead; i++)
	{
		seqIn >> readID[i];
		seqIn >> seqLen[i];
	}
	seqIn.close();
}

int parseTable(std::filesystem::path tableFile, std::filesystem::path outFile)
{
	std::string cmdLine = "cut -f 1,17,18 " + tableFile.string() + " | sed -e '1d' -e 's/,/ /g' > " + outFile.string();
	return system(cmdLine.c_str());
}

int parseFastq(std::string seqkit, std::filesystem::path fastqFile, std::filesystem::path outFile)
{
	std::string cmdLine = seqkit + " fx2tab --name --only-id -l " + fastqFile.string() + " > " + outFile.string();
	return system(cmdLine.c_str());
}

std::vector<int> *sampleSet;
std::vector<int> sampleRange;
std::mt19937 gen(std::random_device{}());
std::set<std::string> *mappedCirc;

void constructSampleSet(int n_thread)
{
	sampleSet = new std::vector<int>[n_thread];
	mappedCirc = new std::set<std::string>[n_thread];
	for (int i = 0; i < n_thread; i++)
	{
		sampleSet[i].resize(cntRead);
	}
	sampleRange.resize(cntRead);
	std::iota(sampleRange.begin(), sampleRange.end(), 0);
}

void sampling(int threadID, int sampleSize)
{
	// shuffle(sampleSet[threadID].begin(), sampleSet[threadID].end(), gen);
	sampleSet[threadID].clear();
	sample(sampleRange.begin(), sampleRange.end(), std::back_inserter(sampleSet[threadID]), sampleSize, gen);
}

void calculate(int threadID, int sampleSize)
{
	sampling(threadID, sampleSize);
	mappedCirc[threadID].clear();
	unsigned long long int lenSum = 0;
	std::string tmpReadID;
	for (int i = 0; i < sampleSize; i++)
	{
		lenSum += (unsigned long long int)seqLen[sampleSet[threadID][i]];
		tmpReadID = readID[sampleSet[threadID][i]];
		if (haveUmi)
			tmpReadID = readIDumi[tmpReadID];
		if (readIDcirc.count(tmpReadID) != 0)
			mappedCirc[threadID].insert(readIDcirc[tmpReadID]);
	}
	mt.lock();
	res.push(std::make_pair(std::make_pair(sampleSize, lenSum), mappedCirc[threadID].size()));
	mt.unlock();
	mappedCirc[threadID].clear();
}

void parseArgs(int argc, char **argv)
{
	argparse::ArgumentParser program{"circfull_satcurve", PROJECT_VERSION, argparse::default_arguments::all};
	program.add_argument("-i", "--interval")
		.metavar("INT")
		.required()
		.default_value(1000)
		.help("interval of sampling size.")
		.scan<'i', int>();
	program.add_argument("-r", "--repeat")
		.metavar("INT")
		.required()
		.default_value(10)
		.help("repeat of sampling.")
		.scan<'i', int>();
	program.add_argument("-u", "--umi")
		.metavar("FILE")
		.default_value("-")
		.help("umi clustr file (for qcircFL).");
	program.add_argument("-o", "--outdir")
		.metavar("DIR")
		.required()
		.default_value(std::string("."))
		.help("output directory.");
	program.add_argument("-t", "--thread")
		.metavar("INT")
		.required()
		.default_value(4)
		.help("number of threads used.")
		.scan<'i', int>();
	program.add_argument("--seqkit")
		.metavar("PATH")
		.default_value(std::string("seqkit"))
		.help("path to seqkit.");
	program.add_argument("--circReads")
		.metavar("PATH")
		.required()
		.help("reads id with supporting circRNA (the first two columns of cirilong output prefix.reads).");
	program.add_argument("--fastq")
		.metavar("PATH")
		.required()
		.help("fastq file.");
	try
	{
		program.parse_args(argc, argv);
		interval = program.get<int>("-i");
		rep = program.get<int>("-r");
		nthread = program.get<int>("-t");
		oPathString = program.get<std::string>("-o");
		readsInfoFileString = program.get<std::string>("--circReads");
		umiFileString = program.get<std::string>("-u");
		if (umiFileString != "-")
		{
			haveUmi = true;
		}
		fastqString = program.get<std::string>("--fastq");
		seqkit = program.get<std::string>("--seqkit");
	}
	catch (const std::runtime_error &err)
	{
		std::cerr << err.what() << std::endl;
		std::cerr << program;
		std::exit(1);
	}
}

int main(int argc, char **argv)
{
	// argument setting
	parseArgs(argc, argv);
	std::filesystem::path readsInfoFile(readsInfoFileString), fastq(fastqString), umiFile(umiFileString), oPath(oPathString);
	// pre parsing circFL output and fastq
	std::filesystem::path tableOutputFile = oPath / "circInfo.tsv";
	// if (circfullOutput) {
	//	if (parseTable(circfullFileName, tableOutputFile) == -1) {
	//		std::cout << "ERROR: Failed when parsing circFL output." << endl;
	//		return -1;
	//	}
	// }
	std::filesystem::path seqkitOutputFile = oPath / "seqInfo.tsv";
	if (parseFastq(seqkit, fastq, seqkitOutputFile) == -1)
	{
		std::cout << "ERROR: Failed when summarize reads info." << std::endl;
		return -1;
	}
	// reading all info
	// if (circfullOutput)
	//	readinCircfull(tableOutputFile, seqkitOutputFile);
	// else
	if (haveUmi)
		readUmi(umiFile);
	readin(readsInfoFile, seqkitOutputFile);
	std::cout << "circRNA info parsed: " << cntCirc << std::endl
			  << "Sequencing read parsed: " << cntRead << std::endl;
	// construct sampling set for each thread
	constructSampleSet(nthread);
	// construct workling queue
	std::queue<int> workQueue;
	unsigned long long int totalRead = 0;
	for (int i = interval; i < cntRead; i += interval)
		for (int j = 0; j < rep; j++)
		{
			workQueue.push(i);
			totalRead += (unsigned long long int)i;
		}
	// calculating all data
	std::thread *th = new std::thread[nthread];
	int tmp;
	unsigned long long int finished = 0;
	int nowt = 0;
	unsigned long long int nowtmp = 0;
	std::cout << "starting calculating saturation curve." << std::endl;
	// if have -q option, do not show progress indicator
	boost::timer::progress_display *bar;
	if (!quiet)
	{
		bar = new boost::timer::progress_display(totalRead);
	}
	while (!workQueue.empty())
	{
		nowt = 0;
		finished = 0;
		for (int i = 0; i < nthread; i++)
		{
			if (workQueue.empty())
			{
				break;
			}
			tmp = workQueue.front();
			finished += (unsigned long long int)tmp;
			th[i] = std::thread(calculate, i, tmp);
			workQueue.pop();
			nowt++;
		}
		for (int i = 0; i < nowt; i++)
			th[i].join();
		if (!quiet)
		{
			(*bar) += finished;
		}
	}
	// output data
	std::filesystem::path outputFile = oPath / "satcurve.tsv";
	std::ofstream outFile(outputFile);
	std::pair<std::pair<int, unsigned long long int>, int> tmpOut;
	while (!res.empty())
	{
		tmpOut = res.front();
		res.pop();
		outFile << tmpOut.first.first << "\t" << tmpOut.first.second << "\t" << tmpOut.second << std::endl;
	}
	outFile.close();
	return 0;
}
