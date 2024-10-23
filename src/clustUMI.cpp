#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <queue>
#include <iomanip>
#include <cmath>
#include <map>
#include <filesystem>
#include <sstream>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include "circfullUtils.hpp"
#include "clust.hpp"

int nthread = 4;

std::map<std::string, int> seqUMIClust;
std::string *UMISeq;
std::string umiFastaString, strandFastqString, oPathString = ".", oPrefix = "circFL", seqkit = "seqkit";
std::filesystem::path umiFasta, umi2Fastq, strandFastq, oPath;

/// @brief generate consensus sequence from a set of sequences
/// @param name sequence id
/// @param seqs sequences set
/// @param quals quality scores set
/// @param output path to output consensus sequence
void consensusSeq(const std::string &name, const std::vector<std::string> *seqs, const std::vector<std::string> *quals, const std::filesystem::path &output)
{
	unsigned idx{0};
	std::string consensus;
	if (seqs->size() > 1)
		consensus = generateConsensusSeq(*seqs, *quals);
	else
		consensus = (*seqs)[0];
	/*std::string out = ">" + name + "\n";
	int lines = ceil(consensus.length() / 60.0);
	for (int i = 0; i < lines; i++)
		out += consensus.substr(i * 60, 60) + "\n";*/
	std::string out = "@" + name + "\n" + consensus + "\n+\n" + std::string(consensus.length(), 'F') + "\n";
	std::ofstream consensusStream{output, std::ios::app};
	consensusStream << out;
	consensusStream.close();
}

/// @brief read UMI sequence from umiFasta file
/// @param strandFasta path to input stranded fasta file
/// @param consensusFastq path to output consensus fastq file
/// @param nthread number of threads
void getConsensusSeq(const std::filesystem::path &strandFasta, const std::filesystem::path &consensusFastq, const int nthread)
{
	// read sequence and generate consensus seq for each UMI clust
	boost::asio::thread_pool pool(nthread);
	seqan3::sequence_file_input seqIn{strandFasta};
	int nowClust = -1, nowt, clustId;
	std::vector<std::string> *seqs = new std::vector<std::string>();
	std::vector<std::string> *quals = new std::vector<std::string>();
	for (auto &record : seqIn)
	{
		clustId = atoi(circfull::getFirstString(record.id(), '_').c_str() + 1);
		if (nowClust != -1 && clustId != nowClust)
		{
			boost::asio::post(pool,
							  [=]
							  {
								  consensusSeq("clust" + std::to_string(nowClust), seqs, quals, consensusFastq);
								  delete seqs, quals;
							  });
			seqs = new std::vector<std::string>;
			quals = new std::vector<std::string>;
			nowClust = clustId;
		}
		if (nowClust == -1)
			nowClust = clustId;
		std::string seq = circfull::seqan3ToString(record.sequence());
		std::string qual = circfull::seqan3ToString(record.base_qualities());
		seqs->push_back(seq);
		quals->push_back(qual);
	}
	boost::asio::post(pool,
					  [=]
					  {
						  consensusSeq("clust" + std::to_string(nowClust), seqs, quals, consensusFastq);
					  });
	pool.join();
}

/// @brief add UMI cluster information into stranded fastq file
/// @param clustInfo map of UMI cluster information
/// @param strandFastq path to input stranded fastq file
/// @param outputFastq path to output fastq file
void addClustInfo2StrandFile(const std::map<std::string, int> &clustInfo, const std::filesystem::path &strandFastq, const std::filesystem::path &outputFastq)
{
	std::ifstream inputStream{strandFastq};
	std::ofstream outputStream{outputFastq};
	std::string name, seq, tmp, qual, out, id;
	while (inputStream.peek() != EOF)
	{
		getline(inputStream, name);
		getline(inputStream, seq);
		getline(inputStream, tmp);
		getline(inputStream, qual);
		id = circfull::getFirstString(name);
		id.erase(0, 1);
		auto iter = clustInfo.find(id);
		if (iter != clustInfo.end())
		{
			out = "@c" + std::to_string(iter->second) + "_" + id + "\n" + seq + "\n+\n" + qual + "\n";
			outputStream << out;
		}
	}
	outputStream.close();
}

namespace circfull
{
	int circfull_umi_clust(CircfullOption &opt)
	{
		printTimeInfo("clustering UMI sequence:");
		int nUMIClust, nStrandClust;
		// std::vector<tuple<int, string>> clustInfo;
		std::map<std::string, int> clustInfo;
		nUMIClust = clustUMISeq(opt.umiFasta, opt.oPath / (opt.oPrefix + "_umi"), 8, opt.nthread, 0.98, clustInfo);
		if (opt.clustSeq)
			nStrandClust = clustUMISeq(opt.strandFastq, opt.oPath / (opt.oPrefix + "_seq"), 10, opt.nthread, 0.95, clustInfo);
		std::cout << "[RESULT]UMI cluster: " << nUMIClust << std::endl;
		if (opt.clustSeq)
			std::cout << "[RESULT]sequence cluster: " << nStrandClust << std::endl;
		/*printTimeInfo("adding information into raw fastq:");
		std::filesystem::path strandFastq = opt.oPath / (opt.oPrefix + "_strand_tmp.fastq");
		addClustInfo2StrandFile(clustInfo, opt.strandFastq, strandFastq);
		std::filesystem::path strandSortedFastq = opt.oPath / (opt.oPrefix + "_strand_sorted.fastq");
		if (sortFastqById(opt.seqkit, strandFastq, strandSortedFastq, opt.nthread))
			exit(-1);
		printTimeInfo("UMI guided generating consensus sequences:");
		/*if (convertFastq2Fasta(opt.seqkit, opt.strandFastq, strandFasta, opt.nthread))
			exit(-1);*/
		// getConsensusSeq(strandSortedFastq, opt.consensusFasta, opt.nthread);
		printTimeInfo("finished.");
		return 0;
	}

}
