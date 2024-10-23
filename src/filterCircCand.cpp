#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <filesystem>
#include <mutex>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include "circfullUtils.hpp"

/// @brief read UMI cluster information from umiClustTab file
/// @param umiFile path to umiClustTab file
/// @return map<string, string> umiInfo: readId -> clustId
std::map<std::string, std::string> readUMIClustInfo(const std::filesystem::path &umiFile)
{
	std::ifstream umiClustFile(umiFile);
	std::map<std::string, std::string> umiInfo;
	// read umi cluster info from umiClustFile: clust_id, read_num, length, read_id, compare_rate, separated by tab
	std::string umiClustId, readNum, length, readId, compareRate;
	std::string line;
	while (getline(umiClustFile, line))
	{
		std::istringstream iss(line);
		iss >> umiClustId >> readNum >> length >> readId >> compareRate;
		umiInfo[readId] = umiClustId;
	}
	return umiInfo;
}

/// @brief struct to store circRNA candidate information
typedef struct
{
	std::string readId, clustId, circId, strand, isoformId, spliceSite;
} CircCand;

/// @brief read circRNA candidate information from circCandFile
/// @param candFile path to circCandFile
/// @return vector<CircCand> circCandInfo: circRNA candidate information
std::vector<CircCand> readCircCandInfo(const std::filesystem::path &candFile, std::map<std::string, std::string> &umiInfo)
{
	std::ifstream circCandFile(candFile);
	std::vector<CircCand> circCandInfo;
	// two line of .fa file for each circRNA candidate, line1: >read_id\tcirc_id\tstrand\tisoform_id\tsplice_site\tCCS_position, line2: sequence
	for (std::string line; getline(circCandFile, line); getline(circCandFile, line))
	{
		std::string readId, clustId, circId, strand, isoformId, spliceSite, CCSPosition;
		std::istringstream iss(line);
		getline(iss, readId, '>');
		iss >> readId >> circId >> strand >> isoformId >> spliceSite >> CCSPosition;
		clustId = umiInfo[readId];
		circCandInfo.push_back({readId, clustId, circId, strand, isoformId, spliceSite});
	}

	sort(circCandInfo.begin(), circCandInfo.end(), [](const auto &a, const auto &b)
		 { return a.clustId < b.clustId; });
	return circCandInfo;
}

namespace circfull
{
	int filterCircCand(CircfullOption &opt)
	{
		printTimeInfo("Start filtering circRNA candidates:");

		// read UMI cluster result from umiClustTab file, and circ_call results from circCandFile file
		printTimeInfo(std::string("Read UMI info from ") + opt.umiClustTab.string());
		std::map<std::string, std::string> umiInfo = readUMIClustInfo(opt.umiClustTab);
		std::cout << "UMI cluster number: " << umiInfo.size() << std::endl;
		printTimeInfo(std::string("Read circRNA candidate info from ") + opt.rawCircRes.string());
		std::vector<CircCand> circCandInfo = readCircCandInfo(opt.rawCircRes, umiInfo);
		std::cout << "circRNA candidate number: " << circCandInfo.size() << std::endl;

		// process each UMI cluster and its circRNA candidates
		// for each UMI cluster, select the circRNA candidate with the highest read number
		boost::asio::thread_pool pool(opt.nthread);
		std::mutex mt;
		std::string nowClustId = circCandInfo.begin()->clustId;
		std::vector<CircCand> nowClustCand;
		std::vector<CircCand> circCand;

		// correct within each UMI cluster
		printTimeInfo("Correct according to UMI");
		for (auto &iter : circCandInfo)
		{
			if (iter.clustId == nowClustId)
			{
				nowClustCand.push_back(iter);
			}
			else
			{
				// select the circRNA candidate with the highest read number
				boost::asio::post(pool, [&, nowClustCand]()
								  {
					int maxReadNum = 0;
					CircCand maxReadCand;
					std::map<std::pair<std::string, std::string>, int> isoformReadNum;
					for (auto &cand : nowClustCand) {
						if (isoformReadNum.find({cand.circId, cand.isoformId}) == isoformReadNum.end()) {
							isoformReadNum[{cand.circId, cand.isoformId}] = 1;
						} else {
							isoformReadNum[{cand.circId, cand.isoformId}]++;
						}
						if (isoformReadNum[{cand.circId, cand.isoformId}] > maxReadNum) {
							maxReadNum = isoformReadNum[{cand.circId, cand.isoformId}];
							maxReadCand = cand;
						}
					}
					// add maxReadCand to circCand
					mt.lock();
					circCand.push_back(maxReadCand);
					mt.unlock(); });
				// update nowClustId and nowClustCand
				nowClustId = iter.clustId;
				nowClustCand.clear();
				nowClustCand.push_back(iter);
			}
		}
		pool.join();

		// filter circRNA candidates with at least 2 UMI cluster
		printTimeInfo("Filtering circRNA candidates with at least 2 UMI cluster");
		std::map<std::pair<std::string, std::string>, int> circReadNum;
		std::map<std::string, int> bsjReadNum;
		for (auto &cand : circCand)
		{
			if (circReadNum.find({cand.circId, cand.isoformId}) == circReadNum.end())
			{
				circReadNum[{cand.circId, cand.isoformId}] = 1;
			}
			else
			{
				circReadNum[{cand.circId, cand.isoformId}]++;
			}
			if (bsjReadNum.find(cand.circId) == bsjReadNum.end())
			{
				bsjReadNum[cand.circId] = 1;
			}
			else
			{
				bsjReadNum[cand.circId]++;
			}
		}
		// output from circReadNum
		std::ofstream circCandOut(opt.circCandFile);
		for (auto &iter : circReadNum)
		{
			if (iter.second >= opt.circNumThres)
			{
				circCandOut << iter.first.first << "\t" << iter.first.second << "\t" << iter.second << std::endl;
			}
		}
		circCandOut.close();
		std::ofstream bsjCandOut(std::filesystem::path(opt.circCandFile) += ".lowconf.tsv");
		for (auto &iter : bsjReadNum)
		{
			if (iter.second >= opt.circNumThres)
			{
				bsjCandOut << iter.first << "\t" << iter.second << std::endl;
			}
		}
		bsjCandOut.close();

		return 0;
	}
}