#include <iostream>
#include <thread>
#include <filesystem>
#include <fstream>
#include <cstring>
#include <tuple>
#include <optional>
#include <map>
#include <seqan3/io/all.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <boost/timer/progress_display.hpp>
#include <queue>
#include <mutex>
#include "cccs/cccs.hpp"
#include "circCall.hpp"
#include "circfullUtils.hpp"
#include "clust.hpp"

using namespace std;
namespace fs = filesystem;

/*template<typename T>
string seqan3ToString(T &&value) {
	string ret;
	for (auto && c : value) {
		ret += seqan3::to_char(c);
	}
	return ret;
}*/

/// @brief find ccs reads from a single read
/// @param header read id
/// @param seq read sequence
/// @param prefix output prefix
/// @return segment, ccs sequence
std::optional<std::tuple<std::string, std::string>> find_ccs(const std::string &header, const std::string &seq, const std::filesystem::path prefix)
{
	if (seq.length() < 100)
		return nullopt;
	// cout << header << seq << endl;
	cccs::CcsResult *res = cccs::ccs(reinterpret_cast<const unsigned char *>(seq.c_str()), seq.length());
	if (res->segment_ptr && res->ccs_seq_ptr)
	{
		std::string seg = std::string(res->segment_ptr); // segment is left open and right closed range
		std::string ccs = std::string(res->ccs_seq_ptr);
		cccs::free_ccs_result(res);
		std::filesystem::path ccsFasta = prefix, trimmedFasta = prefix;
		ccsFasta += ".ccs.fa";
		trimmedFasta += ".raw.fa";
		std::ofstream ccsStream(ccsFasta, ios::app);
		ccsStream << ">" + header + "\t" + seg + "\t" + std::to_string(ccs.length()) + "\n" + ccs + "\n";
		ccsStream.close();
		std::ofstream trimmedStream(trimmedFasta, std::ios::app);
		trimmedStream << ">" + header + "\n" + seq + "\n";
		return make_tuple(seg, ccs);
	}
	cccs::free_ccs_result(res);
	return nullopt;
}

std::map<std::string, std::tuple<std::string, std::string, std::string>> circfull::findCcsReads(std::filesystem::path inFile, std::filesystem::path outPrefix, int nthread)
{
	// readin from inFile and scan for ccs reads
	// write ccs reads to outPrefix_ccs.fa using find_ccs function
	// return map<string, tuple<string, string, string>>: readId -> (segment, ccs, rawSeq)
	std::map<std::string, std::tuple<std::string, std::string, std::string>> ret;
	seqan3::sequence_file_input seqIn{inFile};
	boost::timer::progress_display dp(circfull::getFileNrow(inFile, nthread) / 4);
	boost::asio::thread_pool pool(nthread);
	std::mutex mt;
	for (auto &record : seqIn)
	{
		std::string shortID = circfull::getId(record.id());
		std::string tmpSeq = circfull::seqan3ToString(record.sequence());
		// multiple threads
		boost::asio::post(pool, [shortID, tmpSeq, outPrefix, &ret, &dp, &mt]()
						  {
			auto retCcs { find_ccs(shortID, tmpSeq, outPrefix) };
			if (retCcs.has_value()) {
				auto [segment, ccs] = retCcs.value();
				ret[shortID] = make_tuple(segment, ccs, tmpSeq);
			}
			mt.lock();
			dp += 1;
			mt.unlock(); });
	}
	pool.join();
	std::cout << std::endl;
	return ret;
}

/// @brief get earch segment from raw sequence using segment string information, return a vector of segments
/// @param segment segment string
/// @param rawSeq raw sequence
/// @return vector<string> segments
std::vector<std::string> getSegments(const std::string &segment, const std::string &rawSeq)
{
	std::vector<std::string> ret;
	std::stringstream ss(segment);
	std::string tmp;
	// segments are left open and right closed range, 0-based
	// segments string format as "start1-end1;start2-end2;...startn-endn"
	while (getline(ss, tmp, ';'))
	{
		std::stringstream ss2(tmp);
		std::string start, end;
		std::getline(ss2, start, '-');
		std::getline(ss2, end, '-');
		ret.push_back(rawSeq.substr(stoi(start), stoi(end) - stoi(start)));
	}
	return ret;
}

std::map<std::string, std::string> circfull::constructCCSConsensus(std::map<std::string, std::tuple<std::string, std::string, std::string>> &ccsSeq, const std::map<std::string, std::string> &umiInfo, int nthread, std::filesystem::path outPrefix)
{
	// construct consensus sequence from ccsSeq
	// write consensus sequence to outPrefix_cons.fa
	// return map<string, string> consensus sequence: clustId -> consensus sequence
	// sort ccsSeq by umiInfo (get<2>)
	struct umiSort
	{
		std::string umi, segment, rawSeq;
		bool operator<(const umiSort &other) const
		{
			return umi < other.umi;
		}
	};
	priority_queue<umiSort> ccsSeqQueue;
	for (auto &ccs : ccsSeq)
	{
		auto umiInfoIter = umiInfo.find(ccs.first);
		if (umiInfoIter != umiInfo.end())
			ccsSeqQueue.push({umiInfoIter->second, get<0>(ccs.second), get<2>(ccs.second)});
	}
	std::filesystem::path ccsFasta = outPrefix;
	ccsFasta += "_cons.ccs.fa";
	std::ofstream ccsStream(ccsFasta);
	std::map<std::string, std::string> ret;
	boost::asio::thread_pool pool(nthread);
	std::string lastUmi = ccsSeqQueue.top().umi;
	std::vector<std::tuple<std::string, std::string>> ccsSeqVec;
	// construct consensus sequence for each UMI cluster
	boost::timer::progress_display dp(ccsSeqQueue.size());
	mutex mt;
	while (!ccsSeqQueue.empty())
	{
		// auto ccs = ccsSeqQueue.top();
		auto [umi, segment, rawSeq] = ccsSeqQueue.top();
		ccsSeqQueue.pop();
		if (umi != lastUmi)
		{
			// construct consensus sequence for lastUmi
			boost::asio::post(pool, [lastUmi, ccsSeqVec, outPrefix, &ret, &ccsStream, &dp, &mt]()
							  {
				// get segments from ccsSeqVec
				std::vector<std::string> segments;
				for (auto &ccs : ccsSeqVec) {
					auto [segment, rawSeq] = ccs;
					std::vector<std::string> tmp = getSegments(segment, rawSeq);
					segments.insert(segments.end(), tmp.begin(), tmp.end());
				}
				// construct consensus sequence
				std::string consensus = generateConsensusSeq(segments);
				// write consensus sequence to file
				mt.lock();
				ccsStream << ">" + lastUmi + "\n" + consensus + "\n";
				ret[lastUmi] = consensus;
				dp += ccsSeqVec.size();
				mt.unlock(); });
			// clear ccsSeqVec
			ccsSeqVec.clear();
			lastUmi = umi;
		}
		ccsSeqVec.push_back({segment, rawSeq});
	}
	// post last bundle
	if (!ccsSeqVec.empty())
	{
		// construct consensus sequence for lastUmi
		boost::asio::post(pool, [lastUmi, ccsSeqVec, outPrefix, &ret, &ccsStream]()
						  {
			// get segments from ccsSeqVec
			std::vector<std::string> segments;
			for (auto &ccs : ccsSeqVec) {
				auto [segment, rawSeq] = ccs;
				std::vector<std::string> tmp = getSegments(segment, rawSeq);
				segments.insert(segments.end(), tmp.begin(), tmp.end());
			}
			// construct consensus sequence
			std::string consensus = generateConsensusSeq(segments);
			// write consensus sequence to file
			ccsStream << ">" + lastUmi + "\n" + consensus + "\n";
			ret[lastUmi] = consensus; });
	}
	pool.join();
	std::cout << std::endl;
	return ret;
}