#include <iostream>
#include <string>
#include <filesystem>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <numeric>
#include <mutex>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <boost/timer/progress_display.hpp>
#include <seqan/seq_io.h>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include "circCall.hpp"
#include "align.hpp"
#include "circfullUtils.hpp"

namespace circfull
{
	namespace RG
	{
		std::map<std::string, int> ids;
		seqan3::dna5_vector *seqs;
		BSMotifIndex BSMotifList;

		MapRecord parseExon(std::filesystem::path output, referenceStorage &genome, const seqan3Record &record)
		{
			int refBegin = record.refBegin;
			bool head = true;
			int currentPos = refBegin, currentExonStart = refBegin, refEnd = refBegin, queryBegin = 0, queryEnd = 0;
			std::string exonStart = std::to_string(refBegin), exonEnd = ""; //, exonLen = "";
			unsigned long long int exonLen = 0;
			for (auto [cigarCount, cigarOperation] : record.cigarSeq)
			{
				if (cigarOperation.to_char() == 'M' || cigarOperation.to_char() == 'I' || cigarOperation.to_char() == '=' || cigarOperation.to_char() == 'X')
					queryEnd += cigarCount;
				if (head && (cigarOperation.to_char() == 'S' || cigarOperation.to_char() == 'H'))
					queryBegin += cigarCount, queryEnd += cigarCount;
				else
					head = false;
				switch (cigarOperation.to_char())
				{
				case 'M':
				case 'D':
					currentPos += cigarCount;
					refEnd += cigarCount;
					break;
				case 'N':
					// exonLen += std::to_string(currentPos - currentExonStart + 1) + ",";
					exonLen += currentPos - currentExonStart;
					exonEnd += std::to_string(currentPos) + ",";
					currentPos += cigarCount;
					refEnd += cigarCount;
					exonStart += "," + std::to_string(currentPos);
					currentExonStart = currentPos;
				}
			}
			exonEnd += std::to_string(currentPos);
			exonLen += currentPos - currentExonStart;
			std::string leftSeq, rightSeq;
			leftSeq = genome.getSeqByPosition(record.refName, refBegin - 2, refBegin);
			rightSeq = genome.getSeqByPosition(record.refName, refEnd, refEnd + 2);
			if (circfull::debug)
			{
				std::ofstream fileOut(output, std::ios::app);
				fileOut << record.queryId + "\t" + record.refName + "\t" +
							   std::to_string(refBegin) + "\t" + std::to_string(refEnd) + "\t" +
							   std::to_string(queryBegin) + "\t" + std::to_string(queryEnd) + "\t" +
							   exonStart + "\t" + exonEnd + "\t" + std::to_string(exonLen) + "\t" + (static_cast<bool>(record.flag & seqan3::sam_flag::on_reverse_strand) ? "-" : "+") + "\t" +
							   leftSeq + "\t" + rightSeq + "\n";
				fileOut.close();
			}
			return MapRecord{record.queryId, record.refName, refBegin, refEnd, queryBegin, queryEnd, exonStart, exonEnd, exonLen, static_cast<bool>(record.flag & seqan3::sam_flag::on_reverse_strand) ? '-' : '+', leftSeq, rightSeq};
		}

		std::vector<MapRecord> parseBam(referenceStorage &genome, std::filesystem::path output, std::filesystem::path bam, int nthread)
		{
			std::vector<MapRecord> ret;
			if (checkFile(output))
				delFile(&output);
			std::ifstream mapBamStream(bam);
			seqan3::sam_file_input mapBam{
				mapBamStream,
				seqan3::format_bam{},
				seqan3::fields<seqan3::field::header_ptr, seqan3::field::id, seqan3::field::flag, seqan3::field::cigar, seqan3::field::ref_id, seqan3::field::ref_offset>{}};
			auto refName = &mapBam.header().ref_ids();
			boost::asio::thread_pool tp(nthread);
			std::mutex mt;
			for (auto &&record : mapBam)
			{
				if (static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped) || static_cast<bool>(record.flag() & seqan3::sam_flag::secondary_alignment))
					continue;
				boost::asio::post(tp, [&output, &genome, &refName, &ret, &mt, record]
								  {
					MapRecord res = parseExon(output, genome, seqan3Record(record.reference_position().value(), record.id(), (*refName)[record.reference_id().value()],	record.cigar_sequence(), record.flag()));
					mt.lock();
					ret.push_back(res);
					mt.unlock(); });
			}
			tp.join();
			return ret;
		}

		std::string getRefSeq(const MapRecord &record, referenceStorage &genome)
		{
			std::string tmpSeq;
			std::string returnSeq = "";
			long int exonNum{count(record.exonStart.begin(), record.exonStart.end(), ',') + 1};
			std::istringstream issExonStart{record.exonStart}, issExonEnd{record.exonEnd};
			int currentExonStart, currentExonEnd;
			std::string tmp;
			while (exonNum)
			{
				std::getline(issExonStart, tmp, ',');
				currentExonStart = stoi(tmp);
				std::getline(issExonEnd, tmp, ',');
				currentExonEnd = stoi(tmp);
				returnSeq += genome.getSeqByPosition(record.chr, currentExonStart, currentExonEnd);
				exonNum--;
			}
			return returnSeq;
		}

		std::tuple<bool, int> readRemapRef(MapRecordIter begin, MapRecordIter end, referenceStorage &genome)
		{
			std::string thisMapRegion;
			std::string readSeq = seqan3ToString(seqs[ids[begin->queryId]]);
			for (MapRecordIter iter = begin; iter != end; iter++)
			{
				thisMapRegion = getRefSeq(*iter, genome);
				int nReg, addTmp;
				std::vector<std::pair<int, int>> alignPos;
				mpAligner remap{thisMapRegion};
				mm_reg1_t *mmRecords = remap.map(readSeq, &nReg);
				// calRead start
				for (int i = 0; i < nReg; i++)
				{
					mm_reg1_t *hit = &mmRecords[i];
					assert(hit->p);
					addTmp = MM_CIGAR_STR[hit->p->cigar[0] & 0xf] == 'H' ? (hit->p->cigar[0] >> 4) : 0;
					alignPos.push_back(std::make_pair(hit->qs + addTmp, hit->qe + addTmp));
					free(hit->p);
				}
				free(mmRecords);
				if (alignPos.size() <= 1)
					continue;
				std::sort(alignPos.begin(), alignPos.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b)
						  { return a.first < b.first; });
				int threadholdShort = readSeq.length() / 5, countShort = 0;
				for (std::vector<std::pair<int, int>>::iterator alignIter = alignPos.begin() + 1; alignIter != alignPos.end(); alignIter++)
					if ((alignIter->first) - ((alignIter - 1)->second) > threadholdShort)
						countShort++;
				if (countShort * 2 > alignPos.size() - 1)
					continue;
				// calRead end
				return {true, iter - begin}; // type 'N'
			}
			return {false, -1}; // return type as input, return(type, -1)
		}

		std::tuple<std::string, int> judgeFusionSameChr(MapRecordIter begin, MapRecordIter end, referenceStorage &genome)
		{
			if (end - begin > 3)
			{ // must have full-length sequence for the two fusion part
				// find the longest exon region and store in r1
				MapRecordIter r1 = begin;
				for (MapRecordIter iter = begin + 1; iter != end; iter++)
					if ((iter->exonLength) > (r1->exonLength))
						r1 = iter;
				// check if mostly records overlap with r1->normal circRNA candidate
				int *overlap = new int[end - begin];
				for (MapRecordIter iter = begin; iter != end; iter++)
					overlap[iter - begin] = ((iter->refStart) < (r1->refEnd) && (r1->refStart) < (iter->refEnd)) ? 1 : 0;
				if (std::accumulate(overlap, overlap + (end - begin), 0) > (end - begin - 2))
					return {"sameChrCandidate", -1}; // type 'N'
				// check if have only two mapped sites and switch for each->sameChrFusion
				for (int i = 1; i < (end - begin); i++)
					if (overlap[i] ^ overlap[i - 1] == 0)
					{ // find no switch
						auto [type, rowId]{readRemapRef(begin, end, genome)};
						if (type)
							return {"sameChrCandidate", rowId}; // type 'N'
						else
							return {"sameChrLowConf", rowId}; // type 'C2'
					}
				if ((((begin + 2)->refStart) - (begin->refEnd)) > 40 &&
					(((begin + 3)->refStart) - ((begin + 1)->refEnd)) > 40 &&
					(((begin + 1)->refEnd) - (begin->refEnd)) > 40 &&
					(((begin + 2)->refEnd) - ((begin + 1)->refEnd)) > 40 &&
					(((begin + 3)->refEnd) - ((begin + 2)->refEnd)) > 40 &&
					(((begin + 1)->refStart) - (begin->refEnd)) > -40 &&
					(((begin + 2)->refStart) - ((begin + 1)->refEnd)) > -40 &&
					(((begin + 3)->refStart) - ((begin + 2)->refEnd)) > -40)
				{ // distance between two fusion part is not too close
					auto [type, rowId]{readRemapRef(begin, end, genome)};
					if (type)
						return {"sameChrCandidate", rowId}; // type 'N'
					else
						return {"sameChrFusion", rowId}; // type 'F2'
				}
				else
				{
					auto [type, rowId]{readRemapRef(begin, end, genome)};
					if (type)
						return {"sameChrCandidate", rowId}; // type 'N'
					else
						return {"sameChrFusionLowConf", rowId}; // type 'FC2'
				}
			}
			if (end - begin == 3)
			{ // do not have enough evidence for fusion circRNA
				// get first and last record which may be the same part (not full) of the fusion circRNA
				MapRecordIter r0 = begin, r1 = begin + 1, r2 = begin + 2;
				// if ((r0 overlap with r1 or r0 overlap with r2) and (r2 overlap with r1 or r2 overlap with r0))
				if ((((r0->refStart) < (r1->refEnd) && (r1->refStart) < (r0->refEnd)) || ((r0->refStart) < (r2->refEnd) && (r2->refStart) < (r0->refEnd))) &&
					(((r2->refStart) < (r1->refEnd) && (r1->refStart) < (r2->refEnd)) || ((r2->refStart) < (r0->refEnd) && (r0->refStart) < (r2->refEnd))))
				{
					if ((r0->refStart) < (r1->refEnd) && (r1->refStart) < (r0->refEnd))
						return {"sameChrCandidate", -1}; // type 'N'
					else
					{
						auto [type, rowId]{readRemapRef(begin, end, genome)};
						if (type)
							return {"sameChrCandidate", rowId}; // type 'N'
						else
							return {"sameChrFusionLowConf", rowId}; // type 'FC2'
					}
				}
				else
				{
					auto [type, rowId]{readRemapRef(begin, end, genome)};
					if (type)
						return {"sameChrCandidate", rowId}; // type 'N'
					else
						return {"sameChrLowConf", rowId}; // type 'C2'
				}
			}
			if (end - begin == 2)
			{
				MapRecordIter r1 = begin, r2 = begin + 1;
				// if r1 overlap with r2
				if ((r1->refStart) < (r2->refEnd) && (r2->refStart) < (r1->refEnd))
				{
					return {"sameChrCandidate", -1}; // type 'N'
				}
				auto [type, rowId]{readRemapRef(begin, end, genome)};
				if (type)
					return {"sameChrCandidate", rowId}; // type 'N'
				else
					return {"sameChrNotCandidate", rowId}; // type 'U'
			}
			return {"sameChrNotCandidate", -1}; // type 'U'
		}

		void outputTypeInfo(const std::filesystem::path &mapInfo, const std::string &queryId, const std::string &type, const int rowId, const int &queryEndStartDis)
		{
			std::ofstream mapInfoFile(mapInfo, std::ios::app);
			// mapInfoFile << queryId + "\t" + type + "\t" + std::to_string(queryEndStartDis) + "\t" + to_string(iter) + "\n";
			mapInfoFile << queryId + "\t" + type + "\t" + std::to_string(rowId) + "\t" + std::to_string(queryEndStartDis) + "\n";
			mapInfoFile.close();
		}

		std::tuple<std::string, int> judgeFusionDiffChr(MapRecordIter begin, MapRecordIter end, referenceStorage &genome)
		{
			auto [judgeType, rowId]{readRemapRef(begin, end, genome)};
			if (judgeType)
				return {"sameChrCandidate", rowId}; // type 'N'
			std::set<std::string> chrSet;
			std::map<std::string, std::set<char>> chrStrandSet;
			for (MapRecordIter iter = begin; iter != end; iter++)
				chrSet.insert(iter->chr), chrStrandSet[iter->chr].insert(iter->strand);
			int chrCount = chrSet.size();
			std::string chrs[2]{*(chrSet.begin()), *(chrSet.rbegin())};
			if (begin->chr != chrs[0])
				swap(chrs[0], chrs[1]);
			// in same chromosome but different strand
			if (chrStrandSet[chrs[0]].size() > 1 || chrStrandSet[chrs[1]].size() > 1)
				return {"twoChrsNotCandidate", -1}; // type 'U'
			unsigned int targetNum = 0;
			if (end - begin >= 4)
			{
				for (MapRecordIter iter = begin; iter != end; iter++, targetNum ^= 1)
					if (iter->chr != chrs[0])
						return {"twoChrsFusionLowConf", -1}; // type "FC1"
				return {"twoChrsFusion", -1};				 // type "F1"
			}
			return {"twoChrsFusionLowConf", -1}; // type "FC1"
		}

		// Type
		// M: multipleChrs
		// N: sameChrCandidate
		// U: NotCandidate
		// C1: twoChrLowConf
		// F1: twoChrsFusion
		// FC1: twoChrsFusionLowConf
		// C2: sameChrLowConf
		// F2: sameChrFusion
		// FC2: sameChrFusionLowConf
		std::tuple<std::string, int, int> isFusionCandidate(MapRecordIter begin, MapRecordIter end, referenceStorage &genome, std::filesystem::path mapInfo)
		{ // getType
			MapRecordIter iter;
			// judge if read mapped to different chr
			std::set<std::string> chrs;
			for (iter = begin; iter != end; iter++)
				chrs.insert(iter->chr);
			int chrCount = chrs.size();
			int queryEndStartDis = __INT_MAX__;
			for (iter = begin; iter + 1 != end; iter++)
				queryEndStartDis = std::min(std::abs(((iter + 1)->queryStart) - (iter->queryEnd)), queryEndStartDis);
			if (chrCount > 2)
			{
				// all mapped site located in more than 2 chr
				if (circfull::debug)
					outputTypeInfo(mapInfo, begin->queryId, "multipleChrs", -1, queryEndStartDis);
				return {"multipleChrs", -1, queryEndStartDis}; // type "M"
			}
			if (chrCount == 2)
			{
				// If there are exactly two unique chromosomes, the function checks if the read spans across the two chromosomes
				std::string chr1 = *(chrs.begin()), chr2 = *(chrs.rbegin());
				unsigned int switchChr = 0;
				// If the read spans across the two chromosomes more than 1 time
				for (iter = begin; iter + 1 != end; iter++)
					if (iter->chr != (iter + 1)->chr)
						switchChr++;
				if (switchChr > 1)
				{
					auto [type, rowId]{judgeFusionDiffChr(begin, end, genome)};
					// output only line when rowId is not -1, or output all lines
					if (circfull::debug)
						outputTypeInfo(mapInfo, begin->queryId, type, rowId, queryEndStartDis);
					return {type, rowId, queryEndStartDis};
				}
				if (circfull::debug)
					outputTypeInfo(mapInfo, begin->queryId, "twoChrLowConf", -1, queryEndStartDis); // type 'C1'
				return {"twoChrLowConf", -1, queryEndStartDis};										// type "C1"
			}
			// all mapped site located in the same chr
			auto [type, rowId]{judgeFusionSameChr(begin, end, genome)};
			if (circfull::debug)
				outputTypeInfo(mapInfo, begin->queryId, type, rowId, queryEndStartDis);
			return {type, rowId, queryEndStartDis};
		}

		MapRecord *mergeExons(MapRecordIter first, MapRecordIter second, referenceStorage &genome)
		{
			if (abs((first->refStart) - (second->refStart)) > hangLen)
			{
				std::vector<std::pair<int, bool>> exonPos;
				// parse exon region and merge it
				std::istringstream issExonStart{first->exonStart}, issExonEnd{first->exonEnd};
				std::string tmp;
				while (getline(issExonStart, tmp, ','))
					exonPos.push_back(std::make_pair(std::stoi(tmp), true));
				while (getline(issExonEnd, tmp, ','))
					exonPos.push_back(std::make_pair(std::stoi(tmp), false));
				issExonStart = std::istringstream{second->exonStart}, issExonEnd = std::istringstream{second->exonEnd};
				while (getline(issExonStart, tmp, ','))
					exonPos.push_back(std::make_pair(std::stoi(tmp), true));
				while (getline(issExonEnd, tmp, ','))
					exonPos.push_back(std::make_pair(std::stoi(tmp), false));
				// sort by first (position) and second (true sorted before false)
				std::sort(exonPos.begin(), exonPos.end(), [](const std::pair<int, bool> &a, const std::pair<int, bool> &b)
						  { return a.first < b.first || (a.first == b.first && a.second); });
				int exonNum = 0, exonStart = 0, exonEnd = 0, nowStatus = 0;
				std::string exonStartStr = "", exonEndStr = "";
				unsigned long long int exonLen = 0;
				for (std::vector<std::pair<int, bool>>::iterator iter = exonPos.begin(); iter != exonPos.end(); iter++)
				{
					if ((nowStatus > 0) ^ (nowStatus + (iter->second ? 1 : -1) > 0))
					{
						if (nowStatus)
						{
							exonEnd = iter->first;
							exonStartStr += (exonNum ? "," : "") + std::to_string(exonStart);
							exonEndStr += (exonNum ? "," : "") + std::to_string(exonEnd);
							exonLen += exonEnd - exonStart;
							exonNum++;
						}
						else
						{
							exonStart = iter->first;
						}
					}
					nowStatus += (iter->second ? 1 : -1);
				}
				int refStart = std::min(first->refStart, second->refStart), refEnd = std::max(first->refEnd, second->refEnd);
				int queryStart = std::min(first->queryStart, second->queryStart), queryEnd = std::max(first->queryEnd, second->queryEnd);
				std::string leftSeq = genome.getSeqByPosition(first->chr, refStart - 2, refStart);
				std::string rightSeq = genome.getSeqByPosition(first->chr, refEnd, refEnd + 2);
				return new MapRecord{first->queryId, first->chr, refStart, refEnd, queryStart, queryEnd, exonStartStr, exonEndStr, exonLen, first->strand, leftSeq, rightSeq};
			}
			return NULL;
		}

		std::tuple<std::vector<MapRecord>, std::vector<MapRecord>, std::vector<MapRecord>> filterCandidateCircReads(/*std::filesystem::path mapInfo*/ std::vector<MapRecord> &records, referenceStorage &genome, std::filesystem::path fastq, std::filesystem::path oPrefix, int nthread)
		{
			// int nMap = getFileNrow(mapInfo);
			int nMap = records.size();
			// std::vector<MapRecord> records {nMap};
			// std::ifstream mapInfoFile(mapInfo);
			//  Reading mapping information file and find duplicated reads
			ids.clear();
			int seqCount = 0;
			for (int i = 0; i < nMap; i++)
			{
				// mapInfoFile >> records[i];
				std::map<std::string, int>::iterator iter = ids.find(records[i].queryId);
				if (iter == ids.end())
					ids[records[i].queryId] = seqCount++;
			}
			// mapInfoFile.close();
			std::sort(records.begin(), records.end());
			seqan3::sequence_file_input fastqIn{fastq};
			seqs = new seqan3::dna5_vector[seqCount];
			std::string tmpId;
			for (auto &&record : fastqIn)
			{
				tmpId = getId(record.id());
				std::map<std::string, int>::iterator iter{ids.find(tmpId)};
				if (iter != ids.end())
					seqs[iter->second] = record.sequence();
			}

			// renew output file
			std::filesystem::path typeInfoFile = oPrefix, sameChrFusionInfoFile = oPrefix, diffChrFusionInfoFile = oPrefix, circCandidateInfoFile = oPrefix;
			typeInfoFile += ".txt";
			sameChrFusionInfoFile += ".sameChrFusion.txt";
			diffChrFusionInfoFile += ".diffChrFusion.txt";
			circCandidateInfoFile += ".circCandidate.txt";
			if (checkFile(typeInfoFile))
				delFile(&typeInfoFile);
			if (checkFile(sameChrFusionInfoFile))
				delFile(&sameChrFusionInfoFile);
			if (checkFile(diffChrFusionInfoFile))
				delFile(&diffChrFusionInfoFile);
			if (checkFile(circCandidateInfoFile))
				delFile(&circCandidateInfoFile);

			MapRecordIter iter{records.begin()};
			MapRecordIter mapBegin{iter};
			std::string nowId{iter->queryId};
			boost::asio::thread_pool tp(nthread);
			std::mutex mt;
			std::vector<MapRecord> circCandidateInfo, sameChrFusionInfo, diffChrFusionInfo;
			while (iter != records.end())
			{
				while (++iter != records.end())
					if (iter->queryId != nowId)
						break;
				if (mapBegin + 1 != iter)
				{
					boost::asio::post(tp, [&typeInfoFile, &sameChrFusionInfoFile, &diffChrFusionInfoFile, &circCandidateInfoFile, &circCandidateInfo, &sameChrFusionInfo, &diffChrFusionInfo, &mt, &genome, mapBegin, iter]
									  {
						auto [type, rowId, queryESDis] { isFusionCandidate(mapBegin, iter, genome, typeInfoFile)};
						// output type F1/F2/N into 3 files
						if (type == "twoChrsFusion") {
							mt.lock();
							for (MapRecordIter addIter = mapBegin; addIter != iter; addIter++)
								diffChrFusionInfo.push_back(*addIter);
							mt.unlock();
							if (circfull::debug) {
								std::ofstream diffChrFusionInfoFileOut(diffChrFusionInfoFile, std::ios::app);
								for (MapRecordIter addIter = mapBegin; addIter != iter; addIter++)
									diffChrFusionInfoFileOut << to_string(&(*addIter));
								diffChrFusionInfoFileOut.close();
							}
						}
						if (type == "sameChrFusion") {
							mt.lock();
							for (MapRecordIter addIter = mapBegin; addIter != iter; addIter++)
								sameChrFusionInfo.push_back(*addIter);
							mt.unlock();
							if (circfull::debug) {
								std::ofstream sameChrFusionInfoFileOut(sameChrFusionInfoFile, std::ios::app);
								for (MapRecordIter addIter = mapBegin; addIter != iter; addIter++)
									sameChrFusionInfoFileOut << to_string(&(*addIter));
								sameChrFusionInfoFileOut.close();
							}
						}
						if (type == "sameChrCandidate") {  //type 'N'
							if (rowId == -1) {
								if (iter - mapBegin == 2) { // to generate as full as possible exon region for 2-record alignment (adjExplainFL)
									MapRecord* newExons = mergeExons(mapBegin, mapBegin + 1, genome);
									if (newExons) { // can merge
										mt.lock();
										circCandidateInfo.push_back(*newExons);
										mt.unlock();
										std::ofstream circCandidateInfoFileOut(circCandidateInfoFile, std::ios::app);
										circCandidateInfoFileOut << to_string(newExons);
										circCandidateInfoFileOut.close();
									} else { // no merge
										mt.lock();
										for (MapRecordIter addIter = mapBegin; addIter != iter; addIter++)
											circCandidateInfo.push_back(*addIter);
										mt.unlock();
										std::ofstream circCandidateInfoFileOut(circCandidateInfoFile, std::ios::app);
										for (MapRecordIter addIter = mapBegin; addIter != iter; addIter++)
											circCandidateInfoFileOut << to_string(&(*addIter));
										circCandidateInfoFileOut.close();
									}
								} else {
									mt.lock();
									for (MapRecordIter addIter = mapBegin; addIter != iter; addIter++)
										circCandidateInfo.push_back(*addIter);
									mt.unlock();
									std::ofstream circCandidateInfoFileOut(circCandidateInfoFile, std::ios::app);
									for (MapRecordIter addIter = mapBegin; addIter != iter; addIter++)
										circCandidateInfoFileOut << to_string(&(*addIter));
									circCandidateInfoFileOut.close();
								}
							} else {
								mt.lock();
								circCandidateInfo.push_back(*(mapBegin + rowId));
								mt.unlock();
								if (circfull::debug) {
									std::ofstream circCandidateInfoFileOut(circCandidateInfoFile, std::ios::app);
									circCandidateInfoFileOut << to_string(&(*(mapBegin + rowId)));
									circCandidateInfoFileOut.close();
								}
							}
						} });
				}
				mapBegin = iter;
				if (iter != records.end())
					nowId = iter->queryId;
			}
			tp.join();
			return {circCandidateInfo, sameChrFusionInfo, diffChrFusionInfo};
		}

		std::optional<BSRecord> getBS(MapRecordIter record, referenceStorage &genome, std::string minimap2)
		{
			BSRecord ret;
			ret.queryId = record->queryId;
			ret.chr = record->chr;
			// get the record seq
			std::string readSeq = seqan3ToString(seqs[ids[record->queryId]]);
			// get the reference region and its upstream and downstream hangLen sequence
			std::string refSeq = genome.getSeqByPosition(record->chr, record->refStart - hangLen, record->refEnd + hangLen);
			int refLen = (record->refEnd) - (record->refStart);
			std::string mapRef = refSeq + refSeq;
			mpAligner aligner{mapRef, 14, true, true};
			int nReg;
			mm_reg1_t *mmRecords = aligner.map(readSeq, &nReg);
			// calRead start
			for (int i = 0; i < nReg; i++)
			{
				mm_reg1_t *hit = &mmRecords[i];
				assert(hit->p);
				std::vector<seqan3::cigar> cigarSeq = mpAlignerCigar2Seqan3Cigar(hit);
				// parse cigar sequence of hit to get exons positions
				std::vector<std::pair<int, int>> exons;
				std::vector<int> exonStartDis, exonEndDis;
				int currentPos = hit->rs;
				// char hitStrand = hit->p->trans_strand == 1 ? '+' : hit->p->trans_strand == 2 ? '-' : '.';
				char hitStrand = "+-"[hit->rev];
				int exonStart = currentPos, exonEnd;
				// ====|=========|====|====|=========|====
				// hang    ref    hang hang    ref    hang
				// a   b         c    d    e         f   g
				// find the closest exon_start to the start of the second ref_seq (position e)
				// find the closest exon_end to the end of the first ref_seq (position c)
				for (auto [cigarCount, cigarOperation] : cigarSeq)
				{
					switch (cigarOperation.to_char())
					{
					case 'M':
					case 'D':
						currentPos += cigarCount;
						break;
					case 'N':
						exonEnd = currentPos;
						exons.push_back(std::make_pair(exonStart, exonEnd));
						// exonsDis.push_back(std::make_pair(exonStart - (3 * hangLen + refLen), exonEnd - (hangLen + refLen)));
						exonStartDis.push_back(std::abs(exonStart - (3 * hangLen + refLen)));
						exonEndDis.push_back(std::abs(exonEnd - (hangLen + refLen)));
						currentPos += cigarCount;
						exonStart = currentPos;
					}
				}
				exonEnd = currentPos;
				exons.push_back(std::make_pair(exonStart, exonEnd));
				exonStartDis.push_back(std::abs(exonStart - (3 * hangLen + refLen)));
				exonEndDis.push_back(std::abs(exonEnd - (hangLen + refLen)));
				free(hit->p);
				int minExonStartDis = *std::min_element(exonStartDis.begin(), exonStartDis.end());
				int minExonEndDis = *std::min_element(exonEndDis.begin(), exonEndDis.end());
				// get indexs of minExonStartDis - 1 and minExonEndDis
				std::set<int> minExonStartIdx, minExonEndIdx;
				for (int i = 0; i < exons.size(); i++)
				{
					if (exonStartDis[i] == minExonStartDis)
						minExonStartIdx.insert(i - 1);
					if (exonEndDis[i] == minExonEndDis)
						minExonEndIdx.insert(i);
				}
				std::set<int> commonIdx; // find the exon map to the end of the first ref_seq and the next exon map to the start of the second ref_seq
				std::set_intersection(minExonStartIdx.begin(), minExonStartIdx.end(), minExonEndIdx.begin(), minExonEndIdx.end(), std::inserter(commonIdx, commonIdx.begin()));
				if (commonIdx.empty())
				{																  // if no common index
					if (minExonStartIdx.size() == 1 && minExonEndIdx.size() == 1) // only have one hit?
						if (*(minExonStartIdx.begin()) + 1 == *(minExonEndIdx.begin()))
						{
							int hitPos = *minExonEndIdx.begin();
							if (exons[hitPos].first < 2 * hangLen + refLen)
							{															// this hit locate at the first ref_seq (from a to d)
								if (hitPos < exons.size() - 1)							// it is not the last exon
									if (exons[hitPos + 1].first > 2 * hangLen + refLen) // the next exon map to the second ref_seq (from d to g)
										commonIdx.insert(hitPos);
							}
							else
							{															 // this hit locate at the second ref_seq (from d to g)
								if (hitPos > 0)											 // it is not the first exon
									if (exons[hitPos - 1].second < 2 * hangLen + refLen) // the previous exon map to the first ref_seq (from a to d)
										commonIdx.insert(hitPos - 1);
							}
						}
				}
				for (auto idx : commonIdx)
				{
					int refLeftPos = (record->refStart) + exons[idx + 1].first - (3 * hangLen + refLen) + 1;
					int refRightPos = (record->refStart) - hangLen + exons[idx].second; // (record->refStart) - hangLen + exons[hit + 1].first - (2 * hangLen + refLen)
					if (refLeftPos > refRightPos)
						continue;
					ret.BSPos.push_back(std::make_pair(refLeftPos, refRightPos));
					ret.BSMotif.push_back(
						std::make_pair(
							SpliceMotif{genome.getSeqByPosition(record->chr, refLeftPos - 3, refLeftPos - 1)},
							SpliceMotif{genome.getSeqByPosition(record->chr, refRightPos, refRightPos + 2)}));
					ret.BSstrand += hitStrand;
				}
			}
			free(mmRecords);
			// clean tmp files
			/*delFile(&readSeqFile);
			delFile(&refBSFile);
			delFile(&outSam);*/
			if (ret.BSPos.empty())
				return std::nullopt;
			return ret;
		}

		std::vector<BSRecord> detectBS(referenceStorage &genome, std::filesystem::path BSListFile, std::vector<MapRecord> &circCandidates, int nthread, std::string oPrefix, std::string minimap2)
		{
			MapRecordIter iter{circCandidates.begin()};
			MapRecordIter mapBegin{iter};
			std::string nowId{iter->queryId};
			boost::asio::thread_pool tp(nthread);
			std::mutex mt;
			// display the boost::timer::progress_display
			boost::timer::progress_display pd(circCandidates.size());
			std::vector<BSRecord> BSList;
			std::ofstream BSListFileStream;
			if (circfull::debug)
			{
				BSListFileStream.open(BSListFile, std::ios::out);
			}
			while (iter != circCandidates.end())
			{
				while (++iter != circCandidates.end())
					if (iter->queryId != nowId)
						break;
				boost::asio::post(tp, [&genome, &mt, &pd, &minimap2, &BSList, &BSListFileStream, mapBegin, iter]
								  {
					// sort the same queryId records by exonLength
					std::sort(mapBegin, iter, [](const MapRecord & a, const MapRecord & b) {return a.exonLength > b.exonLength;});
					// get back-splicing junction from the longest record
					auto ret { getBS(mapBegin, genome, minimap2) };
					if (ret.has_value()) {
						mt.lock();
						BSList.push_back(ret.value());
						if (circfull::debug)
							BSListFileStream << to_string(&ret.value());
						mt.unlock();
					}
					mt.lock();
					pd += iter - mapBegin;
					mt.unlock(); });
				mapBegin = iter;
				if (iter != circCandidates.end())
					nowId = iter->queryId;
			}
			if (circfull::debug)
				BSListFileStream.close();
			tp.join();
			return BSList;
		}

		BSRecord adjustBSPosByExon(const BSRecordIter &origin, const ExonIndex<> &exonIndex, const GtfStorage &gtfRecords, referenceStorage &genome)
		{
			BSRecord ret;
			ret.queryId = origin->queryId;
			ret.chr = origin->chr;
			for (int i = 0; i < origin->BSPos.size(); i++)
			{
				const std::string startKey = origin->chr + "_" + std::to_string(origin->BSPos[i].first);
				const std::string endKey = origin->chr + "_" + std::to_string(origin->BSPos[i].second);
				const std::map<std::string, std::set<std::string>>::const_iterator startIter = exonIndex.exonStartIndex.find(startKey);
				const std::map<std::string, std::set<std::string>>::const_iterator endIter = exonIndex.exonEndIndex.find(endKey);
				std::set<std::string> commonGene{};
				if (startIter != exonIndex.exonStartIndex.end() && endIter != exonIndex.exonEndIndex.end())
					std::set_intersection(startIter->second.begin(), startIter->second.end(), endIter->second.begin(), endIter->second.end(), std::inserter(commonGene, commonGene.begin()));
				for (std::set<std::string>::iterator it = commonGene.begin(); it != commonGene.end();)
				{
					const char tmpStrand = exonIndex.exonStrandIndex.find(*it)->second;
					if (tmpStrand != origin->BSstrand[i] && tmpStrand != '.')
					{
						it = commonGene.erase(it);
					}
					else
					{
						it++;
					}
				}
				if (commonGene.size() == 0)
				{
					ret.BSPos.push_back(origin->BSPos[i]);
					ret.BSMotif.push_back(origin->BSMotif[i]);
					ret.BSstrand += origin->BSstrand[i];
					continue;
				}
				if (commonGene.size() == 1)
				{
					const std::string geneName = *(commonGene.begin());
					const GtfIter gtfBegin = gtfRecords.lower_bound(geneName), gtfEnd = gtfRecords.upper_bound(geneName);
					int closestExonStart = __INT_MAX__, closestExonEnd = __INT_MAX__;
					for (GtfIter it = gtfBegin; it != gtfEnd; it++)
					{
						if (closestExonStart == __INT_MAX__ || abs(it->second.start - origin->BSPos[i].first) < abs(closestExonStart - origin->BSPos[i].first))
							closestExonStart = it->second.start;
						if (closestExonEnd == __INT_MAX__ || abs(it->second.end - origin->BSPos[i].second) < abs(closestExonEnd - origin->BSPos[i].second))
							closestExonEnd = it->second.end;
					}
					ret.BSPos.push_back({closestExonStart, closestExonEnd});
					ret.BSMotif.push_back(
						std::make_pair(
							SpliceMotif{genome.getSeqByPosition(origin->chr, closestExonStart - 3, closestExonStart - 1)},
							SpliceMotif{genome.getSeqByPosition(origin->chr, closestExonEnd, closestExonEnd + 2)}));
					ret.BSstrand += origin->BSstrand[i];
					continue;
				}
				int geneClosestExonStart = __INT_MAX__, geneClosestExonEnd = __INT_MAX__;
				for (std::string geneName : commonGene)
				{
					const GtfIter gtfBegin = gtfRecords.lower_bound(geneName), gtfEnd = gtfRecords.upper_bound(geneName);
					int closestExonStart = __INT_MAX__, closestExonEnd = __INT_MAX__;
					for (GtfIter it = gtfBegin; it != gtfEnd; it++)
					{ // find the closer exon boundary to the BS detected
						if (closestExonStart == __INT_MAX__ || abs(it->second.start - origin->BSPos[i].first) < abs(closestExonStart - origin->BSPos[i].first))
							closestExonStart = it->second.start;
						if (closestExonEnd == __INT_MAX__ || abs(it->second.end - origin->BSPos[i].second) < abs(closestExonEnd - origin->BSPos[i].second))
							closestExonEnd = it->second.end;
					}
					if (geneClosestExonStart == __INT_MAX__ || abs(closestExonStart - origin->BSPos[i].first) + abs(closestExonEnd - origin->BSPos[i].second) < abs(geneClosestExonStart - origin->BSPos[i].first) + abs(geneClosestExonEnd - origin->BSPos[i].second)) // if this gene explain the BS better
						geneClosestExonStart = closestExonStart, geneClosestExonEnd = closestExonEnd;
				}
				ret.BSPos.push_back({geneClosestExonStart, geneClosestExonEnd});
				ret.BSMotif.push_back(
					std::make_pair(
						SpliceMotif{genome.getSeqByPosition(origin->chr, geneClosestExonStart - 3, geneClosestExonStart - 1)},
						SpliceMotif{genome.getSeqByPosition(origin->chr, geneClosestExonEnd, geneClosestExonEnd + 2)}));
				ret.BSstrand += origin->BSstrand[i];
			}
			return ret;
		}

		std::optional<BSRecord> filterBSByMotif(const BSRecordIter &origin, const ExonIndex<> &exonIndex, GtfStorage &gtfRecords, referenceStorage &genome)
		{
			BSRecord adjOrigin = adjustBSPosByExon(origin, exonIndex, gtfRecords, genome);
			adjOrigin.adjust = true;
			// find the best mapped BS sites:
			// 1. this BS position is the most mapped one
			// 2. find the most common BS motif like GU-AG listed in BSMotifList
			// 3. as long as possible
			BSRecord ret;
			ret.adjust = true, ret.chr = adjOrigin.chr, ret.queryId = adjOrigin.queryId;
			// retrive all splice site into vector<SpliceSite>
			std::vector<SpliceSite> spliceSiteList;
			for (int i = 0; i < adjOrigin.BSPos.size(); i++)
				spliceSiteList.push_back(
					{adjOrigin.BSPos[i].first,
					 adjOrigin.BSPos[i].second,
					 adjOrigin.BSstrand[i],
					 adjOrigin.BSMotif[i].first.str(),
					 adjOrigin.BSMotif[i].second.str()});
			// find the most mapped BS sites
			std::vector<SpliceSite> filteredSpliceSiteList = getMostFrequent<SpliceSite, std::string>(
				spliceSiteList,
				[](const SpliceSite &a) -> std::string
				{ return std::to_string(a.start) + a.strand + std::to_string(a.end); });

			// find if record can be explained as pair-end motif
			size_t motifSize = BSMotifList.bothMotif.size();
			std::vector<std::string> motifHits;
			int maxCount = 0;
			// get the most frequent BS motif
			for (size_t i = 0; i < motifSize; i++)
				for (size_t j = 0; j < filteredSpliceSiteList.size(); j++)
				{
					std::string bothMotif = filteredSpliceSiteList[j].leftSeq + filteredSpliceSiteList[j].rightSeq;
					if (hitBSMotif(bothMotif, BSMotifList.bothMotif[i], filteredSpliceSiteList[j].strand))
						motifHits.push_back(bothMotif + filteredSpliceSiteList[j].strand);
				}
			if (motifHits.size() > 0)
			{
				std::vector<std::string> mostCommonMotif = getMostFrequent(motifHits);
				filteredSpliceSiteList = filter<SpliceSite>(
					filteredSpliceSiteList,
					[&mostCommonMotif](const SpliceSite &a) -> bool
					{
						return std::find(mostCommonMotif.begin(), mostCommonMotif.end(), a.leftSeq + a.rightSeq + a.strand) != mostCommonMotif.end();
					});
				// find the max length
				int maxLen = 0;
				for (auto spliceSite : filteredSpliceSiteList)
				{
					const int nowLen = spliceSite.end - spliceSite.start;
					if (nowLen > maxLen)
						maxLen = nowLen;
				}
				// filter all BS record by max length
				filteredSpliceSiteList = filter<SpliceSite>(
					filteredSpliceSiteList,
					[&maxLen](const SpliceSite &a) -> bool
					{ return a.end - a.start == maxLen; });
				ret.BSPos = {{filteredSpliceSiteList[0].start, filteredSpliceSiteList[0].end}};
				ret.BSstrand = {filteredSpliceSiteList[0].strand};
				ret.BSMotif = {{SpliceMotif{genome.getSeqByPosition(ret.chr, ret.BSPos[0].first - 3, ret.BSPos[0].first - 1)},
								SpliceMotif{genome.getSeqByPosition(ret.chr, ret.BSPos[0].second, ret.BSPos[0].second + 2)}}};
				// ret.leftSeqMapStatus = "all_bsNumM_BothMotif_BothMax";
				// ret.rightSeqMapStatus = "all_bsNumM_BothMotif_BothMax";
				ret.leftSeqMapStatus = ret.rightSeqMapStatus = true;
				return ret;
			}
			// scan left motif
			for (size_t i = 0; i < motifSize; i++)
				for (size_t j = 0; j < filteredSpliceSiteList.size(); j++)
				{
					if (hitBSMotif(filteredSpliceSiteList[j].leftSeq, BSMotifList.leftMotif[i], filteredSpliceSiteList[j].strand))
						motifHits.push_back(filteredSpliceSiteList[j].leftSeq + filteredSpliceSiteList[j].strand);
				}
			// scan right motif
			for (size_t i = 0; i < motifSize; i++)
				for (size_t j = 0; j < filteredSpliceSiteList.size(); j++)
				{
					if (hitBSMotif(filteredSpliceSiteList[j].rightSeq, BSMotifList.rightMotif[i], filteredSpliceSiteList[j].strand))
						motifHits.push_back(std::to_string(filteredSpliceSiteList[j].strand) + filteredSpliceSiteList[j].rightSeq);
				}
			// if any boundary map to the splicing motif
			if (motifHits.size() > 0)
			{
				std::vector<std::string> mostCommonMotif = getMostFrequent(motifHits);
				std::vector<SpliceSite> leftFilteredSpliceSiteList = filter<SpliceSite>(
					filteredSpliceSiteList,
					[&mostCommonMotif](const SpliceSite &a) -> bool
					{
						return std::find(mostCommonMotif.begin(), mostCommonMotif.end(), a.leftSeq + a.strand) != mostCommonMotif.end();
					});
				std::vector<SpliceSite> rightFilteredSpliceSiteList = filter<SpliceSite>(
					filteredSpliceSiteList,
					[&mostCommonMotif](const SpliceSite &a) -> bool
					{
						return std::find(mostCommonMotif.begin(), mostCommonMotif.end(), std::to_string(a.strand) + a.rightSeq) != mostCommonMotif.end();
					});
				// find the max length for left and right
				int maxLeftLen = 0, maxRightLen = 0;
				for (auto spliceSite : leftFilteredSpliceSiteList)
				{
					const int nowLen = spliceSite.end - spliceSite.start;
					if (nowLen > maxLeftLen)
						maxLeftLen = nowLen;
				}
				for (auto spliceSite : rightFilteredSpliceSiteList)
				{
					const int nowLen = spliceSite.end - spliceSite.start;
					if (nowLen > maxRightLen)
						maxRightLen = nowLen;
				}
				if (maxLeftLen >= maxRightLen)
					filteredSpliceSiteList = filter<SpliceSite>(
						leftFilteredSpliceSiteList,
						[&maxLeftLen](const SpliceSite &a) -> bool
						{ return a.end - a.start == maxLeftLen; });
				else
					filteredSpliceSiteList = filter<SpliceSite>(
						rightFilteredSpliceSiteList,
						[&maxRightLen](const SpliceSite &a) -> bool
						{ return a.end - a.start == maxRightLen; });
				if (!filteredSpliceSiteList.empty())
				{
					ret.BSPos = {{filteredSpliceSiteList[0].start, filteredSpliceSiteList[0].end}};
					ret.BSstrand = {filteredSpliceSiteList[0].strand};
					ret.BSMotif = {{SpliceMotif{genome.getSeqByPosition(ret.chr, ret.BSPos[0].first - 3, ret.BSPos[0].first - 1)},
									SpliceMotif{genome.getSeqByPosition(ret.chr, ret.BSPos[0].second, ret.BSPos[0].second + 2)}}};
					if (maxLeftLen >= maxRightLen)
					{
						// ret.leftSeqMapStatus = "all_bsNumM_LeftMotif_BothMax";
						// ret.rightSeqMapStatus = "all_bsNumM_LeftMotif_RightnoMotif_BothMax";
						ret.leftSeqMapStatus = true, ret.rightSeqMapStatus = false;
					}
					else
					{
						// ret.leftSeqMapStatus = "all_bsNumM_LeftnoMotif_BothMax";
						// ret.rightSeqMapStatus = "all_bsNumM_LeftnoMotif_RightMotif_BothMax";
						ret.leftSeqMapStatus = false, ret.rightSeqMapStatus = true;
					}
					return ret;
				}
			}
			// no motif found, find the longest one
			int maxLen = 0;
			for (auto spliceSite : filteredSpliceSiteList)
			{
				const int nowLen = spliceSite.end - spliceSite.start;
				if (nowLen > maxLen)
					maxLen = nowLen;
			}
			if (maxLen == 0)
				return std::nullopt;
			filteredSpliceSiteList = filter<SpliceSite>(
				filteredSpliceSiteList,
				[&maxLen](const SpliceSite &a) -> bool
				{ return a.end - a.start == maxLen; });
			ret.BSPos = {{filteredSpliceSiteList[0].start, filteredSpliceSiteList[0].end}};
			ret.BSstrand = {filteredSpliceSiteList[0].strand};
			ret.BSMotif = {{SpliceMotif{genome.getSeqByPosition(ret.chr, ret.BSPos[0].first - 3, ret.BSPos[0].first - 1)},
							SpliceMotif{genome.getSeqByPosition(ret.chr, ret.BSPos[0].second, ret.BSPos[0].second + 2)}}};
			// ret.leftSeqMapStatus = "all_bsNumM_noMotif_BothMax";
			// ret.rightSeqMapStatus = "all_bsNumM_noMotif_BothMax";
			ret.leftSeqMapStatus = ret.rightSeqMapStatus = false;
			return ret;
		}

		std::vector<BSRecord> filterBS(referenceStorage &genome, std::vector<BSRecord> &BSList, GtfStorage gtfRecords, const ExonIndex<> &exonIndex, std::filesystem::path filteredBSListFile, const std::string &spliceSignal, int nthread)
		{
			// build motif index, split by comma
			std::stringstream ss{spliceSignal};
			std::string spliceMotif;
			while (std::getline(ss, spliceMotif, ','))
			{
				BSMotifList.addIndex(spliceMotif);
			}
			boost::asio::thread_pool tp(nthread);
			std::mutex mt;
			boost::timer::progress_display pd(BSList.size());
			std::vector<BSRecord> filteredBSList;
			filteredBSList.reserve(BSList.size());
			std::ofstream filteredBSListFileStream;
			if (circfull::debug)
			{
				filteredBSListFileStream.open(filteredBSListFile, std::ios::out);
			}
			// for (auto && record : BSList) {
			for (BSRecordIter BSIter = BSList.begin(); BSIter != BSList.end(); BSIter++)
			{
				boost::asio::post(tp, [&genome, &mt, &pd, &filteredBSList, &exonIndex, &gtfRecords, &filteredBSListFileStream, BSIter]
								  {
					BSRecord filteredRecord;
					auto ret { filterBSByMotif(BSIter, exonIndex, gtfRecords, genome) };
					if (ret.has_value()) {
							filteredRecord = ret.value();
						std::lock_guard<std::mutex> guard(mt);
						filteredBSList.push_back(filteredRecord);
						if (circfull::debug) {
							filteredBSListFileStream << to_string(&filteredRecord);
						}
					}
					std::lock_guard<std::mutex> guard(mt);
					pd += 1; });
			}
			tp.join();
			return filteredBSList;
		}
	}
}