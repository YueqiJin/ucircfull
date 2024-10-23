#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include "circCall.hpp"
#include "circfullUtils.hpp"

namespace circfull
{
	namespace RG
	{

		std::tuple<std::vector<BSRecord>, std::vector<BSRecord>> adjustBSRecordsByHighConfRecords(std::vector<BSRecord> &highConfBSList, std::vector<BSRecord> &toBeAdjustBSList, int nthread)
		{
			// build BSIndex
			ExonIndex<CircID> BSIndex;
			for (auto &record : highConfBSList)
				BSIndex.addExon(record.chr, record.BSPos.begin()->first, record.BSPos.begin()->second, *(record.BSstrand.begin()), CircID{record.chr, *(record.BSstrand.begin()), record.BSPos.begin()->first, record.BSPos.begin()->second, record.BSMotif.begin()->first.str(), record.BSMotif.begin()->second.str()}, errorBSLen);
			boost::asio::thread_pool tp(nthread);
			std::mutex mt;
			for (size_t i = 0; i < toBeAdjustBSList.size(); i++)
			{
				boost::asio::post(tp, [&BSIndex, &toBeAdjustBSList, &highConfBSList, &mt, i]()
								  {
					// find if there is a high conf BS site in nearby region (can be found in BSIndex)
					const std::string startKey = toBeAdjustBSList[i].chr + "_" + std::to_string(toBeAdjustBSList[i].BSPos[0].first);
					const std::string endKey = toBeAdjustBSList[i].chr + "_" + std::to_string(toBeAdjustBSList[i].BSPos[0].second);
					const std::map<std::string, std::set<CircID>>::const_iterator startIter = BSIndex.exonStartIndex.find(startKey);
					const std::map<std::string, std::set<CircID>>::const_iterator endIter = BSIndex.exonEndIndex.find(endKey);
					std::set<CircID> commonGene {};
					if (startIter != BSIndex.exonStartIndex.end() && endIter != BSIndex.exonEndIndex.end())
						std::set_intersection(startIter->second.begin(), startIter->second.end(), endIter->second.begin(), endIter->second.end(), std::inserter(commonGene, commonGene.begin()));
					if (commonGene.size() == 0)
						return;
					if (commonGene.size() == 1) {
						// if at the same strand
						if (commonGene.begin()->strand == *(toBeAdjustBSList[i].BSstrand.begin())) {
							mt.lock();
							toBeAdjustBSList[i].BSPos[0].first = commonGene.begin()->start;
							toBeAdjustBSList[i].BSPos[0].second = commonGene.begin()->end;
							toBeAdjustBSList[i].BSMotif[0].first = SpliceMotif{ commonGene.begin()->leftSeq };
							toBeAdjustBSList[i].BSMotif[0].second = SpliceMotif{ commonGene.begin()->rightSeq };
							toBeAdjustBSList[i].leftSeqMapStatus = toBeAdjustBSList[i].rightSeqMapStatus = true;
							mt.unlock();
						}
						return;
					}
					if (commonGene.size() > 1) {
						// find the closest one
						CircID closestSite;
						int closestDis = __INT_MAX__;
						for (auto & site : commonGene)
							if (site.strand == *(toBeAdjustBSList[i].BSstrand.begin())) {
								int dis = std::abs(site.start - toBeAdjustBSList[i].BSPos[0].first) + std::abs(site.end - toBeAdjustBSList[i].BSPos[0].second);
								if (dis < closestDis)
									closestSite = site;
							}
						if (closestDis == __INT_MAX__)
							return;
						mt.lock();
						toBeAdjustBSList[i].BSPos[0].first = closestSite.start;
						toBeAdjustBSList[i].BSPos[0].second = closestSite.end;
						toBeAdjustBSList[i].BSMotif[0].first = SpliceMotif{ closestSite.leftSeq };
						toBeAdjustBSList[i].BSMotif[0].second = SpliceMotif{ closestSite.rightSeq };
						toBeAdjustBSList[i].leftSeqMapStatus = toBeAdjustBSList[i].rightSeqMapStatus = true;
						mt.unlock();
					} });
			}
			tp.join();
			// move all leftSeqMapStatus & rightSeqMapStatus == true to highConfBSList
			auto [toMoveBSList, lowConfBSList] = splitByBool<BSRecord>(
				toBeAdjustBSList,
				[](const BSRecord &record) -> bool
				{
					return (record.leftSeqMapStatus && record.rightSeqMapStatus);
				});
			highConfBSList.insert(highConfBSList.end(), toMoveBSList.begin(), toMoveBSList.end());
			toBeAdjustBSList.clear();
			return {highConfBSList, lowConfBSList};
		}

		void clustCirc(std::vector<BSRecord> &BSList, int nthread)
		{
			std::string nowChr = BSList.begin()->chr;
			int nowStart = BSList.begin()->BSPos[0].first;
			int nowEnd = BSList.begin()->BSPos[0].second;
			BSRecordIter nowBegin = BSList.begin();
			boost::asio::thread_pool tp(nthread);
			std::mutex mt;
			for (BSRecordIter iter = BSList.begin(); iter != BSList.end(); iter++)
			{
				BSRecordIter endIter = iter + 1;
				if (endIter != BSList.end())
					if (nowChr == endIter->chr && endIter->BSPos.begin()->first <= nowEnd && endIter->BSPos.begin()->second >= nowStart)
					{ // overlap
						nowEnd = std::max(nowEnd, endIter->BSPos.begin()->second);
						nowStart = std::min(nowStart, endIter->BSPos.begin()->first);
						continue;
					}
				// parse this bundle
				boost::asio::post(tp, [&mt, nowBegin, endIter]
								  {
					// count the frequency of each BS site
					std::map<CircID, int> freq = countFrequencies<BSRecord, CircID>(
						std::vector<BSRecord>(nowBegin, endIter),
						[](const BSRecord & record)->CircID {
							return CircID{record.chr, *(record.BSstrand.begin()), record.BSPos.begin()->first, record.BSPos.begin()->second, record.BSMotif.begin()->first.str(), record.BSMotif.begin()->second.str()};
						}
					);
					// merge BS sites for both start and end distance < errorBSLen by frequency
					std::map<CircID, CircID> target;
					std::vector<CircID> originBS;
					for (auto [circID, count] : freq)
						target[circID] = circID, originBS.push_back(circID);
					sort(originBS.begin(), originBS.end(), [](const CircID & a, const CircID & b)->bool {
						return (a.chr < b.chr || (a.chr == b.chr && a.start < b.start) || (a.chr == b.chr && a.start == b.start && a.end < b.end));
					});
					for (size_t i = 0; i < originBS.size(); i++) {
						int maxFreq = freq[originBS[i]];
						if (maxFreq == 0)
							continue;
						int nowTarget = i;
						for (size_t j = 0; j < originBS.size(); j++)
							if (originBS[i].strand == originBS[j].strand && abs(originBS[j].start - originBS[i].start) <= errorBSLen && abs(originBS[j].end - originBS[i].end) <= errorBSLen)
								if (freq[originBS[j]] >= maxFreq) {
									maxFreq = freq[originBS[j]];
									nowTarget = j;
								}
						if (nowTarget != i) {
							target[originBS[i]] = originBS[nowTarget];
							freq[originBS[nowTarget]] += freq[originBS[i]];
							freq[originBS[i]] = 0;
						}
					}
					// update BS record information
					for (BSRecordIter updateIter = nowBegin; updateIter != endIter; updateIter++) {
						CircID pre{updateIter->chr, *(updateIter->BSstrand.begin()), updateIter->BSPos[0].first, updateIter->BSPos[0].second, updateIter->BSMotif[0].first.str(), updateIter->BSMotif[0].second.str()};
						if (target[pre] == pre)
							continue;
						updateIter->BSPos[0].first = target[pre].start;
						updateIter->BSPos[0].second = target[pre].end;
						updateIter->BSMotif[0].first = SpliceMotif{ target[pre].leftSeq };
						updateIter->BSMotif[0].second = SpliceMotif{ target[pre].rightSeq };
					} });
				// new bundle
				if (endIter != BSList.end())
				{
					nowChr = endIter->chr;
					nowStart = endIter->BSPos.begin()->first;
					nowEnd = endIter->BSPos.begin()->second;
					nowBegin = endIter;
				}
			}

			tp.join();
		}

		std::vector<BSRecord> clusterBS(const std::vector<BSRecord> &BSList, std::filesystem::path BSAdjustFile, int nthread)
		{
			// filter both motif mapped and have > 2 UMI-reads supported BS sites
			std::vector<BSRecord> fullMotifBSList = filter<BSRecord>(
				BSList,
				[](const BSRecord &record) -> bool
				{
					return (record.leftSeqMapStatus && record.rightSeqMapStatus);
				});
			std::vector<CircID> highConfSite = filter<std::map<CircID, int>, std::pair<CircID, int>, CircID>(
				countFrequencies<BSRecord, CircID>(
					fullMotifBSList,
					[](const BSRecord &record) -> CircID
					{
						return {record.chr, *(record.BSstrand.begin()), record.BSPos.begin()->first, record.BSPos.begin()->second, record.BSMotif.begin()->first.str(), record.BSMotif.begin()->second.str()};
					}),
				[](const std::pair<CircID, int> &item) -> bool
				{
					return item.second >= 2;
				},
				[](const std::pair<CircID, int> &item) -> CircID
				{
					return item.first;
				});
			auto [highConfBSList, toEditBSList] = splitByBool<BSRecord>(
				BSList,
				[&highConfSite](const BSRecord &record) -> bool
				{
					return (std::find(highConfSite.begin(), highConfSite.end(), std::make_tuple(record.chr, *(record.BSstrand.begin()), record.BSPos.begin()->first, record.BSPos.begin()->second)) != highConfSite.end());
				});
			// for not high conf sites, adjust BS sites according to high conf sites
			fullMotifBSList.clear();
			auto [adjHighConfBSList, lowConfBSList] = adjustBSRecordsByHighConfRecords(highConfBSList, toEditBSList, nthread);
			// clustering high conf BS sites
			sort(adjHighConfBSList.begin(), adjHighConfBSList.end(), [](const BSRecord &a, const BSRecord &b) -> bool
				 { return (a.chr < b.chr || (a.chr == b.chr && a.BSPos[0].first < b.BSPos[0].first) || (a.chr == b.chr && a.BSPos[0].first == b.BSPos[0].first && a.BSPos[0].second < b.BSPos[0].second)); });
			clustCirc(adjHighConfBSList, nthread);
			// clustering low conf BS sites
			sort(lowConfBSList.begin(), lowConfBSList.end(), [](const BSRecord &a, const BSRecord &b) -> bool
				 { return (a.chr < b.chr || (a.chr == b.chr && a.BSPos[0].first < b.BSPos[0].first) || (a.chr == b.chr && a.BSPos[0].first == b.BSPos[0].first && a.BSPos[0].second < b.BSPos[0].second)); });
			if (!lowConfBSList.empty())
				clustCirc(lowConfBSList, nthread);
			// merge highConfBSList and lowConfBSList
			adjHighConfBSList.insert(adjHighConfBSList.end(), lowConfBSList.begin(), lowConfBSList.end());
			lowConfBSList.clear();
			sort(adjHighConfBSList.begin(), adjHighConfBSList.end(), [](const BSRecord &a, const BSRecord &b) -> bool
				 { return (a.chr < b.chr || (a.chr == b.chr && a.BSPos[0].first < b.BSPos[0].first) || (a.chr == b.chr && a.BSPos[0].first == b.BSPos[0].first && a.BSPos[0].second < b.BSPos[0].second)); });
			// cluster BS sites
			clustCirc(adjHighConfBSList, nthread);
			// output
			if (circfull::debug)
			{
				std::ofstream fout(BSAdjustFile);
				for (auto &record : adjHighConfBSList)
					fout << to_string(&record);
				fout.close();
			}
			return adjHighConfBSList;
		}

		std::vector<std::pair<int, int>> adjustExonByRef(const std::string &chr, const char strand, const std::vector<std::pair<int, int>> &origin, const ExonIndex<> &exonIndex, const GtfStorage &gtfRecords)
		{
			std::vector<std::pair<int, int>> ret;
			for (int i = 0; i < origin.size(); i++)
			{
				const std::string startKey = chr + "_" + std::to_string(origin[i].first);
				const std::string endKey = chr + "_" + std::to_string(origin[i].second);
				const std::map<std::string, std::set<std::string>>::const_iterator startIter = exonIndex.exonStartIndex.find(startKey);
				const std::map<std::string, std::set<std::string>>::const_iterator endIter = exonIndex.exonEndIndex.find(endKey);
				std::set<std::string> commonGene{};
				if (startIter != exonIndex.exonStartIndex.end() && endIter != exonIndex.exonEndIndex.end())
					std::set_intersection(startIter->second.begin(), startIter->second.end(), endIter->second.begin(), endIter->second.end(), std::inserter(commonGene, commonGene.begin()));
				for (std::set<std::string>::iterator it = commonGene.begin(); it != commonGene.end();)
				{
					const char tmpStrand = exonIndex.exonStrandIndex.find(*it)->second;
					if (tmpStrand != strand && tmpStrand != '.')
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
					ret.push_back(origin[i]);
					continue;
				}
				if (commonGene.size() == 1)
				{
					const std::string geneName = *(commonGene.begin());
					const GtfIter gtfBegin = gtfRecords.lower_bound(geneName), gtfEnd = gtfRecords.upper_bound(geneName);
					int closestExonStart = __INT_MAX__, closestExonEnd = __INT_MAX__;
					for (GtfIter it = gtfBegin; it != gtfEnd; it++)
					{
						if (closestExonStart == __INT_MAX__ || abs(it->second.start - origin[i].first) < abs(closestExonStart - origin[i].first))
							closestExonStart = it->second.start;
						if (closestExonEnd == __INT_MAX__ || abs(it->second.end - origin[i].second) < abs(closestExonEnd - origin[i].second))
							closestExonEnd = it->second.end;
					}
					ret.push_back({closestExonStart, closestExonEnd});
					continue;
				}
				int geneClosestExonStart = __INT_MAX__, geneClosestExonEnd = __INT_MAX__;
				for (std::string geneName : commonGene)
				{
					const GtfIter gtfBegin = gtfRecords.lower_bound(geneName), gtfEnd = gtfRecords.upper_bound(geneName);
					int closestExonStart = __INT_MAX__, closestExonEnd = __INT_MAX__;
					for (GtfIter it = gtfBegin; it != gtfEnd; it++)
					{ // find the closer exon boundary to the BS detected
						if (closestExonStart == __INT_MAX__ || abs(it->second.start - origin[i].first) < abs(closestExonStart - origin[i].first))
							closestExonStart = it->second.start;
						if (closestExonEnd == __INT_MAX__ || abs(it->second.end - origin[i].second) < abs(closestExonEnd - origin[i].second))
							closestExonEnd = it->second.end;
					}
					if (geneClosestExonStart == __INT_MAX__ || abs(closestExonStart - origin[i].first) + abs(closestExonEnd - origin[i].second) < abs(geneClosestExonStart - origin[i].first) + abs(geneClosestExonEnd - origin[i].second)) // if this gene explain the exon better
						geneClosestExonStart = closestExonStart, geneClosestExonEnd = closestExonEnd;
				}
				ret.push_back({geneClosestExonStart, geneClosestExonEnd});
			}
			ret.begin()->first = origin.begin()->first;
			ret.rbegin()->second = origin.rbegin()->second;
			return ret;
		}

		std::vector<CircRecord> constructFullStruct(std::vector<BSRecord> &BSList, std::vector<MapRecord> &mapList, const referenceStorage &genome, const ExonIndex<> &exonIndex, GtfStorage &gtfRecords, std::filesystem::path fullStructFile, int nthread)
		{
			// sort by queryID
			sort(BSList.begin(), BSList.end(), [](const BSRecord &a, const BSRecord &b) -> bool
				 { return (a.queryId < b.queryId); });
			sort(mapList.begin(), mapList.end(), [](const MapRecord &a, const MapRecord &b) -> bool
				 { return (a.queryId < b.queryId); });
			MapRecordIter mapIter = mapList.begin();
			MapRecordIter nowBegin, nowEnd;
			BSRecordIter thisBSIter;
			boost::asio::thread_pool tp(nthread);
			std::ofstream fout;
			if (circfull::debug)
				fout.open(fullStructFile, std::ios::out);
			std::vector<CircRecord> circList;
			std::mutex mt;
			for (size_t i = 0; i < BSList.size(); i++)
			{
				std::string nowQueryId = BSList[i].queryId;
				while (mapIter != mapList.end() && mapIter->queryId < nowQueryId)
					mapIter++;
				nowBegin = mapIter;
				while (mapIter != mapList.end() && mapIter->queryId == nowQueryId)
					mapIter++;
				nowEnd = mapIter;
				thisBSIter = BSList.begin() + i;
				boost::asio::post(tp, [&circList, &exonIndex, &gtfRecords, &fout, &mt, nowQueryId, nowBegin, nowEnd, thisBSIter]
								  {
					std::vector<MapRecord> candidateRecord = filter<MapRecord>(std::vector<MapRecord>(nowBegin, nowEnd), [&thisBSIter](const MapRecord & record)->bool {
						return (abs(record.refStart - thisBSIter->BSPos[0].first) <= hangLen && abs(record.refEnd - thisBSIter->BSPos[0].second) <= hangLen);
					});
					if (candidateRecord.size() == 0)
						return;
					// select the record with max exon number and longest exon length
					candidateRecord = maxElements<MapRecord, int>(
						candidateRecord,
						[](const MapRecord & record)->int {
							return std::count(record.exonStart.begin(), record.exonStart.end(), ',') + 1;
						}
					);
					if (candidateRecord.size() > 1)
						candidateRecord = maxElements<MapRecord, unsigned long long int> (
							candidateRecord,
							[](const MapRecord & record)->unsigned long long int {
								return record.exonLength;
							}
						);
					CircRecord ret;
					ret.queryId = nowQueryId;
					ret.chr = thisBSIter->chr;
					ret.strand = *(thisBSIter->BSstrand.begin());
					ret.start = thisBSIter->BSPos[0].first;
					ret.end = thisBSIter->BSPos[0].second;
					// parse exons
					int exonNum = std::count(candidateRecord[0].exonStart.begin(), candidateRecord[0].exonStart.end(), ',') + 1;
					if (exonNum == 1) {
						ret.exons.push_back({ret.start, ret.end});
						mt.lock();
						circList.push_back(ret);
						mt.unlock();
						return;
					}
					std::vector<std::vector<std::pair<int, int>>> exonMerge(exonNum);
					// init
					for (auto record : candidateRecord) {
						std::vector<int> exonStart = split(record.exonStart, ',');
						std::vector<int> exonEnd = split(record.exonEnd, ',');
						for (size_t i = 0; i < exonNum; i++)
							exonMerge[i].push_back({exonStart[i], exonEnd[i]});
					}
					std::pair<int, int> firstEnd = getMostFrequent<std::pair<int, int>, int>(exonMerge[0], [](const std::pair<int, int> & item)->int {
						return item.second;
					})[0];
					std::vector<std::pair<int, int>> tmpExon;
					tmpExon.push_back({ret.start, firstEnd.second});
					for (size_t i = 1; i + 1 < exonMerge.size(); i++) {
						std::pair<int, int> bestExon = getMostFrequent<std::pair<int, int>>(exonMerge[i])[0];
						bestExon.first++;
						tmpExon.push_back(bestExon);
					}
					std::pair<int, int> lastStart = getMostFrequent<std::pair<int, int>, int>(exonMerge[exonMerge.size() - 1], [](const std::pair<int, int> & item)->int {
						return item.first;
					})[0];
					tmpExon.push_back({lastStart.first + 1, ret.end});
					// remove exon with len <= 10
					for (size_t i = 0; i < tmpExon.size(); i++)
						if (tmpExon[i].second - tmpExon[i].first <= 10) {
							tmpExon.erase(tmpExon.begin() + i);
							i--;
						}
					if (tmpExon.size() == 0) {
						return;
					}
					tmpExon = adjustExonByRef(ret.chr, ret.strand, tmpExon, exonIndex, gtfRecords);
					ret.len = 0;
					for (auto exon : tmpExon)
						ret.exons.push_back(exon);//, ret.len += exon.second - exon.first + 1;
					ret.updateLength();
					mt.lock();
					circList.push_back(ret);
					if (circfull::debug)
						fout << to_string(&ret);
					mt.unlock(); });
			}
			tp.join();
			fout.close();
			return circList;
		}

		std::vector<CircRecord> clustFullStruct(std::vector<CircRecord> &circList, std::filesystem::path adjFullStructFile, int nthread)
		{
			// filter circRNA transcript with > 1 reads
			std::map<std::string, int> transcriptCount = countFrequencies<CircRecord, std::string>(
				circList,
				[](const CircRecord &record) -> std::string
				{
					std::string ret = record.chr + ":";
					std::string exonStart = std::to_string(record.exons.begin()->first);
					std::string exonEnd = std::to_string(record.exons.begin()->second);
					for (auto iter = record.exons.begin() + 1; iter != record.exons.end(); iter++)
					{
						exonStart += "," + std::to_string(iter->first);
						exonEnd += "," + std::to_string(iter->second);
					}
					ret += exonStart + "|" + exonEnd + ":" + record.strand;
					return ret;
				});
			std::vector<std::string> highConfTransIdList = filter<std::map<std::string, int>, std::pair<std::string, int>, std::string>(
				transcriptCount,
				[](const std::pair<std::string, int> &item) -> bool
				{
					return item.second > 1;
				},
				[](const std::pair<std::string, int> &item) -> std::string
				{
					return item.first;
				});
			// split into high conf transcript list and low conf list
			// high conf transcript list: transcript with > 1 UCCS
			auto [highConfCircList, lowConfCircList] = splitByBool<CircRecord>(
				circList,
				[&highConfTransIdList](const CircRecord &record) -> bool
				{
					std::string ret = record.chr + ":";
					std::string exonStart = std::to_string(record.exons.begin()->first);
					std::string exonEnd = std::to_string(record.exons.begin()->second);
					for (auto iter = record.exons.begin() + 1; iter != record.exons.end(); iter++)
					{
						exonStart += "," + std::to_string(iter->first);
						exonEnd += "," + std::to_string(iter->second);
					}
					ret += exonStart + "|" + exonEnd + ":" + record.strand;
					return (std::find(highConfTransIdList.begin(), highConfTransIdList.end(), ret) != highConfTransIdList.end());
				});
			// adjust low conf circRNA transcript by high conf circRNA transcript
			// sort by chr:start|end:exonNum
			auto circKeyComp = [](const CircRecord &a, const CircRecord &b) -> bool
			{
				std::string aStr = a.chr + ":" + std::to_string(a.start) + "|" + std::to_string(a.end) + ":" + std::to_string(a.exons.size());
				std::string bStr = b.chr + ":" + std::to_string(b.start) + "|" + std::to_string(b.end) + ":" + std::to_string(b.exons.size());
				return aStr < bStr;
			};
			std::sort(highConfCircList.begin(), highConfCircList.end(), circKeyComp);
			std::sort(lowConfCircList.begin(), lowConfCircList.end(), circKeyComp);
			std::vector<CircRecord>::iterator nowHighConf = highConfCircList.begin();
			std::vector<CircRecord>::iterator nowBegin, nowEnd;
			boost::asio::thread_pool tp(nthread);
			for (size_t i = 0; i < lowConfCircList.size(); i++)
			{
				while (nowHighConf != highConfCircList.end() && circKeyComp(*nowHighConf, lowConfCircList[i]))
					nowHighConf++;
				nowBegin = nowHighConf;
				while (nowHighConf != highConfCircList.end() && !circKeyComp(lowConfCircList[i], *nowHighConf))
					nowHighConf++;
				nowEnd = nowHighConf;
				if (nowBegin == nowEnd)
					continue;
				boost::asio::post(tp, [&lowConfCircList, nowBegin, nowEnd, i]
								  {
					size_t exonNum = lowConfCircList[i].exons.size();
					for (size_t j = 0; j < exonNum; j++) {
						int minStartDis = errorFSLen, minEndDis = errorFSLen;
						int startIndex = -1, endIndex = -1;
						for (auto iter = nowBegin; iter != nowEnd; iter++) {
							int startDis = abs(iter->exons[j].first - lowConfCircList[i].exons[j].first);
							int endDis = abs(iter->exons[j].second - lowConfCircList[i].exons[j].second);
							if (startDis < minStartDis)
								minStartDis = startDis, startIndex = iter - nowBegin;
							if (endDis < minEndDis)
								minEndDis = endDis, endIndex = iter - nowBegin;
						}
						if (startIndex != -1)
							lowConfCircList[i].exons[j].first = (nowBegin + startIndex)->exons[j].first;
						if (endIndex != -1)
							lowConfCircList[i].exons[j].second = (nowBegin + endIndex)->exons[j].second;
					} });
			}
			tp.join();
			// merge highConfCircList and lowConfCircList
			highConfCircList.insert(highConfCircList.end(), lowConfCircList.begin(), lowConfCircList.end());
			// lowConfCircList.clear();
			//  sort record by BSJ position
			std::sort(highConfCircList.begin(), highConfCircList.end(), [](const CircRecord &a, const CircRecord &b) -> bool
					  { return (a.chr < b.chr || (a.chr == b.chr && a.start < b.start) || (a.chr == b.chr && a.start == b.start && a.end < b.end)); });
			// rebuild exon list for each BSJ
			std::vector<CircRecord>::iterator endCircIter = highConfCircList.begin(), nowBeginCirc = highConfCircList.begin();
			while (endCircIter != highConfCircList.end())
			{
				while (endCircIter != highConfCircList.end() && endCircIter->chr == nowBeginCirc->chr && endCircIter->start == nowBeginCirc->start && endCircIter->end == nowBeginCirc->end)
					endCircIter++;
				// parse this bundle
				boost::asio::post(tp, [nowBeginCirc, endCircIter]
								  {
					std::vector<std::pair<int, int>> exonLists;
					for (auto iter = nowBeginCirc; iter != endCircIter; iter++)
						exonLists.insert(exonLists.end(), iter->exons.begin(), iter->exons.end());
					// count exonStart freq and link origin exonStart to new array
					std::map<int, int> exonStartFreq = countFrequencies<std::pair<int, int>, int>(
						exonLists,
						[](const std::pair<int, int> & item)->int {
							return item.first;
						}
					);
					std::map<int, int> exonStartMap;
					for (auto i : exonStartFreq)
						exonStartMap[i.first] = i.first;
					// count exonEnd freq and link origin exonEnd to new array
					std::map<int, int> exonEndFreq = countFrequencies<std::pair<int, int>, int>(
						exonLists,
						[](const std::pair<int, int> & item)->int {
							return item.second;
						}
					);
					std::map<int, int> exonEndMap;
					for (auto i : exonEndFreq)
						exonEndMap[i.first] = i.first;
					// find best exonStarts
					for (auto iterI = exonStartFreq.begin(); iterI != exonStartFreq.end(); iterI++) {
						int maxFreq = iterI->second;
						if (maxFreq == 0)
							continue;
						int nowTarget = iterI->first;
						for (auto iterJ = exonStartFreq.begin(); iterJ != exonStartFreq.end(); iterJ++)
							if (abs(iterJ->first - iterI->first) <= errorFSLen)
								if (iterJ->second >= maxFreq) {
									maxFreq = iterJ->second;
									nowTarget = iterJ->first;
								}
						if (nowTarget != iterI->first) {
							exonStartMap[iterI->first] = nowTarget;
							exonStartFreq[nowTarget] += iterI->second;
							iterI->second = 0;
						}
					}
					// find best exonEnds
					for (auto iterI = exonEndFreq.begin(); iterI != exonEndFreq.end(); iterI++) {
						int maxFreq = iterI->second;
						if (maxFreq == 0)
							continue;
						int nowTarget = iterI->first;
						for (auto iterJ = exonEndFreq.begin(); iterJ != exonEndFreq.end(); iterJ++)
							if (abs(iterJ->first - iterI->first) <= errorFSLen)
								if (iterJ->second >= maxFreq) {
									maxFreq = iterJ->second;
									nowTarget = iterJ->first;
								}
						if (nowTarget != iterI->first) {
							exonEndMap[iterI->first] = nowTarget;
							exonEndFreq[nowTarget] += iterI->second;
							iterI->second = 0;
						}
					}
					// update exons
					for (auto iter = nowBeginCirc; iter != endCircIter; iter++) {
						for (auto & exon : iter->exons) {
							exon.first = exonStartMap[exon.first];
							exon.second = exonEndMap[exon.second];
						}
						iter->updateLength();
					} });
				nowBeginCirc = endCircIter;
			}
			tp.join();
			std::ofstream adjFullStructFileOut(adjFullStructFile);
			for (auto &record : highConfCircList)
				adjFullStructFileOut << to_string(&record);
			adjFullStructFileOut.close();
			return highConfCircList;
		}

	}
}