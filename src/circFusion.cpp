#include <algorithm>
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>
#include "circCall.hpp"
#include "circfullUtils.hpp"

namespace circfull
{
	namespace RG
	{
		namespace
		{
			using LocusSplit = std::pair<std::vector<MapRecord>, std::vector<MapRecord>>;

			std::vector<std::pair<int, int>> parseExons(const MapRecord &record)
			{
				std::vector<int> exonStart = split(record.exonStart, ',');
				std::vector<int> exonEnd = split(record.exonEnd, ',');
				std::vector<std::pair<int, int>> exons;
				exons.reserve(exonStart.size());
				for (size_t i = 0; i < exonStart.size() && i < exonEnd.size(); i++)
					exons.push_back({exonStart[i], exonEnd[i]});
				return exons;
			}

			std::string joinInts(const std::vector<int> &values)
			{
				if (values.empty())
					return "";
				std::string ret = std::to_string(values.front());
				for (size_t i = 1; i < values.size(); i++)
					ret += "," + std::to_string(values[i]);
				return ret;
			}

			std::string joinStrings(const std::vector<std::string> &values)
			{
				if (values.empty())
					return "";
				std::string ret = values.front();
				for (size_t i = 1; i < values.size(); i++)
					ret += "," + values[i];
				return ret;
			}

			bool overlaps(const MapRecord &a, const MapRecord &b)
			{
				return a.refStart < b.refEnd && b.refStart < a.refEnd;
			}

			char getMajorStrand(const std::vector<MapRecord> &records)
			{
				std::map<char, int> strandCount;
				for (const auto &record : records)
					strandCount[record.strand]++;
				if (strandCount.empty())
					return '.';
				return std::max_element(strandCount.begin(), strandCount.end(), [](const auto &a, const auto &b)
									 { return a.second < b.second; })
					->first;
			}

			std::optional<LocusSplit> splitDiffChrFusionRecords(std::vector<MapRecord> records)
			{
				std::stable_sort(records.begin(), records.end(), [](const MapRecord &a, const MapRecord &b)
								 { return a.queryStart < b.queryStart; });
				std::vector<MapRecord> first, second;
				std::map<std::string, std::vector<MapRecord>> byChr;
				std::vector<std::string> chrOrder;
				for (const auto &record : records)
				{
					if (!byChr.count(record.chr))
						chrOrder.push_back(record.chr);
					byChr[record.chr].push_back(record);
				}
				if (chrOrder.size() != 2)
					return std::nullopt;
				first = byChr[chrOrder[0]];
				second = byChr[chrOrder[1]];
				if (first.empty() || second.empty())
					return std::nullopt;
				return LocusSplit{first, second};
			}

			std::optional<LocusSplit> splitSameChrFusionRecords(std::vector<MapRecord> records)
			{
				std::stable_sort(records.begin(), records.end(), [](const MapRecord &a, const MapRecord &b)
								 { return a.queryStart < b.queryStart; });
				if (records.size() < 2)
					return std::nullopt;
				auto anchor = std::max_element(records.begin(), records.end(), [](const MapRecord &a, const MapRecord &b)
									   { return a.exonLength < b.exonLength; });
				std::vector<MapRecord> overlapGroup, otherGroup;
				for (const auto &record : records)
				{
					if (overlaps(record, *anchor))
						overlapGroup.push_back(record);
					else
						otherGroup.push_back(record);
				}
				if (overlapGroup.empty() || otherGroup.empty())
				{
					overlapGroup.clear();
					otherGroup.clear();
					for (size_t i = 0; i < records.size(); i++)
						(i % 2 == 0 ? overlapGroup : otherGroup).push_back(records[i]);
				}
				if (overlapGroup.empty() || otherGroup.empty())
					return std::nullopt;
				if (overlapGroup.front().queryStart > otherGroup.front().queryStart)
					return LocusSplit{otherGroup, overlapGroup};
				return LocusSplit{overlapGroup, otherGroup};
			}

			std::optional<CircRecord> constructFusionLocus(const std::string &queryId, const std::vector<MapRecord> &locusRecords, const ExonIndex<> &exonIndex, GtfStorage &gtfRecords)
			{
				if (locusRecords.empty())
					return std::nullopt;
				std::vector<MapRecord> candidateRecord = maxElements<MapRecord, int>(
					locusRecords,
					[](const MapRecord &record) -> int
					{
						return std::count(record.exonStart.begin(), record.exonStart.end(), ',') + 1;
					});
				if (candidateRecord.size() > 1)
					candidateRecord = maxElements<MapRecord, unsigned long long int>(
						candidateRecord,
						[](const MapRecord &record) -> unsigned long long int
						{
							return record.exonLength;
						});
				if (candidateRecord.empty())
					return std::nullopt;
				std::vector<std::pair<int, int>> firstExons = parseExons(candidateRecord.front());
				if (firstExons.empty())
					return std::nullopt;
				const size_t exonNum = firstExons.size();
				std::vector<std::vector<std::pair<int, int>>> exonMerge(exonNum);
				for (const auto &record : candidateRecord)
				{
					auto exons = parseExons(record);
					if (exons.size() != exonNum)
						continue;
					for (size_t i = 0; i < exonNum; i++)
						exonMerge[i].push_back(exons[i]);
				}
				std::vector<std::pair<int, int>> consensusExons;
				consensusExons.reserve(exonNum);
				for (size_t i = 0; i < exonMerge.size(); i++)
				{
					if (exonMerge[i].empty())
						continue;
					consensusExons.push_back(getMostFrequent<std::pair<int, int>>(exonMerge[i])[0]);
				}
				if (consensusExons.empty())
					return std::nullopt;
				consensusExons.erase(std::remove_if(consensusExons.begin(), consensusExons.end(), [](const std::pair<int, int> &exon)
												   { return exon.second - exon.first <= 10; }),
					consensusExons.end());
				if (consensusExons.empty())
					return std::nullopt;
				char strand = getMajorStrand(candidateRecord);
				consensusExons = adjustExonByRef(candidateRecord.front().chr, strand, consensusExons, exonIndex, gtfRecords);
				std::sort(consensusExons.begin(), consensusExons.end());
				CircRecord locus;
				locus.queryId = queryId;
				locus.chr = candidateRecord.front().chr;
				locus.strand = strand;
				locus.exons = consensusExons;
				locus.start = consensusExons.front().first;
				locus.end = consensusExons.back().second;
				locus.updateLength();
				return locus;
			}

			std::string getGeneName(const CircRecord &record, const ExonIndex<> &exonIndex, GtfStorage &gtfRecords)
			{
				std::set<std::string> genes;
				for (const auto &exon : record.exons)
				{
					const std::string startKey = record.chr + "_" + std::to_string(exon.first);
					const std::string endKey = record.chr + "_" + std::to_string(exon.second);
					auto startIter = exonIndex.exonStartIndex.find(startKey);
					auto endIter = exonIndex.exonEndIndex.find(endKey);
					if (startIter == exonIndex.exonStartIndex.end() || endIter == exonIndex.exonEndIndex.end())
						continue;
					std::set<std::string> commonGenes;
					std::set_intersection(startIter->second.begin(), startIter->second.end(), endIter->second.begin(), endIter->second.end(), std::inserter(commonGenes, commonGenes.begin()));
					for (const auto &geneId : commonGenes)
					{
						auto strandIter = exonIndex.exonStrandIndex.find(geneId);
						if (strandIter != exonIndex.exonStrandIndex.end() && strandIter->second != '.' && strandIter->second != record.strand)
							continue;
						genes.insert(geneId);
					}
				}
				if (genes.empty())
					return "intergenic";
				std::vector<std::string> geneNames;
				for (const auto &geneId : genes)
				{
					auto range = gtfRecords.equal_range(geneId);
					std::string geneName = geneId;
					for (auto iter = range.first; iter != range.second; iter++)
					{
						if (!iter->second.geneName.empty())
						{
							geneName = iter->second.geneName;
							break;
						}
					}
					geneNames.push_back(geneName);
				}
				std::sort(geneNames.begin(), geneNames.end());
				geneNames.erase(std::unique(geneNames.begin(), geneNames.end()), geneNames.end());
				return joinStrings(geneNames);
			}

			std::tuple<std::string, std::string, std::string, std::string> getExonStrings(const CircRecord &record, const referenceStorage &genome)
			{
				std::vector<int> exonStart;
				std::vector<int> exonEnd;
				std::vector<std::string> leftSeq;
				std::vector<std::string> rightSeq;
				for (const auto &exon : record.exons)
				{
					exonStart.push_back(exon.first);
					exonEnd.push_back(exon.second);
					leftSeq.push_back(genome.getSeqByPosition(record.chr, exon.first - 2, exon.first));
					rightSeq.push_back(genome.getSeqByPosition(record.chr, exon.second, exon.second + 2));
				}
				return {joinInts(exonStart), joinInts(exonEnd), joinStrings(leftSeq), joinStrings(rightSeq)};
			}
		}

		std::string to_string(const FusionCircRecord *record)
		{
			if (record == NULL)
				return "";
			return record->queryId + "\t" + record->getCircId() + "\t" + record->getIsoformId() + "\n";
		}

		std::vector<FusionCircRecord> identifyFusionCirc(std::vector<MapRecord> &sameChrFusionInfo, std::vector<MapRecord> &diffChrFusionInfo, const referenceStorage &genome, const ExonIndex<> &exonIndex, GtfStorage &gtfRecords, int nthread)
		{
			(void)genome;
			(void)nthread;
			std::vector<FusionCircRecord> fusionList;
			auto buildFusionRecords = [&fusionList, &exonIndex, &gtfRecords](std::vector<MapRecord> &records, bool sameChr)
			{
				if (records.empty())
					return;
				std::sort(records.begin(), records.end());
				auto iter = records.begin();
				while (iter != records.end())
				{
					auto begin = iter;
					std::string queryId = iter->queryId;
					while (iter != records.end() && iter->queryId == queryId)
						iter++;
					std::vector<MapRecord> group(begin, iter);
					auto split = sameChr ? splitSameChrFusionRecords(group) : splitDiffChrFusionRecords(group);
					if (!split.has_value())
						continue;
					auto first = constructFusionLocus(queryId, split->first, exonIndex, gtfRecords);
					auto second = constructFusionLocus(queryId, split->second, exonIndex, gtfRecords);
					if (!first.has_value() || !second.has_value())
						continue;
					FusionCircRecord fusionRecord;
					fusionRecord.queryId = queryId;
					fusionRecord.first = first.value();
					fusionRecord.second = second.value();
					if (fusionRecord.first.start > fusionRecord.second.start && fusionRecord.first.chr == fusionRecord.second.chr)
						std::swap(fusionRecord.first, fusionRecord.second);
					fusionList.push_back(fusionRecord);
				}
			};
			buildFusionRecords(diffChrFusionInfo, false);
			buildFusionRecords(sameChrFusionInfo, true);
			return fusionList;
		}

		void outputFusionCirc(std::vector<FusionCircRecord> &fusionList, const referenceStorage &genome, const ExonIndex<> &exonIndex, GtfStorage &gtfRecords, std::filesystem::path output)
		{
			std::ofstream fout(output, std::ios::out);
			fout << "circID\tisoID\tchr_first\tstart_first\tend_first\tlen_first\texonNum_first\texon_start_first\texon_end_first\texon_leftSeq_first\texon_rightSeq_first\tstrand_first\tgeneName_first\tchr_second\tstart_second\tend_second\tlen_second\texonNum_second\texon_start_second\texon_end_second\texon_leftSeq_second\texon_rightSeq_second\tstrand_second\tgeneName_second\treadCount\treadID\n";
			std::map<std::string, std::vector<FusionCircRecord>> isoformGroups;
			for (auto &record : fusionList)
				isoformGroups[record.getIsoformId()].push_back(record);
			std::vector<std::string> isoIds;
			isoIds.reserve(isoformGroups.size());
			for (const auto &[isoId, _] : isoformGroups)
				isoIds.push_back(isoId);
			std::sort(isoIds.begin(), isoIds.end());
			for (const auto &isoId : isoIds)
			{
				auto &records = isoformGroups[isoId];
				auto &record = records.front();
				auto [exonStartFirst, exonEndFirst, exonLeftFirst, exonRightFirst] = getExonStrings(record.first, genome);
				auto [exonStartSecond, exonEndSecond, exonLeftSecond, exonRightSecond] = getExonStrings(record.second, genome);
				std::vector<std::string> readIds;
				readIds.reserve(records.size());
				for (const auto &item : records)
					readIds.push_back(item.queryId);
				std::sort(readIds.begin(), readIds.end());
				fout << record.getCircId() << "\t"
					 << isoId << "\t"
					 << record.first.chr << "\t"
					 << record.first.start << "\t"
					 << record.first.end << "\t"
					 << record.first.len << "\t"
					 << record.first.exons.size() << "\t"
					 << exonStartFirst << "\t"
					 << exonEndFirst << "\t"
					 << exonLeftFirst << "\t"
					 << exonRightFirst << "\t"
					 << record.first.strand << "\t"
					 << getGeneName(record.first, exonIndex, gtfRecords) << "\t"
					 << record.second.chr << "\t"
					 << record.second.start << "\t"
					 << record.second.end << "\t"
					 << record.second.len << "\t"
					 << record.second.exons.size() << "\t"
					 << exonStartSecond << "\t"
					 << exonEndSecond << "\t"
					 << exonLeftSecond << "\t"
					 << exonRightSecond << "\t"
					 << record.second.strand << "\t"
					 << getGeneName(record.second, exonIndex, gtfRecords) << "\t"
					 << readIds.size() << "\t"
					 << joinStrings(readIds) << "\n";
			}
			fout.close();
		}
	}
}