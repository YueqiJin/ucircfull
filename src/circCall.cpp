#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include "circfullUtils.hpp"
#include "circCall.hpp"
#include "ucircfull.hpp"
#include "align.hpp"

/// @brief read UMI cluster information from umiClustTab file
/// @param umiFile path to umiClustTab file
/// @return map<string, string> umiInfo: readId->clustId
std::map<std::string, std::string> readUMIClustInfo(const std::filesystem::path &umiFile)
{
	std::ifstream umiClustFile(umiFile);
	std::map<std::string, std::string> umiInfo;
	// read umi cluster info from umiClustFile: clust_id, read_num, length, read_id, compare_rate, separated by tab
	std::string umiClustId, readNum, length, readId, compareRate;
	std::string line;
	while (std::getline(umiClustFile, line))
	{
		std::istringstream iss(line);
		iss >> umiClustId >> readNum >> length >> readId >> compareRate;
		umiInfo[readId] = umiClustId;
	}
	return umiInfo;
}

namespace circfull
{

	void outputCirc(std::vector<circfull::RG::CircRecord> &circLists, std::filesystem::path output, int nthread)
	{
		std::ofstream circGtf(output, std::ios::out);
		// sort by circId and isoforms
		sort(
			circLists.begin(),
			circLists.end(),
			[](circfull::RG::CircRecord &a, circfull::RG::CircRecord &b)
			{
				std::string circIdA = a.getCircId(), circIdB = b.getCircId();
				if (circIdA == circIdB)
					return a.getExonsString() < b.getExonsString();
				return circIdA < circIdB;
			});
		// output circRNA gtf
		std::vector<circfull::RG::CircRecord>::iterator nowEnd = circLists.begin(), nowBegin = circLists.begin();
		std::string nowCircId = nowBegin->getCircId();
		boost::asio::thread_pool tp(nthread);
		std::mutex mt;
		while (nowEnd != circLists.end())
		{
			// for circRNA BSJ
			while (nowEnd != circLists.end() && nowEnd->getCircId() == nowCircId)
				nowEnd++;
			boost::asio::post(tp, [&mt, &circGtf, nowBegin, nowEnd, nowCircId]()
							  {
				std::string nowTranscriptId = nowBegin->getTranscriptId();
				std::vector<circfull::RG::CircRecord>::iterator it = nowBegin, transcriptBegin = nowBegin;
				int totalCount = 0;
				std::string bsjAnnoAttr = " circ_type \"" + nowBegin->circ_type + "\"; host_gene_id \"" + nowBegin->host_gene_id + "\"; host_gene_name \"" + nowBegin->host_gene_name + "\";";
				std::string transcriptList = "";
				while (it != nowEnd) {
					while (it != nowEnd && it->getTranscriptId() == nowTranscriptId)
						it++;
					if (it - transcriptBegin > 1) {
						transcriptList += transcriptBegin->chr + "\tcircfull\ttranscript\t" + std::to_string(transcriptBegin->start) + "\t" + std::to_string(transcriptBegin->end) + "\t.\t" + transcriptBegin->strand + "\t.\tgene_id \"" + transcriptBegin->getCircId(false) + "\"; transcript_id \"" + nowTranscriptId + "\"; uniform_id \"" + transcriptBegin->uniform_id + "\";" + bsjAnnoAttr + " bsj \"" + std::to_string(it - transcriptBegin) + "\";\n";
						auto exons = transcriptBegin->exons;
						for (size_t i = 0; i < exons.size(); i++) {
							transcriptList += transcriptBegin->chr + "\tcircfull\texon\t" + std::to_string(exons[i].first) + "\t" + std::to_string(exons[i].second) + "\t.\t" + transcriptBegin->strand + "\t.\tgene_id \"" + transcriptBegin->getCircId(false) + "\"; transcript_id \"" + nowTranscriptId + "\"; uniform_id \"" + transcriptBegin->uniform_id + "\"; exon_number \"" + std::to_string(i + 1) + "\";\n";
						}
						totalCount += it - transcriptBegin;
					}
					if (it != nowEnd) {
						transcriptBegin = it;
						nowTranscriptId = it->getTranscriptId();
					}
				}
				std::string circBSJLine = nowBegin->chr + "\tcircfull\tBSJ\t" + std::to_string(nowBegin->start) + "\t" + std::to_string(nowBegin->end) + "\t.\t" + nowBegin->strand + "\t.\tgene_id \"" + nowBegin->getCircId(false) + "\";" + bsjAnnoAttr + " bsj \"" + std::to_string(totalCount) + "\";\n";
				if (totalCount > 1) {
					mt.lock();
					circGtf << circBSJLine + transcriptList;
					mt.unlock();
				} });
			if (nowEnd != circLists.end())
			{
				nowBegin = nowEnd;
				nowCircId = nowBegin->getCircId();
			}
		}
		tp.join();
		circGtf.close();
	}

	int circ_call_cpp(circfull::CircfullOption &opt)
	{
		// Output running options
		std::cout << "[INFO] Running mod: ";
		switch (opt.modCircCall)
		{
		case circfull::CircfullOption::CircCallMod::ucRG:
			std::cout << "ucRG" << std::endl;
			break;
		case circfull::CircfullOption::CircCallMod::RG:
			std::cout << "RG" << std::endl;
			break;
		case circfull::CircfullOption::CircCallMod::cRG:
			std::cout << "cRG" << std::endl;
		}
		std::cout << "[INFO] UMI cluster file: " << opt.umiClustTab << std::endl;
		std::cout << "[INFO] Input stranded fastq file: " << opt.strandFastq << std::endl;
		std::cout << "[INFO] Reference fasta file: " << opt.referenceFasta << std::endl;
		std::cout << "[INFO] Reference annotation file: " << opt.annotationGTF << std::endl;
		if (opt.modCircCall == CircfullOption::CircCallMod::CIRI) {
			std::cout << "[INFO] Output path: " << opt.circCallOutputPath << std::endl;
		} else {
			std::cout << "[INFO] Output path: " << opt.oPath << std::endl;
		}
		if (!opt.useBam)
		{
			// Checking input files exist
			if (!checkFile(opt.strandFastq))
			{
				std::cerr << "Error: " << opt.strandFastq << " does not exist" << std::endl;
				return 1;
			}
			if (!checkFile(opt.referenceFasta))
			{
				std::cerr << "Error: " << opt.referenceFasta << " does not exist" << std::endl;
				return 1;
			}

			if (opt.modCircCall == CircfullOption::CircCallMod::ucRG) {
				if (!checkFile(opt.umiClustTab))
				{
					std::cerr << "Error: " << opt.umiClustTab << " does not exist" << std::endl;
					return 1;
				}

				// read in UMI cluster file
				printTimeInfo("Scaning candidate circRNA reads:");
				std::map<std::string, std::string> umiCluster = readUMIClustInfo(opt.umiClustTab);

				// scan for ccs reads
				std::map<std::string, std::tuple<std::string, std::string, std::string>> ccsSeq = circfull::findCcsReads(opt.strandFastq, opt.oPath / (opt.oPrefix + "_ccs"), opt.nthread);

				// UMI guided ccs calling
				printTimeInfo("UMI guided ccs calling:");
				std::map<std::string, std::string> ccsConsensus = circfull::constructUCCSConsensus(ccsSeq, umiCluster, opt.nthread, opt.oPath);

				// write consensus sequence (doubled) to file
				opt.ccsFasta = opt.oPath / (opt.oPrefix + "_RG_pseudo.fa");
				std::ofstream consensusFasta(opt.ccsFasta, std::ios::out);
				for (auto &it : ccsConsensus)
				{
					consensusFasta << ">" << it.first << std::endl;
					consensusFasta << it.second << it.second << it.second << std::endl;
				}
				consensusFasta.close();
			}
			if (opt.modCircCall ==  CircfullOption::CircCallMod::cRG) {
				// read in UMI cluster file
				printTimeInfo("Scaning candidate circRNA reads:");
				// scan for ccs reads
				std::map<std::string, std::tuple<std::string, std::string, std::string>> ccsSeq = circfull::findCcsReads(opt.strandFastq, opt.oPath / (opt.oPrefix + "_ccs"), opt.nthread);
				// write consensus sequence (doubled) to file
				opt.ccsFasta = opt.oPath / (opt.oPrefix + "_RG_pseudo.fa");
				std::ofstream consensusFasta(opt.ccsFasta, std::ios::out);
				for (auto &it : ccsSeq)
				{
					consensusFasta << ">" << it.first << std::endl;
					consensusFasta << get<1>(it.second) << get<1>(it.second) << get<1>(it.second) << std::endl;
				}
				consensusFasta.close();
			}
			if (opt.modCircCall == CircfullOption::CircCallMod::RG) {
					opt.ccsFasta = opt.strandFastq;
			}
		}
		else
		{
			std::cout << "[INFO] Use provided BAM file: " << opt.bamFile << std::endl;
			opt.ccsFasta = opt.strandFastq;
		}
		// RG
		printTimeInfo("RG:");
		std::vector<circfull::RG::CircRecord> circLists = circfull::RG::RG(opt);

		// output circRNA gtf
		printTimeInfo("Output circRNA gtf.");
		std::filesystem::path outputGtfPath{opt.oPath / (opt.oPrefix + ".circ.gtf")};
		outputCirc(circLists, outputGtfPath, opt.nthread);
		printTimeInfo("Done.");
		return 0;
	}

	namespace RG
	{
		const int ANNO_BOUNDARY_TOL = 10;

		AnnotationIndex::AnnotationIndex(const GtfStorage &gtfRecords, const ExonIndex<> &_exonIndex) : exonIndex(_exonIndex)
		{
			for (const auto &[geneId, record] : gtfRecords)
			{
				if (record.type != GtfRecord::FeatureType::exon)
					continue;
				auto gi = geneInfo.find(record.geneId);
				if (gi == geneInfo.end())
				{
					geneInfo[record.geneId] = {record.chr, record.geneName, record.strand, record.start, record.end};
					chrToGenes[record.chr].push_back(record.geneId);
				}
				else
				{
					gi->second.geneStart = std::min(gi->second.geneStart, record.start);
					gi->second.geneEnd = std::max(gi->second.geneEnd, record.end);
					if (gi->second.geneName.empty() && !record.geneName.empty())
						gi->second.geneName = record.geneName;
				}
				auto &transVec = geneTranscripts[record.geneId];
				auto it = std::find_if(transVec.begin(), transVec.end(), [&](const TranscriptExonInfo &t)
									   { return t.transcriptId == record.transcriptId; });
				if (it == transVec.end())
					transVec.push_back({record.transcriptId, {{record.start, record.end}}});
				else
					it->exons.push_back({record.start, record.end});
			}
			for (auto &[geneId, transVec] : geneTranscripts)
				for (auto &t : transVec)
					std::sort(t.exons.begin(), t.exons.end());
			for (auto &[chr, genes] : chrToGenes)
			{
				std::sort(genes.begin(), genes.end());
				genes.erase(std::unique(genes.begin(), genes.end()), genes.end());
			}
		}

		namespace
		{
			int getExonNum5p3p(const TranscriptExonInfo &trans, int idx, char strand)
			{
				if (strand == '-')
					return static_cast<int>(trans.exons.size()) - idx;
				return idx + 1;
			}

			bool isRetainedIntron(const std::pair<int, int> &circExon, const TranscriptExonInfo &trans)
			{
				if (trans.exons.size() < 2)
					return false;
				for (size_t i = 0; i < trans.exons.size() - 1; i++)
				{
					int intronStart = trans.exons[i].second + 1;
					int intronEnd = trans.exons[i + 1].first - 1;
					if (intronStart > intronEnd)
						continue;
					if (circExon.first >= intronStart - ANNO_BOUNDARY_TOL && circExon.second <= intronEnd + ANNO_BOUNDARY_TOL &&
						std::abs(circExon.first - intronStart) <= ANNO_BOUNDARY_TOL &&
						std::abs(circExon.second - intronEnd) <= ANNO_BOUNDARY_TOL)
						return true;
				}
				return false;
			}

			int findMatchingRefExon(const std::pair<int, int> &circExon, const TranscriptExonInfo &trans)
			{
				int bestIdx = -1, bestOverlap = 0;
				for (size_t i = 0; i < trans.exons.size(); i++)
				{
					int overlap = std::min(circExon.second, trans.exons[i].second) - std::max(circExon.first, trans.exons[i].first);
					if (overlap > bestOverlap && overlap > 0)
					{
						bestOverlap = overlap;
						bestIdx = static_cast<int>(i);
					}
				}
				return bestIdx;
			}

			std::string generateExonToken(const std::pair<int, int> &circExon, const std::pair<int, int> &refExon,
										 int exonNum, char strand)
			{
				bool threePrimeLonger = false;
				bool fivePrimeShorter = false;
				if (strand == '+')
				{
					if (circExon.second > refExon.second + ANNO_BOUNDARY_TOL)
						threePrimeLonger = true;
					if (circExon.first > refExon.first + ANNO_BOUNDARY_TOL)
						fivePrimeShorter = true;
				}
				else
				{
					if (circExon.first < refExon.first - ANNO_BOUNDARY_TOL)
						threePrimeLonger = true;
					if (circExon.second < refExon.second - ANNO_BOUNDARY_TOL)
						fivePrimeShorter = true;
				}
				std::string token = std::to_string(exonNum);
				if (threePrimeLonger)
					token = "L" + token;
				if (fivePrimeShorter)
					token = token + "S";
				return token;
			}
		}

		AnnotationResult annotateOneCirc(const CircRecord &circ, const AnnotationIndex &idx)
		{
			AnnotationResult result;
			const int searchRange = errorBSLen;

			auto chrIt = idx.chrToGenes.find(circ.chr);
			if (chrIt == idx.chrToGenes.end())
			{
				result.circ_type = "intergenic_region";
				return result;
			}

			const std::string bsjStartKey = circ.chr + "_" + std::to_string(circ.start);
			const std::string bsjEndKey = circ.chr + "_" + std::to_string(circ.end);
			auto bsjSStart = idx.exonIndex.exonStartIndex.find(bsjStartKey);
			auto bsjEEnd = idx.exonIndex.exonEndIndex.find(bsjEndKey);
			auto bsjSEnd = idx.exonIndex.exonEndIndex.find(bsjStartKey);
			auto bsjEStart = idx.exonIndex.exonStartIndex.find(bsjEndKey);

			struct ExonLK
			{
				const std::set<std::string> *s1, *e1, *s2, *e2;
			};
			std::vector<ExonLK> exonLK;
			for (const auto &exon : circ.exons)
			{
				std::string sk = circ.chr + "_" + std::to_string(exon.first);
				std::string ek = circ.chr + "_" + std::to_string(exon.second);
				auto ss1 = idx.exonIndex.exonStartIndex.find(sk);
				auto ee1 = idx.exonIndex.exonEndIndex.find(ek);
				auto se2 = idx.exonIndex.exonEndIndex.find(sk);
				auto es2 = idx.exonIndex.exonStartIndex.find(ek);
				exonLK.push_back({
					ss1 != idx.exonIndex.exonStartIndex.end() ? &ss1->second : nullptr,
					ee1 != idx.exonIndex.exonEndIndex.end() ? &ee1->second : nullptr,
					se2 != idx.exonIndex.exonEndIndex.end() ? &se2->second : nullptr,
					es2 != idx.exonIndex.exonStartIndex.end() ? &es2->second : nullptr});
			}

			auto checkBSJ = [&](const auto *a, const auto *b, const std::string &gid) -> bool
			{
				if (!a || !b)
					return false;
				return a->find(gid) != a->end() && b->find(gid) != b->end();
			};

			bool foundExon = false, foundIntron = false;
			int bestGeneScore = -1;

			for (const auto &geneId : chrIt->second)
			{
				auto gi = idx.geneInfo.find(geneId);
				if (gi == idx.geneInfo.end())
					continue;
				const auto &gene = gi->second;

				if (circ.strand != '.' && gene.strand != '.' && gene.strand != circ.strand)
					continue;
				if (circ.end < gene.geneStart - searchRange || circ.start > gene.geneEnd + searchRange)
					continue;

				if (!foundExon)
				{
					bool bs1 = (bsjSStart != idx.exonIndex.exonStartIndex.end() && bsjEEnd != idx.exonIndex.exonEndIndex.end() &&
							   checkBSJ(&bsjSStart->second, &bsjEEnd->second, geneId));
					bool bs2 = (bsjSEnd != idx.exonIndex.exonEndIndex.end() && bsjEStart != idx.exonIndex.exonStartIndex.end() &&
							   checkBSJ(&bsjSEnd->second, &bsjEStart->second, geneId));
					if (bs1 || bs2)
					{
						foundExon = true;
					}
					else
					{
						bool startIn = false, endIn = false;
						auto gIt = idx.geneTranscripts.find(geneId);
						if (gIt != idx.geneTranscripts.end())
						{
							for (const auto &trans : gIt->second)
							{
								for (const auto &ref : trans.exons)
								{
									if (!startIn && circ.start >= ref.first - ANNO_BOUNDARY_TOL && circ.start <= ref.second + ANNO_BOUNDARY_TOL)
										startIn = true;
									if (!endIn && circ.end >= ref.first - ANNO_BOUNDARY_TOL && circ.end <= ref.second + ANNO_BOUNDARY_TOL)
										endIn = true;
								}
								if (startIn && endIn)
									break;
							}
						}
						if (startIn && endIn)
							foundExon = true;
					}
				}
				if (!foundExon && circ.start >= gene.geneStart - searchRange && circ.end <= gene.geneEnd + searchRange)
					foundIntron = true;

				int score = 0;
				if (gene.strand == circ.strand)
					score += 1000;
				for (const auto &lk : exonLK)
				{
					if (lk.s1 && lk.s1->find(geneId) != lk.s1->end())
						score += 10;
					if (lk.e1 && lk.e1->find(geneId) != lk.e1->end())
						score += 10;
					if (lk.s2 && lk.s2->find(geneId) != lk.s2->end())
						score += 10;
					if (lk.e2 && lk.e2->find(geneId) != lk.e2->end())
						score += 10;
				}
				auto gIt = idx.geneTranscripts.find(geneId);
				if (gIt != idx.geneTranscripts.end())
				{
					for (const auto &trans : gIt->second)
						for (const auto &ref : trans.exons)
						{
							for (const auto &exon : circ.exons)
							{
								int overlap = std::min(exon.second, ref.second) - std::max(exon.first, ref.first);
								if (overlap > 0)
									score += 1 + overlap;
							}
						}
				}
				if (score > bestGeneScore)
				{
					bestGeneScore = score;
					result.host_gene_id = geneId;
					result.host_gene_name = gene.geneName;
					if (gIt != idx.geneTranscripts.end())
					{
						int bestTS = -1;
						for (const auto &trans : gIt->second)
						{
							int ts = 0;
							for (const auto &exon : circ.exons)
							{
								for (const auto &ref : trans.exons)
								{
									if (std::abs(exon.first - ref.first) <= ANNO_BOUNDARY_TOL && std::abs(exon.second - ref.second) <= ANNO_BOUNDARY_TOL)
									{
										ts += 20;
										break;
									}
									int overlap = std::min(exon.second, ref.second) - std::max(exon.first, ref.first);
									if (overlap > 0)
										ts += 1 + overlap / 10;
								}
							}
							if (ts > bestTS)
							{
								bestTS = ts;
								result.bestTranscriptId = trans.transcriptId;
							}
						}
					}
				}
			}

			result.circ_type = foundExon ? "exon" : (foundIntron ? "intron" : "intergenic_region");
			return result;
		}

		std::string generateUniformId(const CircRecord &circ, const AnnotationIndex &idx, const std::string &bestTranscriptId)
		{
			if (circ.host_gene_id.empty())
				return "circIntergenic(" + circ.chr + ":" + std::to_string(circ.start) + "-" + std::to_string(circ.end) + ":" + circ.strand + ")";
			std::string geneName = circ.host_gene_name;
			if (geneName.empty())
				geneName = circ.host_gene_id;
			const TranscriptExonInfo *bestTrans = nullptr;
			auto gIt = idx.geneTranscripts.find(circ.host_gene_id);
			if (gIt != idx.geneTranscripts.end())
			{
				for (const auto &trans : gIt->second)
				{
					if (trans.transcriptId == bestTranscriptId)
					{
						bestTrans = &trans;
						break;
					}
				}
			}
			std::vector<std::string> tokens;
			auto exonList = circ.exons;
			if (circ.strand == '-')
				std::reverse(exonList.begin(), exonList.end());
			for (const auto &exon : exonList)
			{
				if (!bestTrans)
				{
					tokens.push_back("NE");
					continue;
				}
				int matchIdx = findMatchingRefExon(exon, *bestTrans);
				if (matchIdx < 0)
				{
					if (isRetainedIntron(exon, *bestTrans))
						tokens.push_back("RI");
					else
						tokens.push_back("NE");
				}
				else
				{
					int exonNum = getExonNum5p3p(*bestTrans, matchIdx, circ.strand);
					tokens.push_back(generateExonToken(exon, bestTrans->exons[matchIdx], exonNum, circ.strand));
				}
			}
			std::string id = "circ" + geneName + "(";
			for (size_t i = 0; i < tokens.size(); i++)
			{
				if (i > 0)
					id += ",";
				id += tokens[i];
			}
			id += ")";
			return id;
		}

		void annotateCircRNAs(std::vector<CircRecord> &circRecords, const GtfStorage &gtfRecords,
							  const ExonIndex<> &exonIndex, int nthread)
		{
			AnnotationIndex idx(gtfRecords, exonIndex);
			boost::asio::thread_pool tp(nthread);
			for (size_t i = 0; i < circRecords.size(); i++)
				boost::asio::post(tp, [&, i]()
								  {
					CircRecord &circ = circRecords[i];
					auto result = annotateOneCirc(circ, idx);
					circ.circ_type = result.circ_type;
					circ.host_gene_id = result.host_gene_id;
					circ.host_gene_name = result.host_gene_name;
					circ.uniform_id = generateUniformId(circ, idx, result.bestTranscriptId); });
			tp.join();
		}

		std::vector<circfull::RG::CircRecord> RG(CircfullOption &opt)
		{
			// align ccs pseudoSeq to the reference
			printTimeInfo("Align pseudo Fa to reference genome:");
			if (!opt.useBam)
			{
				makeDir(opt.oPath / "RG");
				std::filesystem::path samFile{opt.oPath / "RG" / (opt.oPrefix + ".minimap2.sam")};
				opt.bamFile = opt.oPath / "RG" / (opt.oPrefix + ".minimap2.bam");
				minimap2Align(opt.ccsFasta, opt.referenceFasta, samFile, opt.minimap2, opt.nthread, opt.strand);
				printTimeInfo("Transform SAM to BAM.\n");
				sam2bam(samFile, opt.bamFile, opt.samtools, opt.nthread);
			}
			printTimeInfo("Reading sequences.");
			referenceStorage genome{readReference(opt.referenceFasta)};
			printTimeInfo("Analysis bam file.");
			std::filesystem::path mapInfo{opt.oPath / "RG" / (opt.oPrefix + ".info.txt")};
			std::vector<MapRecord> mapRecords = parseBam(genome, mapInfo, opt.bamFile, opt.nthread);
			printTimeInfo("Filter candidate fusion circRNAs.");
			std::filesystem::path typeInfoPrefix{opt.oPath / "RG" / (opt.oPrefix + ".type")};
			auto [circCandidateInfo, sameChrFusionInfo, diffChrFusionInfo]{filterCandidateCircReads(mapRecords, genome, opt.ccsFasta, typeInfoPrefix, opt.nthread)};
			printTimeInfo("Detect back-splicing junction.");
			std::filesystem::path BSListFile{opt.oPath / "RG" / (opt.oPrefix + ".BS.txt")};
			std::vector<BSRecord> BSList = detectBS(genome, BSListFile, circCandidateInfo, opt.nthread, opt.oPrefix, opt.minimap2);
			printTimeInfo("Reading annotation file.");
			GtfStorage gtfRecords{readGtf(opt.annotationGTF)};
			ExonIndex<> exonBSIndex = getExonIndex(gtfRecords, errorBSLen);
			ExonIndex<> exonFSIndex = getExonIndex(gtfRecords, errorFSLen);
			printTimeInfo("Identify fusion circRNAs.");
			std::vector<FusionCircRecord> fusionList = identifyFusionCirc(sameChrFusionInfo, diffChrFusionInfo, genome, exonFSIndex, gtfRecords, opt.nthread);
			std::filesystem::path fusionOutputFile{opt.oPath / (opt.oPrefix + ".fusion.txt")};
			outputFusionCirc(fusionList, genome, exonFSIndex, gtfRecords, fusionOutputFile);
			printTimeInfo("filter BS signal according to annotation file.");
			std::filesystem::path BSFilteredListFile{opt.oPath / "RG" / (opt.oPrefix + ".BS.filtered.txt")};
			BSList = filterBS(genome, BSList, gtfRecords, exonBSIndex, BSFilteredListFile, opt.spliceSignal, opt.nthread);
			printTimeInfo("Clustering BS sites");
			std::filesystem::path BSAdjustFile{opt.oPath / "RG" / (opt.oPrefix + ".BS.adj.txt")};
			BSList = clusterBS(BSList, BSAdjustFile, opt.nthread);
			printTimeInfo("Inferring internal structure");
			std::filesystem::path fullStructFile{opt.oPath / "RG" / (opt.oPrefix + ".full_struct.txt")};
			std::vector<CircRecord> circList{constructFullStruct(BSList, circCandidateInfo, genome, exonFSIndex, gtfRecords, fullStructFile, opt.nthread)};
			std::filesystem::path adjFullStructFile{opt.oPath / "RG" / (opt.oPrefix + ".adj_full_struct.txt")};
			circList = clustFullStruct(circList, adjFullStructFile, opt.nthread);
			if (opt.circNumThres > 1) {
				std::map<std::string, int> tcounts;
				for (auto &circ : circList)
					tcounts[circ.getTranscriptId()]++;
				circList.erase(std::remove_if(circList.begin(), circList.end(),
					[&](const CircRecord &c) { return tcounts[c.getTranscriptId()] < opt.circNumThres; }),
					circList.end());
			}
			printTimeInfo("Annotating circRNAs.");
			annotateCircRNAs(circList, gtfRecords, exonFSIndex, opt.nthread);
			{
				std::map<std::pair<std::string, std::string>, int> bsjsNum;
				std::map<std::string, std::vector<std::string>> uniGroups;
				for (auto &circ : circList)
					uniGroups[circ.uniform_id].push_back(circ.getCircId());
				for (auto &[base, bsjs] : uniGroups)
				{
					std::sort(bsjs.begin(), bsjs.end());
					bsjs.erase(std::unique(bsjs.begin(), bsjs.end()), bsjs.end());
					for (size_t n = 0; n < bsjs.size(); n++)
						bsjsNum[{base, bsjs[n]}] = static_cast<int>(n + 1);
				}
				for (auto &circ : circList)
				{
					auto it = bsjsNum.find({circ.uniform_id, circ.getCircId()});
					if (it != bsjsNum.end())
						circ.uniform_id = circ.uniform_id + "." + std::to_string(it->second);
				}
			}
			return circList;
		}
	}
}
