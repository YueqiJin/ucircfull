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
			// for circRNA BSJ (gene)
			while (nowEnd != circLists.end() && nowEnd->getCircId() == nowCircId)
				nowEnd++;
			boost::asio::post(tp, [&mt, &circGtf, nowBegin, nowEnd, nowCircId]()
							  {
				std::string nowTranscriptId = nowBegin->getTranscriptId();
				std::vector<circfull::RG::CircRecord>::iterator it = nowBegin, transcriptBegin = nowBegin;
				int totalCount = 0;
				// transcript list
				std::string transcriptList = "";
				//TODO: add circRNA type: exonic, intronic, intergenic etc.
				while (it != nowEnd) {
					while (it != nowEnd && it->getTranscriptId() == nowTranscriptId)
						it++;
					if (it - transcriptBegin > 1) {
						transcriptList += transcriptBegin->chr + "\tcircfull\ttranscript\t" + std::to_string(transcriptBegin->start) + "\t" + std::to_string(transcriptBegin->end) + "\t.\t" + transcriptBegin->strand + "\t.\tgene_id \"" + nowCircId + "\"; transcript_id \"" + nowTranscriptId + "\"; bsj \"" + std::to_string(it - transcriptBegin) + "\";\n";
						// exon list
						auto exons = transcriptBegin->exons;
						for (size_t i = 0; i < exons.size(); i++) {
							transcriptList += transcriptBegin->chr + "\tcircfull\texon\t" + std::to_string(exons[i].first) + "\t" + std::to_string(exons[i].second) + "\t.\t" + transcriptBegin->strand + "\t.\tgene_id \"" + nowCircId + "\"; transcript_id \"" + nowTranscriptId + "\"; exon_number \"" + std::to_string(i + 1) + "\";\n";
						}
						totalCount += it - transcriptBegin;
					}
					if (it != nowEnd) {
						transcriptBegin = it;
						nowTranscriptId = it->getTranscriptId();
					}
				}
				std::string circBSJLine = nowBegin->chr + "\tcircfull\tgene\t" + std::to_string(nowBegin->start) + "\t" + std::to_string(nowBegin->end) + "\t.\t" + nowBegin->strand + "\t.\tgene_id \"" + nowBegin->getCircId() + "\"; bsj \"" + std::to_string(totalCount) + "\";\n";
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
		std::cout << "[INFO] Output path: " << opt.circCallOutputPath << std::endl;
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

			if (opt.modCircCall == CircfullOption::CircCallMod::ucRG)
			{
				// read in UMI cluster file
				printTimeInfo("Scaning candidate circRNA reads:");
				std::map<std::string, std::string> umiCluster = readUMIClustInfo(opt.umiClustTab);

				// scan for ccs reads
				std::map<std::string, std::tuple<std::string, std::string, std::string>> ccsSeq = circfull::findCcsReads(opt.strandFastq, opt.circCallOutputPath, opt.nthread);

				// UMI guided ccs calling
				printTimeInfo("UMI guided ccs calling:");
				std::map<std::string, std::string> ccsConsensus = circfull::constructCCSConsensus(ccsSeq, umiCluster, opt.nthread, opt.circCallOutputPath);

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
			else
			{
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
			printTimeInfo("filter BS signal according to annotation file.");
			std::filesystem::path BSFilteredListFile{opt.oPath / "RG" / (opt.oPrefix + ".BS.filtered.txt")};
			BSList = filterBS(genome, BSList, gtfRecords, exonBSIndex, BSFilteredListFile, opt.spliceSignal, opt.nthread);
			printTimeInfo("Clustering BS sites");
			std::filesystem::path BSAdjustFile{opt.oPath / "RG" / (opt.oPrefix + ".BS.adj.txt")};
			BSList = clusterBS(BSList, BSAdjustFile, opt.nthread);
			printTimeInfo("Inferring internal structure");
			std::filesystem::path fullStructFile{opt.oPath / "RG" / (opt.oPrefix + ".full_struct.txt")};
			ExonIndex<> exonFSIndex = getExonIndex(gtfRecords, errorFSLen);
			std::vector<CircRecord> circList{constructFullStruct(BSList, circCandidateInfo, genome, exonFSIndex, gtfRecords, fullStructFile, opt.nthread)};
			std::filesystem::path adjFullStructFile{opt.oPath / "RG" / (opt.oPrefix + ".adj_full_struct.txt")};
			circList = clustFullStruct(circList, adjFullStructFile, opt.nthread);
			return circList;
		}
	}
}
