#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <set>
#include <thread>
#include <mutex>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <string>
#include <cstring>
#include <map>
#include <queue>
#include <filesystem>
#include <utility>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alignment/all.hpp>
#include <seqan3/alphabet/views/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include "circfullUtils.hpp"
#include "align.hpp"
#include "ucircfull.hpp"

std::string umi{"CTCNNNYRNNNYRNNNYRNNNGAG"};
std::string polyA{"AAAAAAAAAAAAAAAAAAAAAAAACTC"};
seqan3::dna4_vector anchorx, anchorxRev;

int umiLength;
int keepCountSum{0}, discardCountSum{0};
std::mutex mt;

// extract umi sequence from the string seq
// return umi, start of umi, end of umi, score
std::optional<std::tuple<std::string, int, int, int>> extract(const DegenerateBaseVector &seq, const DegenerateBaseVector &umiSeq)
{
	int l1 = seq.size();
	int alignMatrix[40][40];
	short path[40][40];
	degenerateAlign(seq, umiSeq, l1, umiLength, *alignMatrix, *path, defaultAlignOption, false, true);
	AlignPatternResult aln = matchSeq(seq, umiSeq, *alignMatrix, *path, l1, umiLength, false, true);
	std::size_t l1Mapped = aln.mappedString.length();
	if (l1Mapped > static_cast<int>(round(umiLength * 1.2)) || l1Mapped < static_cast<int>(round(umiLength * 0.8)))
	{
		return std::nullopt;
	}
	return std::tuple(aln.mappedString, aln.seq1Site.first, aln.seq1Site.second, aln.score);
}

seqan3::nucleotide_scoring_scheme scoreMat{seqan3::match_score{3}, seqan3::mismatch_score{-6}};
auto aConfig = seqan3::align_cfg::method_local{} | seqan3::align_cfg::scoring_scheme{scoreMat} | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-5}, seqan3::align_cfg::extension_score{-2}};

// identified strand from the string seq, also extract umi sequence
// return anchor mapping info, strand, umi, seqence removed umi, quality removed umi
// +: cDNA strand, -: RNA strand
// output RNA strand (-) as circFL_strand.fastq
std::optional<std::tuple<std::string, char, std::string, std::string, std::string>> identifyStrandWithUMI(
	std::string &id, const seqan3::dna5_vector &seq, std::string &&qual,
	seqan3::dna4_vector &anchorx, // seqan3::dna4_vector &anchory,
	// int anchorxThreshold, int anchoryThreshold, int endSize,
	DegenerateBaseVector &umiSeq)
{
	// discard too short reads less than 200 bp
	if (seq.size() < 200)
		return std::nullopt;
	using namespace seqan3;

	// align anchorx to the seq
	auto alignForwardStrandAnchorX = align_pairwise(std::tie(seq, anchorx), aConfig);
	auto alignReverseStrandAnchorX = align_pairwise(std::tie(seq, anchorxRev), aConfig);

	size_t beginForwardStrandAnchorX, endForwardStrandAnchorX, beginReverseStrandAnchorX, endReverseStrandAnchorX;
	int scoreForwardStrandAnchorX, scoreReverseStrandAnchorX;
	auto &res = *alignForwardStrandAnchorX.begin();
	beginForwardStrandAnchorX = res.sequence1_begin_position();
	endForwardStrandAnchorX = res.sequence1_end_position();
	scoreForwardStrandAnchorX = res.score();
	res = *alignReverseStrandAnchorX.begin();
	beginReverseStrandAnchorX = res.sequence1_begin_position();
	endReverseStrandAnchorX = res.sequence1_end_position();
	scoreReverseStrandAnchorX = res.score();

	// align anchory to the seq
	/*auto alignForwardStrandAnchorY = align_pairwise(std::tie(seq, anchory), aConfig);
	auto alignReverseStrandAnchorY = align_pairwise(std::tie(seq, anchoryRev), aConfig);

	size_t beginForwardStrandAnchorY, endForwardStrandAnchorY, beginReverseStrandAnchorY, endReverseStrandAnchorY;
	int scoreForwardStrandAnchorY, scoreReverseStrandAnchorY;
	res = *alignForwardStrandAnchorY.begin();
	beginForwardStrandAnchorY = res.sequence1_begin_position();
	endForwardStrandAnchorY = res.sequence1_end_position();
	scoreForwardStrandAnchorY = res.score();
	res = *alignReverseStrandAnchorY.begin();
	beginReverseStrandAnchorY = res.sequence1_begin_position();
	endReverseStrandAnchorY = res.sequence1_end_position();
	scoreReverseStrandAnchorY = res.score();*/

	size_t puStartPos, puEndPos;
	std::string strandOut = id + "\t+\t" + std::to_string(beginForwardStrandAnchorX) + "\t" + std::to_string(endForwardStrandAnchorX) + "\t" + std::to_string(scoreForwardStrandAnchorX) + "\t-\t" +
							std::to_string(beginReverseStrandAnchorX) + "\t" + std::to_string(endReverseStrandAnchorX) + "\t" + std::to_string(scoreReverseStrandAnchorX) + "\t" +
							(scoreForwardStrandAnchorX > scoreReverseStrandAnchorX ? "+" : "-");
	if (scoreForwardStrandAnchorX > scoreReverseStrandAnchorX)
	{
		// writeRecord(strandFqFileOut, shortID, seq, qual);
		puStartPos = endForwardStrandAnchorX < 5 ? 0 : endForwardStrandAnchorX - 5;
		puEndPos = std::min(endForwardStrandAnchorX + umiSeq.size() + 5, seq.size());
		if (puStartPos > puEndPos)
			return std::nullopt;
		auto umiret{extract(
			convertBaseArray(circfull::seqan3ToString(seq | seqan3::views::slice(puStartPos, puEndPos))),
			umiSeq)};
		if (umiret.has_value())
		{
			auto [umi, uStart, uEnd, uScore] = umiret.value();
			std::string umiDepleteSequence = circfull::seqan3ToString(seq | seqan3::views::slice(puStartPos + uEnd, seq.size()) | std::views::reverse | seqan3::views::complement);
			std::string qualDepleteSequence = circfull::seqan3ToString(qual | seqan3::views::slice(puStartPos + uEnd, seq.size()) | std::views::reverse);
			strandOut = strandOut + "\tumiscore=" + std::to_string(uScore) + "_" + std::to_string(uStart) + "_" + std::to_string(uEnd) + "_" + std::to_string(uEnd - uStart) + "_" + std::to_string(uEnd - uStart - umi.length()) + "\n";
			return std::make_tuple(strandOut, '+', umi, umiDepleteSequence, qualDepleteSequence);
		}
		else
		{
			return std::nullopt;
		}
	}
	else
	{
		std::string seqRev = circfull::seqan3ToString(seq | seqan3::views::complement | std::views::reverse);
		std::string qualRevString = circfull::seqan3ToString(qual | std::views::reverse);
		puStartPos = seq.size() - beginReverseStrandAnchorX < 5 ? 0 : seq.size() - beginReverseStrandAnchorX - 5;
		puEndPos = std::min(seq.size() - beginReverseStrandAnchorX + umiSeq.size() + 5, seq.size());
		if (puStartPos > puEndPos)
			return std::nullopt;
		auto umiret{extract(
			convertBaseArray(circfull::seqan3ToString(seqRev | seqan3::views::slice(puStartPos, puEndPos))),
			umiSeq)};
		if (umiret.has_value())
		{
			auto [umi, uStart, uEnd, uScore] = umiret.value();
			std::string umiDepleteSequence = circfull::seqan3ToString(seqRev | seqan3::views::char_to<seqan3::dna5> | seqan3::views::slice(puStartPos + uEnd, seqRev.length()) | std::views::reverse | seqan3::views::complement);
			std::string qualDepleteSequence = circfull::seqan3ToString(qualRevString | seqan3::views::slice(puStartPos + uEnd, seqRev.length()) | std::views::reverse);
			strandOut = strandOut + "\tumiscore=" + std::to_string(uScore) + "_" + std::to_string(uStart) + "_" + std::to_string(uEnd) + "_" + std::to_string(uEnd - uStart) + "_" + std::to_string(uEnd - uStart - umi.length()) + "\n";
			return std::make_tuple(strandOut, '-', umi, umiDepleteSequence, qualDepleteSequence);
		}
		else
		{
			return std::nullopt;
		}
	}
}

int checkFiles(std::filesystem::path fastq, std::filesystem::path oPath, std::string seqkit)
{
	if (!circfull::checkFile(fastq))
	{
		std::cout << "[ERROR] fastq file not existed." << std::endl;
		return -1;
	}
	if (std::filesystem::exists(oPath))
	{
		if (!oPath.empty())
		{
			std::cout << "[WARNING] output dir existed." << std::endl;
		}
	}
	else
	{
		circfull::makeDir(oPath);
	}
	return 0;
}

namespace circfull
{

	int circfull_umi_extract(circfull::CircfullOption &opt)
	{
		umiLength = opt.umiLength;
		if (checkFiles(opt.rawFastq, opt.oPath, opt.seqkit) != 0)
			return -1;
		boost::asio::thread_pool pool{opt.nthread};
		std::filesystem::path tmpFile, tmpOut, strandFastq, umiFasta, umiFasta2, preumiFasta, preumiFasta2;
		std::filesystem::path tmpDir;
		printTimeInfo("identify strand:");
		// read input sequence from opt.rawFastq
		seqan3::sequence_file_input seqFileIn{opt.rawFastq};
		std::ofstream outStream{opt.oPath / (opt.oPrefix + "_barcode.tsv")};
		std::ofstream umiFileOut{opt.umiFasta};
		std::ofstream strandFileOut{opt.strandNotTrimFastq};
		std::ranges::copy(
			opt.anchorx | seqan3::views::char_to<seqan3::dna4>,
			std::back_inserter(anchorx));
		std::ranges::copy(
			anchorx | seqan3::views::complement | std::views::reverse,
			std::back_inserter(anchorxRev));
		seqan3::debug_stream << "[INFO] anchorx: " << anchorx << std::endl;
		seqan3::debug_stream << "[INFO] anchorxRev: " << anchorxRev << std::endl;
		std::cout << "[INFO] umiLength: " << umiLength << std::endl;
		// identify strand and extract umi sequence using boost thread pool
		for (auto &record : seqFileIn)
		{
			boost::asio::post(pool, [&, record]()
							  {
				// get short id as splitted with space
				std::string shortID = getId(record.id());
				auto ret {identifyStrandWithUMI(
					shortID,
					//to_std::string(record.sequence()),
					record.sequence(),
					circfull::seqan3ToString(record.base_qualities()),
					anchorx,
					opt.umiSeq
				)};
				if (ret.has_value()) {
					mt.lock();
					outStream << std::get<0>(ret.value());
					umiFileOut << outputAsFasta(shortID, std::get<2>(ret.value()));
					strandFileOut << "@" + shortID + "\n" + std::get<3>(ret.value()) + "\n+\n" + std::get<4>(ret.value()) + "\n";
					keepCountSum++;
					mt.unlock();
				} else {
					mt.lock();
					discardCountSum++;
					mt.unlock();
				} });
		}
		pool.join();
		std::cout << "discard " << discardCountSum << " reads in " << keepCountSum + discardCountSum << " reads" << std::endl;

		// merge result files
		printTimeInfo("trimming adapters from stranded fastq:");
		if (trimAdaptor(opt.porechop, opt.strandNotTrimFastq, opt.strandFastq, opt.nthread))
			return -1;
		printTimeInfo("finished.");
		return 0;
	}

}