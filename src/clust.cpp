#include <fstream>
#include <string>
#include <iostream>
#include <queue>
#include <filesystem>
#include <thread>
#include "circfullUtils.hpp"
#include "clust.hpp"
#include "cdhit/cdhit-common.h"
#include "spoa/spoa.hpp"

// over-write some defs in cd-hi.h for est version
#undef MAX_UAA
#define MAX_UAA 4

// over-write some defs in cd-hi-init.h for est version

void setaa_to_na();
void make_comp_short_word_index(int NAA, int *NAAN_array, Vector<int> &Comp_AAN_idx);
void make_comp_iseq(int len, char *iseq_comp, char *iseq);

/// @brief print clustering information to output file
/// @param clustId cluster id
/// @param seq sequence
/// @param id seq id
/// @param fout output file
/// @param options cd-hit options
/// @param buf output buffer
void clustPrintInfo(int clustId, Sequence *seq, int id, FILE *fout, const Options &options, char *buf)
{
	fprintf(fout, "Cluster_%i\t", clustId);
	bool print = options.print != 0;
	bool strand = options.isEST;
	fprintf(fout, "%i\t%int\t%s", id, seq->size, seq->identifier + 1);
	if (seq->identity)
	{
		int *c = seq->coverage;
		fprintf(fout, "\tat ");
		if (print)
			fprintf(fout, "%i:%i:%i:%i/", c[0], c[1], c[2], c[3]);
		if (strand)
			fprintf(fout, "%c/", (seq->state & IS_MINUS_STRAND) ? '-' : '+');
		fprintf(fout, "%.2f%%", seq->identity * 100);
		if (options.useDistance)
			fprintf(fout, "/%.2f%%", seq->distance * 100);
		fprintf(fout, "\n");
	}
	else
	{
		fprintf(fout, "\t*\n");
	}
}

/// @brief write clustering information to output file
/// @param seqDb cd-hit sequence database
/// @param options cd-hit options
/// @param output path to output file
/// @param clustInfo map of UMI cluster information
/// @return number of clusters
int WriteClustInfo(const SequenceDB &seqDb, const Options &options, std::filesystem::path &output, std::map<std::string, int> &clustInfo)
{
	int i, i0, k, N = seqDb.sequences.size();
	std::vector<long long> sorting(N);
	for (i = 0; i < N; i++)
		sorting[i] = (static_cast<long long>(seqDb.sequences[i]->index) << 32) | i;
	sort(sorting.begin(), sorting.end());

	FILE *fout;
	char *buf = new char[MAX_DES + 1];

	std::cout << "writing clustering information" << std::endl;
	int M = seqDb.rep_seqs.size();
	std::vector<Vector<int>> clusters(M);
	for (i = 0; i < N; i++)
	{
		int k = sorting[i] & 0xffffffff;
		int id = seqDb.sequences[k]->cluster_id;
		clusters[id].Append(k);
	}

	fout = fopen(output.c_str(), "w+");
	for (i = 0; i < M; i++)
	{
		for (k = 0; k < static_cast<int>(clusters[i].size()); k++)
		{
			// seqDb.sequences[clusters[i][k]] -> PrintInfo(k, fout, options, buf);
			clustPrintInfo(i, seqDb.sequences[clusters[i][k]], k, fout, options, buf);
			// clustInfo.push_back(make_tuple(i, string(seqDb.sequences[clusters[i][k]] -> identifier + 1)));
			clustInfo[std::string(seqDb.sequences[clusters[i][k]]->identifier + 1)] = i;
		}
	}
	fclose(fout);
	delete[] buf;
	return M;
}

int clustUMISeq(const std::filesystem::path &fastqFile, const std::filesystem::path &oPrefix, int wordSize, int nTread, double identity, std::map<std::string, int> &clustInfo)
{
	Options options;
	SequenceDB seqDb;
	std::filesystem::path outputFile = std::filesystem::path(oPrefix) += ".clstr";
	// set cd-hit options
	options.cluster_thd = identity;
	options.NAA = 10;
	options.NAAN = NAA8;
	seqDb.NAAN = NAA8;
	options.NAA_top_limit = 12;
	setaa_to_na();
	mat.set_to_na();
	options.SetOption("-i", fastqFile.c_str());
	options.SetOption("-c", std::to_string(identity).c_str());
	options.SetOption("-n", std::to_string(wordSize).c_str());
	options.SetOption("-T", std::to_string(nTread).c_str());
	options.SetOption("-bak", "");
	options.SetOption("-d", "100");
	options.SetOption("-sf", "0");
	options.SetOption("-M", "0");
	options.SetOption("-g", "1");
	options.SetOption("-o", (std::filesystem::path(oPrefix) += ".clstr").c_str());
	options.Validate();
	InitNAA(MAX_UAA);
	options.NAAN = NAAN_array[options.NAA];
	seqDb.NAAN = NAAN_array[options.NAA];
	// read fastq
	seqDb.Read(fastqFile.c_str(), options);
	std::cout << "total seq: " << seqDb.sequences.size() << std::endl;
	seqDb.SortDivide(options);
	seqDb.DoClustering(options);
	return WriteClustInfo(seqDb, options, outputFile, clustInfo);
}

/// @brief generate consensus sequence from a set of sequences
/// @param seqs vector of sequences
/// @param quals vector of quality scores
/// @return consensus sequence
std::string generateConsensusSeq(const std::vector<std::string> &seqs, const std::vector<std::string> &quals)
{
	auto alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);
	spoa::Graph graph{};
	for (size_t i = 0; i < seqs.size(); i++)
	{
		const auto &it = seqs[i];
		const auto &qu = quals[i];
		auto alignment = alignmentEngine->Align(it, graph);
		graph.AddAlignment(alignment, it, qu);
	}
	std::string consensus = graph.GenerateConsensus();
	return consensus;
}

/// @brief generate consensus sequence from a set of sequences
/// @param seqs vector of sequences
/// @return consensus sequence
std::string generateConsensusSeq(const std::vector<std::string> &seqs)
{
	if (seqs.size() == 0)
		return "";
	if (seqs.size() == 1)
		return seqs[0];
	// auto alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3); //Alignment Score for UMI guided consensus generation
	auto alignmentEngine = spoa::AlignmentEngine::Create(static_cast<spoa::AlignmentType>(2), 10, -4, -8, -2, -24, -1); // Alignment Score for CCS consensus generation
	spoa::Graph graph{};
	for (size_t i = 0; i < seqs.size(); i++)
	{
		const auto &it = seqs[i];
		auto alignment = alignmentEngine->Align(it, graph);
		graph.AddAlignment(alignment, it);
	}
	std::string consensus = graph.GenerateConsensus();
	return consensus;
}
