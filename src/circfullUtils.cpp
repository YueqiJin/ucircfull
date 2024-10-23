#include <fstream>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <filesystem>
#include <thread>
#include <mutex>
#include <seqan3/alphabet/all.hpp>
#include "argparse/argparse.hpp"
#include "circfullUtils.hpp"
#include "ucircfull.hpp"
#include "version.hpp"
#include "degenerateBase.hpp"

using namespace std;
namespace fs = filesystem;

namespace circfull
{

	bool debug = false;

	std::istream &operator>>(std::istream &is, const literal &lit)
	{
		for (const char *p = lit.str; *p; ++p)
		{
			if (std::isspace(*p))
				std::istream::sentry s(is); // discard whitespace
			else if (is.get() != *p)
				is.setstate(std::ios::failbit); // stop extracting
		}
		return is;
	}

	std::vector<int> getIntFromCommaSepStr(const std::string &str)
	{
		std::vector<int> ret;
		std::stringstream ss(str);
		std::string temp;
		while (getline(ss, temp, ','))
			ret.push_back(stoi(temp));
		return ret;
	}

	bool checkFile(const fs::path &file)
	{
		return fs::exists(file);
	}

	int getFileNrow(fs::path filename, size_t nthread)
	{
		if (!fs::exists(filename))
			return -1;
		filename = fs::canonical(filename);
		ifstream fin;
		fin.open(filename, ios::in);
		if (nthread == 1)
		{
			int cnt = 0;
			string temp;
			while (fin.peek() != EOF)
				getline(fin, temp), cnt++;
			return cnt;
		}
		// multi-thread version of counting '\n'
		fin.seekg(0, std::ios::end);
		unsigned long long int fileSize = fin.tellg();
		fin.seekg(0, std::ios::beg);
		unsigned long long int chunkSize = fileSize / nthread;
		thread *threads = new thread[nthread];
		int totalLine = 0;
		// Launch threads to count '\n' in parallel
		mutex mt;
		for (int i = 0; i < nthread; i++)
		{
			unsigned long long int start = i * chunkSize;
			unsigned long long int end = (i == nthread - 1) ? fileSize : (i + 1) * chunkSize;
			threads[i] = thread([filename, start, end, &totalLine, &mt]
								{
				std::ifstream file(filename, ios::in);
				std::string line;
				int count = 0;
				// Seek to the start of the portion
				file.seekg(start);
				while (file.tellg() < end && std::getline(file, line)) {
					count++;
				}
				mt.lock();
				totalLine += count;
				mt.unlock(); });
		}
		for (int i = 0; i < nthread; i++)
		{
			threads[i].join();
		}
		return totalLine;
	}

	int makeDir(const fs::path &targetDir)
	{
		if (fs::exists(targetDir))
			delDir(targetDir);
		return !fs::create_directory(targetDir);
	}

	int delDir(fs::path targetDir)
	{
		targetDir = fs::canonical(targetDir);
		return !fs::remove_all(targetDir);
	}

	int checkSeqkit(const string &seqkit)
	{
		string cmdLine = seqkit + " --help >/dev/null 2>&1";
		return system(cmdLine.c_str());
	}

	int splitFastq(const string &seqkit, const fs::path &fastqFile, const fs::path &outputDir, const string &outputPrefix, int nSplit)
	{
		string cmdLine = seqkit + " split2 " + fastqFile.string() + " -j " + to_string(nSplit) + " -p " + to_string(nSplit) + " -O " + outputDir.string() + " -o " + outputPrefix + " --force --quiet";
		return system(cmdLine.c_str());
	}

	template <class T>
	int getArrayLength(T &array)
	{
		return sizeof(array) / sizeof(array[0]);
	}

	int mergeFile(fs::path *files, const fs::path &outputFile, int nFiles)
	{
		ifstream fin;
		ofstream fout;
		fout.open(outputFile, ios::out);
		for (int i = 0; i < nFiles; i++)
		{
			fin.open(files[i], ios::in);
			fout << fin.rdbuf();
			fin.close();
		}
		fout.close();
		return 0;
	}

	int delFile(fs::path *files, int nFiles)
	{
		for (int i = 0; i < nFiles; i++)
		{
			if (fs::remove(files[i]) == false)
				return -1;
		}
		return 0;
	}

	void printTimeInfo(const string &text)
	{
		time_t nowtime;
		time(&nowtime);
		cout << "[" << std::put_time(std::localtime(&nowtime), "%Y-%m-%d %H:%M:%S") << "] " << text << endl;
		cout.flush();
	}

	string getFirstString(const string &text, const char &delim)
	{
		stringstream ss(text);
		string ret;
		getline(ss, ret, delim);
		return ret;
	}

	string outputAsFasta(const string &id, const string &seq)
	{
		// output id and seq as fasta format with 80 chars per line
		string ret = ">" + id + "\n";
		int n = seq.length();
		for (int i = 0; i < n; i += 80)
		{
			ret += seq.substr(i, 80) + "\n";
		}
		return ret;
	}

	int sortFastqById(const string &seqkit, const fs::path &inputFastq, const fs::path &outputFastq, size_t nthread)
	{
		string cmdLine = seqkit + " sort -n " + inputFastq.string() + " -o " + outputFastq.string() + " -j " + to_string(nthread);
		return system(cmdLine.c_str());
	}

	int convertFastq2Fasta(const string &seqkit, const fs::path &fastqFile, const fs::path &fastaFile, size_t nthread)
	{
		string cmdLine = seqkit + " fq2fa " + fastqFile.string() + " -o " + fastaFile.string() + " -j " + to_string(nthread);
		return system(cmdLine.c_str());
	}

	int trimAdaptor(const string &porechop, const fs::path &fastqFile, const fs::path &outFile, size_t nthread)
	{
		string cmdLine = porechop + " -i " + fastqFile.string() + " -o " + outFile.string() + " --barcode_threshold 95 --check_reads 1000 -v 1 --no_split -t " + to_string(nthread);
		return system(cmdLine.c_str());
	}

	std::string getId(const std::string &idString)
	{
		std::istringstream iss(idString);
		std::string ret;
		iss >> ret;
		return ret;
	}

	void applyExtractOption(CircfullOption &opt, const argparse::ArgumentParser &args)
	{
		opt.call = CircfullOption::mod::UMI_EXTRACT;
		// input info
		opt.rawFastq = args.get<string>("-i");
		opt.anchorx = args.get<string>("-x");
		opt.umiSeq = convertBaseArray(args.get<string>("-u"));
		opt.noUMI = args.get<bool>("-n");
		opt.umiLength = opt.umiSeq.size();
		// output info
		opt.oPath = args.get<string>("-o");
		opt.oPrefix = args.get<string>("-p");
		opt.strandNotTrimFastq = opt.oPath / (opt.oPrefix + "_strand_ntrim.fastq");
		opt.strandFastq = opt.oPath / (opt.oPrefix + "_strand.fastq");
		opt.umiFasta = opt.oPath / (opt.oPrefix + "_umi.fasta");
		// running options
		opt.nthread = args.get<int>("-t");
		opt.porechop = args.get<string>("--porechop");
	}

	void applyClustOption(CircfullOption &opt, const argparse::ArgumentParser &args)
	{
		opt.call = CircfullOption::mod::UMI_CLUST;
		// input info
		opt.strandFastq = args.get<string>("-i");
		opt.umiFasta = args.get<string>("-u");
		// output info
		opt.oPath = args.get<string>("-o");
		opt.oPrefix = args.get<string>("-p");
		opt.umiClustTab = opt.oPath / (opt.oPrefix + "umi_clust_info.clstr");
		opt.consensusFasta = opt.oPath / (opt.oPrefix + "_cons.fastq");
		// running options
		opt.nthread = args.get<int>("-t");
		opt.seqkit = args.get<string>("--seqkit");
		opt.clustSeq = args.get<bool>("-s");
	}

	void applyCircCallOption(CircfullOption &opt, const argparse::ArgumentParser &args)
	{
		opt.call = CircfullOption::mod::CIRC_CALL;
		// input info
		opt.umiClustTab = args.get<string>("-u");
		// opt.consensusFasta = args.get<string>("-c");
		opt.referenceFasta = args.get<string>("-r");
		opt.annotationGTF = args.get<string>("-a");
		// output info
		opt.oPath = args.get<string>("-o");
		opt.oPrefix = args.get<string>("-p");
		opt.circCallOutputPath = opt.oPath / (opt.oPrefix + "_ciri_call");
		opt.circCollapsePath = opt.oPath / (opt.oPrefix + "_ciri_collapse");
		opt.spliceSignal = args.get<string>("--splice");
		// running options
		opt.nthread = args.get<int>("-t");
		// opt.ciriLong = args.get<string>("--ciri_long");
		//  software path
		opt.minimap2 = args.get<string>("--minimap2");
		opt.samtools = args.get<string>("--samtools");
		opt.strandFastq = args.get<string>("-i");
		debug = args.get<bool>("--debug");
		// circ call options
		// opt.strand = !(args.get<bool>("-sn"));
		opt.strand = true; // default to find splice motif in forward strand
		if (opt.useBam = args.is_used("--bam"))
			opt.bamFile = args.get<string>("--bam");
		std::string circCallMod = args.get<string>("-m");
		if (circCallMod == "RG")
			opt.modCircCall = CircfullOption::CircCallMod::RG;
		else if (circCallMod == "cRG")
			opt.modCircCall = CircfullOption::CircCallMod::cRG;
		else if (circCallMod == "ucRG")
			opt.modCircCall = CircfullOption::CircCallMod::ucRG;
		else if (circCallMod == "CIRI")
			opt.modCircCall = CircfullOption::CircCallMod::CIRI;
		// RG mod options
		if (opt.modCircCall == CircfullOption::CircCallMod::RG || opt.modCircCall == CircfullOption::CircCallMod::cRG)
		{
			opt.minimap2 = args.get<string>("--minimap2");
			opt.samtools = args.get<string>("--samtools");
		}
		// CIRI mod options
		if (opt.modCircCall == CircfullOption::CircCallMod::CIRI)
		{
			opt.ciriLong = args.get<string>("--ciri_long");
		}
	}
	/*
		void applyCircFiltOption(CircfullOption &opt, const argparse::ArgumentParser &args)
		{
			opt.call = CircfullOption::mod::CIRC_FILT;
			// input info
			opt.rawCircRes = args.get<string>("-i");
			opt.umiClustTab = args.get<string>("-u");
			opt.consCircRes = args.get<string>("-c");
			// output info
			opt.oPath = args.get<string>("-o");
			opt.oPrefix = args.get<string>("-p");
			opt.circCallOutputPath = opt.oPath / (opt.oPrefix + "_ciri_call");
			opt.circCollapsePath = opt.oPath / (opt.oPrefix + "_ciri_collapse");
			opt.circCandFile = opt.oPath / (opt.oPrefix + "_circ_cand.tsv");
			// filter options
			opt.circNumThres = args.get<int>("-n");
			// running options
			opt.nthread = args.get<int>("-t");
		}
	*/

	std::vector<int> split(const std::string &str, const char &delim)
	{
		std::vector<int> ret;
		std::istringstream iss(str);
		std::string token;
		while (std::getline(iss, token, delim))
		{
			ret.push_back(std::stoi(token));
		}
		return ret;
	}

	void parseArgs(CircfullOption &opt, int argc, char **argv)
	{
		argparse::ArgumentParser program{"ucircfull", PROJECT_VERSION, argparse::default_arguments::all};

		argparse::ArgumentParser extract_umi_command{"extract_umi"};
		extract_umi_command.add_description("Extract UMI sequence and identify strand from ucircFL-seq raw fastq");
		extract_umi_command.add_argument("-i", "--input")
			.metavar("FQ")
			.required()
			.help("ucircFL-seq raw fastq file.");
		extract_umi_command.add_argument("-x", "--anchorx")
			.metavar("SEQ")
			.required()
			.help("anchor sequence used in 1st strand cDNA synthesis.");
		extract_umi_command.add_argument("-u", "--umi")
			.metavar("SEQ")
			.required()
			.default_value(string("CTCNNNYRNNNYRNNNYRNNNGAG"))
			.help("umi pattern.");
		extract_umi_command.add_argument("-n", "--noumi")
			.default_value(false)
			.help("no UMIs were added to 1st strand anchor.")
			.implicit_value(true);
		extract_umi_command.add_argument("-o", "--outdir")
			.metavar("DIR")
			.required()
			.default_value(string("."))
			.help("output directory.");
		extract_umi_command.add_argument("-p", "--prefix")
			.metavar("PREFIX")
			.required()
			.default_value(string("circFL"))
			.help("output prefix.");
		extract_umi_command.add_argument("-t", "--thread")
			.metavar("INT")
			.required()
			.default_value(4)
			.help("number of threads used.")
			.scan<'i', int>();
		extract_umi_command.add_argument("--seqkit")
			.metavar("PATH")
			.default_value(string("seqkit"))
			.help("path to seqkit.");
		extract_umi_command.add_argument("--porechop")
			.metavar("PATH")
			.default_value(string("porechop"))
			.help("path to porechop.");

		argparse::ArgumentParser clust_umi_command{"clust_umi"};
		clust_umi_command.add_description("UMI clusting guided consensus generation.");
		clust_umi_command.add_argument("-i", "--input")
			.metavar("FQ")
			.required()
			.help("stranded fastq file.");
		clust_umi_command.add_argument("-u", "--umi")
			.metavar("FA")
			.required()
			.help("umi fasta file.");
		clust_umi_command.add_argument("-s", "--seq")
			.default_value(false)
			.help("generate sequence cluster.");
		clust_umi_command.add_argument("-o", "--outdir")
			.metavar("DIR")
			.required()
			.default_value(string("."))
			.help("output directory.");
		clust_umi_command.add_argument("-p", "--prefix")
			.metavar("PREFIX")
			.required()
			.default_value(string("circFL"))
			.help("output prefix.");
		clust_umi_command.add_argument("-t", "--thread")
			.metavar("INT")
			.required()
			.default_value(4)
			.help("number of threads used.")
			.scan<'i', int>();
		clust_umi_command.add_argument("--seqkit")
			.metavar("PATH")
			.default_value(string("seqkit"))
			.help("path to seqkit.");

		argparse::ArgumentParser circ_call_command{"circ_call"};
		circ_call_command.add_description("circRNA identification.");
		circ_call_command.add_argument("-m", "--mode")
			.metavar("STR")
			.required()
			//.help("circRNA calling mode. (RG, ucRG, CIRI)");
			.help("circRNA calling mode. (RG, ucRG)");
		/*
			circ_call_command.add_argument("-sn", "--notstranded")
				.default_value(false)
				.help("stranded mode.");
		*/
		circ_call_command.add_argument("-i", "--input")
			.metavar("FQ")
			.required()
			.help("stranded fastq file.");
		circ_call_command.add_argument("--bam")
			.metavar("BAM")
			.help("mapped bam file as input.");
		circ_call_command.add_argument("-r", "--ref")
			.metavar("REF")
			.required()
			.help("CIRI-long reference directory.");
		circ_call_command.add_argument("-a", "--anno")
			.metavar("GTF")
			.required()
			.help("reference annotation GTF file.");
		// circ_call_command.add_argument("-c", "--consensus")
		//	.metavar("FA")
		//	.default_value("-")
		//	.help("consensus fasta file.");
		circ_call_command.add_argument("-u", "--umi")
			.metavar("CLSTR")
			.default_value("-")
			.help("umi clust results file.");
		circ_call_command.add_argument("--splice")
			.metavar("MOTIF")
			.required()
			.default_value(string("AGGT,AGGC,ACAT,ACGT,AGAT"))
			.help("output directory.");
		circ_call_command.add_argument("-o", "--outdir")
			.metavar("DIR")
			.required()
			.default_value(string("."))
			.help("output directory.");
		circ_call_command.add_argument("-p", "--prefix")
			.metavar("PREFIX")
			.required()
			.default_value(string("circFL"))
			.help("output prefix.");
		circ_call_command.add_argument("-t", "--thread")
			.metavar("INT")
			.required()
			.default_value(4)
			.help("number of threads used.")
			.scan<'i', int>();
		/*
			circ_call_command.add_argument("--ciri_long")
				.metavar("PATH")
				.default_value(string("CIRI-long"))
				.help("path to CIRI-long.");
		*/
		circ_call_command.add_argument("--minimap2")
			.metavar("PATH")
			.default_value(string("minimap2"))
			.help("path to minimap2.");
		circ_call_command.add_argument("--samtools")
			.metavar("PATH")
			.default_value(string("samtools"))
			.help("path to samtools.");
		circ_call_command.add_argument("--debug")
			.default_value(false)
			.help("enable debug output.")
			.implicit_value(true);

		/*
			argparse::ArgumentParser circ_filt_command{"circ_filt"};
			circ_filt_command.add_description("filter CIRI-long high confidence circRNA candidate by using UMI cluster.");
			circ_filt_command.add_argument("-i", "--input")
				.metavar("FA")
				.required()
				.help("CIRI-long circ_call result for raw sequence.");
			circ_filt_command.add_argument("-u", "--umi")
				.metavar("CLSTR")
				.required()
				.help("UMI cluster result.");
			circ_filt_command.add_argument("-a", "--anno")
				.metavar("GTF")
				.required()
				.help("reference annotation GTF file.");
			circ_filt_command.add_argument("-c", "--consensus")
				.metavar("FA")
				.default_value("-")
				.help("CIRI-long circ_call result for consensus sequence file.");
			circ_filt_command.add_argument("-o", "--outdir")
				.metavar("DIR")
				.required()
				.default_value(string("."))
				.help("output directory.");
			circ_filt_command.add_argument("-p", "--prefix")
				.metavar("PREFIX")
				.required()
				.default_value(string("circFL"))
				.help("output prefix.");
			circ_filt_command.add_argument("-n", "--circ_num_thresh")
				.metavar("INT")
				.required()
				.default_value(3)
				.help("circRNA number threshold.");
			circ_filt_command.add_argument("-t", "--thread")
				.metavar("INT")
				.required()
				.default_value(4)
				.help("number of threads used.")
				.scan<'i', int>();
		*/
		program.add_subparser(extract_umi_command);
		program.add_subparser(clust_umi_command);
		program.add_subparser(circ_call_command);
		// program.add_subparser(circ_filt_command);
		try
		{
			program.parse_args(argc, argv);
			if (program.is_subcommand_used("extract_umi"))
				applyExtractOption(opt, program.at<argparse::ArgumentParser>("extract_umi"));
			if (program.is_subcommand_used("clust_umi"))
				applyClustOption(opt, program.at<argparse::ArgumentParser>("clust_umi"));
			if (program.is_subcommand_used("circ_call"))
				applyCircCallOption(opt, program.at<argparse::ArgumentParser>("circ_call"));
			/*
				if (program.is_subcommand_used("circ_filt"))
					applyCircFiltOption(opt, program.at<argparse::ArgumentParser>("circ_filt"));
			*/
		}
		catch (const std::runtime_error &err)
		{
			std::cerr << err.what() << std::endl;
			std::cerr << program;
			std::exit(1);
		}
	}
}
