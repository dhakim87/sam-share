#include "main.h"

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <unordered_set>

//Daniel Hakim
//Goals:
// 1) See how fast woltka could be if it was in C++
// 2) Create a confusion matrix tracking reads that map to multiple locations
// 3) Examine reads that map to >1 genome
// 4) Examine reads that are <X% coverage in a sample

// Global scratch space
static std::map<std::string, int> private_reads;
static std::map<std::string, float> split_reads;  //This could easily use fixed point, eh.
static std::map<std::string, std::map<std::string, float> > shared2_reads;

static std::map<std::string, std::map<std::string, std::vector<IndexRange> > > shared2_row_coverage;
static std::map<std::string, std::map<std::string, std::vector<IndexRange> > > shared_powerset_row_coverage;

static std::vector<std::string> akkermansia_muciniphila_list;
static std::unordered_set<std::string> akkermansia_muciniphila;

static std::vector<std::string> yersinia_list;
static std::unordered_set<std::string> yersinia;

static std::vector<std::string> escherichia_list;
static std::unordered_set<std::string> escherichia;

static std::vector<std::string> klebsiella_list;
static std::unordered_set<std::string> klebsiella;

//Bacillus lists
static std::vector<std::string> b_anthracis_list;
static std::unordered_set<std::string> b_anthracis;
static std::vector<std::string> b_cereus_list;
static std::unordered_set<std::string> b_cereus;
static std::vector<std::string> b_thuringiensis_list;
static std::unordered_set<std::string> b_thuringiensis;

static long total_yersinia_reads = 0;
static long yersinia_escherichia_reads = 0;
static long yersinia_klebsiella_reads = 0;
static long yersinia_only_reads = 0;
static long yersinia_escherichia_and_klebsiella_reads = 0;

static long total_anthracis_reads = 0;
static long anthracis_cereus_reads = 0;
static long anthracis_thuringiensis_reads = 0;
static long anthracis_only_reads = 0;
static long anthracis_cereus_and_thuringiensis_reads = 0;

static void init_akkermansia_lists()
{
    akkermansia_muciniphila_list.push_back("G001683795");
    akkermansia_muciniphila_list.push_back("G900097105");
    akkermansia_muciniphila_list.push_back("G000020225");
    akkermansia_muciniphila_list.push_back("G000723745");
    akkermansia_muciniphila_list.push_back("G000436395");
    akkermansia_muciniphila_list.push_back("G001917295");
    akkermansia_muciniphila_list.push_back("G000437075");
    akkermansia_muciniphila_list.push_back("G001647615");
    akkermansia_muciniphila_list.push_back("G001578645");
    akkermansia_muciniphila_list.push_back("G001580195");
    akkermansia_muciniphila_list.push_back("G001940945");
    akkermansia_muciniphila_list.push_back("G000980515");
    for (int i = 0; i < akkermansia_muciniphila_list.size(); i++)
        akkermansia_muciniphila.insert(akkermansia_muciniphila_list[i]);
}

static void init_yersinia_lists()
{
    yersinia_list.push_back("G000009065");
    yersinia_list.push_back("G000009345");
    yersinia_list.push_back("G000168835");
    yersinia_list.push_back("G000169655");
    yersinia_list.push_back("G000253175");
    yersinia_list.push_back("G000320465");
    yersinia_list.push_back("G000323285");
    yersinia_list.push_back("G000323565");
    yersinia_list.push_back("G000323585");
    yersinia_list.push_back("G000323645");
    yersinia_list.push_back("G000323725");
    yersinia_list.push_back("G000323785");
    yersinia_list.push_back("G000323825");
    yersinia_list.push_back("G000324005");
    yersinia_list.push_back("G000324145");
    yersinia_list.push_back("G000324265");
    yersinia_list.push_back("G000324305");
    yersinia_list.push_back("G000324445");
    yersinia_list.push_back("G000324485");
    yersinia_list.push_back("G000324505");
    yersinia_list.push_back("G000324625");
    yersinia_list.push_back("G000324685");
    yersinia_list.push_back("G000324725");
    yersinia_list.push_back("G000324905");
    yersinia_list.push_back("G000324985");
    yersinia_list.push_back("G000325085");
    yersinia_list.push_back("G000325145");
    yersinia_list.push_back("G000325165");
    yersinia_list.push_back("G000325245");
    yersinia_list.push_back("G000325265");
    yersinia_list.push_back("G000325285");
    yersinia_list.push_back("G000964565");
    yersinia_list.push_back("G001122605");
    for (int i = 0; i < yersinia_list.size(); i++)
        yersinia.insert(yersinia_list[i]);

}

static void init_ecoli_lists()
{
    escherichia_list.push_back("G000008865");
    escherichia_list.push_back("G000026325");
    escherichia_list.push_back("G000026345");
    escherichia_list.push_back("G000183345");
    escherichia_list.push_back("G000299455");
    escherichia_list.push_back("G000759795");
    escherichia_list.push_back("G001283625");
    for (int i = 0; i < escherichia_list.size(); i++)
        escherichia.insert(escherichia_list[i]);
}

static void init_klebsiella_lists()
{
    klebsiella_list.push_back("G000215745");
    klebsiella_list.push_back("G000240185");
    klebsiella_list.push_back("G001022195");
    klebsiella_list.push_back("G001187865");
    klebsiella_list.push_back("G001887595");
    for (int i = 0; i < klebsiella_list.size(); i++)
        klebsiella.insert(klebsiella_list[i]);
}

static void init_bacillus_lists()
{
    b_anthracis_list.push_back("G000007845");
    b_anthracis_list.push_back("G000008165");
    for (int i = 0; i < b_anthracis_list.size(); i++)
        b_anthracis.insert(b_anthracis_list[i]);

    b_cereus_list.push_back("G000007825");
    b_cereus_list.push_back("G000291295");
    for (int i = 0; i < b_cereus_list.size(); i++)
        b_cereus.insert(b_cereus_list[i]);

    b_thuringiensis_list.push_back("G000008505");
    b_thuringiensis_list.push_back("G001420855");
    for (int i = 0; i < b_thuringiensis_list.size(); i++)
        b_thuringiensis.insert(b_thuringiensis_list[i]);
}

int main(int argc, char** argv)
{
  init_akkermansia_lists();
  init_ecoli_lists();
  init_yersinia_lists();
  init_klebsiella_lists();
  init_bacillus_lists();
  std::vector<SAMLine> batch;
  std::string last_qname;
  std::string line;
  while (std::getline(std::cin, line))
  {
    if (line[0] == '@')
      continue;  // Header info

    std::istringstream parseline(line);
    SAMLine samline;

    bool success = parseline
      >> samline.qname
      >> samline.flag
      >> samline.rname
      >> samline.pos
      >> samline.mapq
      >> samline.cigar
      >> samline.rnext
      >> samline.pnext
      >> samline.tlen
      >> samline.seq
      >> samline.qual;

    if (!success)
    {
      std::cout << "Oops, parsing failure";
      std::cout << "Last line\n" << line << std::endl;
      exit(-1);
    }

    int read_len;
    int reference_len;
    parse_cigar(samline.cigar, read_len, reference_len);

    samline.reference_cover.start = samline.pos;
    samline.reference_cover.end = samline.pos + reference_len - 1;

    if (last_qname != samline.qname)
      batch_process(batch);

    batch.push_back(samline);

    last_qname = samline.qname;
  }
  batch_process(batch);


  // Print out scratch space for visualization
  for (std::map<std::string, int>::iterator private_it = private_reads.begin(); private_it != private_reads.end(); private_it++)
  {
    std::cout << private_it->first << ',' << private_it->second << std::endl;
  }

  std::cout << "---" << std::endl;

  for (std::map<std::string, float>::iterator split_it = split_reads.begin(); split_it != split_reads.end(); split_it++)
  {
    std::cout << split_it->first << ',' << split_it->second << std::endl;
  }

  std::cout << "---" << std::endl;

  for (std::map<std::string, std::map<std::string, float> >::iterator shared2_it = shared2_reads.begin(); shared2_it != shared2_reads.end(); shared2_it++)
  {
    const std::string& rname1 = shared2_it->first;
    for (std::map<std::string, float>::iterator shared2_it2 = shared2_it->second.begin(); shared2_it2 != shared2_it->second.end(); shared2_it2++)
    {
      const std::string& rname2 = shared2_it2->first;
      float count = shared2_it2->second;

      std::cout << rname1 << ',' << rname2 << ',' << count << std::endl;
    }
  }

//  std::map<std::string, std::map<std::string, std::vector<IndexRange> > >& to_write = shared2_row_coverage;
  std::map<std::string, std::map<std::string, std::vector<IndexRange> > >& to_write = shared_powerset_row_coverage;

  std::cout << "---" << std::endl;
  for (std::map<std::string, std::map<std::string, std::vector<IndexRange> > >::iterator shared2_it = to_write.begin(); shared2_it != to_write.end(); shared2_it++)
  {
    const std::string& rname1 = shared2_it->first;
    for (std::map<std::string, std::vector<IndexRange> >::iterator shared2_it2 = shared2_it->second.begin(); shared2_it2 != shared2_it->second.end(); shared2_it2++)
    {
      const std::string& rname2 = shared2_it2->first;
      std::vector<IndexRange>& row_reads = shared2_it2->second;

      std::vector<IndexRange> compressed = compress_ranges(row_reads);

      std::cout << rname1 << ',' << rname2 << ',' << compressed.size() << std::endl;
      for (int i = 0; i < compressed.size(); i++)
        std::cout << compressed[i].start << '-' << compressed[i].end << std::endl;
    }
  }

  std::cout << "---" << std::endl;
  std::cout << "TotalYersinia" << "," << "PrivateYersinia" << "," << "YersiniaEscherichia" << "," << "YersiniaKlebsiella" << "," << "YersiniaEscherichiaKlebsiella" << std::endl;
  std::cout << total_yersinia_reads << "," << yersinia_only_reads << "," << yersinia_escherichia_reads << "," << yersinia_klebsiella_reads << "," << yersinia_escherichia_and_klebsiella_reads << std::endl;
  std::cout << "---" << std::endl;
  std::cout << "TotalAnthracis" << "," << "PrivateAnthracis" << "," << "AnthracisCereus" << "," << "AnthracisThuringiensis" << "," << "AnthracisCereusThuringiensis" << std::endl;
  std::cout << total_anthracis_reads << "," << anthracis_only_reads << "," << anthracis_cereus_reads << "," << anthracis_thuringiensis_reads << "," << anthracis_cereus_and_thuringiensis_reads << std::endl;
}

//https://en.wikipedia.org/wiki/Sequence_alignment#CIGAR_Format
void parse_cigar(std::string cigar, int& read_len, int& reference_len)
{
  read_len = 0;
  reference_len = 0;

  std::istringstream cigar_parser(cigar);
  while (cigar_parser.rdbuf()->in_avail() != 0)
  {
    int len;
    char type;
    bool success = cigar_parser >> len >> type;

    if (!success)
    {
      std::cout << "Couldn't parse cigar string: " << cigar << std::endl;
      exit(-1);
    }

    switch(type)
    {
      case 'M':
      case 'I':
      case 'S':
      case '=':
      case 'X':
        read_len += len;
        break;
    }
    switch(type)
    {
      case 'M':
      case 'D':
      case 'N':
      case '=':
      case 'X':
        reference_len += len;
        break;
    }
  }
}

void batch_process(std::vector<SAMLine>& batch)
{
  if (batch.size() == 0)
    return;

  std::vector<batch_processor*> procs;
  // procs.push_back(debug_batch);
  // procs.push_back(count_private_and_split_reads);
  // procs.push_back(count_shared_reads);
  // procs.push_back(calc_akkermansia_coverage);
  // procs.push_back(calc_akkermansia_coverage_powerset);
  // procs.push_back(track_yersinia_reads);
  procs.push_back(track_bacillus_reads);

  for (int i = 0; i < procs.size(); i++)
    procs[i](batch);

  batch.clear();
}

void debug_batch(std::vector<SAMLine>& batch)
{
  std::string& qname = batch[0].qname;
  int match_count = batch.size();

  std::cout << qname << " Matches:" << match_count << " [";
  for (int i = 0; i < batch.size(); i++)
  {
    SAMLine& samline = batch[i];
    std::cout << samline.rname;
    if (i != batch.size() - 1)
      std::cout << ", ";
  }
  std::cout << ']' << std::endl;
}


void count_private_and_split_reads(std::vector<SAMLine>& batch)
{
  std::unordered_set<std::string> s;
  for (int i = 0; i < batch.size(); i++)
    s.insert(batch[i].rname);

  if (s.size() == 1)
  {
    //Read matches only to this genome (possibly to multiple spots in this genome)
    private_reads[batch[0].rname]++;
  }
  else
  {
    // Need to ask exact specifics of read splitting
    // Does it split based on how many copies exist in a single genome?
    // Or does it split based on total number of genomes?
    for (int i = 0; i < batch.size(); i++)
    {
      split_reads[batch[i].rname] += 1.0f / batch.size();
    }
  }
}

void count_shared_reads(std::vector<SAMLine>& batch)
{
  std::set<std::string> set;
  for (int i = 0; i < batch.size(); i++)
    set.insert(batch[i].rname);

  if (set.size() > 1)
  {
    for (int i = 0; i < batch.size(); i++)
      for (int j = 0; j < batch.size(); j++)
        shared2_reads[batch[i].rname][batch[j].rname] += 1.0f / (batch.size() * batch.size());
  }
}

void calc_akkermansia_coverage(std::vector<SAMLine>& batch)
{
  for (int i = 0 ; i < batch.size(); i++)
  {
    for (int j = 0; j < batch.size(); j++)
    {
      SAMLine& samline_row = batch[i];
      SAMLine& samline_col = batch[j];

      //If its not akkermansia, we don't care for now
      if (akkermansia_muciniphila.count(samline_row.rname) == 0 ||
         akkermansia_muciniphila.count(samline_col.rname) == 0)
         continue;

      shared2_row_coverage[samline_row.rname][samline_col.rname].push_back(samline_row.reference_cover);
    }
  }
}

void calc_akkermansia_coverage_powerset(std::vector<SAMLine>& batch)
{
  std::set<std::string> set;
  for (int i = 0; i < batch.size(); i++)
  {
    SAMLine& samline_row = batch[i];
    //If its not akkermansia, we don't care for now
    if (akkermansia_muciniphila.count(samline_row.rname) == 0)
      continue;
    set.insert(batch[i].rname);
  }

  std::vector<std::string> sorted_set(set.begin(), set.end());
  std::sort(sorted_set.begin(), sorted_set.end());

  std::ostringstream outstream;
  for (int i = 0; i < sorted_set.size(); i++)
  {
    outstream << sorted_set[i];
    if (i < sorted_set.size() - 1)
      outstream << ";";
  }
  std::string setname = outstream.str();

  for (int i = 0; i < batch.size(); i++)
  {
    SAMLine& samline_row = batch[i];
    //If its not akkermansia, we don't care for now
    if (akkermansia_muciniphila.count(samline_row.rname) == 0)
      continue;

    shared_powerset_row_coverage[samline_row.rname][setname].push_back(samline_row.reference_cover);
  }
}

void track_yersinia_reads(std::vector<SAMLine>& batch)
{
  bool found_yersinia = false;
  for (int i = 0 ; i < batch.size(); i++)
  {
    SAMLine& samline_row = batch[i];

    //If its not yersinia, we don't care
    if (yersinia.count(samline_row.rname) == 0)
       continue;
    found_yersinia = true;
  }

  if (!found_yersinia)
    return;

  bool found_ecoli = false;
  bool found_klebsiella = false;

  for (int i = 0 ; i < batch.size(); i++)
  {
    SAMLine& samline_row = batch[i];

    if (escherichia.count(samline_row.rname) != 0)
        found_ecoli = true;
    if (klebsiella.count(samline_row.rname) != 0)
        found_klebsiella = true;
  }

  total_yersinia_reads++;
  if (found_ecoli && found_klebsiella)
      yersinia_escherichia_and_klebsiella_reads++;
  else if (found_ecoli && !found_klebsiella)
      yersinia_escherichia_reads++;
  else if (!found_ecoli && found_klebsiella)
      yersinia_klebsiella_reads++;
  else if (!found_ecoli && !found_klebsiella)
      yersinia_only_reads++;
}

void track_bacillus_reads(std::vector<SAMLine>& batch)
{
  bool found_anthracis = false;
  for (int i = 0 ; i < batch.size(); i++)
  {
    SAMLine& samline_row = batch[i];

    //If its not anthrax, we don't care
    if (b_anthracis.count(samline_row.rname) == 0)
       continue;
    found_anthracis = true;
  }

  if (!found_anthracis)
    return;

  bool found_cereus = false;
  bool found_thuringiensis = false;

  for (int i = 0 ; i < batch.size(); i++)
  {
    SAMLine& samline_row = batch[i];

    if (b_cereus.count(samline_row.rname) != 0)
        found_cereus = true;
    if (b_thuringiensis.count(samline_row.rname) != 0)
        found_thuringiensis = true;
  }

  total_anthracis_reads++;
  if (found_cereus && found_thuringiensis)
      anthracis_cereus_and_thuringiensis_reads++;
  else if (found_cereus && !found_thuringiensis)
      anthracis_cereus_reads++;
  else if (!found_cereus && found_thuringiensis)
      anthracis_thuringiensis_reads++;
  else if (!found_cereus && !found_thuringiensis)
      anthracis_only_reads++;
}
