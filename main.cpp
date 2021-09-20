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


int main(int argc, char** argv)
{
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

    // int read_len;
    // int reference_len;
    // parse_cigar(samline.cigar, read_len, reference_len);

    // pos is 1 based index into the genome
    // for an inclusive 1 based start and end range in the reference, use:
    // start = pos
    // end = pos + reference_len - 1;

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
  procs.push_back(count_private_and_split_reads);
  procs.push_back(count_shared_reads);

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
