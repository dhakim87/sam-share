#pragma once

#include <iostream>
#include <vector>

struct SAMLine;

void parse_cigar(std::string, int&, int&);
void batch_process(std::vector<SAMLine>&);

typedef void (batch_processor)(std::vector<SAMLine>&);
void debug_batch(std::vector<SAMLine>&);
void count_private_and_split_reads(std::vector<SAMLine>&);
void count_shared_reads(std::vector<SAMLine>&);


struct SAMLine
{
  std::string qname;
  int flag;
  std::string rname;
  int pos;
  int mapq;
  std::string cigar;
  std::string rnext;
  int pnext;
  int tlen;
  std::string seq;
  std::string qual;
};
