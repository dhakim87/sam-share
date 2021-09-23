#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>


struct SAMLine;
struct IndexRange
{
  int start;
  int end;
};

void parse_cigar(std::string, int&, int&);
void batch_process(std::vector<SAMLine>&);

typedef void (batch_processor)(std::vector<SAMLine>&);
void debug_batch(std::vector<SAMLine>&);
void count_private_and_split_reads(std::vector<SAMLine>&);
void count_shared_reads(std::vector<SAMLine>&);
void calc_akkermansia_coverage(std::vector<SAMLine>&);

bool compare_range_starts(IndexRange i,IndexRange j) { return (i.start<j.start); }
//https://github.com/ucsd-cmi/zebra_filter/blob/master/cover.py
std::vector<IndexRange> compress_ranges(std::vector<IndexRange>& interleaved_ranges)
{
  // Sort ranges by start index
  std::sort(interleaved_ranges.begin(), interleaved_ranges.end(), compare_range_starts);

  std::vector<IndexRange> new_ranges;
  int start_val = -1;
  int end_val = -1;

  IndexRange new_range;

  for (int i = 0; i < interleaved_ranges.size(); i++)
  {
    IndexRange& r = interleaved_ranges[i];
    if (end_val == -1)
    {
      //case 1: no active range, start active range.
      start_val = r.start;
      end_val = r.end;
    }
    else if (end_val >= r.start - 1)
    {
      //case 2: active range continues through this range
      //extend active range
      if (end_val < r.end)
        end_val = r.end;
    }
    else {
      // if end_val < r.start - 1:
      //case 3: active range ends before this range begins
      //write new range out, then start new active range
      new_range.start = start_val;
      new_range.end = end_val;

      new_ranges.push_back(new_range);
      start_val = r.start;
      end_val = r.end;
    }
  }

  if (end_val != -1)
  {
    new_range.start = start_val;
    new_range.end = end_val;
    new_ranges.push_back(new_range);
  }

  return new_ranges;
}

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
  IndexRange reference_cover;
};
