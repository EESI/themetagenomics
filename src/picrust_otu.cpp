#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <Rcpp.h>

#include "gzstream.h"

using namespace std;
using namespace Rcpp;

vector<string> split_line(const string &source, const char *delimiter = " ") {

  vector<string> results;

  size_t prev = 0;
  size_t next = 0;

  while ((next = source.find_first_of(delimiter,prev)) != string::npos) {

    if (next - prev != 0) {
      results.push_back(source.substr(prev,next - prev));
    }

    prev = next + 1;

  }

  if (prev < source.size()) {
    results.push_back(source.substr(prev));
  }

  return results;

}

List sweep_picrust(string file_path, StringVector otu_id_targets) {

  int a=0;
  int o=0;
  List out;
  string rec;
  StringVector lines;
  StringVector matches;
  StringVector genemeta;

  ifstream file(file_path.c_str()); // load in file path

  getline(file,rec,'\n'); // start with line 1, column names
  int id_pos = rec.find('\t');
  int id_pos_end = rec.find("metadata");
  int pimeta_pos_end = rec.find("\n");

  // extract from #OTU to picrust metadata_*
  string header_left = rec.substr(id_pos + 1,id_pos_end - 10);
  vector<string> gene_ids = split_line(header_left,"\t");

  // extract from picrust metadata_* to end
  string header_right = rec.substr(id_pos_end,pimeta_pos_end);
  vector<string> pimeta_ids = split_line(header_right,"\t");

  int gene_n=gene_ids.size();
  int pimeta_n=pimeta_ids.size();

  IntegerMatrix genome_table(otu_id_targets.size(),gene_n);
  NumericMatrix pimeta_table(otu_id_targets.size(),pimeta_n);

  while(getline(file,rec,'\n')) { // starting at line 2, data row 1

    if (a % 5000 == 0){
      checkUserInterrupt();
    }
    a += 1;

    // check if at gene metadata lines (final 2)
    if (rec.substr(0,8) == "metadata"){

      int line_end = rec.find('\n');
      genemeta.push_back(rec.substr(id_pos + 1,line_end));

      // otherwise, extract counts
    }else{

      // get row name
      int id_pos = rec.find('\t');
      string otu_id = rec.substr(0,id_pos);

      // look over all otu ids from otu table to check
      // if row name match otu id from otu table
      for (int j=0;j<otu_id_targets.size();j++){

        if (otu_id_targets[j] == otu_id){

          int line_end = rec.find('\n');
          string counts_start = rec.substr(id_pos + 1,line_end);
          vector<string> counts = split_line(counts_start,"\t");

          for (int k=0;k<gene_n;k++){
            genome_table(o,k) = atoi(counts[k].c_str()); // str to int
          }

          int c=0;
          for (std::size_t m=gene_n;m<counts.size();m++){
            pimeta_table(o,c) = atof(counts[m].c_str()); // str to dbl
            c += c;
          }

          o += 1;

          // for (int s=0; s<otu_id_targets.size();s++){
          //   Rcout << otu_id_targets(s) << " ";
          // }

          matches.push_back(otu_id); // record otu id of match
          otu_id_targets.erase(j); // remove otu id from targets

          // printf("\niter %d, found a match: %s\n",o,otu_id.c_str());

          break;
        }

      }

    }

  }
  file.close();

  if (o == 0){

    return(R_NilValue);

  }else{

    SubMatrix<INTSXP> genome_table_out = genome_table(Range(0,o - 1),_);
    SubMatrix<REALSXP> pimeta_table_out = pimeta_table(Range(0,o - 1),_);

    out["gene_ids"] = gene_ids;
    out["pimeta_ids"] = pimeta_ids;
    out["matches"] = matches;
    out["genemeta"] = genemeta;
    out["genome_table_out"] = genome_table_out;
    out["pimeta_table_out"] = pimeta_table_out;

    return(out);

  }

}

List sweep_picrust_gz(string file_path, StringVector otu_id_targets) {

  int a=0;
  int o=0;
  List out;
  string rec;
  StringVector lines;
  StringVector matches;
  StringVector genemeta;

  igzstream file;
  file.open(file_path.c_str());

  getline(file,rec,'\n'); // start with line 1, column names
  int id_pos = rec.find('\t');
  int id_pos_end = rec.find("metadata");
  int pimeta_pos_end = rec.find("\n");

  // extract from #OTU to picrust metadata_*
  string header_left = rec.substr(id_pos + 1,id_pos_end - 10);
  vector<string> gene_ids = split_line(header_left,"\t");

  // extract from picrust metadata_* to end
  string header_right = rec.substr(id_pos_end,pimeta_pos_end);
  vector<string> pimeta_ids = split_line(header_right,"\t");

  int gene_n=gene_ids.size();
  int pimeta_n=pimeta_ids.size();

  IntegerMatrix genome_table(otu_id_targets.size(),gene_n);
  NumericMatrix pimeta_table(otu_id_targets.size(),pimeta_n);

  while(getline(file,rec,'\n')) { // starting at line 2, data row 1

    if (a % 5000 == 0){
      checkUserInterrupt();
    }
    a += 1;

    // check if at gene metadata lines (final 2)
    if (rec.substr(0,8) == "metadata"){

      int line_end = rec.find('\n');
      genemeta.push_back(rec.substr(id_pos + 1,line_end));

      // otherwise, extract counts
    }else{

      // get row name
      int id_pos = rec.find('\t');
      string otu_id = rec.substr(0,id_pos);

      // look over all otu ids from otu table to check
      // if row name match otu id from otu table
      for (int j=0;j<otu_id_targets.size();j++){

        if (otu_id_targets[j] == otu_id){

          int line_end = rec.find('\n');
          string counts_start = rec.substr(id_pos + 1,line_end);
          vector<string> counts = split_line(counts_start,"\t");

          for (int k=0;k<gene_n;k++){
            genome_table(o,k) = atoi(counts[k].c_str()); // str to int
          }

          int c=0;
          for (std::size_t m=gene_n;m<counts.size();m++){
            pimeta_table(o,c) = atof(counts[m].c_str()); // str to dbl
            c += c;
          }

          o += 1;

          matches.push_back(otu_id); // record otu id of match
          otu_id_targets.erase(j); // remove otu id from targets

          break;
        }

      }

    }

  }
  file.close();

  if (o == 0){

    return(R_NilValue);

  }else{

    SubMatrix<INTSXP> genome_table_out = genome_table(Range(0,o - 1),_);
    SubMatrix<REALSXP> pimeta_table_out = pimeta_table(Range(0,o - 1),_);

    out["gene_ids"] = gene_ids;
    out["pimeta_ids"] = pimeta_ids;
    out["matches"] = matches;
    out["genemeta"] = genemeta;
    out["genome_table_out"] = genome_table_out;
    out["pimeta_table_out"] = pimeta_table_out;

    return(out);

  }

}

//' Predict OTU functional potential via PICRUSt
//'
//' Given an OTU abundance table prepared with the GreenGenes reference database,
//' this function predicts the functional content using either COG or KO
//' precalculated mapping tables that map the taxonomic abundance for a given OTU
//' to functional abundance content across a set of functional genes.
//'
//' @param file_path Path to the precalculated table
//' @param otu_id_targets Character vector of OTU IDs to predict
//'
//' @return A list containing
//' \describe{
//' \item{gene_ids}{String vector of KO IDs, the column names in genome_table_out.}
//' \item{pimeta_ids}{String vector of names for the PICRUSt metadata categories,
//' the column names of pimeta_table_out.}
//' \item{matches}{String vector of OTU IDs from otu_id_targets that were present
//' in the mapping file.}
//' \item{genemeta}{String vector of functional metadata corresponding to gene_ids}
//' \item{genome_table_out}{Integer matrix of gene counts across topics}
//' \item{pimeta_table_out}{Numeric matrix of method specific metadata (NSTI)}
//' }
// [[Rcpp::export]]

List picrust_otu(std::string file_path, StringVector otu_id_targets) {

  string file_ext = file_path.substr(file_path.rfind('.'),file_path.size());

  if (file_ext == ".gz"){

    return(sweep_picrust_gz(file_path,otu_id_targets));

  }else{

    return(sweep_picrust(file_path,otu_id_targets));

  }

}
