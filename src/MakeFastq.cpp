#include "CLI/CLI.hpp"
#include "CanonicalKmerIterator.hpp"
#include "FastxParser.hpp"
#include "string_view.hpp"
#include "FileOps.hpp"
#include "ScopedTimer.hpp"
#include <cmath>
#include <iterator>
#include <iostream>
#include <vector>


void sampleKmers(std::string fastaFile,
                      uint32_t k){

  CanonicalKmer::k(k);

  {
    //ScopedTimer st ;
    std::vector<std::string> read_file = {fastaFile} ;
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1) ;
    parser.start() ;

    CanonicalKmerIterator kit_end;
    std::cerr<< "\n Num of line is fastq "<< "\n" ;
    size_t globalId{0};
    size_t rn{0} ;
    size_t numOfKmers{0} ;
    auto rg = parser.getReadGroup() ;
    while(parser.refill(rg)){
      for(auto& rp: rg){

        auto& r1 = rp.seq ;
        CanonicalKmerIterator kit1(r1);

        while(kit1 != kit_end){
            std::cout << kit1->first.to_str() << "\n";
            numOfKmers++ ;
            if(numOfKmers == 2){
                std::cerr << "skipping 10 kmers" << "\n" ;
                kit1 += 3 ;
            }else{
                ++kit1 ;
            }
        }
        std::cerr << "\nTotal kmers " << r1.length() - k + 1 << "\n" ;
        std::cerr << "\nnum of seen kmers " << numOfKmers << "\n" ;

        ++rn ;


      }
    }
  }
}


int main(int argc, char* argv[]){
  CLI::App app{"Playing tools"} ;
  auto fastaqApp = app.add_subcommand("sample", "sample CanonicalKmers ") ;

  //options for creating the required files
  uint32_t k{31} ;
  std::string fastaFile ;
  //end

  fastaqApp
    ->add_option("-f,--fasta", fastaFile,"the full path to the fasta file")
    ->required() ;


  //->add_flag("-p,--library", isPairedEnd, "paired end or not") ;

  fastaqApp
    ->add_option("-k,--k-mer", k,"read length") ;

  try{
    app.parse(argc,argv) ;
  } catch (const CLI::ParseError& e){
    std::cerr << "\n Not proper arguments\n " ;
    app.exit(e) ;
  }

  if(app.got_subcommand(fastaqApp)){
    sampleKmers(fastaFile, k) ;
  }else{
    std::cerr << "wrong subcommand\n" ;
  }


}
