//
//  Diskutil.hpp
//  NB_C++
//
//  Created by Alexandru Cristian on 05/05/2017.
//
//

#ifndef Diskutil_hpp
#define Diskutil_hpp

#include <vector>
#include <utility>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;
using namespace std;

class Diskutil {
private:
public:
  static vector<path> getItemsInDir(path parent);
  static bool isFolderValid(path dir);
  static bool hasFileExtension(path dir, string ext);
  static vector<tuple<string, path, path> >
    getTrainingGenomePaths(path source_folder, string extension);

  static const string SEQUENCE_FILE_EXT; // = ".fasta"
  static const string SAVE_FILE_EXT; // = ".dat"
protected:
};

#endif /* Class_hpp */
