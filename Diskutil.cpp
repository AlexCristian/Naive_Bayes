//
//  Diskutil.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 05/05/2017.
//
//

#include "Diskutil.hpp"

const string Diskutil::SEQUENCE_FILE_EXT = ".fasta";
const string Diskutil::SAVE_FILE_EXT = ".dat";
const int BTOKB = 1024;

vector<path> Diskutil::getItemsInDir(path parent){
  vector<path> result;

  if(!exists(parent) || !is_directory(parent)){
    cout<<"Could not open directory path\n";
    exit(1);
  }
  directory_iterator iter(parent);
  for(;iter != directory_iterator(); iter++){
    result.push_back(iter->path());
  }
  return result;
}

int Diskutil::getFileSize(path file){
  return file_size(file) / BTOKB;
}

bool Diskutil::isFolderValid(path dir){
  if(dir.native().compare(".") == 0 || dir.native().compare("..") == 0){
    return false;
  }
  return exists(dir) && is_directory(dir);
}

bool Diskutil::hasFileExtension(path dir, string ext){
  if(exists(dir) && is_regular_file(dir)){
    return dir.has_extension() && (dir.extension().compare(ext) == 0);
  }else{
    return false;
  }
}

vector<tuple<string, path, path> >
  Diskutil::getTrainingGenomePaths(path source_folder, string extension){
    vector<tuple<string, path, path> > result;

    vector<path> class_paths = getItemsInDir(source_folder);
    for(vector<path>::iterator iter=class_paths.begin();
        iter != class_paths.end(); iter++){
      if(!isFolderValid(*iter)){
        continue;
      }

      string cls(iter->native());
      string cls_s;
      try{
        cls_s = cls.substr(cls.rfind(path::preferred_separator)+1);
      }catch(invalid_argument){
        continue;
      }

      vector<path> genomes = getItemsInDir(*iter);
      for(vector<path>::iterator gen=genomes.begin();
          gen != genomes.end(); gen++){
        if(!hasFileExtension(*gen, extension)){
          continue;
        }

        string seqFile = gen->native();
        seqFile = seqFile.substr(0, seqFile.rfind("."));
        seqFile += SEQUENCE_FILE_EXT;

        result.push_back(make_tuple(cls_s, *gen, path(seqFile)));
      }
    }

    return result;
  }
