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
#include <iostream>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;
using namespace std;

class Diskutil {
private:
public:

  /**
   * Given the path to a directory, fetches all items (files, folders) in it.
   * @param  dir Path to the target directory.
   * @return     An array with all the items contained in the directory.
   */
  static vector<path> getItemsInDir(path parent);

  /**
   * Checks if a folder is valid for use.
   * @param  dir Path to the target directory.
   * @return     True if the directory is not a file or a wildcard (".", "..").
   */
  static bool isFolderValid(path dir);

  /**
   * Checks if a file has the given extension.
   * @param  dir Path to the target file.
   * @param  ext Desired extension.
   * @return     True if the file has the desired extension.
   */
  static bool hasFileExtension(path dir, string ext);

  /**
   * Get paths for all genomes present in the target directory.
   * Expects target directory to be filled with separate subdirectories for
   * each class, containing the genomes pertaining to it.
   * @param  source_folder  The target directory.
   * @param  extension      The extension used for kmer counts.
   * @return                An array containing tuples representing the
   *                        genome's class, kmer count path and sequence path,
   *                        in that order.
   */
  static vector<tuple<string, path, path> >
    getTrainingGenomePaths(path source_folder, string extension);

  /**
   * Gets the size of the file, in kilobytes.
   * @param  file Path to the file on disk.
   * @return      The filesize in kilobytes.
   */
  static int getFileSize(path file);

  static const string SEQUENCE_FILE_EXT; // = ".fasta"
  static const string SAVE_FILE_EXT; // = ".dat"
protected:
};

#endif /* Class_hpp */
