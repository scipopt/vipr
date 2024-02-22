/*
*
*   Copyright (c) 2024 Zuse Institute Berlin
*
*   Permission is hereby granted, free of charge, to any person obtaining a
*   copy of this software and associated documentation files (the "Software"),
*   to deal in the Software without restriction, including without limitation
*   the rights to use, copy, modify, merge, publish, distribute, sublicense,
*   and/or sell copies of the Software, and to permit persons to whom the
*   Software is furnished to do so, subject to the following conditions:
*
*   The above copyright notice and this permission notice shall be included in
*   all copies or substantial portions of the Software.
*
*   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
*   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
*   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
*   DEALINGS IN THE SOFTWARE.
*
*/

// Includes
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <limits>
#include "soplex.h"
#include <memory>
#include <soplex/dsvector.h>
#include <soplex/lprow.h>
#include <boost/bimap.hpp>
#include "CMakeConfig.hpp"

// Timing
#include <sys/time.h>

// Parallelization
#include <tbb/tbb.h>

// Namespaces
using namespace std;
using namespace soplex;

// sparse boolean vectors
typedef map<int, bool> SVectorBool;
typedef shared_ptr<DSVectorRational> DSVectorPointer;
class SVectorRat : public map<int, Rational>
{
   public:
      void compactify() {  if( !_compact )
                           {
                              auto it = this->begin();
                              while( it != this->end() )
                              {
                                 if( it->second == 0 ) this->erase(it++);
                                 else ++it;
                              }
                              _compact = true;
                           }
                        }
      bool operator!=(SVectorRat &other);
      bool operator==(SVectorRat &other) { return !(*this != other);}

   private:
      bool _compact = false;
};


// Input and output files
ifstream certificateFile;
ofstream completedFile;

// Settings
bool debugmode = false;
bool usesoplex = true;

// Global variables
DSVectorRational dummycol(0); // SoPlex placeholder column
vector<string> variableNames; // variable names
DSVectorPointer ObjCoeff(make_shared<DSVectorRational>()); // sparse vector of objective coefficients
size_t numberOfVariables;

unsigned int nthreads = thread::hardware_concurrency();

struct Constraint {DSVectorPointer vec; Rational side; int sense; string line;}; // @todo: take care of deep copies here!
vector<Constraint> constraints;

typedef boost::bimap<int, long> bimap; // maps rows of the LP to corresponding indices in the certificate used for updating the local LPs
map<int, long> origConsCertIndex;

vector<tuple<Rational,Rational,long>> lowerBounds; // rational boundval, multiplier, long certindex
vector<tuple<Rational,Rational,long>> upperBounds; // rational boundval, multiplier, long certindex


struct passData { SoPlex passLp; bimap LProwCertificateMap;};
struct parallelData { passData* bufData; stringstream linestream; bool needscompletion; std::vector<size_t> conidx; string line;};

// circular buffer (also known as ring buffer, circular queue, cyclic queue), references:
// https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Using_Circular_Buffers.html
// https://stackoverflow.com/a/15167828/15777342
class circBuf
{
   int tail, head, size;
   std::vector<passData*> arr;
   public:
   circBuf( int maxtokens )
   {
      head = tail = -1;
      size = maxtokens;
      arr.resize(maxtokens);
   }
   /// assignment operator
   void resize( int maxtokens )
   {
      head = tail = -1;
      size = maxtokens;
      arr.resize(maxtokens);
   }
   void enqueue(passData* warmStartData);
   passData* dequeue();
   bool isEmpty() { return head == -1; }
};

void circBuf::enqueue( passData* warmStartData )
{
   if( head == -1 )
   {
      head = tail = 0;
      arr[tail] = warmStartData;
   }
   else if( (tail == size-1) && (head != 0) )
   {
      tail = 0;
      arr[tail] = warmStartData;
   }
   else
   {
      tail++;
      arr[tail] = warmStartData;
   }
}

passData* circBuf::dequeue()
{
   passData* data = arr[head];
   if( head == tail )
   {
      head = -1;
      tail = -1;
   }
   else if( head == size-1 )
   {
      head = 0;
   }
   else
   {
      head++;
   }
   return data;
}



// Forward declaration
void modifyFileName( string &path, const string &newExtension );
bool checkversion( string ver );
bool processVER();
bool processVAR( SoPlex &workinglp );
bool processINT();
bool processOBJ( SoPlex &workinglp );
bool processCON( SoPlex &workinglp );
bool processRTP();
bool processSOL();
bool processDER( SoPlex workinglp );
bool getConstraints( SoPlex &workinglp, string &consense, Rational &rhs, int &activeConstraint, size_t currentDerivation);
std::string completelin( SoPlex &workinglp, bimap& LProwCertificateMap, Constraint &constraint);

static bool completeWeakDomination(DSVectorRational &row, int consense, Rational &rhs, stringstream& completedLine,
                                    stringstream &initialLine);
static bool readLinComb( int &sense, Rational &rhs, SVectorRat& coefficients, SVectorRat& mult,
                  int currentConstraintIndex, SVectorBool &assumptionList, stringstream &baseLine);
static bool readMultipliers( int &sense, SVectorRat &mult, stringstream &initialLine );
bool completeIncomplete( SoPlex &workinglp, bimap& LProwCertificateMap, vector<long> &activeDerivations, string label, stringstream& completedLine,
                        stringstream &baseLine );
bool printReasoningToLine(DVectorRational &dualmultipliers, DVectorRational &reducedcosts, stringstream& completedLine, bimap& LProwCertificateMap);

// Global bound changes for completing weak
static void processGlobalBoundChange(Rational rhs, Rational boundmult, int varindex,
                                       long boundindex, int sense)
{
   char restline[256];
   string tmpstring;
   string lastword;

   if( certificateFile.peek() != '\n' )
   {
      certificateFile.getline(restline, 256);
      tmpstring = restline;
      istringstream tmpstream(tmpstring);

      while(tmpstream >> lastword)
      {
         if( lastword == "global" )
         {
            if( sense <= 0 ) // eq or ub
               upperBounds[varindex] = make_tuple(rhs, boundmult, boundindex);
            if( sense >= 0 ) // eq or lb
               lowerBounds[varindex] = make_tuple(rhs, boundmult, boundindex);
         }
      }
   }
}

static double getTimeSecs(timeval start, timeval end)
{
   auto seconds = end.tv_sec - start.tv_sec;
   auto microseconds = end.tv_usec - start.tv_usec;
   auto clock_dur = seconds + microseconds*1e-6;
   return clock_dur;
}

// Set usage options
static void printUsage(const char* const argv[], int idx)
{
   const char* usage =
      "general options:\n"
      "  --soplex=on/off       use soplex to complete derivations? needs to be on if incomplete derivations in certificate file;\
      \n                        turn off to boost performance if only weak derivations are present.\n"
      "  --debugmode=on/off    enable extra debug output from viprcomp\n"
      "  --verbosity=<level>   set verbosity level inside SoPlex\n"
      "  --threads=<number>    maximal number of threads to use \n"
      "\n";
   if(idx <= 0)
      cerr << "missing input file\n\n";
   else
      cerr << "invalid option \"" << argv[idx] << "\"\n\n";

   cerr << "usage: " << argv[0] << " " << "[options] <certificateFile>\n"
             << "  <certificateFile>               .vipr file to be completed\n\n"
             << usage;
}

// Main function
int main(int argc, char *argv[])
{
   int returnStatement = -1;
   int optidx;
   const char* certificateFileName;
   int verbosity = 0;
   string path = "";

   if( argc == 0 )
   {
      printUsage(argv, -1);
      return 1;
   }

   // read arguments from command line
   for(optidx = 1; optidx < argc; optidx++)
   {
      char* option = argv[optidx];

      // we reached <certificateFile>
      if(option[0] != '-')
      {
         certificateFileName = argv[optidx];
         continue;
      }

      // option string must start with '-', must contain at least two characters, and exactly two characters if and
      // only if it is -x, -y, -q, or -c
      if(option[0] != '-' || option[1] == '\0'
            || ((option[2] == '\0') != (option[1] == 'x' || option[1] == 'X' || option[1] == 'y'
                                        || option[1] == 'Y' || option[1] == 'q' || option[1] == 'c')))
      {
         printUsage(argv, optidx);
         return 1;
      }

      switch(option[1])
      {
      case '-' :
         option = &option[2];

         // set verbosity of SoPlex
         if(strncmp(option, "verbosity=", 10) == 0)
         {
            char* str = &option[10];
            if( isdigit(option[10]))
            {
               verbosity = atoi(str);
               if( verbosity < 0 || verbosity > 5 )
               {
                  cerr << "Verbosity level outside range 0 to 5. Read " << verbosity << " instead." << endl;
                  printUsage(argv, optidx);
                  return 1;
               }
            }
         }
         // enable additional debug output
         else if(strncmp(option, "debugmode=", 10) == 0)
         {
            char* str = &option[10];
            // Set Debugmode
            if( string(str) == "on")
            {
               debugmode = true;
               cout << "Debugmode turned on." << endl;
            }
            else if( string(str) == "off")
            {
               debugmode = false;
               cout << "Debugmode turned off." << endl;
            }
            else
            {
               cout << "Unknown input for debug settings (on/off expected). Read "
               << string(str) << " instead" << endl;
               cout << "Continue with default setings (debugmode off)" << endl;
            }
         }
         // enable/disable soplex -> faster if no incomplete, necessary if incomplete
         else if(strncmp(option, "soplex=", 7) == 0)
         {
            char* str = &option[7];
            // Set Debugmode
            if( string(str) == "on")
            {
               usesoplex = true;
               cout << "SoPlex turned on." << endl;
            }
            else if( string(str) == "off")
            {
               usesoplex = false;
               cout << "SoPlex turned off." << endl;
            }
            else
            {
               cout << "Unknown input for SoPlex suopport (on/off expected). Read "
               << string(str) << " instead" << endl;
               cout << "Continue with default setings (SoPlex on)" << endl;
            }
         }
         // set maximal number of threads that should be used
         else if(strncmp(option, "threads=", 8) == 0)
         {
            char* str = &option[8];
            if( isdigit(option[8]))
            {
               int maxthreads = atoi(str);
               if( maxthreads < 0 || maxthreads > thread::hardware_concurrency() )
               {
                  cerr << "threads outside of valid range: specified " << maxthreads << " but maximal number is " << thread::hardware_concurrency() << endl;
                  printUsage(argv, optidx);
                  return 1;
               }
               else
                  nthreads = maxthreads;
            }
         }
         else if(strncmp(option, "outfile=", 8) == 0)
         {
            path = string(&option[8]);
         }
         else
         {
            printUsage(argv, optidx);
            return 1;
         }
      }
   }

   certificateFile.open(certificateFileName);

   if( certificateFile.fail() )
   {
      cerr << "Failed to open file " << argv[1] << endl;
      printUsage(argv, -1);
      return returnStatement;
   }

   if(path == "")
   {
      path = certificateFileName;
      modifyFileName(path, "_complete.vipr");
   }

   completedFile.open( path.c_str(), ios::out );

   if( completedFile.fail() )
   {
      cerr << "Failed to open file " << path << endl;
      printUsage(argv, -1);
      return returnStatement;
   }

   // Set parameters for exact solving
   SoPlex baselp;
   baselp.setIntParam(SoPlex::READMODE, SoPlex::READMODE_RATIONAL);
   baselp.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
   baselp.setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_RATIONAL);
   baselp.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
   baselp.setRealParam(SoPlex::FEASTOL, 0.0);
   baselp.setRealParam(SoPlex::OPTTOL, 0.0);

   baselp.setIntParam(SoPlex::VERBOSITY, verbosity); // verbosity var prob. obsolete for parallelization

   struct timeval start, end;
   gettimeofday( &start, 0);
   if( processVER() )
      if( processVAR(baselp) )
         if( processINT() )
            if( processOBJ(baselp) )
               if( processCON(baselp) )
                  if( processRTP() )
                     if( processSOL() )
                     {
                        gettimeofday( &end, 0 );

                        cout << endl << "reading took " << getTimeSecs(start, end)
                             << " seconds (Wall Clock)" << endl;
                        gettimeofday( &start, 0 );
                        if( processDER(baselp) )
                        {
                              cout << "Completion of File successful!" <<endl;
                              returnStatement = 0;
                              gettimeofday( &end, 0 );

                              cout << endl << "completing and printing took " << getTimeSecs(start, end)
                                   << " seconds (Wall Clock)" << endl;
                        }
                     }
   return returnStatement;
}

// append "_complete.vipr" to filename
void modifyFileName(string &path, const string& newExtension)
{
   string::size_type position = path.find_last_of('.');
   if( position == string::npos )
      position = path.size();
   path.replace(position, newExtension.length(), newExtension);
}

// Version control for .vipr input file. Backward compatibility possible for minor versions
bool checkVersion(string version)
{
   bool returnStatement = false;

   size_t position = version.find(".");
   completedFile << " " + version;
   int major = atoi(version.substr(0, position).c_str());
   int minor = atoi(version.substr(position+1, version.length()-position).c_str());

   cout << "Certificate format version " << major << "." << minor << endl;

   if( (major ==VIPR_VERSION_MAJOR) && (minor <=VIPR_VERSION_MINOR) )
   {
      returnStatement = true;
   }
   else
   {
      cerr << "Version unsupported" << endl;
   }

   return returnStatement;
}

// Check version and correct file format allows next process to start
// Error if the version is incompatible or not specified
bool processVER()
{
   bool returnStatement = false;
   string tmpStr;

   for( ;; )
   {
      certificateFile >> tmpStr;
      if( tmpStr == "VER" )
      {
         completedFile << tmpStr;
         certificateFile >> tmpStr;

         returnStatement = checkVersion(tmpStr);
break;
      }
      else if( tmpStr == "%" )
      {
         getline( certificateFile, tmpStr );
      }
      else
      {
         cerr << endl << "Comment or VER expected. Read instead "
                << tmpStr << endl;
break;
      }
   }

   return returnStatement;
}

// Processes number of variables
// Puts placeholder variables in SoPlex LP
// Error if nr of variables invalid or number of variables > specified variables or section missing
bool processVAR(SoPlex &workinglp)
{
   cout << endl << "Processing VAR section..." << endl;

   auto returnStatement = true;

   string section;

   certificateFile >> section;



   // Check section
   if( section != "VAR" )
   {
      cerr << "VAR expected.   Read instead " << section << endl;
   }
   else
   {
      completedFile << "\r\n" + section;
      certificateFile >> numberOfVariables; // number of variables

      if( certificateFile.fail() || numberOfVariables < 0 )
      {
         cerr << "Invalid number after VAR" << endl;
         returnStatement = false;
      }

      // Store variables
      else
      {
         completedFile << " " << numberOfVariables;

         upperBounds.resize( numberOfVariables );
         lowerBounds.resize( numberOfVariables );
         fill(upperBounds.begin(), upperBounds.end(), make_tuple(infinity, 1, -1));
         fill(lowerBounds.begin(), lowerBounds.end(), make_tuple(-infinity, 1, -1));
         for( int i = 0; i < numberOfVariables; ++i )
         {
            string tmp;
            certificateFile >> tmp;
            if( certificateFile.fail() )
            {
               cerr << "Error reading variable for index " << i << endl;
               returnStatement = false;
         break;
            }
            completedFile <<"\r\n" + tmp;
            if( usesoplex )
            {
               workinglp.addColRational( LPColRational( 1, dummycol, infinity, -infinity ) );
            }
            variableNames.push_back( tmp );
         }
      }
   }

   return returnStatement;
}

// Processes number of integer variables and specifies their indices
// Produces vector of indices of integer variables
// Error if nr of integers invalid or nr of integers > specified integers or section missing
bool processINT()
{
   size_t numberOfIntegers;

   cout << endl << "Processing INT section..." << endl;

   bool returnStatement = false;

   string section;
   certificateFile >> section;

   if( section!= "INT" )
   {
      cerr << "INT expected. Read instead: " << section << endl;
      return returnStatement;
   }

   completedFile << "\r\n" + section;
   certificateFile >> numberOfIntegers;


   if( certificateFile.fail() )
   {
      cerr << "Failed to read number after INT" << endl;
      return returnStatement;
   }

   completedFile << " " <<numberOfIntegers << "\r\n";
   int index;

   for( int i = 0; i < numberOfIntegers; ++i )
   {
      certificateFile >> index;
      if( certificateFile.fail() )
      {
         cerr << "Error reading integer index " << i << endl;
         return returnStatement;
      }
      completedFile << index << " ";
   }
   returnStatement = true;

   return returnStatement;
}

// Processes the sense of the objective function and coefficients for variables
// Stores sense and coefficients for possible future use
// Error if objective sense invalid (other than -1, 0, 1 for min, equality, max) or subroutine fails
bool processOBJ(SoPlex &workinglp)
{
   bool returnStatement = false;
   string section, objectiveSense;
   int numberOfObjCoeff, idx;
   vector<Rational> values;
   vector<int> indices;
   Rational val;
   VectorRational Objective(0); // full objective Vector

   cout << endl << "Processing OBJ section..." << endl;

   certificateFile >> section;

   if( section != "OBJ" )
   {
      cerr << "OBJ expected. Read instead: " << section << endl;
      return returnStatement;
   }
   completedFile << "\r\n" + section;
   certificateFile >> objectiveSense;
   completedFile << " " + objectiveSense;
   if( objectiveSense == "min" )
   {
      workinglp.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
   }
   else if( objectiveSense == "max" )
   {
      workinglp.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
   }
   else
   {
      cerr << "Invalid objective sense: " << objectiveSense << endl;
      return returnStatement;
   }

   Objective.reSize(numberOfVariables);
   Objective.reDim(numberOfVariables);

   certificateFile >> numberOfObjCoeff;
   completedFile << "\r\n" << numberOfObjCoeff << " ";

   values.reserve(numberOfObjCoeff);
   indices.reserve(numberOfObjCoeff);

   for( int j = 0; j < numberOfObjCoeff; ++j)
   {
      certificateFile >> idx >> val;
      completedFile << idx << " " << val << " ";
      values.push_back(val);
      indices.push_back(idx);
   }

   ObjCoeff->add(numberOfObjCoeff, indices.data(), values.data());
   Objective.assign(*ObjCoeff);

   if( usesoplex )
      workinglp.changeObjRational(Objective);

   returnStatement = true;
   return returnStatement;
}

// Processes the constraint section of the file
// Stores constraints for future access
bool processCON(SoPlex &workinglp)
{
   bool returnStatement = false;
   string section;
   string label, consense;
   Rational rhs;
   size_t numberOfConstraints, numberOfBoundedCons, currentDerivation;

   cout << endl << "Processing CON section..." << endl;

   certificateFile >> section;
   if( section!= "CON" )
   {
      cerr << "CON expected. Read instead: " << section << endl;
      return returnStatement;
   }
   completedFile << "\n" + section;

   certificateFile >> numberOfConstraints >> numberOfBoundedCons;
   completedFile << " " << numberOfConstraints << " " << numberOfBoundedCons;

   currentDerivation = 0;

   for( int i = 0; i < numberOfConstraints; ++i)
   {
      certificateFile >> label >> consense >> rhs;
      completedFile << "\r\n" + label + "  " + consense + " " << rhs << "  ";
      currentDerivation++;
      returnStatement = getConstraints(workinglp, consense, rhs, i, currentDerivation);
   }

   return returnStatement;
}

// Processes the relation to prove in order to complete file
bool processRTP()
{
   cout << endl << "Processing RTP section..." << endl;

   bool returnStatement = false;

   string section, relationToProveTypeStr, upperString, lowerString;

   certificateFile >> section;

   if( section != "RTP" )
   {
      cerr << "RTP expected.  Read instead " << section << endl;
   }
   else
   {
      completedFile << "\r\n" + section;
      certificateFile >> relationToProveTypeStr;
      if( relationToProveTypeStr == "infeas")
      {
         completedFile << " " + relationToProveTypeStr;
         returnStatement = true;
      }
      else if( relationToProveTypeStr != "range" )
      {
         cerr << "RTP: unrecognized verification type: " << relationToProveTypeStr << endl;
         return returnStatement;
      }
      else
      {
         completedFile << " " + relationToProveTypeStr;
         certificateFile >> lowerString >> upperString;
         completedFile << " " + lowerString + " " + upperString;
         returnStatement = true;
      }

   }

   return returnStatement;
}

// Processes the solutions to check in order to complete file
bool processSOL()
{
   bool returnStatement = false;
   string section, label;
   int numberOfVarSol = 0;
   int idx;
   Rational val;
   size_t numberOfSolutions;

   cout << endl << "Processing SOL section... " << endl;

   certificateFile >> section;

   if( section != "SOL" )
   {
      cerr << "SOL expected.   Read instead " << section << endl;
      return returnStatement;
   }

   completedFile << "\r\n" + section;
   certificateFile >> numberOfSolutions;

   if( certificateFile.fail() )
   {
      cerr << "Failed to read number after SOL" << endl;
      return returnStatement;
   }
   else if( numberOfSolutions < 0 )
   {
      cerr << "Invalid number after SOL: " << numberOfSolutions << endl;
      return returnStatement;
   }
   else
   {
      completedFile << " " << numberOfSolutions;

      for( int i = 0; i < numberOfSolutions; ++i )
      {
         certificateFile >> label >> numberOfVarSol;
         completedFile << "\r\n" << label << " " << numberOfVarSol;
         for( int i = 0; i < numberOfVarSol; ++i )
         {
            certificateFile >> idx >> val;

            completedFile << "  " << idx << " " << val;
         }
      }

      returnStatement = true;
   }
   return returnStatement;
}

static bool isEqual(DSVectorRational row1, DSVectorRational row2)
{
   if( row1.size() != row2.size() )
      return false;
   else
   {
      for( int i = 0; i < row1.size(); i++ )
      {
         if( row1.index(i) != row2.index(i))
            return false;
         if( row1[i] != row2[i] )
            return false;
      }
   }
   return true;
}

static size_t pushLineToConstraints(string& line, size_t conidx)
{
   stringstream linestream(line);

   string label, numberOfCoefficients, consense;
   int sense, intOfCoefficients = 0;
   long idx;
   Rational val, rhs;

   DSVectorPointer row;
   vector<Rational> values;
   vector<int> indices;
   bool isobjective, global;

   linestream >> label >> consense >> rhs;
   linestream >> numberOfCoefficients;

   if( numberOfCoefficients == "OBJ" )
   {
      isobjective = true;
      row = ObjCoeff;
      intOfCoefficients = row->size();
   }
   else
   {
      isobjective = false;
      intOfCoefficients = atoi(numberOfCoefficients.c_str());
      values.reserve(intOfCoefficients);
      indices.reserve(intOfCoefficients);
      for( int j = 0; j < intOfCoefficients; ++j )
      {
         linestream >> idx >> val;
         values.push_back(val);
         indices.push_back(idx);
      }

      row = make_shared<DSVectorRational>();

      (*row).add(intOfCoefficients, indices.data(), values.data());
   }

   switch(consense[0])
   {
      case 'E':
         sense = 0; break;
      case 'L':
         sense = -1; break;
      case 'G':
         sense = 1; break;
      default:
         cerr << "wrong sense for constraints " << consense << endl;
         break;
   }
   constraints[conidx] = {row, rhs, sense, line};
   return constraints.size() - 1;
}

// read data from certificateFile and buffer it correctly in the stringstream of returnData
static void processSequentialInputFilter(parallelData& returnData, size_t& lineindex, tbb::flow_control& fc)
{
   string line;
   bool stop = false;
   int nbufferedlines = 0;

   stringstream passstream;
   bool lineneedscompletion = false;
   size_t filepos = certificateFile.tellg();

   while( !stop && getline(certificateFile, line) )
   {
      lineneedscompletion = (line.find("weak") != string::npos) || (line.find("incomplete") != string::npos);
      // line needs to be completed -> only one line at a time
      if( lineneedscompletion )
      {
         if( nbufferedlines == 0 )
         {
            returnData.conidx.push_back(lineindex);
            pushLineToConstraints(line, lineindex);
            lineindex++;
            nbufferedlines++;
            passstream << line;
            stop = true;
         }
         else
         {
            // reset to beginning of line
            certificateFile.seekg(filepos);
            stop = true;
            lineneedscompletion = false;
         }
      }
      else
      {
         returnData.conidx.push_back(lineindex);
         pushLineToConstraints(line, lineindex);
         lineindex++;
         passstream << line << std::endl;
         nbufferedlines++;
         if( nbufferedlines >= 10 )
            stop = true;
      }
      filepos = certificateFile.tellg();
   }

   returnData.needscompletion = lineneedscompletion;

   if( nbufferedlines > 0 )
      returnData.linestream << passstream.rdbuf();
   else
      fc.stop();
}

// read data from certificateFile and buffer it correctly in the stringstream of returnData
static void readFileIntoMemory(size_t initialLine, std::vector<size_t>& toCompleteLines)
{
   string line;
   size_t lineindex = initialLine;
   bool lineneedscompletion = false;

   while( getline(certificateFile, line) )
   {
      lineneedscompletion = (line.find("weak") != string::npos) || (line.find("incomplete") != string::npos);
      // line needs to be completed -> only one line at a time
      if( lineneedscompletion )
         toCompleteLines.push_back(lineindex);

      pushLineToConstraints(line, lineindex);
      lineindex++;
   }
}

bool processDER(SoPlex workinglp)
{
   bool returnStatement = false;
   string section;
   size_t numberOfDerivations;
   size_t numberOfConstraints = constraints.size();
   size_t completedLines = 0;
   size_t lineindex;
   timeval start, end;
   std::vector<size_t> toCompleteLines;
   std::vector<string> completedlines;

   cout << endl << "Processing DER section... " << endl;
   certificateFile >> section;

   gettimeofday( &start, 0 );

   if( section!= "DER" )
   {
      cerr << "DER expected.  Read instead " << section << endl;
      return false;
   }

   completedFile << "\r\n" + section;

   certificateFile >> numberOfDerivations;
   completedFile << " " << numberOfDerivations;
   constraints.resize(constraints.size() + numberOfDerivations);
   cout << "Number of Derivations is " << numberOfDerivations << endl;

   if( numberOfDerivations == 0 )
   {
      cout << "Number of derivations = 0. Nothing to complete." << endl;
      return true;
   }

   // skip to next line
   certificateFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

   cout << "Available threads: " << nthreads << endl;

   lineindex = numberOfConstraints;

   readFileIntoMemory(lineindex, toCompleteLines);

   circBuf circQueue(0);

   // generate and queue the LPs that are needed in parallel
   if( usesoplex )
   {
      circQueue.resize(2 * nthreads);
      for(int i = 0; i < 2 * nthreads; ++i)
      {
         passData* queueData = new passData[1];
         //queueData->passLp();
         queueData->passLp.setIntParam(SoPlex::READMODE, SoPlex::READMODE_RATIONAL);
         queueData->passLp.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
         queueData->passLp.setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_RATIONAL);
         queueData->passLp.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
         queueData->passLp.setRealParam(SoPlex::FEASTOL, 0.0);
         queueData->passLp.setRealParam(SoPlex::OPTTOL, 0.0);
         queueData->passLp = workinglp;
         circQueue.enqueue(queueData);
      }
   }

   // reset the line index
   lineindex = 0;

   gettimeofday( &end, 0 );
   cout << endl << "processing file into memory took " << getTimeSecs(start, end)
        << " seconds (Wall Clock)" << endl;

   gettimeofday( &start, 0 );


   // pipeline to complete the derivations
   tbb::parallel_pipeline(nthreads,
      // sequential filter to manage the circular buffer and to ensure output is in correct order
      tbb::make_filter<void, parallelData>( tbb::filter_mode::serial_in_order,
         [&]( tbb::flow_control& fc) {
            parallelData returnData;
            while(lineindex != toCompleteLines.size())
            {
               if( usesoplex )
                  returnData.bufData = circQueue.dequeue();
               returnData.conidx.push_back(toCompleteLines[lineindex]);
               lineindex++;
               return returnData;
            }
            fc.stop();
            return returnData;
         }
      ) &
      // parallel processing and completion of derivations -> passes completed line to last filter
      tbb::make_filter<parallelData, parallelData>( tbb::filter_mode::parallel,
         [&]( parallelData returnData ) {
            returnData.line = completelin(returnData.bufData->passLp, returnData.bufData->LProwCertificateMap, constraints[returnData.conidx[0]]);
            return returnData;
         }
      ) &
      //  manages the queueing of the circular buffer and pushes the completed lines to a vector in the correct order
      tbb::make_filter<parallelData, void>( tbb::filter_mode::serial_in_order,
            [&]( parallelData returnData ) {

            completedlines.push_back(returnData.line);
            if( usesoplex )
               circQueue.enqueue(returnData.bufData);
            }
         )
      );

   if( usesoplex )
   {
      while(!circQueue.isEmpty())
      {
         passData* data = circQueue.dequeue();
         delete[] data;
      }
   }


   gettimeofday( &end, 0 );
   cout << endl << "processing completion pipeline took " << getTimeSecs(start, end)
        << " seconds (Wall Clock)" << endl;

   gettimeofday(&start, 0 );

   // output the lines to the certificate
   int c = 0;
   completedFile << endl;
   for( auto i = numberOfConstraints; i < constraints.size(); ++i )
   {
      if( toCompleteLines.size() == 0 || c >= toCompleteLines.size() || i != toCompleteLines[c] )
         completedFile << constraints[i].line << endl;
      else
      {
         completedFile << completedlines[c] << endl;
         c++;
      }
   }

   gettimeofday( &end, 0 );
   cout << endl << "writing to file took " << getTimeSecs(start, end)
        << " seconds (Wall Clock)" << endl << endl;

   std::cout << "Completed " << toCompleteLines.size() << " out of " << numberOfDerivations << endl;
   return true;
}



// Read Coefficients for any constraint
// Modify main LP
bool getConstraints(SoPlex &workinglp, string &consense, Rational &rhs, int &activeConstraint, size_t currentDerivation)
{
   bool returnStatement = true;
   string numberOfCoefficients, normalizedSense, actualsense;
   int intOfCoefficients = 0, sense;
   DSVectorPointer row;
   vector<Rational> values;
   vector<int> indices;
   long idx, lastrow;
   Rational val, normalizedRhs;

   certificateFile >> numberOfCoefficients;
   completedFile << " " + numberOfCoefficients;


   if( numberOfCoefficients == "OBJ" )
   {
      row = ObjCoeff;
   }
   else
   {
      row = make_shared<DSVectorRational>();
      intOfCoefficients = atoi(numberOfCoefficients.c_str());
      values.reserve(intOfCoefficients);
      indices.reserve(intOfCoefficients);
      for( int j = 0; j < intOfCoefficients; ++j )
      {
         certificateFile >> idx >> val;
         completedFile << " " << idx << " " << val;
         values.push_back(val);
         indices.push_back(idx);
      }
      if( intOfCoefficients == 1)
      {

         // normalize bound constraints
         if( consense == "E")
            normalizedSense = "E";
         else if( consense == "L" )
            if( sign(val) <= 0 )
               normalizedSense = "G";
            else normalizedSense = "L";
         else
            if( sign(val) <= 0 )
               normalizedSense = "L";
            else normalizedSense = "G";

         normalizedRhs = rhs/val;

         // check if normalized bound is an improvement over existing bound

         if( normalizedSense == "E" )
         {
            if( get<0>(upperBounds[idx]) > normalizedRhs )
            {
               upperBounds[idx] = make_tuple(normalizedRhs, Rational(val), currentDerivation -1);
            }
            if( get<0>(lowerBounds[idx]) < normalizedRhs )
            {
               lowerBounds[idx] = make_tuple(normalizedRhs, Rational(val), currentDerivation -1);
            }
         }
         else if( normalizedSense == "L" && get<0>(upperBounds[idx]) > normalizedRhs )
         {
            upperBounds[idx] = make_tuple(normalizedRhs, Rational(val), currentDerivation -1);
         }
         else if( normalizedSense == "G" && get<0>(lowerBounds[idx]) < normalizedRhs )
         {
            lowerBounds[idx] = make_tuple(normalizedRhs, Rational(val), currentDerivation -1);
         }
      }

      (*row).add(intOfCoefficients, indices.data(), values.data());
   }

   /* only populate soplex LP if soplex is actually run */
   if( usesoplex )
   {
      if( consense == "E")
      {
         workinglp.addRowRational( LPRowRational( rhs, *row, rhs) );
         sense = 0;
      }
      else if( consense == "L" )
      {
         workinglp.addRowRational( LPRowRational( -infinity, *row, rhs) );
         sense = -1;
      }
      else if( consense == "G" )
      {
         workinglp.addRowRational( LPRowRational( rhs, *row, infinity) );
         sense = 1;
      }
      else
         returnStatement = false;

      lastrow = workinglp.numRows();
      origConsCertIndex[lastrow-1] = activeConstraint;
   }
   else
   {
      /* code */
      if( consense == "E")
         sense = 0;
      else if( consense == "L" )
         sense = -1;
      else if( consense == "G" )
         sense = 1;
      else
         returnStatement = false;
   }

   constraints.push_back({row, Rational(rhs), sense, ""});

   return returnStatement;
}

string completelin( SoPlex &workinglp,  bimap& LProwCertificateMap, Constraint& constraint)
{
   string& line = constraint.line;
   int consense = constraint.sense;
   Rational& rhs = constraint.side;
   DSVectorPointer row = constraint.vec;

   auto derivationstart = line.find("lin");
   stringstream linestream(line.substr(derivationstart + 3));
   stringstream completedDerivation;
   string numberOfCoefficients;
   vector<long> activeDerivations;

   bool retval = true;

   linestream >> numberOfCoefficients;

   if( numberOfCoefficients == "incomplete" )
   {
      string tmp;

      assert(usesoplex);
      if( !usesoplex )
      {
         cerr << "Soplex support must be enabled to process incomplete constraint type. Rerun with parameter soplex=ON." << endl;
         return "";
      }
      VectorRational newObjective(0);
      newObjective.reSize(workinglp.numColsRational());
      newObjective.reDim(workinglp.numColsRational());
      newObjective = *row;
      workinglp.changeObjRational(newObjective);

      // workinglp.setRealParam(SoPlex::TIMELIMIT, 10.0);

      assert( consense == 0 || consense == -1 || consense == 1);

      if( consense >= 0 )
         workinglp.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
      else
         workinglp.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);

      linestream >> tmp;
      while( tmp != "}" )
      {
         activeDerivations.push_back(stol(tmp));
         linestream >> tmp;
      }
      linestream.seekg(0);
      retval = completeIncomplete( workinglp, LProwCertificateMap, activeDerivations, "", completedDerivation, linestream );
#ifndef NDEBUG
      cout << "Completed derivation: " << line.substr(0,10) << endl;
#endif
   }
   else if( numberOfCoefficients == "weak" )
      retval = completeWeakDomination( *row, consense, rhs, completedDerivation, linestream );
   else
   {
      cerr << "Wrong type of derivation: " << numberOfCoefficients << endl;
   }

   completedDerivation << " -1 ";

   string completedstring = line.substr(0,derivationstart + 3) + completedDerivation.str();
   return completedstring;

}

// Reads multipliers for completing weak domination
static bool readMultipliers( int &sense, SVectorRat &mult, stringstream &baseLine )
{

   int k;
   bool returnStatement = true;

   mult.clear();

   baseLine >> k;

   for( auto j = 0; j < k; ++j )
   {
      Rational a;
      int index;

      baseLine >> index >> a;

      if( a == 0 ) continue; // ignore 0 multiplier

      mult[index] = a;

      if( sense == 0 )
      {
         sense = constraints[index].sense * a.sign();
      }
      else
      {
         int tmp = constraints[index].sense * a.sign();
         if( tmp != 0 && sense != tmp )
         {
            cerr << "Coefficient has wrong sign for index " << index << endl;
            returnStatement = false;
            goto TERMINATE;
         }
      }
   }

TERMINATE:
   return returnStatement;
}


// Reads linear combinations for completing weak domination
static bool readLinComb( int &sense, Rational &rhs, SVectorRat& coefficients, SVectorRat& mult,
                  int currentConstraintIndex, SVectorBool &assumptionList, stringstream &baseLine)
{
   bool returnStatement = true;

   if( !readMultipliers(sense, mult, baseLine) )
   {
      returnStatement = false;
   }
   else
   {
      rhs = 0;
      coefficients.clear();
      assumptionList.clear();

      for( auto it = mult.begin(); it != mult.end(); ++it )
      {
         auto index = it->first;
         auto a = it->second;

         // auto myassumptionList = constraint[index].getassumptionList();

         // for( auto it2 = myassumptionList.begin(); it2 != myassumptionList.end(); ++it2 )
         //    assumptionList[it2->first] = true;

         auto& con = constraints[index];

         auto& c = con.vec;

         for( auto i = 0; i < c->size(); ++i )
         {
            (coefficients)[c->index(i)] += a * (*c)[c->index(i)];
         }

         rhs += a * con.side;
      }
   }

   return returnStatement;
}

// Complete "lin"-type derivations marked "weak"
static bool completeWeakDomination(DSVectorRational &row, int consense, Rational &rhs, stringstream& completedLine,
                                    stringstream &baseLine)
{
   SVectorRat coefDer;
   SVectorRat multDer;
   Rational rhsDer;
   Rational correctedSide;
   SVectorBool asmlist;
   int nbounds;
   bool success;
   tuple<Rational,Rational,long> tup;
   string bracket;
   map<int,tuple<long, Rational>> localLowerBoundsToUse;
   map<int,tuple<long, Rational>> localUpperBoundsToUse;
   Rational boundval;
   Rational boundfactor;
   long long boundindex;

   baseLine >> bracket >> nbounds;
   for( size_t i = 0; i < nbounds; i++ )
   {
      int varIndex;
      long boundIndex;
      Rational val;
      string type;

      baseLine >> type >> varIndex >> boundIndex >> val;
      if( type == "L" )
      {
         localLowerBoundsToUse[varIndex] = make_tuple(boundIndex, val);
      }
      else if( type == "U" )
      {
         localUpperBoundsToUse[varIndex] = make_tuple(boundIndex, val);
      }
      else
      {
         cerr << "type does not match L/U, but is instead " << type << endl;
         abort();
         return false;
      }
   }

   baseLine >> bracket;

   if( !readLinComb(consense, rhsDer, coefDer, multDer, 0, asmlist, baseLine) )
   {
      baseLine.ignore(numeric_limits<streamsize>::max(), '\n');
      return false;
   }

   correctedSide = rhsDer;

   for( auto it = coefDer.begin() ; it != coefDer.end() ; ++it )
   {
      int idx = it->first;
      Rational derivedVal = it->second;
      Rational valToDerive = (row)[idx];

      coefDer[idx] = valToDerive;

      if( derivedVal == valToDerive )
         continue;
      else
      {
         // multiplier for the bound
         Rational boundmult = valToDerive - derivedVal;
         bool islower;

         if(consense == -1)
            islower = (boundmult <= 0);
         // con is >= -> need to use ub for negative, lb for positive boundmult
         else if(consense == 1)
            islower = (boundmult >= 0);
         // cons is == -> this can currently not be handled, would need to be split in two parts
         else if(consense == 0)
         {
            cerr << "  cannot complete weak dominated equality constraints" << endl;
            baseLine.ignore(numeric_limits<streamsize>::max(), '\n');
            return false;
         }

         // get the correct bound index and value
         if( islower && (localLowerBoundsToUse.find(idx) != localLowerBoundsToUse.end()) )
         {
            boundindex = get<0>(localLowerBoundsToUse[idx]);
            boundval = get<1>(localLowerBoundsToUse[idx]);
            boundfactor = 1;
         }
         else if( (!islower) && (localUpperBoundsToUse.find(idx) != localUpperBoundsToUse.end()) )
         {
            boundindex = get<0>(localUpperBoundsToUse[idx]);
            boundval = get<1>(localUpperBoundsToUse[idx]);
            boundfactor = 1;
         }
         else
         {
            tup = islower ? lowerBounds[idx] : upperBounds[idx];
            boundfactor = get<1>(tup);
            boundval = get<0>(tup);
            boundindex = get<2>(tup);
         }

         multDer[boundindex] += boundmult / boundfactor;
         correctedSide += boundmult * boundval;

         if( debugmode == true)
         {
            if( boundmult != 0 )
            {
               cout << "    correcting variable " << variableNames[idx] << " (idx " << idx << ") by " <<
                        boundmult << " (" << static_cast<double>(boundmult) << ") using";
               if( islower )
               {
                  if( localLowerBoundsToUse.count(idx) )
                     cout << " local ";
                  cout << " lower bound ";
               }
               else
               {
                  if( localUpperBoundsToUse.count(idx) )
                     cout << " local ";
                  cout << " upper bound ";
               }
               cout << boundval << " (" << static_cast<double>(boundval) << ")" << endl;
            }
         }
      }
   }
   // now go the other way
   for( int i = 0; i < row.size(); ++i )
   {
      int idx = row.index(i);
      Rational derivedVal = coefDer[idx];
      Rational valToDerive = (row)[idx];

      if( derivedVal == valToDerive )
         continue;
      else
      {
         // multiplier for the bound
         Rational boundmult = valToDerive - derivedVal;
         bool islower;

         // con is <= -> need to use ub for positive, lb for negative boundmult
         if(consense == -1)
            islower = (boundmult <= 0);
         // con is >= -> need to use ub for negative, lb for positive boundmult
         else if(consense == 1)
            islower = (boundmult >= 0);
         // cons is == -> this can currently not be handled, would need to be split in two parts
         else if(consense == 0)
         {
            cerr << "  cannot complete weak dominated equality constraints" << endl;
            baseLine.ignore(numeric_limits<streamsize>::max(), '\n');
            return false;
         }

         // get the correct bound index and value
         if( islower && (localLowerBoundsToUse.find(idx) != localLowerBoundsToUse.end()) )
         {
            boundindex = get<0>(localLowerBoundsToUse[idx]);
            boundval = get<1>(localLowerBoundsToUse[idx]);
            boundfactor = 1;
         }
         else if( (!islower) && (localUpperBoundsToUse.find(idx) != localUpperBoundsToUse.end()) )
         {
            boundindex = get<0>(localUpperBoundsToUse[idx]);
            boundval = get<1>(localUpperBoundsToUse[idx]);
            boundfactor = 1;
         }
         else
         {
            tup = islower ? lowerBounds[idx] : upperBounds[idx];
            boundfactor = get<1>(tup);
            boundval = get<0>(tup);
            boundindex = get<2>(tup);
         }

         multDer[boundindex] += boundmult / boundfactor;
         correctedSide += boundmult * boundval;

         if( debugmode == true )
         {
            if( boundmult != 0 )
            {
               cout << "    correcting variable " << variableNames[idx] << " (idx " << idx << ") by "
                    << boundmult << " (" << static_cast<double>(boundmult) << ") using";
               if( islower )
               {
                  if( localLowerBoundsToUse.count(idx) )
                     cout << " local ";
                  cout << " lower bound ";
               }
               else
               {
                  if( localUpperBoundsToUse.count(idx) )
                     cout << " local ";
                  cout << " upper bound ";
               }
               cout << boundval << " (" << static_cast<double>(boundval) << ")" << endl;
            }
         }
      }
   }

   if( debugmode == true )
   {
      cout.precision(numeric_limits<double>::max_digits10);
      cout << "  exact rhs before correction: " << rhsDer << " ("
           <<    static_cast<double>(rhsDer) << ")" << endl;

      cout << "  exact rhs after  correction: " << correctedSide << " ("
           <<    static_cast<double>(correctedSide) << ")" << endl;

      cout << "  rhs that was printed       : " << rhs << " ("
           <<    static_cast<double>(rhs) << ")" << endl;
   }

   // first case: < and the side is larger, second case: > and side is smaller
   if( (consense == -1 && correctedSide > rhs) || (consense == 1 && correctedSide < rhs) )
   {
      if( row.size() == 0 )
      {
         if( (consense == -1 && correctedSide < 0) || consense == 1 && correctedSide > 0)
            success = true;
         else
         {
            cerr << "invalid claim of infeasibility " << endl;
            success = false;
         }
      }
      else
      {
         cerr.precision(numeric_limits<double>::max_digits10);
         cerr << "Constraint does not dominate original one." << endl << "  Corrected Side is "
            << correctedSide << "(" << static_cast<double>(correctedSide) << ")" << endl
            <<  "  Original rhs is " << rhs << "(" << static_cast<double>(rhs) << ")" << endl;

         cerr << "  difference: " << static_cast<double>(correctedSide - rhs) << endl;
         success = false;
      }
   }
   else
      success = true;

   completedLine << " " << multDer.size();
   for(auto it = multDer.begin(); it != multDer.end(); it++)
   {
      completedLine << " " << it->first << " " << it->second;
   }

   completedLine << " }";

   return success;
}

bool completeIncomplete( SoPlex &localLP,  bimap& LProwCertificateMap, vector<long> &newActiveDerivations, string label, stringstream& completedLine,
                        stringstream &baseLine )
{
   string tmp;
   long numrows, derHierarchy;
   int normalizedSense;
   int lpIndex;
   Rational normalizedRhs;

   DVectorRational dualmultipliers(0);
   DVectorRational reducedcosts(0);

   vector<long>::iterator derIterator;

   vector<long> previousActiveDerivations;
   vector<long> toDeleteDerivations;
   vector<long> toAddDerivations;

   const int lpSize = localLP.numRows();

   // Create array for SoPlex removeRowsRational in order to delete Rows
   int idxShift[lpSize];
   for (size_t i = 0; i < lpSize; i++)
   {
      idxShift[i] = 0;
   }

   // retrieve all active cert indices from map
   for( auto it = LProwCertificateMap.left.begin(); it != LProwCertificateMap.left.end(); it++ )
   {
      previousActiveDerivations.push_back(it->second);
   }

   // sort both vectors
   sort(previousActiveDerivations.begin(), previousActiveDerivations.end(), less<int>());
   sort(newActiveDerivations.begin(), newActiveDerivations.end(), less<int>());

   // compute the set difference
   set_difference( previousActiveDerivations.begin(), previousActiveDerivations.end(),
                        newActiveDerivations.begin(), newActiveDerivations.end(),
                        inserter( toDeleteDerivations, toDeleteDerivations.begin() ) );

   set_difference( newActiveDerivations.begin(), newActiveDerivations.end(),
                     previousActiveDerivations.begin(), previousActiveDerivations.end(),
                     inserter( toAddDerivations, toAddDerivations.begin() ) );

   // Delete rows and reset bounds from LP which are not used in the current completion attempt
   // This is done by setting the corresponding entries in idxShift to -1
   for( derIterator = toDeleteDerivations.begin(); derIterator != toDeleteDerivations.end(); derIterator++)
   {
      lpIndex = LProwCertificateMap.right.at(*derIterator);
      idxShift[lpIndex] = -1;
   }

   // Use RemoveRowsRational from SoPlex
   // gets an array of length numRows with indices -1 indicating that row is to be deleted
   // Returns array with new position of row and -1 if row is deleted
   localLP.removeRowsRational( idxShift );

   // sizeof(idxShift)/sizeof(int) is necessary to avoid usage of additional includes
   for( int i = 0; i < sizeof(idxShift)/sizeof(int); ++i)
   {
      if( i != idxShift[i] )
      {
         if( idxShift[i] == -1 )
         {
            LProwCertificateMap.left.erase(i);
         }
         else
         {
            long updatedCertIndex = LProwCertificateMap.left.at(i);
            LProwCertificateMap.left.erase(i);
            LProwCertificateMap.insert( bimap::value_type( idxShift[i], updatedCertIndex ) );
         }
      }
   }

   // Update LP with new derivations to be used
   for( derIterator = toAddDerivations.begin(); derIterator != toAddDerivations.end(); derIterator++ )
   {
      assert(localLP.numRows() > 0);
      assert(*derIterator < constraints.size());
      auto& missingCon = constraints[*derIterator];
      auto ncurrent = localLP.numRowsRational();
      DSVectorRational row = *missingCon.vec;
      Rational rhs = missingCon.side;
      int consense = missingCon.sense;

      if( consense == 0 )
            localLP.addRowRational( LPRowRational( rhs, row, rhs ) );
      else if( consense == -1 )
            localLP.addRowRational( LPRowRational( -infinity, row, rhs ) );
      else if( consense == 1 )
            localLP.addRowRational( LPRowRational( rhs, row, infinity ) );
      else
         return false;

      LProwCertificateMap.insert( bimap::value_type( localLP.numRows()-1, *derIterator ));
      assert(ncurrent + 1 == localLP.numRowsRational());
   }

   SPxSolver::Status stat;
   numrows = localLP.numRows();
   dualmultipliers.reDim(numrows);
   reducedcosts.reDim(numberOfVariables);

   stat = localLP.optimize();

   if( stat == SPxSolver::OPTIMAL )
   {
      localLP.getDualRational(dualmultipliers);
      localLP.getRedCostRational(reducedcosts);

      printReasoningToLine(dualmultipliers, reducedcosts, completedLine, LProwCertificateMap);
      completedLine<< " }";
   }
   else if( stat == SPxSolver::INFEASIBLE )
   {
      localLP.getDualFarkasRational(dualmultipliers);
      localLP.getRedCostRational(reducedcosts);

      printReasoningToLine(dualmultipliers, reducedcosts, completedLine, LProwCertificateMap);
      completedLine<< "}";

   }
   else
   {
      cerr << "Warning: Completion attempt of Derivation "<< label <<" returned with status " << stat << ".\n";
      cerr << "Skip and continue completion of certificate.\n";
      completedLine << " incomplete";
      for( auto i: newActiveDerivations )
         completedLine << " " << i;

   }
   return true;
}

bool printReasoningToLine(DVectorRational &dualmultipliers, DVectorRational &reducedcosts, stringstream& completedLine, bimap& LProwCertificateMap)
{
   long certIndex;

   completedLine << " " << reducedcosts.dim() + dualmultipliers.dim();

   for( int i= 0; i < reducedcosts.dim(); ++i )
   {
      completedLine << " " << i << " " << reducedcosts[i].str();
   }

   for( int i = 0; i < dualmultipliers.dim(); ++i )
   {
      auto res =  LProwCertificateMap.left.find(i);
      if( res != LProwCertificateMap.left.end() )
         certIndex = res->second;
      else
         certIndex = origConsCertIndex[i];

      completedLine << " " << certIndex << " " << dualmultipliers[i].str();
   }

   return true;
}
