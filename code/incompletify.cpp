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

#include <iostream>
#include <fstream>
#include <random>
#include "soplex.h"
#include "CMakeConfig.hpp"

using namespace std;
using namespace soplex;

ifstream certificateFile;
ofstream incompleteFile;
std::string incomptype = "incomplete";
std::string incompobj = "all";

long numberOfConstraints = 0;
double percentageIncomplete = 100;

void modifyFileName(string &path, const string &newExtension);
bool checkversion(string ver);
bool processVER();
bool processVAR();
bool processINT();
bool processOBJ();
bool processCON();
bool processRTP();
bool processSOL();
bool processDER();

int main(int argc, char *argv[])
{
   int returnStatement = -1;
   std::string modifyName;

   std::cout << "Usage <filename> <percentage> (0-100) <type> (incomplete/weak) <all?> (all/noobj)" << endl;

   certificateFile.open(argv[1]);

   if (certificateFile.fail())
   {
      cerr << "Failed to open file " << argv[1] << endl;
      return returnStatement;
   }

   if(argc > 2)
   {
      percentageIncomplete = atof(argv[2]);
      modifyName = argv[2];
   }

   if(argc > 3)
   {
      incomptype = argv[3];
   }

   if(argc > 4)
   {
      incompobj = argv[4];
   }


   string path = argv[1];
   modifyName.append("_");
   modifyName.append(incomptype);
   modifyName.append("_");
   modifyName.append(incompobj);
   modifyName.append(".vipr");
   modifyFileName(path, modifyName);

   incompleteFile.open(path.c_str(), std::ios::out);

   if (incompleteFile.fail())
   {
      cerr << "Failed to open file " << path << endl;
      return returnStatement;
   }

   double start_cpu_tm = clock();
   if (processVER())
      if (processVAR())
         if (processINT())
            if (processOBJ())
               if (processCON())
                  if (processRTP())
                     if (processSOL())
                        if (processDER())
                        {
                           cout << "Incompletion of File successful!" << endl;

                           returnStatement = 0;
                           double cpu_dur = (clock() - start_cpu_tm) / (double)CLOCKS_PER_SEC;

                           cout << endl
                                << "Completed in " << cpu_dur
                                << " seconds (CPU)" << endl;
                        }
   return returnStatement;
}

void modifyFileName(string &path, const string &newExtension)
{
   string::size_type position = path.find_last_of('.');
   if (position == string::npos)
      position = path.size();
   path.replace(position, newExtension.length(), newExtension);
}

bool checkVersion(string version)
{
   bool returnStatement = false;

   size_t position = version.find(".");
   incompleteFile << " " + version;
   int major = atoi(version.substr(0, position).c_str());
   int minor = atoi(version.substr(position + 1, version.length() - position).c_str());

   cout << "Certificate format version " << major << "." << minor << " mode: ";
#ifndef NDEBUG
   cout << "[debug] [incompletify]"<< endl;
#else
   cout << "[optimized] [incompletify]"<< endl;
#endif

   if ((major == VIPR_VERSION_MAJOR) && (minor <= VIPR_VERSION_MINOR))
   {
      returnStatement = true;
   }
   else
   {
      cerr << "Version unsupported" << endl;
   }

   return returnStatement;
}

bool processVER()
{
   bool returnStatement = false;
   string tmpStr;

   for (;;)
   {
      certificateFile >> tmpStr;
      if (tmpStr == "VER")
      {
         incompleteFile << tmpStr;
         certificateFile >> tmpStr;

         returnStatement = checkVersion(tmpStr);
         break;
      }
      else if (tmpStr == "%")
      {
         getline(certificateFile, tmpStr);
      }
      else
      {
         cerr << endl
              << "Comment or VER expected. Read instead "
              << tmpStr << endl;
         break;
      }
   }

   return returnStatement;
}

bool processVAR()
{
   cout << endl
        << "Processing VAR section ..." << endl;
   auto returnStatement = true;
   string section;
   int numberOfVariables;
   certificateFile >> section;

   if (section != "VAR")
   {
      cerr << "VAR expected.   Read instead " << section << endl;
   }
   else
   {
      incompleteFile << "\r\n" + section;
      certificateFile >> numberOfVariables; // number of variables

      if (certificateFile.fail() || numberOfVariables < 0)
      {
         cerr << "Invalid number after VAR" << endl;
         returnStatement = false;
      }

      // Store variables
      else
      {
         incompleteFile << " " + to_string(numberOfVariables);
         for (int i = 0; i < numberOfVariables; ++i)
         {
            string tmp;
            certificateFile >> tmp;
            if (certificateFile.fail())
            {
               cerr << "Error reading variable for index " << i << endl;
               returnStatement = false;
               break;
            }
            incompleteFile << "\r\n" + tmp;
         }
      }
   }
   return returnStatement;
}

bool processINT()
{
   int numberOfIntegers;
   cout << endl
        << "Processing INT section..." << endl;

   bool returnStatement = false;

   string section;

   certificateFile >> section;
   if (section != "INT")
   {
      cerr << "INT expected. Read instead: " << section << endl;
      return returnStatement;
   }
   incompleteFile << "\r\n" + section;
   certificateFile >> numberOfIntegers;

   if (certificateFile.fail())
   {
      cerr << "Failed to read number after INT" << endl;
      return returnStatement;
   }

   incompleteFile << " " + to_string(numberOfIntegers) + "\r\n";
   int index;

   for (int i = 0; i < numberOfIntegers; ++i)
   {
      certificateFile >> index;
      if (certificateFile.fail())
      {
         cerr << "Error reading integer index " << i << endl;
         return returnStatement;
      }
      incompleteFile << to_string(index) + " ";
   }
   returnStatement = true;

   return returnStatement;
}

bool processOBJ()
{
   cout << endl
        << "Processing OBJ section..." << endl;

   bool returnStatement = false;
   string section, objectiveSense;
   int numberOfObjCoeff, idx;
   Rational val;

   certificateFile >> section;

   if (section != "OBJ")
   {
      cerr << "OBJ expected. Read instead: " << section << endl;
      return returnStatement;
   }
   incompleteFile << "\r\n" + section;
   certificateFile >> objectiveSense;
   incompleteFile << " " + objectiveSense;

   certificateFile >> numberOfObjCoeff;
   incompleteFile << "\r\n" + to_string(numberOfObjCoeff) + " ";

   for (int j = 0; j < numberOfObjCoeff; ++j)
   {
      certificateFile >> idx >> val;
      incompleteFile << idx << " " << val << " ";
   }

   returnStatement = true;
   return returnStatement;
}

bool processCON()
{
   cout << endl
        << "Processing CON section..." << endl;

   bool returnStatement = false;
   int numberOfBoundedCons;
   long numberOfCoefficients;
   int idx;
   Rational val;
   string section;
   string label, consense;
   Rational rhs;

   certificateFile >> section;
   if (section != "CON")
   {
      cerr << "CON expected. Read instead: " << section << endl;
      return returnStatement;
   }
   incompleteFile << "\n" + section;

   certificateFile >> numberOfConstraints >> numberOfBoundedCons;
   incompleteFile << " " + to_string(numberOfConstraints) << " " + to_string(numberOfBoundedCons);

   for (int i = 0; i < numberOfConstraints; ++i)
   {
      // correspondingCertRow[i] = i; // initially each constraint is active on itself
      certificateFile >> label >> consense >> rhs >> numberOfCoefficients;
      incompleteFile << "\r\n"
                     << label + "  " + consense + " " << rhs << "  " << numberOfCoefficients << " ";

      for (int j = 0; j < numberOfCoefficients; ++j)
      {
         certificateFile >> idx >> val;
         incompleteFile << idx << " " << val << " ";
      }
   }
   returnStatement = true;

   return returnStatement;
}

bool processRTP()
{
   cout << endl
        << "Processing RTP section..." << endl;

   bool returnStatement = false;

   string section, relationToProveTypeStr, upperString, lowerString;

   certificateFile >> section;

   if (section != "RTP")
   {
      cerr << "RTP expected.  Read instead " << section << endl;
   }
   else
   {
      incompleteFile << "\r\n" + section;
      certificateFile >> relationToProveTypeStr;
      if (relationToProveTypeStr == "infeas")
      {
         incompleteFile << " " + relationToProveTypeStr;
         returnStatement = true;
      }
      else if (relationToProveTypeStr != "range")
      {
         cerr << "RTP: unrecognized verification type: " << relationToProveTypeStr << endl;
         return returnStatement;
      }
      else
      {
         incompleteFile << " " + relationToProveTypeStr;
         certificateFile >> lowerString >> upperString;
         incompleteFile << " " + lowerString + " " + upperString;
         returnStatement = true;
      }
   }

   return returnStatement;
}

bool processSOL()
{
   cout << endl
        << "Processing SOL section... " << endl;
   bool returnStatement = false;
   int numberOfSolutions;
   string section, label;
   int numberOfVarSol = 0;
   int idx;
   Rational val;

   certificateFile >> section;

   if (section != "SOL")
   {
      cerr << "SOL expected.   Read instead " << section << endl;
      return returnStatement;
   }

   incompleteFile << "\r\n" + section;
   certificateFile >> numberOfSolutions;

   if (certificateFile.fail())
   {
      cerr << "Failed to read number after SOL" << endl;
      return returnStatement;
   }
   else if (numberOfSolutions < 0)
   {
      cerr << "Invalid number after SOL: " << numberOfSolutions << endl;
      return returnStatement;
   }
   else
   {
      incompleteFile << " " + to_string(numberOfSolutions);

      for (int i = 0; i < numberOfSolutions; ++i)
      {
         certificateFile >> label >> numberOfVarSol;
         incompleteFile << "\r\n" + label + " " + to_string(numberOfVarSol);
         for (int i = 0; i < numberOfVarSol; ++i)
         {
            certificateFile >> idx >> val;

            incompleteFile << "  " + to_string(idx) + " " << val;
         }
      }

      returnStatement = true;
   }
   return returnStatement;
}

bool processDER()
{
   // random number generation
   std::random_device rd;  // a seed source for the random number engine
   std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
   std::uniform_int_distribution<> distrib(0, 100);

   bool returnStatement = false;
   string section, numberOfCoefficients, label, consense, bracket, kind;
   int intOfCoefficients = 0;
   long idx, sense, numberOfDerivations;
   long derHir;
   string val, rhs;
   bool incomplete = (incomptype == "incomplete");
   bool incompleteObj = (incompobj == "all");

   std::cout << endl
        << "Processing DER section... " << endl;
   certificateFile >> section;

   if (section != "DER")
   {
      cerr << "DER expected.   Read instead " << section << endl;
      return false;
   }

   incompleteFile << "\r\n" + section;

   certificateFile >> numberOfDerivations;
   incompleteFile << " " + to_string(numberOfDerivations);

   if (numberOfDerivations == 0)
   {
      cout << "Number of derivations = 0. Nothing to complete." << endl;
      return true;
   }

   for (int i = 0; i < numberOfDerivations; ++i)
   {

      certificateFile >> label >> consense >> rhs;
      incompleteFile << "\r\n"
                     << label << " " << consense << " " << rhs;

      // get derived constraints
      certificateFile >> numberOfCoefficients;

      if (numberOfCoefficients == "OBJ")
      {
         incompleteFile << " " << numberOfCoefficients;
      }
      else
      {
         intOfCoefficients = atoi(numberOfCoefficients.c_str());
         incompleteFile << " " << numberOfCoefficients;
         for (long j = 0; j < intOfCoefficients; ++j)
         {
            certificateFile >> idx >> val;
            incompleteFile << " " << idx << " " << val;
         }
      }

      // obtain derivation kind
      certificateFile >> bracket >> kind;

      if (bracket != "{")
      {
         cerr << "Expecting { but read instead " << bracket << endl;
         returnStatement = false;
      }
      else
      {
         if (kind == "asm")
         {

            incompleteFile << " " + bracket + " " + kind;
            string bracket;
            long derHierarchy;
            certificateFile >> bracket;
            if (bracket != "}")
            {
               cerr << "Expecting } but read instead" << bracket << endl;
               returnStatement = false;
            }
            else
            {
               incompleteFile << " " + bracket;
               certificateFile >> derHierarchy;
               incompleteFile << " " << derHierarchy;
               returnStatement = true;
            }
            returnStatement = true;
         }

         else if (kind == "lin")
         {
            long numberOfDer;
            long a_i;
            string todelete;
            bool makeincomplete = distrib(gen) <= percentageIncomplete;

            incompleteFile << " " + bracket + " " + kind;

            if( makeincomplete && (numberOfCoefficients != "OBJ" || incompleteObj) )
            {
               if( incomplete )
                  incompleteFile << " incomplete";
               else
               {
                  incompleteFile << " weak { 0 } ";
                  makeincomplete = false;
               }
            }

            certificateFile >> numberOfDer;
            incompleteFile << " " << numberOfDer;

            if( numberOfCoefficients == "OBJ" && !incompleteObj )
            {
               for (long i = 0; i < numberOfDer; ++i)
               {
                  certificateFile >> a_i >> todelete;
                  incompleteFile << " " << a_i << " " << todelete;
               }
            }
            else
            {
               for (long i = 0; i < numberOfDer; ++i)
               {
                  certificateFile >> a_i >> todelete;
                  if (makeincomplete)
                  {
                     if (a_i >= numberOfConstraints)
                        incompleteFile << " " << a_i;
                  }
                  else
                     incompleteFile << " " << a_i << " " << todelete;
               }
            }
            certificateFile >> bracket >> derHir;
            incompleteFile << " " << bracket << " " << derHir;
            returnStatement = true;
         }
         else if (kind == "rnd")
         {

            incompleteFile << " " + bracket + " " + kind;

            long numberOfCoefficients, idx, dernr;
            Rational val;
            string bracket;

            certificateFile >> numberOfCoefficients;
            incompleteFile << " " << numberOfCoefficients;

            for (int i = 0; i < numberOfCoefficients; ++i)
            {
               certificateFile >> idx >> val;
               incompleteFile << " " + to_string(idx) + " " << val;
            }

            certificateFile >> bracket;
            if (bracket != "}")
            {
               cerr << "Expecting } but read instead" << bracket << endl;
               returnStatement = false;
            }
            else
            {
               incompleteFile << " " + bracket;
               certificateFile >> dernr;
               incompleteFile << " " + to_string(dernr);
               returnStatement = true;
            }
            returnStatement = true;
         }

         else if (kind == "uns")
         {

            incompleteFile << " " + bracket + " " + kind;

            long i1, l1, i2, l2, dernr;
            string bracket;

            certificateFile >> i1 >> l1 >> i2 >> l2;
            incompleteFile << " " + to_string(i1) + " " + to_string(l1) +
                                  " " + to_string(i2) + " " + to_string(l2);

            certificateFile >> bracket;
            if (bracket != "}")
            {
               cerr << "Expecting } but read instead" << bracket << endl;
               returnStatement = false;
            }
            else
            {
               incompleteFile << " " + bracket;
               certificateFile >> dernr;
               incompleteFile << " " + to_string(dernr);
               returnStatement = true;
            }
            returnStatement = true;
         }
         else
         {
            cerr << "Unknown reason. Nothing to complete." << endl;
            returnStatement = false;
         }
      }
   }
   return returnStatement;
}
