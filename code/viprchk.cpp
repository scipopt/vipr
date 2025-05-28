/*
*
*   Copyright (c) 2016 Kevin K. H. Cheung
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
#include <string>
#include <map>
#include <vector>
#include <cassert>
#include <gmpxx.h>
#include <gmp.h>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <memory>
#include "CMakeConfig.hpp"


// Avoid using namespace std to avoid non-obvious complications (ambiguities)
using std::map;
using std::string;
using std::shared_ptr;
using std::vector;
using std::ifstream;
using std::make_shared;
using std::cerr;
using std::endl;
using std::cout;

// Types
typedef map<int, bool> SVectorBool;

// The type of derivation used to derive a constraint
enum DerivationType
{
   ASM,   // assumption
   LIN,   // simple implication
   RND,   // simple implication with integer rounding; i.e. a CG cut
   UNS,   // unsplit operation
   SOL,   // cutoff bound from primal solution
   UNKNOWN
};

// The type of relation to prove
enum RelationToProveType
{
   INFEAS,   // infeasible
   RANGE   // lower bound (-inf if none) and upper bound (inf if none) to be verified
};


// Classes
// Sparse vectors of rational numbers as maps
class SVectorGMP : public map<int, mpq_class>
{
   public:
      void compactify() {
                           auto it = this->begin();
                           while( it != this->end() )
                           {
                              if( it->second == 0 )
                                 it = this->erase(it);
                              else ++it;
                           }
#ifndef NDEBUG
                           cout << "compactified vector with entries ";
                           for( it = this->begin(); it != this->end(); ++it )
                           {
                              cout << it->second << " ";
                              assert(it->second != 0);
                           }
                           cout << "and size " << this->size() << endl;
#endif
                        }
      void canonicalize() {
                           {
                              auto it = this->begin();
                              while( it != this->end() )
                              {
                                 it->second.canonicalize();
                                 ++it;
                              }
                           }
                        }
      bool operator!=(SVectorGMP &other);
      bool operator==(SVectorGMP &other) { return !(*this != other);}
      SVectorGMP operator-(const SVectorGMP &other)
      {
         SVectorGMP returnsvec(*this);
         for(auto it = other.begin(); it != other.end(); ++it)
         {
            returnsvec[it->first] -= it->second;
         }
         return returnsvec;
      }
};

// Constraint format
class Constraint
{
   public:
      Constraint() {}

      Constraint( const string label, const int sense, const mpq_class rhs,
                  shared_ptr<SVectorGMP> coefficients, const bool isAssumptionCon,
                  const SVectorBool assumptionList):

                  _label(label), _sense(sense), _rhs(rhs), _coefficients(coefficients),
                  _isAssumption(isAssumptionCon), _assumptionList(assumptionList)
                  {
                     _coefficients->compactify();
                     _trashed = false;
                     _falsehood = _isFalsehood();
                     _coefsEqualObj = false;
                  }

      void canonicalize() { _coefficients->canonicalize(); }
      void compactify() { _coefficients->compactify(); }
      bool round();

      mpq_class getRhs() const { return _rhs; }
      mpq_class getCoef(const int index) { if( _coefficients->find(index) != _coefficients->end() )
                                              return (*_coefficients)[index];
                                           else return mpq_class(0); }

      shared_ptr<SVectorGMP> coefSVec() const { return _coefficients; }

      int getSense() const { return _sense; }

      bool isAssumption() { return _isAssumption; }

      bool isFalsehood() const { return _falsehood; }
                  // true iff the constraint is a contradiction like 0 >= 1

      bool isTautology();
                  // true iff the constraint is a tautology like 0 <= 1

      bool hasAsm(const int index) {
               return (_assumptionList.find(index) != _assumptionList.end());
            }

      bool hasObjectiveCoefficients() const { return _coefsEqualObj; }
                  // true iff the coefficient vector equals the objective coefficient vector

      void markObjectiveCoefficients() { _coefsEqualObj = true; }
                  // marks the coefficient vector to be equal to the objective coefficient vector

      void setassumptionList(const SVectorBool assumptionList) { _assumptionList = assumptionList; }
      SVectorBool getassumptionList() const { return _assumptionList; }

      bool dominates(Constraint &other) const;
      void print();

      void trash() { _trashed = true; _falsehood = false; _coefficients = nullptr;
                     _rhs = 0; _assumptionList.clear(); }
      bool isTrashed() const { return _trashed; }

      string label() const { return _label; }

      void setMaxRefIdx(int refIdx) { _refIdx = refIdx; }
      int getMaxRefIdx() { return _refIdx; }

      Constraint operator-(const Constraint& other)
      {
         Constraint returncons(*this);
         for(auto it = other._coefficients->begin(); it != other._coefficients->end(); ++it)
         {
            (*returncons._coefficients)[it->first] -= it->second;
         }

         returncons._rhs -= other._rhs;
         return returncons;
      }

   private:
      string _label;
      int _sense;
      mpq_class _rhs;
      shared_ptr<SVectorGMP> _coefficients;
      int _refIdx = -1;
      bool _isAssumption;
      SVectorBool _assumptionList; // constraint index list that are assumptions
      bool _falsehood;
      bool _coefsEqualObj;

      bool _isFalsehood();
      bool _trashed;
};


// Globals
const SVectorBool emptyList;

int numberOfVariables = 0; // number of variables
int numberOfConstraints = 0; // number of constraints
int numberOfBounds = 0; // number of bounds
int numberOfDerivations = 0; // number od derivations
int numberOfSolutions = 0; // number of solutions
vector<bool> isInt; // integer variable indices
vector<string> variable; // variable names
vector<Constraint> constraint; // all the constraints, including derived ones
vector<SVectorGMP> solution; // all the solutions for checking feasibility
ifstream certificateFile;   // certificate file stream

RelationToProveType relationToProveType;
mpq_class bestObjectiveValue; // best objective function value of specified solutions
mpq_class lowerBound; // lower bound for optimal value to be checked
mpq_class upperBound; // upper bound for optimal value to be checked
string lowerStr, upperStr;
bool isMin; // is minimization problem
bool checkLower; // true iff need to verify lower bound
bool checkUpper; // true iff need to verify upper bound
Constraint relationToProve; // constraint to be derived in the case of bound checking
shared_ptr<SVectorGMP> objectiveCoefficients(make_shared<SVectorGMP>()); // obj coefficients
bool objectiveIntegral;


// Forward declaration
bool checkVersion(string ver);
bool processVER();
bool processVAR();
bool processINT();
bool processOBJ();
bool processCON();
bool processRTP();
bool processSOL();
bool processDER();

bool readMultipliers(int &sense, SVectorGMP &mult);
bool readConstraintCoefficients(shared_ptr<SVectorGMP> &v, bool& coefEqualsObj);
bool readConstraint( string &label, int &sense, mpq_class &rhs,
                     shared_ptr<SVectorGMP> &coef, bool& coefEqualsObj);

inline mpq_class floor(const mpq_class &q); // rounding down
inline mpq_class ceil(const mpq_class &q); // rounding up
bool isInteger(const mpq_class &q); // check if variable is integer

mpq_class scalarProduct(shared_ptr<SVectorGMP> u, shared_ptr<SVectorGMP> v);

bool canUnsplit(  Constraint &toDer, const int con1, const int a1, const int con2,
                  const int a2, SVectorBool &assumptionList);

bool readLinComb( int &sense, mpq_class &rhs, shared_ptr<SVectorGMP> coef,
                  int currConIdx, SVectorBool &amsList);

// Main function
int main(int argc, char *argv[])
{

   int returnStatement = -1;

   if( argc != 2 )
   {
      cerr << "Usage: " << argv[0] << " <certificate filename>\n";
      return returnStatement;
   }

   certificateFile.open(argv[1]);

   if( certificateFile.fail() )
   {
      cerr << "Failed to open file " << argv[1] << endl;
      return returnStatement;
   }

   double start_cpu_tm = clock();
   if( processVER() )
      if( processVAR() )
         if( processINT() )
            if( processOBJ() )
               if( processCON() )
                  if( processRTP() )
                     if( processSOL() )
                        if( processDER() ) {
                           returnStatement = 0;
                           double cpu_dur = (clock() - start_cpu_tm)
                                            / (double)CLOCKS_PER_SEC;

                           cout << endl << "Completed in " << cpu_dur
                                << " seconds (CPU)" << endl;
                        }

   return returnStatement;
}


// Processes in order of appearance

// Version control for .vipr input file. Backward compatibility possible for minor versions
bool checkVersion(string version)
{
   bool returnStatement = false;

   size_t position = version.find(".");

   int major = atoi(version.substr(0, position).c_str());
   int minor = atoi(version.substr(position+1, version.length()-position).c_str());

   cout << "Certificate format version " << major << "." << minor << " mode: ";
#ifndef NDEBUG
   cout << "[debug]"<< endl;
#else
   cout << "[optimized]"<< endl;
#endif

   if( (major == VIPR_VERSION_MAJOR) && (minor <= VIPR_VERSION_MINOR) )
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
         certificateFile >> tmpStr;
         returnStatement = checkVersion(tmpStr);
break;
      }
      else if( tmpStr == "%" )
      {
         getline(certificateFile, tmpStr);
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


// Processes number of variables then indexes them in order from 0 to n-1
// Produces vector of variables
// Error if nr of variables invalid or number of variables > specified variables or section missing
bool processVAR()
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
      certificateFile >> numberOfVariables; // number of variables

      if( certificateFile.fail() || numberOfVariables < 0 )
      {
         cerr << "Invalid number after VAR" << endl;
         returnStatement = false;
      }


      // Store variables
      else
      {
         for( int i = 0; i < numberOfVariables; i++ )
         {
            string tmp;
            certificateFile >> tmp;
            if( certificateFile.fail() )
            {
               cerr << "Error reading variable for index " << i << endl;
               returnStatement = false;
         break;
            }
            variable.push_back(tmp);
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
   cout << endl << "Processing INT section..." << endl;

   bool returnStatement = false;

   string section;

   certificateFile >> section;


   // Check section
   if( section != "INT" )
   {
      cerr << "INT expected.   Read instead " << section << endl;
   }
   else
   {
      auto numberOfIntegers = 0;

      certificateFile >> numberOfIntegers; // number of integer variables

      if( certificateFile.fail() || numberOfIntegers < 0 )
      {
         cerr << "Invalid number after INT" << endl;
      }


      // Store integer variables
      else
      {
         isInt.resize(variable.size());

         for( auto it = isInt.begin(); it != isInt.end(); ++it )
         {
            *it = false;
         }

         if( numberOfIntegers > 0 ) {

            for( int i = 0; i < numberOfIntegers; ++i )
            {
               int index;

               certificateFile >> index;
               if( certificateFile.fail() )
               {
                  cerr << "Error reading integer index " << i << endl;
            goto TERMINATE;
               }
               isInt[index] = true;
            }
         }
         returnStatement = true;
      }
   }

TERMINATE:
   return returnStatement;
}


// Processes the sense of the objective function and coefficients for variables
// Stores sense and runs subroutine to store coefficients
// Error if objective sense invalid (other than -1, 0, 1 for min, equality, max) or subroutine fails
bool processOBJ()
{
   cout << endl << "Processing OBJ section..." << endl;

   bool returnStatement = false;

   string section;

   certificateFile >> section;


   // Check section
   if( section != "OBJ" )
   {
      cerr << "OBJ expected.   Read instead " << section << endl;
   }
   else
   {
      string objectiveSense;
      bool dummy;

      certificateFile >> objectiveSense;

      if( objectiveSense == "min" )
      {
          isMin = true;
      }
      else if( objectiveSense == "max" )
      {
          isMin = false;
      }
      else
      {
          cerr << "Invalid objective sense: " << objectiveSense << endl;
          goto TERMINATE;
      }

      returnStatement = readConstraintCoefficients(objectiveCoefficients, dummy);
      assert(!dummy);
      objectiveIntegral = true;

      for( auto it = objectiveCoefficients->begin(); it != objectiveCoefficients->end(); ++it )
      {
         if ( !isInteger(it->second) || !isInt[it->first] )
            objectiveIntegral = false;
      }

      if( !returnStatement )
      {
         cerr << "Failed to read objective coefficients" << endl;
      }
   }

TERMINATE:
   return returnStatement;
}


// Processes constraints
// Produces vector of constraints and calls subroutine
// Error if number of constraints or bounds smaller 0
bool processCON()
{
   cout << endl << "Processing CON section..." << endl;

   bool returnStatement = false;

   string section;

   certificateFile >> section;


   // Check section
   if( section != "CON" )
   {
      cerr << "CON expected.   Read instead " << section << endl;
   }
   else
   {
      certificateFile >> numberOfConstraints >> numberOfBounds;
      // numberOfBounds not used in verification but useful for debugging

      if( certificateFile.fail() || numberOfConstraints < 0 || numberOfBounds < 0 )
      {
         cerr << "Invalid number(s) after CON" << endl;
      }


      // Store constraints
      else
      {
         string label;
         int sense;
         mpq_class rhs;

         for( int i = 0; i < numberOfConstraints; i++ )
         {
            shared_ptr<SVectorGMP> coef(make_shared<SVectorGMP>());
            bool dummy;

            returnStatement = readConstraint(label, sense, rhs, coef, dummy);
            if( label[0] == '%' )
            {
               i--;
//             TODO: this is not working if the line with the comment is not containing a whitespace 
//               if( label[ label.size( ) - 1 ] != '\n' )
               certificateFile.ignore(std::numeric_limits<std::streamsize>::max( ), '\n');
               continue;
            }
            if( !returnStatement ) break;

            constraint.push_back(Constraint(label, sense, rhs, coef, false, emptyList));
         }
      }
   }

   return returnStatement;
}


// Processes the relation to prove - either infeasibility or given range
// Stores type of relation
// Error if invalid verification type or bounds
bool processRTP()
{

   cout << endl << "Processing RTP section..." << endl;

   bool returnStatement = false;

   string section;

   certificateFile >> section;

   while( section[0] == '%' )
      certificateFile >> section;

   // Checking section
   if( section != "RTP" )
   {
      cerr << "RTP expected.   Read instead " << section << endl;
   }
   else
   {
      string relationToProveTypeStr;

      certificateFile >> relationToProveTypeStr;


      // Check verification type
      if( relationToProveTypeStr == "infeas" )
      {
         relationToProveType = RelationToProveType::INFEAS;

         cout << endl << "Need to verify infeasibility. " << endl;
      }
      else if( relationToProveTypeStr != "range" )
      {
         cerr << "RTP: unrecognized verification type: " << relationToProveTypeStr << endl;
         goto TERMINATE;
      }
      else
      {
         relationToProveType = RelationToProveType::RANGE;

         checkLower = checkUpper = false;

         certificateFile >> lowerStr >> upperStr;

         if( lowerStr != "-inf" )
         {
            checkLower = true;
            lowerBound = mpq_class(lowerStr);
         }

         if( upperStr != "inf" )
         {
            checkUpper = true;
            upperBound = mpq_class(upperStr);
         }


         // Check bounds
         if( checkLower && checkUpper && (lowerBound > upperBound) )
         {
            cerr << "RTP: invalid bounds. " << endl;
            goto TERMINATE;
         }


         // Stores relation as Constraint variable
         if( isMin && checkLower )
         {
            relationToProve = Constraint("rtp", 1, lowerBound, objectiveCoefficients, false, emptyList);
         }
         else if( !isMin && checkUpper )
         {
            relationToProve = Constraint("rtp", -1, upperBound, objectiveCoefficients, false, emptyList);
         }
         else
         {
            assert(relationToProve.isTautology());
         }

         cout << "Need to verify optimal value range "
                << (lowerStr == "-inf" ? "(" : "[")
                << lowerStr << ", " << upperStr
                << (upperStr == "inf" ? ")" : "]")
                << "." << endl;

      }
      returnStatement = true;
   }

TERMINATE:
   return returnStatement;

}


// Processes solutions to be verified
// Checks constraints and bounds
// Error if wrong format, type or if bounds violated by solution
bool processSOL()
{
   cout << endl << "Processing SOL section..." << endl;

   bool returnStatement = false;
   mpq_class value;

   string section, label;

   certificateFile >> section;

   // Check format
   if( section != "SOL" )
   {
      cerr << "SOL expected.   Read instead " << section << endl;
      return returnStatement;
   }

   certificateFile >> numberOfSolutions;

   if( certificateFile.fail() )
   {
      cerr << "Failed to read number after SOL" << endl;
   }
   else if( numberOfSolutions < 0 )
   {
      cerr << "Invalid number after SOL: " << numberOfSolutions << endl;
   }
   else
   {
      auto satisfies = [] (Constraint &con, shared_ptr<SVectorGMP> &x)
      {
         bool returnStat = false;

         mpq_class prod = scalarProduct(con.coefSVec(), x);

         if( con.getSense() < 0 )
         {
            returnStat = (prod <= con.getRhs());
         }
         else if( con.getSense() > 0 )
         {
            returnStat = (prod >= con.getRhs());
         }
         else
         {
            returnStat = (con.getRhs() == prod);
         }

         return returnStat;
      };

      shared_ptr<SVectorGMP> solutionSpecified(make_shared<SVectorGMP>());
      vector<mpq_class> sol(numberOfVariables);

      for( int i = 0; i < numberOfSolutions; ++i )
      {
         bool dummy;

         certificateFile >> label;
         cout << "checking solution " << label << endl;

         if( !readConstraintCoefficients(solutionSpecified, dummy) )
         {
            cerr << "Failed to read solution." << endl;
            goto TERMINATE;
         }
         else
         {
            for( int j = 0; j < numberOfVariables; ++j )
            {
              sol[j] = 0;
            }

            // Check integrality constraints
            for( auto it = solutionSpecified->begin(); it != solutionSpecified->end(); ++it )
            {
               if( isInt[it->first] && !isInteger(it->second) )
               {
                  cerr << "Noninteger value for integer variable "
                       << it->first << endl;
                  goto TERMINATE;
               }
              sol[it->first] = it->second;
            }

            for( int j = 0; j < numberOfConstraints; ++j )
            {
               if( !satisfies(constraint[j], solutionSpecified) )
               {
                  cerr << "Constraint " << j << " not satisfied." << endl;
                  goto TERMINATE;
               }
            }
         }

         value = scalarProduct(objectiveCoefficients , solutionSpecified);

         cout << "   objval = " << value << endl;

         // Update best value for the objective function
         if( i )
         {
            if( isMin && value < bestObjectiveValue )
            {
               bestObjectiveValue = value;
            }
            else if( !isMin && value > bestObjectiveValue )
            {
               bestObjectiveValue = value;
            }
         }
         else
         {
            bestObjectiveValue = value;
         }
      }

      if( numberOfSolutions )
      {
         cout << "Best objval: " << bestObjectiveValue << endl;

         // Check if bounds are already violated
         if( isMin && checkUpper && bestObjectiveValue > upperBound )
         {
            cerr << "Best objective values (" << bestObjectiveValue<< ")  exceeds upper bound (" << upperBound << ")." << endl;
            goto TERMINATE;
         }
         else if( !isMin && checkLower && bestObjectiveValue < lowerBound )
         {
            cerr << "Best objective values (" << bestObjectiveValue<< ")  exceeds lower bound (" << lowerBound << ")." << endl;
            goto TERMINATE;
         }
      }

      returnStatement = true;
   }

TERMINATE:
   return returnStatement;
}

// Processes derived constraints
// Checks derivation types and derived constraints
// Finally confirms or rejects Solution and/or relation to prove
// Error if wrong format, derived constraints differ from given
bool processDER()
{

   cout << endl << "Processing DER section..." << endl;

   if( relationToProve.isTautology() )
   {
      cout << "RTP is a tautology: Successfully verified." << endl;
      return true;
   }

   bool returnStatement = false;

   string section;

   certificateFile >> section;

   if( section != "DER" )
   {
      cerr << "DER expected.   Read instead " << section << endl;
      return false;
   }

   certificateFile >> numberOfDerivations;

   cout << "numberOfDerivations = " << numberOfDerivations << endl;


   // No lower bound to check and no deriviations -> nothing to do
   if( numberOfDerivations == 0 && !checkLower )
   {
      cout << "Successfully checked solution for feasibility" << endl;
      return true;
   }

   // no upper bound to check and no deriviations -> nothing to do
   if( numberOfDerivations == 0 && !checkUpper )
   {
      cout << "Successfully checked solution for feasibility" << endl;
      return true;
   }

   string label;
   int sense;
   mpq_class rhs;

   for( int i = 0; i < numberOfDerivations; ++i )
   {
      shared_ptr<SVectorGMP> coef(make_shared<SVectorGMP>());
      bool coefEqualsObj = false;

      if( !readConstraint(label, sense, rhs, coef, coefEqualsObj) )
         return false;

      if( label[0] == '%' )
      {
         i--;
//       TODO: this is not working if the line with the comment is not containing a whitespace
//         if( label[ label.size( ) - 1 ] != '\n' )
         certificateFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
         continue;
      }

      // Obtain derivation method and info
      string bracket, kind;
      int refIdx;

      certificateFile >> bracket >> kind;

      if( bracket != "{" )
      {
         cerr << "Expecting { but read instead " << bracket << endl;
         return false;
      }

      DerivationType derivationType = DerivationType::UNKNOWN;

      if( kind == "asm" )
         derivationType = DerivationType::ASM;
      else if( kind == "sol" )
         derivationType = DerivationType::SOL;
      else if( kind == "lin" )
         derivationType = DerivationType::LIN;
      else if( kind == "rnd" )
         derivationType = DerivationType::RND;
      else if( kind == "uns" )
         derivationType = DerivationType::UNS;

      // The constraint to be derived
      Constraint toDer(label, sense, rhs, coef, (derivationType == DerivationType::ASM), emptyList);
      if( coefEqualsObj )
         toDer.markObjectiveCoefficients();

#ifdef MORE_DEBUG_OUTPUT
      cout << numberOfConstraints + i << " - deriving..." << label << endl;
#endif

      SVectorBool assumptionList;

      int newConIdx = constraint.size();

      switch( derivationType )
      {

         // Assumption, i.e. set of assumptions only contains index of constraint
         case DerivationType::ASM:
            assumptionList[ newConIdx ] = true;
            certificateFile >> bracket;

            if( bracket != "}" )
            {
               cerr << "Syntax Error in " << label << ": Expecting } but read instead " << bracket << endl;
               return false;
            }
            break;
         // Linear combination or rounding
         case DerivationType::LIN:
         case DerivationType::RND:
            {
               shared_ptr<SVectorGMP> coefDer(make_shared<SVectorGMP>());
               mpq_class rhsDer;
               int senseDer;

              if( !readLinComb(senseDer, rhsDer, coefDer, newConIdx, assumptionList) )
                 return false;

               certificateFile >> bracket;

               if( bracket != "}" )
               {
                  cerr << "Syntax Error in " << label << ": Expecting } but read instead " << bracket << endl;
                  return false;
               }

               Constraint derived("", senseDer, rhsDer, coefDer, toDer.isAssumption(),
                                           toDer.getassumptionList());


               if( derivationType == DerivationType::RND )   // round the coefficients
                  if( !derived.round() )
                     return false;


               // check the from reason derived constraint against the given
               // very rarely it happens that the equality check fails because values are not canonicalized or contain
               // zeros, so and we perform that lazily
               if( !derived.dominates(toDer) )
               {
                  derived.canonicalize();
                  toDer.canonicalize();
               }

               if( !derived.dominates(toDer) )
               {
                  cout << "Failed to derive constraint " << label << endl;
                  toDer.print();

                  cout << "Derived instead " << endl;
                  derived.print();

                  cout << "difference: " << endl;
                  (derived - toDer).print();

                  return false;
               }
            }
            break;

            // Unsplit
         case DerivationType::UNS:
            {
               int con1, asm1, con2, asm2;

               certificateFile >> con1 >> asm1 >> con2 >> asm2;

               if( certificateFile.fail() )
               {
                  cerr << "Error reading con1 asm1 con2 asm2" << endl;
                  return false;
               }

               if( (con1 < 0) || (con1 >= newConIdx) )
               {
                  cerr << "con1 out of bounds: " << con1 << endl;
                  return false;
               }

               if( (con2 < 0) || (con2 >= newConIdx) )
               {
                  cerr << "con2 out of bounds: " << con2 << endl;
                  return false;
               }

               if( !canUnsplit(toDer, con1, asm1, con2, asm2, assumptionList) )
               {
                  cerr << label << ": unsplit failed" << endl;
                  return false;
               }

               certificateFile >> bracket;
               if( bracket != "}" )
               {
                  cerr << "Expecting } but read instead " << bracket << endl;
                  return false;
               }
            }
            break;
         case DerivationType::SOL:
         {
            mpq_class cutoffbound = bestObjectiveValue;
            if (objectiveIntegral)
            {
               cutoffbound -= 1;
            }
            certificateFile >> bracket;
            if (coef != objectiveCoefficients)
            {
               cerr << "Cutoff bound can only be applied to objective value " << endl;
               return false;
            }
            else if (sense != -1)
            {
               cerr << "Cutoff bound should have sense 'L'" << endl;
               return false;
            }
            else if (rhs < cutoffbound )
            {
               cerr << "No solution known with objective at most " << rhs << ", best solution is " << bestObjectiveValue << endl;
               return false;
            }
            else if( bracket != "}" )
            {
               cerr << "Expecting } but read instead " << bracket << endl;
               return false;
            }
            break;
         }
         default:
            cout << label << ": unknown derivation type " << kind << endl;
            return false;
            break;
      }

      // Set the list of assumptions
      toDer.setassumptionList(assumptionList);

      // Constraint hierarchy handling (??)
      certificateFile >> refIdx;
      toDer.setMaxRefIdx(refIdx);
      constraint.push_back(toDer);

      // Check whether we have globally proven the dual side of RTP
      if( assumptionList == emptyList )
      {
         if( relationToProveType == RelationToProveType::INFEAS && constraint.back().isFalsehood() )
         {
            cout << "Infeasibility verified." << endl;
            return true;
         }
         else if( relationToProveType == RelationToProveType::RANGE && constraint.back().hasObjectiveCoefficients() && constraint.back().dominates(relationToProve) )
         {
            cout << endl << "Terminated after " << i << " derivations." << endl;

            if( numberOfSolutions ) {
               cout << "Best objval over all solutions: " << bestObjectiveValue << endl;
            }

            cout << "Successfully verified optimal value range "
                 << (lowerStr == "-inf" ? "(" : "[")
                 << lowerStr << ", " << upperStr
                 << (upperStr == "inf" ? ")" : "]")
                 << "." << endl;

            return true;
         }
      }

      // skip to next line
      certificateFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      if( i < numberOfDerivations - 1 ) // Never trash last constraint
         if( (refIdx >= 0) && (refIdx < int(constraint.size())) )
         {
             constraint.back().trash();
         }

#ifdef MORE_DEBUG_OUTPUT
      toDer.print();
#endif
   }

   cout << endl;

   auto assumptionList = constraint.back().getassumptionList();

   // Print potential reasons for errors
   if( assumptionList != emptyList )
   {
      cout << "Final derived constraint contains undischarged assumptions:" << endl;
      for(auto & it : assumptionList)
         cout << it.first << ": " << constraint[ it.first ].label() << endl;
   }
   else
   {
      if( relationToProveType == RelationToProveType::INFEAS )
         cout << "Failed to verify infeasibility." << endl;
      else if( isMin && checkLower )
         cout << "Failed to derive lower bound." << endl;
      else if( !isMin && checkUpper )
            cout << "Failed to derive upper bound." << endl;

      cout << "Proved: " << endl;
      constraint.back().print();
      cout << "Instead of: " << endl;
      relationToProve.print();
   }

   return false;
}


// Classes and Functions
inline mpq_class floor(const mpq_class &q)
{
   mpz_t z;
   mpq_class result;
   mpz_init (z);
   mpz_fdiv_q(z, q.get_num_mpz_t(), q.get_den_mpz_t()); // Divide numerator by denominator and floor the result
   result = mpz_class(z);
   mpz_clear (z);
   return result;
}


inline mpq_class ceil(const mpq_class &q)
{
   mpz_t z;
   mpq_class result;
   mpz_init (z);
   mpz_cdiv_q(z, q.get_num_mpz_t(), q.get_den_mpz_t()); // Divide numerator by denominator and ceil the result
   result = mpz_class(z);
   mpz_clear (z);
   return result;
}

bool isInteger(const mpq_class &q)
{
   return q == floor(q);
}


bool readLinComb( int &sense, mpq_class &rhs, shared_ptr<SVectorGMP> coefficients,
                  int currentConstraintIndex,SVectorBool &assumptionList)
{
   bool returnStatement = true;

   SVectorGMP mult;

#ifndef NDEBUG
   std::cout << "reading linear combination" << std::endl;
#endif

   if( !readMultipliers(sense, mult) )
   {
      returnStatement = false;
   }
   else
   {
      rhs = 0;
      coefficients->clear();
      assumptionList.clear();
      mpq_class t;

      for(auto & it : mult)
      {
         auto index = it.first;
         auto a = it.second;

         auto myassumptionList = constraint[index].getassumptionList();

         for(auto & it2 : myassumptionList)
            assumptionList[it2.first] = true;

         const Constraint &con = constraint[index];

         if( con.isTrashed() )
         {
            cerr << "Accessing trashed constraint: " << con.label() << endl;
            returnStatement = false;
         }
         else
         {
            shared_ptr<SVectorGMP> c = constraint[index].coefSVec();

            for(auto & itr : *c)
               (*coefficients)[ itr.first ] += a * itr.second;

            rhs += a * constraint[index].getRhs();

            if( (constraint[index].getMaxRefIdx() <= currentConstraintIndex) &&
                (constraint[index].getMaxRefIdx() >= 0) )
               constraint[index].trash();
         }
      }
   }

   return returnStatement;
}


bool readMultipliers(int &sense, SVectorGMP &mult)
{

   int k;
   bool returnStatement = true;

   sense = 0;
   mult.clear();

   certificateFile >> k;

   for( auto j = 0; j < k; ++j )
   {
      mpq_class a;
      int index;

      certificateFile >> index >> a;

      if( index < 0 )
      {
         cerr << "Index is out of bounds " << index << endl;
         returnStatement = false;
         goto TERMINATE;
      }

      if( a == 0 ) continue; // ignore 0 multiplier

      mult[index] = a;

#ifndef NDEBUG
      if( index < 0 || index >= constraint.size( ) )
      {
         cerr << "Index out of range " << index << " (0," << constraint.size() << ")" << endl;
         returnStatement = false;
         goto TERMINATE;
      }
      std::cout << "multiplying " << a << " * " << constraint[index].label() << "( ";
      constraint[index].print();
      std::cout << " ) " << std::endl;
#endif

      if( sense == 0 )
      {
         sense = constraint[index].getSense() * sgn(a);
      }
      else
      {
         int tmp = constraint[index].getSense() * sgn(a);
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

// Read and store constraints
bool readConstraintCoefficients(shared_ptr<SVectorGMP> &coefficients, bool& coefEqualsObj)
{
   auto returnStatement = false;
   int k = 0;
   string tmp;
   coefEqualsObj = false;

   coefficients->clear();
   certificateFile >> tmp;

   if( tmp == "OBJ" ) // case that constraint = objective function
   {
      coefficients = objectiveCoefficients;
      coefEqualsObj = true;
      returnStatement = true;
   }
   else
   {
      k = atoi(tmp.c_str());

      if( certificateFile.fail() )
      {
         cerr << "Error reading number of elements " << endl;
         goto TERMINATE;
      }
      else
      {
         for( int j = 0; j < k; j++ )
         {
            int index;
            mpq_class a;

            certificateFile >> index >> a;
            if( certificateFile.fail() )
            {
               cerr << "Error reading integer-rational pair " << endl;
               goto TERMINATE;
            }
            else if( index < 0 || index >= numberOfVariables )
            {
               cerr << "Index out of bounds: " << index << endl;
               goto TERMINATE;
            }
            (*coefficients)[index] = a;
         }
         returnStatement = true;
      }
   }

   coefficients->compactify();

TERMINATE:
   return returnStatement;
}


bool readConstraint(string &label, int &sense, mpq_class &rhs, shared_ptr<SVectorGMP> &coefficients, bool& coefEqualsObj)
{

   auto returnStatement = false;
   char senseChar;

   coefEqualsObj = false;

   certificateFile >> label >> senseChar;

   if( !certificateFile.fail() )
   {
      if( senseChar == 'E' )
         sense = 0;
      else if( senseChar == 'L' )
         sense = -1;
      else if( senseChar == 'G' )
         sense = 1;
      else if( label[0] == '%' )
         return true;
      else
      {
        cerr << "Unknown sense for " << label << ": " << senseChar << endl;
        goto TERMINATE;
      }

      certificateFile >> rhs;

      if( !certificateFile.fail() )
         returnStatement = readConstraintCoefficients(coefficients, coefEqualsObj);

      if( !returnStatement ) cerr << label <<   ": Error reading constraint " << endl;
   }

TERMINATE:
   return returnStatement;
}


// con1 and con2 must be inequalities for an integer disjunction
// e.g. mx <= d and mx >= d+1 such that the variables indexed by
// the support of m are integers.   The function checks this.
// a1 and a2 are assumptions.
bool canUnsplit(  Constraint &toDer, const int con1, const int a1,
                  const int con2, const int a2, SVectorBool &assumptionList)
{

   bool returnStatement = false;

   const Constraint &c1 = constraint[con1];
   const Constraint &c2 = constraint[con2];

   const Constraint &branchAsm1 = constraint[a1];
   const Constraint &branchAsm2 = constraint[a2];

   if( c1.isTrashed() )
   {
      cerr << "unsplitting trashed constraint: " << c1.label() << endl;
      goto TERMINATE;
   }
   else if( c2.isTrashed() )
   {
      cerr << "unsplitting trashed constraint: " << c2.label() << endl;
      goto TERMINATE;
   }

   if( c1.dominates(toDer) && c2.dominates(toDer) )
   {
      SVectorGMP asm1Coef, asm2Coef;
      mpq_class asm1Rhs, asm2Rhs;

      SVectorBool asm1 = c1.getassumptionList();
      SVectorBool asm2 = c2.getassumptionList();

      // remove the indices involved in unsplitting
#ifdef MORE_DEBUG_OUTPUT
      if (asm1.find(a1) == asm1.end())
         cout << "Warning: " << a1 << " not present in unsplit" << endl;
      if (asm2.find(a2) == asm2.end())
         cout << "Warning: " << a2 << " not present in unsplit" << endl;
#endif

      asm1.erase(a1);
      asm2.erase(a2);

      assumptionList.clear();
      assumptionList = asm1;

#ifdef MORE_DEBUG_OUTPUT
      cout << "asm1: ";
      for( auto it = asm1.begin(); it != asm1.end(); ++it ) {
         cout << it->first << " ";
      }
      cout << endl;

      cout << "asm2: ";
      for( auto it = asm2.begin(); it != asm2.end(); ++it ) {
         cout << it->first << " ";
      }
      cout << endl;
#endif

      for(auto & it : asm2) {
         assumptionList[it.first] = true;
      }

      if( branchAsm1.isTrashed() )
      {
         cerr << "accessing trashed constraint: " << branchAsm1.label() << endl;
         goto TERMINATE;
      }
      else if( c2.isTrashed() )
      {
         cerr << "accessing trashed constraint: " << c2.label() << endl;
         goto TERMINATE;
      }


      // the constraints must have opposite senses
      if( -1 != branchAsm1.getSense() * branchAsm2.getSense() )
      {
         cerr << "canUnsplit: Failed sense requirement for assumptions" << endl;
         cerr << "branchAsm1 sense:: " << branchAsm1.getSense() << endl;
         cerr << "branchAsm2 sense:: " << branchAsm2.getSense() << endl;
         goto TERMINATE;
      }
      else
      {

         // check if disjunction gives a tautology with respect to the variable
         // integrality requirements
         bool stat = true;
         if( branchAsm1.getSense() < 0 )
            stat = ((branchAsm1.getRhs() + 1) == branchAsm2.getRhs());
         else // must be > 0
            stat = (branchAsm1.getRhs() == (branchAsm2.getRhs() + 1));

         if( !stat ) {
            cerr << branchAsm1.label() << " and " << branchAsm2.label()
                 << " do not form a tautology" << endl;
            goto TERMINATE;
         };

         shared_ptr<SVectorGMP> c1ptr = branchAsm1.coefSVec();
         shared_ptr<SVectorGMP> c2ptr = branchAsm2.coefSVec();
         if( (c1ptr == c2ptr) || (*c1ptr == *c2ptr) ) // coefSVec can both point to objectiveCoefficients
         {

            for(auto & it : *c1ptr)
            {
               if( !isInt[it.first] )
               {
                  cerr << "canUnsplit: noninteger variable index " << it.first
                         << endl;
                  goto TERMINATE;
               }
               else if( !isInteger(it.second) )
               {
                  cerr << "canUnsplit: noninteger coefficient for index "
                         << it.first << endl;
                  goto TERMINATE;
               }
            }
         }
         else
         {
            cerr << "canUnsplit: coefs of asm constraints differ" << endl;
            goto TERMINATE;
         }
         returnStatement = true;
      }
   }

TERMINATE:
   return returnStatement;
}


// SVectorGMP methods
bool SVectorGMP::operator!=(SVectorGMP &other)
{
   SVectorGMP &cf1 = *this;
   SVectorGMP &cf2 = other;

   // get rid of all zero entries
   if( cf1.size() != cf2.size() )
   {
      cf1.compactify();
      cf2.compactify();
   }

   if( cf1.size() != cf2.size() )
   {
#ifndef NDEBUG
      cout << "vectors found different due to size: " << cf1.size() << " vs. " << cf2.size() << endl;
#endif
      return true;
   }
   else
   {
      for(auto & it1 : cf1)
      {
         auto it2 = cf2.find(it1.first);
         if( ( it2 == cf2.end() && it1.second != 0 ) )
         {
#ifndef NDEBUG
            cout << "vectors found different in coefficient comparison: "
                 << " variable " << it1.first << " not found" << endl;
#endif
            return true;
         }
         else if( it2 != cf2.end() && it2->second != it1.second )
         {
#ifndef NDEBUG
            cout << "vectors found different in coefficient comparison: "
                 << " variable " << it1.first << " has coefficients "
                 << it1.second << " vs. " << it2->second << endl;
#endif
            return true;
         }
      }
   }

   return false;
}


// use non-sparse vector for v to reduce lookup time
mpq_class scalarProduct(shared_ptr<SVectorGMP> u, shared_ptr<SVectorGMP> v)
{
   mpq_class product = 0;

   for( auto it = u->begin(); it != u->end(); ++it )
   {
      auto it2 = v->find(it->first);
      if( it2 != v->end() )
         product += it->second * it2->second;
   }

   return product;
}


// Constraint methods
bool Constraint::round()
{
   bool returnStatement = true;

   for(auto & it : *_coefficients)
   {
      auto j = it.first;
      auto a = it.second;

      if( isInt[j] )
      {   // needs to be an integer variable
         if( !isInteger(a) )
         {
            cerr << "Coefficient of integer variable with index "
                 << j << " is not an integer" << endl;
            returnStatement = false;
            goto TERMINATE;
         }
      }
   }

   if( getSense() < 0 ) // round down
      _rhs = floor(_rhs);
   else if( getSense() > 0 ) // round up
      _rhs = ceil(_rhs);

TERMINATE:
   return returnStatement;
}


bool Constraint::_isFalsehood()
{
   bool returnStatement = false;

   if( _coefficients->empty() )
   {
      if( ((getSense() <= 0) && (_rhs < 0)) || ((getSense() >= 0) && (_rhs > 0)) )
         returnStatement = true;
   }

   return returnStatement;
}


bool Constraint::dominates(Constraint &other) const
{
   bool returnStatement = false;

   if( this->isFalsehood() )
   {
      returnStatement = true;
   }
   else if( *(this->_coefficients) == *(other._coefficients) )  // force object comparison
   {
      if( (other.getSense() > 0 && this->getSense() >= 0 &&
             this->_rhs >= other._rhs)
         || (other.getSense() < 0 && this->getSense() <= 0 &&
             this->_rhs <= other._rhs)
         || (other.getSense() == 0 && this->getSense() == 0 &&
             this->_rhs == other._rhs) )
      {
         returnStatement = true;
      }
   }

   return returnStatement;
}


bool Constraint::isTautology() {
   bool returnStatement = false;

   if( _coefficients->empty() )
   {
      if( ((getSense() == 0) && (0 == _rhs))
         || ((getSense() < 0) && (_rhs >= 0))
         || ((getSense() > 0) && (_rhs <= 0)) )
      {
         returnStatement = true;
      }
   }
   return returnStatement;
}


void Constraint::print() {
   bool first = true;
   mpq_class myCoefficient;
   cout.precision(std::numeric_limits<double>::max_digits10);

   int count = 0;

   if( _isAssumption )
      cout << "Is assumption: ";

   for(auto & it : *_coefficients)
   {
      auto index = it.first;
      auto a = it.second;

      myCoefficient = abs(a);

      if( a > 0 )
      {
         if( !first ) cout << " + ";
         cout << (myCoefficient == 1 ? string("") : string(myCoefficient.get_str()) + " ")
              << "( " << a.get_d() << " ) "
              << variable[index];
         ++count;
         first = false;
      }
      else if( a < 0 )
      {
         cout << " - "
                << (myCoefficient == 1 ? string("") : string(myCoefficient.get_str()) + " ")
                << "( " << a.get_d() << " ) "
                << variable[index];
         first = false;
         ++count;
      }

      if( (count+1) % 4 == 0 ) cout << endl;
   }

   if( first ) // coefficients are all zero
   {
      cout << "0";
   }

   switch( _sense )
   {
      case -1: cout << " <= "; break;
      case 1: cout << " >= "; break;
      case 0: cout << " = "; break;
      default : assert(false);
   }
   cout << _rhs << " ( " << _rhs.get_d() << " )" << endl;

#ifndef NDEBUG
   if( !_isAssumption && !_assumptionList.empty() )
   {
      cout << " -- assumptions: " << endl;
      for( auto it = _assumptionList.begin(); it != _assumptionList.end(); ++it )
      {
         auto index = it->first;
         cout << "   "<< it->first << ": " << constraint[index].label() << endl;
      }
      cout << endl;
   }
#endif
}
