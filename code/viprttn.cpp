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
#include <vector>
#include <functional>
#include "CMakeConfig.hpp"

using namespace std;

enum Mark {
   NONE,
   TEMP,
   PERM
};

class Node {

 public:
   vector<int> neededBy;
   vector<int> needs;
   streampos fpos;
   Mark mark = NONE;
   int newIdx = -1;
};

bool firstPass( ifstream &pf, int &numCon, vector<Node> &nodes, streampos &fposDer );
bool writeReorderedDER( ifstream &pf, ofstream &optF, streampos fposDer, int &numCon, vector<Node> &nodes, vector<int> &L );

int main(int argc, char *argv[])
{

   bool stat = false;
   int numCon;

   ifstream pf; // input vipr file
   streampos fposDer = -1;
   ofstream optF; // optimized vipr file

   vector<Node> nodes; // node list for derived constraints


   int rs = -1;
   int farg = 1;

   if( argc != 2 )
   {
      cerr << "Usage: " << argv[0] << " filename\n" << endl;
      return rs;
   }


   pf.open( argv[farg] );

   if( pf.fail() )
   {
      cerr << "Failed to open file " << argv[farg] << endl;
      return rs;
   }


   string optFname = string(argv[farg]) + ".opt";

   optF.open( optFname.c_str());

   if ( optF.fail() )
   {
      cerr << "Failed to open file " << optFname << endl;
      return rs;
   }

   if( !firstPass( pf, numCon, nodes, fposDer ) ) goto TERMINATE;


#ifndef NDEBUG
   for( size_t i = 0; i < nodes.size(); ++i )
   {
      cout << "Node " << i << endl;
      cout << "  fpos = " << nodes[i].fpos << endl;
      cout << "  needed by: ";
      for( auto it : nodes[i].neededBy )
          cout << it << " ";
      cout << endl;
      cout << "  needs: ";
      for( auto it : nodes[i].needs )
          cout << it << " ";
      cout << endl;
   }
#endif

   // Topological sort using DFS and using only the last constraint as root.
   // Nodes not connected to the root are discarded

   {

      vector<int> L;

      std::function<bool (int)> _visit = [&nodes, &_visit, &L] (int n)
      {
          bool rval = false;
          if( nodes[n].mark == PERM) rval = true;
          else if( nodes[n].mark == NONE)
          {
             rval = true;
             nodes[n].mark = TEMP;
             for( auto m : nodes[n].needs )
             {
                rval = _visit(m);
                if( !rval ) break;
             }
             nodes[n].mark = PERM;
             L.push_back( n );
          }
          return rval;
      };

      stat = _visit( nodes.size() - 1 );

      if( stat )
      {
         for( size_t i = 0; i < L.size(); ++i )
         {
             nodes[ L[i] ].newIdx = i;
         }

#ifndef NDEBUG
         cout << "Nodes: " << endl;
         for( size_t i = 0; i < nodes.size(); ++i )
         {
             if( nodes[ i ].newIdx < 0) continue;
             cout << i <<  ": " << nodes[ i ].newIdx << " " << endl;
         }
         cout << endl;
#endif

         if( !writeReorderedDER( pf, optF, fposDer, numCon, nodes, L ) )
         {
            goto TERMINATE;
         }
      }

   }


   stat = true;

TERMINATE:
   if( !stat ) {
      cerr << "Error encountered while processing file" << endl;
   }

   optF.close();
   pf.close();

   return rs;
}


// reads the entire file and construct the digraph for reordering derived
// constraints and outputs the vipr file up to right before DER.
// returns the file position right after numDer.
// returns -1 if an error has occurred.
bool firstPass( ifstream &pf, int &numCon, vector<Node> &nodes, streampos &fposDer )
{
   string section, tmp, label;
   char sense;
   int con1, asm1, con2, asm2; // for reading unsplitting indices
   int numBnd, numSol, numDer, idx, num;
   bool stat = false;

   auto _checkVersion = [](string ver)
   {
      bool rstat = false;
      int major, minor;

      size_t pos = ver.find( "." );

      major = atoi( ver.substr( 0, pos ).c_str() );
      minor = atoi( ver.substr( pos+1, ver.length()-pos ).c_str() );

      if ( (major ==VIPR_VERSION_MAJOR) && (minor <=VIPR_VERSION_MINOR ) )
      {
         rstat = true;
      }
      else
      {
         cerr << "Version " << ver << " unsupported" << endl;
      }

      return rstat;
   };


   auto _processSparseVec = [ &pf ]()
   {
      bool rval = true;

      string input, val;
      int k, index;

      pf >> input;
      if (input != "OBJ")
      {
         k = atoi(input.c_str());

         for( int i = 0; i < k; ++i )
         {
            pf >> index >> val;
            if( pf.fail() )
            {
               cerr << "Failed reading coefficient " << i << endl;
               rval = false;
         break;
            }
         }
      }
      return rval;
   };

   auto _processLinCombSparseVec = [ &pf, &numCon, &nodes ]( int derConIdx )
   {
      bool rval = true;

      string input, val, tmp;
      int k, index;

      pf >> tmp;

      if( tmp == "weak" )
      {
         // forward until we hit linear derivation
         while( tmp != "}" )
         {
            pf >> tmp;
         }
         pf >> k;
         if( pf.fail() )
         {
            cerr << "Failed to read number of coefficients" << endl;
            rval = false;
         }
         else
         {
            for( int i = 0; i < k; ++i )
            {
               pf >> index >> val;
               if( pf.fail() )
               {
                  cerr << "Failed reading coefficient " << i << endl;
                  rval = false;
                  break;
               }
               else if( index >= numCon )
               {
                  index -= numCon;
                  nodes[ derConIdx ].needs.push_back( index );
                  nodes[ index ].neededBy.push_back( derConIdx );
               }
            }
         }
      }
      else if( tmp == "incomplete")
      {
         pf >> tmp;
         while( tmp != "}" )
         {
            index = atoi(tmp.c_str());
            if( index >= numCon )
            {
               index -= numCon;
               nodes[ derConIdx ].needs.push_back( index );
               nodes[ index ].neededBy.push_back( derConIdx );
            }
            pf >> tmp;
         }
         pf.putback( '}' );
      }
      else
      {
         k = atoi(tmp.c_str());
         if( pf.fail() )
         {
            cerr << "Failed to read number of coefficients" << endl;
            rval = false;
         }
         else
         {
            for( int i = 0; i < k; ++i )
            {
               pf >> index >> val;
               if( pf.fail() )
               {
                  cerr << "Failed reading coefficient " << i << endl;
                  rval = false;
            break;
               }
               else if( index >= numCon )
               {
                  index -= numCon;
                  nodes[ derConIdx ].needs.push_back( index );
                  nodes[ index ].neededBy.push_back( derConIdx );
               }
            }
         }
      }
      return rval;
   };


   // Eat up comment lines, if any, until hitting VER
   for(;;)
   {
      pf >> section;

      if( pf.fail() ) goto TERMINATE;

      if (section == "VER")
      {
         pf >> tmp;
         if( _checkVersion( tmp ))
         {
   break;
         }
         else
         {
            goto TERMINATE;
         }
      }
      else if (section == "%")
      {
         getline( pf, tmp );
      }
      else
      {
         cerr << endl << "\% or VER expected. Read instead "
                << section << endl;
         goto TERMINATE;
      }
   }


   pf >> section;
   if( section != "VAR" )
   {
      cerr << "VAR expected. Read instead: " << section << endl;
      goto TERMINATE;
   }

   pf >> num;

   if( pf.fail() ) goto TERMINATE;

   for( int i = 0; i < num; ++i )
   {
      pf >> tmp; // read variable name
   }


   pf >> section;
   if( section != "INT" )
   {
      cerr << "INT expected. Read instead: " << section << endl;
      goto TERMINATE;
   }

   pf >> num;

   if( pf.fail() )
   {
      cerr << "Failed to read number after INT" << endl;
      goto TERMINATE;
   }

   for( int i = 0; i < num; ++i )
   {
      pf >> idx;
      if( pf.fail() ) goto TERMINATE;

   }

   pf >> section;
   if( section != "OBJ" )
   {
      cerr << "OBJ expected. Read instead: " << section << endl;
      goto TERMINATE;
   }

   pf >> tmp;
   if( (tmp != "min") && ( tmp != "max" ) )
   {
      cerr << "Unrecognized string after OBJ: " << tmp << endl;
      goto TERMINATE;
   }

   if( ! _processSparseVec() ) goto TERMINATE;


   pf >> section;
   if( section != "CON" )
   {
      cerr << "CON expected. Read instead: " << section << endl;
      goto TERMINATE;
   }

   pf >> numCon >> numBnd;

   for( int i = 0; i < numCon; ++i )
  {
      pf >> label >> sense >> tmp;

      if( !_processSparseVec() ) goto TERMINATE;

   }


   pf >> section;
   if( section != "RTP" )
   {
      cerr << "RTP expected. Read instead: " << section << endl;
      goto TERMINATE;
   }

   pf >> tmp;

   if( tmp == "range" )
   {
      pf >> tmp; // lower bound

      pf >> tmp; // upper bound
   }
   else if ( tmp != "infeas" )
   {
      cerr << "Unrecognized string after RTP: " << tmp << endl;
      goto TERMINATE;
   }

   pf >> section;
   if( section != "SOL" )
   {
      cerr << "SOL expected. Read instead: " << section << endl;
      goto TERMINATE;
   }

   pf >> numSol;

   if (numSol)
   {
      for( int i = 0; i < numSol; ++i ) {
         pf >> label;
         if( !_processSparseVec() ) goto TERMINATE;

      }
   }


   pf >> section;
   if( section != "DER" )
   {
      cerr << "DER expected. Read instead: " << section << endl;
      goto TERMINATE;
   }

   fposDer = pf.tellg(); // remember where DER begins for second pass to
                         // write derived constraints to file

   pf >> numDer;

   cout << "numCon = " << numCon << endl;
   cout << "numDer = " << numDer << endl;
   cout << "fposDer = " << fposDer << endl;

   nodes.resize( numDer );

   for(auto i = 0; i < numDer; ++i )
   {

      nodes[i].fpos = pf.tellg();

      pf >> label >> sense >> tmp;

      if( pf.fail() )
      {
         cerr << "Error reading " << i << " " << label << endl;
         goto TERMINATE;
      }

      stat = _processSparseVec(); // just eat up the derived constraint
      if( !stat )
      {
         cerr << "Error processing " << label << endl;
         goto TERMINATE;
      }

      pf >> tmp;

      if( pf.fail() )
      {
         cerr << "Error reading reason for " << label << endl;
         goto TERMINATE;
      }

      if( tmp != "{" )
      {
         cerr << "'{' expected.   Reading instead: " << tmp << " in "
              << label << endl;
         goto TERMINATE;
      }
      else
      {
         pf >> tmp;
         if( pf.fail() )
         {
            cerr << "Error reading reason type for " << label << endl;
            goto TERMINATE;
         }
         if( tmp == "asm" || tmp == "sol" )
         {
            pf >> tmp;
            if( tmp != "}")
            {
               cerr << "'}' expected. Read instead: " << tmp << endl;
               goto TERMINATE;
            }
         }
         else if( tmp == "lin" )
         {
            stat = _processLinCombSparseVec( i );
            if( stat )
            {
               pf >> tmp;
               if( tmp != "}")
               {
                  cerr << "'}' expected. Read instead: " << tmp << endl;
                  goto TERMINATE;
               }
            }
         }
         else if( tmp == "rnd" )
         {
            stat = _processLinCombSparseVec( i );
            if( stat )
            {
               pf >> tmp;
               if( tmp != "}")
               {
                  cerr << "'}' expected. Read instead: " << tmp << endl;
                  goto TERMINATE;
               }
            }
         }
         else if( tmp == "uns" )
         {
            pf >> con1 >> asm1 >> con2 >> asm2;
            if( pf.fail() )
            {
               cerr << "Error reading unsplit indices for " << label << endl;
               goto TERMINATE;
            }
            else
            {

               auto _insertArc = [ &nodes ](int tail, int head)
               {
                   nodes[ tail ].neededBy.push_back( head );
                   nodes[ head ].needs.push_back( tail );
               };
               if( con1 >= numCon ) _insertArc( con1 - numCon, i );
               if( con2 >= numCon ) _insertArc( con2 - numCon, i );
               if( asm1 >= numCon ) _insertArc( asm1 - numCon, i );
               if( asm2 >= numCon ) _insertArc( asm2 - numCon, i );

               pf >> tmp;
               if( tmp != "}")
               {
                  cerr << "'}' expected. Read instead: " << tmp << endl;
                  goto TERMINATE;
               }
            }
         }
         else
         {
            cerr << "Unrecognized reason type: " << tmp << endl;
            goto TERMINATE;
         }

      }

      pf >> idx; // read off current max con index and ignore it

   }

   stat = true;


TERMINATE:

   return stat;
}

bool writeReorderedDER( ifstream &pf, ofstream &optF, streampos fposDer, int &numCon, vector<Node> &nodes, vector<int> &L )
{
   string section, tmp, label;
   char sense;
   int con1, asm1, con2, asm2; // for reading unsplitting indices
   bool stat = false;

   auto _processSparseVec = [ &pf, &optF, &numCon, &nodes ]( bool useNewIdx )
   {
      bool rval = true;

      string input, val;
      int k, index;

      pf >> input;
      if (input == "OBJ")
      {
         optF << " OBJ ";
      }
      else
      {
         if(input == "weak")
         {
            // pass over the weak derivation
            optF << " weak ";
            while( input != "}" )
            {
               pf >> input;
               optF << " " << input;
            }

            // now read the linear derivation
            pf >> k;

            optF << " " << k;

            for( int i = 0; i < k; ++i )
            {
               pf >> index >> val;
               if( pf.fail() )
               {
                  cerr << "Failed reading coefficient " << i << endl;
                  rval = false;
            break;
               }
               else
               {
                  // if ( (i+1) % 20 == 0) optF << endl;

                  if( useNewIdx && index >= numCon )
                  {
                     index -= numCon;
                     index = nodes[ index ].newIdx + numCon;
                  }

                  optF << "  " << index << " " << val;
               }
            }
         }
         else if(input == "incomplete")
         {
            optF << " incomplete ";
            while(input != "}")
            {
               pf >> input;
               if(input != "}")
               {
                  index = atoi(input.c_str());
                  if( useNewIdx && index >= numCon )
                  {
                     index -= numCon;
                     index = nodes[ index ].newIdx + numCon;
                  }

                  optF << "  " << index;
               }
            }
            pf.putback( '}' );
         }
         else
         {
            k = atoi(input.c_str());

            optF << " " << k;

            for( int i = 0; i < k; ++i )
            {
               pf >> index >> val;
               if( pf.fail() )
               {
                  cerr << "Failed reading coefficient " << i << endl;
                  rval = false;
            break;
               }
               else
               {
                  // if ( (i+1) % 20 == 0) optF << endl;

                  if( useNewIdx && index >= numCon )
                  {
                     index -= numCon;
                     index = nodes[ index ].newIdx + numCon;
                  }

                  optF << "  " << index << " " << val;
               }
            }
         }
      }
      return rval;
   };


   pf.seekg( 0 );
   // copy up to fposDer
   auto inStream = istreambuf_iterator<char>(pf);
   auto outStream = ostreambuf_iterator<char>(optF);

   for(int i = 0; i < fposDer; ++i)
   {
      *outStream++ = *inStream++;
   }

   // copy_n(istreambuf_iterator<char>(pf), fposDer, ostreambuf_iterator<char>(optF));  This doesn't compile in Ubuntu!!!

   cout << "NDER = " << L.size() << endl;

   optF << " " << L.size() << endl;

   for(auto i : L )
   {

      pf.seekg( nodes[ i ].fpos );

      pf >> label >> sense >> tmp;
      if( pf.fail() )
      {
         cerr << "Error reading " << label << " " << i << endl;
         goto TERMINATE;
      }

      optF << label << " " << sense << " " << tmp;

      stat = _processSparseVec( false );
      if( !stat )
      {
         cerr << "Error processing " << label << endl;
         goto TERMINATE;
      }

      pf >> tmp;
      if( pf.fail() )
      {
         cerr << "Error reading reason for " << label << endl;
         goto TERMINATE;
      }

      optF << " " << tmp;

      if( tmp != "{" )
      {
         cerr << "'{' expected.   Reading instead: " << tmp << " in "
              << label << endl;
         goto TERMINATE;
      }
      else
      {
         pf >> tmp;
         if( pf.fail() )
         {
            cerr << "Error reading reason type for " << label << endl;
            goto TERMINATE;
         }

         optF << " " << tmp;
         if( tmp == "asm" || tmp == "sol" )
         {
            pf >> tmp;
            if( tmp != "}")
            {
               cerr << "'}' expected. Read instead: " << tmp << endl;
               goto TERMINATE;
            }
            optF << " " << tmp;
         }
         else if( tmp == "lin" )
         {
            stat = _processSparseVec( true );
            if( stat )
            {
               pf >> tmp;
               if( tmp != "}")
               {
                  cerr << "'}' expected. Read instead: " << tmp << endl;
                  goto TERMINATE;
               }
               optF << " " << tmp;
            }
         }
         else if( tmp == "rnd" )
         {
            stat = _processSparseVec( true );
            if( stat )
            {
               pf >> tmp;
               if( tmp != "}")
               {
                  cerr << "'}' expected. Read instead: " << tmp << endl;
                  goto TERMINATE;
               }
               optF << " " << tmp;
            }
         }
         else if( tmp == "uns" )
         {
            pf >> con1 >> asm1 >> con2 >> asm2;
            if( pf.fail() )
            {
               cerr << "Error reading unsplit indices for " << label << endl;
               goto TERMINATE;
            }
            else
            {

               auto _printNewIdx = [ &nodes, &optF, &numCon ](int idx)
               {
                   if( idx >= numCon ) idx = nodes[idx - numCon].newIdx + numCon;
                   optF << " " << idx;
               };
               _printNewIdx( con1 );
               _printNewIdx( asm1 );
               _printNewIdx( con2 );
               _printNewIdx( asm2 );

               pf >> tmp;
               if( tmp != "}")
               {
                  cerr << "'}' expected. Read instead: " << tmp << endl;
                  goto TERMINATE;
               }
               optF << " " << tmp;
            }
         }
         else
         {
            cerr << "Unrecognized reason type: " << tmp << endl;
            goto TERMINATE;
         }

      }

      pf >> tmp; // read off current max con index and ignore it

      auto _maxIdx = [ &nodes, &numCon ](int m)
      {
         int maxIdx = -1;
         for( auto n : nodes[m].neededBy )
         {
            int newIdx = nodes[n].newIdx;
            if( newIdx > maxIdx ) maxIdx = newIdx;

         }
         if (maxIdx != -1) maxIdx += numCon;
         return maxIdx;
      };

      optF << " " << _maxIdx( i ) << endl;

   }

   stat = true;


TERMINATE:

   return stat;
}

