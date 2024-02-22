/*
*
*   Copyright (c) 2016 Kevin K. H. Cheung
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
#include <set>

#define VERSION_MAJOR 1
#define VERSION_MINOR 1 

using namespace std;



int main(int argc, char *argv[])
{

   int num, idx, numCon, numSol, numDer, numBnd;
   int maxWidth = 0;

   ifstream pf; // input vipr file

   vector<int> width;


   int rs = -1;
   int farg = 1;

   if( (argc == 1) || (argc > 2) )
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
   else
   {

      auto _checkVersion = [](string ver)
      {
         bool rstat = false;
         int major, minor;
      
         size_t pos = ver.find( "." );
      
         major = atoi( ver.substr( 0, pos ).c_str() );
         minor = atoi( ver.substr( pos+1, ver.length()-pos ).c_str() );
      
         if ( (major == VERSION_MAJOR) && (minor <= VERSION_MINOR ) )
         {
            rstat = true;
         }
         else
         {
            cerr << "Version " << ver << " unsupported" << endl;
         }
      
         return rstat;
      };

      auto _updateWidth = [ &width ] (int beginIdx, int endIndex)
      {
         for( auto j = beginIdx ; j < endIndex; ++j )
         {
             ++( width[ j ] );
         }
      };

      auto _processSparseVec = [ &pf, &numCon, &_updateWidth ]( bool updateWidth, int curDerIdx )
      {
         bool rval = true;

         string input, val;
         int k, index;

         pf >> input;
         if( !(input == "OBJ") )
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
               else if( updateWidth )
               {
                  if( index >= numCon)
                  {
                     _updateWidth( index - numCon, curDerIdx );
 
                  }
               }
            }
         }
         return rval;
      };


      string section, tmp, label;
      char sense;
      int con1, asm1, con2, asm2; // for reading unsplitting indices
      bool stat = false;

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


      stat = _processSparseVec( false, 0 );

      if( !stat ) goto TERMINATE;


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

         stat = _processSparseVec( false, i );
         if( !stat ) break;

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
            stat = _processSparseVec( false, 0 );
         }
      }


      pf >> section;
      if( section != "DER" )
      {
         cerr << "DER expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> numDer;

      width.resize( numDer );

      // Go through derived constraints and set max constraint indices
      for( int i = 0; i < numDer; ++i )
      {
         pf >> label >> sense >> tmp;

         stat = _processSparseVec( false, i );
         if( !stat ) 
         {
            cerr << "Error processing " << label << endl;
            goto TERMINATE;
         }

         pf >> tmp;
         if( tmp != "{" )
         {
            cerr << "'{' expected.   Reading instead: " << tmp << " in " 
                 << label << endl;
            goto TERMINATE;
         }
         else
         {
            pf >> tmp;
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
               stat = _processSparseVec( true, i );
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
               stat = _processSparseVec( true, i );
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
                  goto TERMINATE;
               }
               else
               {
                  if( con1 >= numCon ) _updateWidth( con1 - numCon, i );
                  if( con2 >= numCon ) _updateWidth( con2 - numCon, i );
                  if( asm1 >= numCon ) _updateWidth( asm1 - numCon, i );
                  if( asm2 >= numCon ) _updateWidth( asm2 - numCon, i );

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

      idx = 0;
      for( auto it : width )
      {
         cout << "width[" << idx++ << "] = " << it << endl;
         if( it > maxWidth ) maxWidth = it;
      }

      cout << "Cutwidth on derived constraints is: " << maxWidth << endl;


TERMINATE:
      if( !stat ) {
         cerr << "Error encountered while processing file" << endl;
      }
   }

   pf.close();

   return rs;
}
