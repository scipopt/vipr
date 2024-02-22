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

#define VERSION_MAJOR 1
#define VERSION_MINOR 1

using namespace std;

ifstream pf;
ofstream html;
vector<string> colName;
vector<string> rowName;
vector<bool> isInt; // boolean indicators for integer variable indices



int main(int argc, char *argv[])
{
   int rs = -1;

   if( argc != 2 )
   {
      cerr << "Usage: " << argv[0] << " filename\n";
      return rs;
   }

   pf.open( argv[1] );

   if( pf.fail() )
   {
      cerr << "Failed to open file " << argv[1] << endl;
      return rs;
   }

   string htmlFname = string(argv[1]) + ".html";

   html.open( htmlFname.c_str());

   if ( html.fail() )
   {
      cerr << "Failed to open file " << htmlFname << endl;
   }
   else
   {

      auto _checkVersion = [](string ver) {
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

      auto _processSparseVec = [](vector<string> &name, bool isSol)
      {
         bool rval = true;

         string input, val;
         int k, index;

         pf >> input;
         if (input == "OBJ")
         {
            html << " OBJ ";
         }
         else
         {
            k = atoi(input.c_str());                

            if( k == 0 )
            {
               html << "0";
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
                  else
                  {
                     if( isSol ) // print sparse vec as a solution
                     {
                        if( i ) html << ", ";
                        html << name[index] << " = " << val;
                     }
                     else
                     {
                        if( val[0] == '-')
                        {
                           html << " - ";
                           if (val != "-1") 
                              html << val.substr(1, val.length()-1);
                           html << " " << name[index];
                        }
                        else
                        {
                           if( i ) html << " + ";
                           if( val != "1" ) html << val;
                           html << " " << name[index];
                        }
                     }
                  }
               }
            }
         }
         return rval;
      };

      string section, tmp, label;
      char sense;
      int num, idx, numCon, numSol, numDer, numBnd;
      int con1, asm1, con2, asm2; // for reading unsplitting indices
      bool stat = false;

      html << "<HTML>" << endl;
      html << "<HEAD>" << endl;
      html << "<STYLE>" << endl;
      html << "TABLE, TH, TD {" << endl;
      html << "   border-collapse: collapse;" << endl;
      html << "   border: 1px solid black;" << endl;
      html << "}" << endl;
      html << "</STYLE>" << endl;
      html << "</HEAD>" << endl;
      html << "<BODY>" << endl;

      // Eat up comment lines, if any, until hitting VER
      for(;;)
      {
         pf >> section;

         if( pf.fail() ) goto TERMINATE;

         if (section == "VER")
         {
            pf >> tmp;
            html << "<P> Certificate version " << tmp << "</P>" << endl;
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
            getline( pf, tmp ); // throw away the rest of the line
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

      colName.resize( num );
      isInt.resize( num );

      for( int i = 0; i < num; ++i )
      {
         pf >> colName[i];
         isInt[i] = false;
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
         else isInt[idx] = true;
      }


      pf >> section;
      if( section != "OBJ" )
      {
         cerr << "OBJ expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> tmp;
      if( tmp == "min" )
      {
         html << "<P><B>Minimize:</B></P>";
      }
      else if( tmp == "max" )
      {
         html << "<P><B>Maximize:</B></P> ";
      }
      else
      {
         cerr << "Unrecognized string after OBJ: " << tmp << endl;
         goto TERMINATE;
      }

      html << "<TABLE cellpadding='8'>" << endl;
      html << "<TR><TD>OBJ</TD>" << endl;
      html << "<TD>" << endl;

      stat = _processSparseVec( colName, false );

      html << "</TD>" << endl;
      html << "</TABLE>" << endl;


      if( !stat ) goto TERMINATE;


      pf >> section;
      if( section != "CON" )
      {
         cerr << "CON expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> numCon >> numBnd;

      html << "<P><B>Subject To:</B></P>" << endl;
      html << "<TABLE cellpadding='8'>" << endl;

      for( int i = 0; i < numCon; ++i )
      {
         html << "<TR>" << endl;
         pf >> label >> sense >> tmp;

         rowName.push_back( label );

         html << "<TD> " << i << " </TD>" << endl;
         html << "<TD> " << label << " </TD>" << endl;
         html << "<TD> ";

         stat = _processSparseVec( colName, false );
         if( stat )
         {
             if (sense == 'E') html << " = ";
             else if (sense == 'L') html << " &le; ";
             else if (sense == 'G') html << " &ge; ";
             else stat = false;
             html << tmp;
         }
         html << " </TD>" << endl;

         if( i >= numCon - numBnd ) html << "<TD> bound </TD>" << endl;

         html << "</TR>" << endl;

         if( !stat ) break;
      }

      html << "</TABLE>" << endl;


      pf >> section;
      if( section != "RTP" )
      {
         cerr << "RTP expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> tmp;

      html << "<P><B>Check:</B></P>" << endl;
      if( tmp == "infeas" )
      {
         html << "<P>infeasible</P>" << endl;
      }
      else if( tmp == "range" )
      {
         html << "<TABLE cellpadding='8'>" << endl;
         pf >> tmp;
         if( tmp != "-inf" )
         {
            html << "<TR><TD>lower bound</TD><TD>"
                   << tmp << "</TD></TR>" << endl;
         }
         pf >> tmp;
         if( tmp != "inf" )
         {
            html << "<TR><TD>upper bound</TD><TD>"
                   << tmp << "</TD></TR>" << endl;
         }
         html << "</TABLE>" << endl;
      }
      else
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
         html << "<P><B>Solutions:</B></P>" << endl;
         html << "<TABLE cellpadding='8'>" << endl;
         for( int i = 0; i < numSol; ++i ) {
            pf >> label;
            html << "<TR><TD>" << label << "</TD><TD>";
            stat = _processSparseVec( colName, true );
            html << "</TD></TR>" << endl;
         }
         html << "</TABLE>" << endl;
      }


      pf >> section;
      if( section != "DER" )
      {
         cerr << "DER expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> numDer;

      html << "<P><B>Derivations:</B></P>" << endl;
      html << "<TABLE cellpadding='8'>" << endl;

      for( int i = 0; i < numDer; ++i )
      {
         html << "<TR>" << endl;
         pf >> label >> sense >> tmp;

         rowName.push_back( label );

         html << "<TD> " << numCon + i << " </TD>" << endl;
         html << "<TD> " << label << " </TD>" << endl;
         html << "<TD> ";
         stat = _processSparseVec( colName, false );
         if( stat )
         {
            if (sense == 'E') html << " = ";
            else if (sense == 'L') html << " &le; ";
            else if (sense == 'G') html << " &ge; ";
            else stat = false;
            html << tmp;
         }
         html << " </TD>" << endl;

         html << "<TD> ";
         pf >> tmp;
         if( tmp != "{" )
         {
            cerr << "'{' expected.   Reading instead: " << tmp << endl;
            stat = false;
         }
         else
         {
            pf >> tmp;
            if( tmp == "asm" || tmp == "sol" )
            {
               html << tmp;
               pf >> tmp;
               if( tmp != "}")
               {
                  cerr << "'}' expected. Read instead: " << tmp << endl;
                  stat = false;
               }
            }
            else if( tmp == "lin" )
            {
               html << "lin ";
               stat = _processSparseVec( rowName, false );
               if( stat )
               {
                  pf >> tmp;
                  if( tmp != "}")
                  {
                     cerr << "'}' expected. Read instead: " << tmp << endl;
                     stat = false;
                  }
               }
            }
            else if( tmp == "rnd" )
            {
               html << "rnd ";
               stat = _processSparseVec( rowName, false );
               if( stat )
               {
                  pf >> tmp;
                  if( tmp != "}")
                  {
                     cerr << "'}' expected. Read instead: " << tmp << endl;
                     stat = false;
                  }
               }
            }
            else if( tmp == "uns" )
            {
               html << "uns ";
               pf >> con1 >> asm1 >> con2 >> asm2;
               if( pf.fail() )
               {
                  stat = false;
               }
               else
               {
                  html << rowName[con1] << ", " << rowName[con2];
                  html << " on ";
                  html << rowName[asm1] << ", " << rowName[asm2];

                  pf >> tmp;
                  if( tmp != "}")
                  {
                     cerr << "'}' expected. Read instead: " << tmp << endl;
                     stat = false;
                  }
               }
            }
            else
            {
               cerr << "Unrecognized reason type: " << tmp << endl;
               stat = false;
               for(;;)
               {
                  pf >> tmp;
                  if (pf.fail()) break;
                  if( tmp == "}" ) break;
               }
            }

         }

         html << " </TD>" << endl;

         if( stat )
         {
            html << "<TD> ";
            pf >> idx;
            html << idx;
            html << "</TD>";
         }

         html << "</TR>" << endl;

         if( !stat ) break;
      }

      html << "</TABLE>" << endl;
      if( stat ) rs = 0;

TERMINATE:
      if( !stat ) {
         html << "<P style='color:red;'>";
         html << "Error encountered while processing file";
         html << "</P>" << endl;
      }
      html << "</BODY>" << endl;
      html << "</HTML>" << endl;
      html.close();
   }

   pf.close();

   return rs;
}
