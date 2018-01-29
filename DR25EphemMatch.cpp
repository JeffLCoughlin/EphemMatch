/*
DR25 Kepler Ephemeris Matching
Written by Jeffrey L. Coughlin
jcoughlin@seti.org

See https://github.com/JeffLCoughlin/EphemMatch for instructions on how to compile and run the code, an example, and how to cite the code if used.

*/

// Include libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <vector>

// Declare constants
const double RACENT=19.3777777777777;  // Center of Kepler FoV in RA
const double DECCENT=44.5;  // Center of Kepler FoV in Dec
const double PI=3.14159265359;
const double PF1 = 1.08228446E-05;  // Period matching tolerance fit variable 1
const double PF2 = 1.37838040;  // Period matching tolerance fit variable 2
const double TF1 = 0.00969916;  // Epoch matching tolerance fit variable 1
const double TF2 = 0.36253537;  // Epoch matching tolerance fit variable 2

// Standard namespace
using namespace std;

// Declare struct to hold all the data used for matching
struct datastruct {int kic; string id; double kepmag; double period; double tzero; string disp; double duration; double depth; double ra; double dec; int chan[4]; int mod[4]; int out[4]; int row[4]; int col[4];};

// Functions
void READTCEFILE();  // Read in TCE information
void READKOIFILE();  // Read in KOI information
void READKEBFILE();  // Read in Kepler EB information
void READGEBFILE();  // Read in Ground-based EB information
void READCCPARAM();  // Read in CCD parameter file
void DOMATCH(int, int, vector <datastruct> A1, vector <datastruct> A2, string);  // Compare ephemerides of two objects to see if match is found
double STARDIST(double, double, double, double);  // Input decimal RA and DEC, both in DEGREES, output distance between them in arseconds
double PERR(double), TERR(double);  // Functions that return the error associated with a given period or epoch (based on fit to injected/recovered transits).
double INVERFC(double); // Inverse complementary error function

// Variables - Sorry they are not all documented.
int ntce,nkoi,nkeb,ngeb;  // Total number of TCEs, KOIs, Kepler EBs, and Ground-based EBs
double ptol,ttol,dtol,postol,rowtol,coltol,reflpostol,perratmax;
double pdiff1,tdiff1,pdiff2,tdiff2,ddiff1,ddiff2,posdiff,reflfacra,reflfacdec;
string tmpstr,tmpstr3[25];
stringstream iss;
ofstream outfile;
ifstream infile;
ofstream plotfile,posfile;
datastruct tmpdata;
double prat,trat;
vector <datastruct> tcedat,koidat,kebdat,gebdat;


// Main
int main() {

cout << "Reading Input Files..." << endl;

READTCEFILE();

READKOIFILE();

READKEBFILE();

READGEBFILE();

cout << "Combining KIC data with input data..." << endl;

READCCPARAM();

cout << "Performing Matching..." << endl;

ptol = 5.0;  // Sigma tolerance for Period.
ttol = 5.0;  // Sigma tolerance for T0
dtol = -1.0;  // Sigma tolerance for duration. Set to -1 because have not found duration to be adequetly measured to be used in matching (so -1 means no duration matching done.)
reflpostol = 100; // Reflection position tolerance in arcseconds
rowtol = 5.0;   // Row difference tolerance in pixels
coltol = 5.0;   // Column difference tolerance pixels
perratmax = 50.5;  // Maximum period ratio - Set to 50.5 so no period ratios higher than 50:1 are found.

// Open GNUPlot plotting file
plotfile.open("plotfile.gnu");
plotfile << "set term postscript" << endl << "set output 'MatchGraph.ps'" << endl << "set xlabel 'RA (Hours)'" << endl << "set ylabel 'DEC (Deg)'" << endl;
posfile.open("posfile.dat");


// Do the matching - I'm matching every TCE to every other TCE, KOI, KEB, and GEB
cout << "TCE-TCE Matches..." << endl;
DOMATCH(ntce,ntce,tcedat,tcedat,"TCETCEMatches.txt");

cout << "TCE-KOI Matches..." << endl;
DOMATCH(ntce,nkoi,tcedat,koidat,"TCEKOIMatches.txt");

cout << "TCE-KEB Matches..." << endl;
DOMATCH(ntce,nkeb,tcedat,kebdat,"TCEKEBMatches.txt");

cout << "TCE-GEB Matches..." << endl;
DOMATCH(ntce,ngeb,tcedat,gebdat,"TCEGEBMatches.txt");

// Make plots
cout << "Making plotting file..." << endl;
posfile << "-1 " << RACENT << " " << DECCENT << endl;
plotfile << "plot 'posfile.dat' u ($1==1 ? $2 : 1/0):3 pt 7 lc -1 ps 0.5 title '1:1', 'posfile.dat' u ($1==2 ? $2 : 1/0):3 pt 7 lc  1 ps 0.5 title '2:1', 'posfile.dat' u ($1>2 ? $2 : 1/0):3 pt 7 lc  3 ps 0.5 title '3:1+', 'posfile.dat' u ($1==-1 ? $2 : 1/0):3 pt 3 lc 2 ps 2 title 'FoV Center'" << endl;
plotfile.close();
posfile.close();

// Bunch of system calls to make nicely formatted output files using unix tools
system("tail -n +3 TCETCEMatches.txt > tmp1");
system("tail -n +3 TCEKOIMatches.txt >> tmp1");
system("tail -n +3 TCEKEBMatches.txt >> tmp1");
system("tail -n +3 TCEGEBMatches.txt >> tmp1");

system("cat tmp1 | paste - - - | sort -k 3n | sed 's/\\t/\\n/g' - > tmp2");
system("echo '#             ID       KIC      Period  T0-2454900  Duration    Depth DISP KepMag   P1/P2   dT0/P  Pmatch(σ)  Tmatch(σ)  Dmatch(σ)  Sep(\") Chan0 Mod0 Out0 Row0 Col0 Chan1 Mod1 Out1 Row1 Col1 Chan2 Mod2 Out2 Row2 Col2 Chan3 Mod3 Out3 Row3 Col3\n' | cat - tmp2 > AllMatchesSorted.txt");

cout << "Making Best Matches..." << endl;
system("awk '{if($1!=\"\") print($0)}' AllMatchesSorted.txt > tmp1");
system("awk '{if(NR==1 || NR%2==1) print($0); else printf($0\" NEWLINE \")}' tmp1 > tmp2");
system("sort -k 1,1 tmp2 > tmp3");

// Column 41 is depth - choose most likely parent as one with highest depth
system("awk 'BEGIN{a=0;b=0;d=0} {if(NR==1) {print($0); printf(\"\\n\")} else {if(a!=$1) {a=$1; d=0; if(NR>2) print(b)}; if(a==$1 && $41>d) {b=$0; d=$41}}} END{print(b)}' tmp3 > tmp4");
system("sed 's/ NEWLINE /\\n/g' tmp4 > tmp5");
system("awk '{print($0); if(NR>2 && NR%2==0) printf(\"\\n\")}' tmp5 > BestMatches.txt");

// Make BestMatchesSorted.txt
system("tail -n +3 BestMatches.txt > tmp1");
system("cat tmp1 | paste - - - | sort -k 3n | sed 's/\\t/\\n/g' - > tmp2");
system("echo '#             ID       KIC      Period  T0-2454900  Duration    Depth DISP KepMag   P1/P2   dT0/P  Pmatch(σ)  Tmatch(σ)  Dmatch(σ)  Sep(\") Chan0 Mod0 Out0 Row0 Col0 Chan1 Mod1 Out1 Row1 Col1 Chan2 Mod2 Out2 Row2 Col2 Chan3 Mod3 Out3 Row3 Col3\n' | cat - tmp2 > BestMatchesSorted.txt");

// This takes BestMatchesSorted.txt and groups systems by period, eliminating duplicates
cout << "Making Period Groups..." << endl;
system("awk '{if($1!=\"\") print($0)}' BestMatchesSorted.txt > tmp1");
system("awk '{if(NR==1 || NR%2==1) print($0); else printf($0\" NEWLINE \")}' tmp1 > tmp2");
system("awk 'BEGIN{a=1} {if(($3/a<0.9995 || $3/a>1.0005) && NR>1) printf(\"\\n\"); print($0); if(NR>1) a=$3}' tmp2 > tmp3");
system("sed 's/ NEWLINE /\\n/g' tmp3 > tmp4");
system("awk '{for(j=0;j<i;j++) if($1==a[j]) sw=1; if(sw==0) print($0); a[i]=$1; sw=0; i++; if($1==\"\") i=0}' tmp4 > PeriodGroups.txt");


cout << "Making Results File..." << endl;
system("awk -f /home/jeff/KeplerJob/DR25/Code/BestMatch-To-Results.awk < BestMatches.txt | sort -k 1n | awk -f /home/jeff/KeplerJob/DR25/Code/RemoveRedundant.awk > Results.txt");

// Make the plots
system("gnuplot plotfile.gnu");

// Clean up temp files
system("rm tmp*");

exit(0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PERR (double per) {

return PF1*pow(per,PF2);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double TERR (double per) {

return TF1*pow(per,TF2);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double STARDIST(double ra1, double dec1, double ra2, double dec2) {  // Input decimal ra and dec, both in degrees!

ra1*=PI/180.0;
dec1*=PI/180.0;
ra2*=PI/180.0;
dec2*=PI/180.0;

return 3600*(180.0/PI)*acos(sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra1-ra2));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DOMATCH(int N1, int N2, vector <datastruct> A1, vector <datastruct> A2, string outname)  {

int i,j,k,l;

outfile.open(outname.c_str());
outfile << "#             ID       KIC      Period  T0-2454900  Duration    Depth DISP KepMag   P1/P2   dT0/P  Pmatch(σ)  Tmatch(σ)  Dmatch(σ)  Sep(\") Chan0 Mod0 Out0 Row0 Col0 Chan1 Mod1 Out1 Row1 Col1 Chan2 Mod2 Out2 Row2 Col2 Chan3 Mod3 Out3 Row3 Col3" << endl << endl;
for(i=0;i<N1;i++)
  for(j=0;j<N2;j++)
    if(A1[i].id!=A2[j].id)
      {
      if(A1[i].kepmag < A2[j].kepmag && A1[i].kepmag > 0)  // Choose the brighter star for maximum interaction distance.
        postol = 55*sqrt(pow(10,-0.4*A1[i].kepmag+6.0)+1); //    40+pow(10,7.25)*pow(10,-0.4*A1[i].kepmag);
      else
        {
        if(A2[j].kepmag > 0)
          postol = 55*sqrt(pow(10,-0.4*A2[j].kepmag+6.0)+1); //40+pow(10,7.25)*pow(10,-0.4*A2[j].kepmag);
        else
          postol=55;
        }

      posdiff = STARDIST((360.0/24.0)*A1[i].ra,A1[i].dec,(360.0/24.0)*A2[j].ra,A2[j].dec); // 3600*sqrt(pow((360.0/24.0)*(A1[i].ra-A2[j].ra),2)+pow(A1[i].dec-A2[j].dec,2));  // Position difference in arcseconds
      reflfacra  = 3600*(360.0/24.0)*fabs((A1[i].ra+A2[j].ra)/2.0 - RACENT);  // Center of line between match deviation from FoV RA center
      reflfacdec = 3600*fabs((A1[i].dec+A2[j].dec)/2.0 - DECCENT);  // Center of line between match deviation from FoV Dec center

      if(  // Impose constraints on what matches are allowed, based on relative position of objects and other factors
        (A1[i].kic != A2[j].kic || A1[i].kic==0 || A2[j].kic==0)                                          &&  // Can't be the identical KIC number, unless at least one is zero, AND
        (A1[i].period/A2[j].period < perratmax && A2[j].period/A1[i].period < perratmax)                  &&  // Require period ratios of 50:1 or less, AND
          (
          posdiff < postol                                                                                ||  // The stars can be close to each other - Direct PRF
          (reflfacra < reflpostol && reflfacdec < reflpostol)                                             ||  // Or reflect through the FoV center - Antidopodal Reflection
          (A1[i].mod[0]==A2[j].mod[0] && A1[i].col[0] > 0 && A2[j].col[0] > 0 && fabs(A1[i].col[0]-A2[j].col[0])<=coltol) ||  // Or can be same column on same module - col anomaly
          (A1[i].mod[1]==A2[j].mod[1] && A1[i].col[1] > 0 && A2[j].col[1] > 0 && fabs(A1[i].col[1]-A2[j].col[1])<=coltol) ||  // Or can be same column on same module - col anomaly
          (A1[i].mod[2]==A2[j].mod[2] && A1[i].col[2] > 0 && A2[j].col[2] > 0 && fabs(A1[i].col[2]-A2[j].col[2])<=coltol) ||  // Or can be same column on same module - col anomaly
          (A1[i].mod[3]==A2[j].mod[3] && A1[i].col[3] > 0 && A2[j].col[3] > 0 && fabs(A1[i].col[3]-A2[j].col[3])<=coltol) ||  // Or can be same column on same module - col anomaly
          (A1[i].mod[0]==A2[j].mod[0] && A1[i].row[0] > 0 && A2[j].row[0] > 0 && A1[i].col[0] > 0 && A2[j].col[0] > 0 && fabs(A1[i].row[0]-A2[j].row[0])<=rowtol && fabs(A1[i].col[0]-A2[j].col[0])<=coltol) ||  // Or can be same row AND column on same mod - cross-talk
          (A1[i].mod[1]==A2[j].mod[1] && A1[i].row[1] > 0 && A2[j].row[1] > 0 && A1[i].col[1] > 0 && A2[j].col[1] > 0 && fabs(A1[i].row[1]-A2[j].row[1])<=rowtol && fabs(A1[i].col[1]-A2[j].col[1])<=coltol) ||  // Or can be same row AND column on same mod - cross-talk
          (A1[i].mod[2]==A2[j].mod[2] && A1[i].row[2] > 0 && A2[j].row[2] > 0 && A1[i].col[2] > 0 && A2[j].col[2] > 0 && fabs(A1[i].row[2]-A2[j].row[2])<=rowtol && fabs(A1[i].col[2]-A2[j].col[2])<=coltol) ||  // Or can be same row AND column on same mod - cross-talk
          (A1[i].mod[3]==A2[j].mod[3] && A1[i].row[3] > 0 && A2[j].row[3] > 0 && A1[i].col[3] > 0 && A2[j].col[3] > 0 && fabs(A1[i].row[3]-A2[j].row[3])<=rowtol && fabs(A1[i].col[3]-A2[j].col[3])<=coltol)     // Or can be same row AND column on same mod - cross-talk
          )
        )
        {

      // Do match - have to know which object has shorter period to do matching correctly and figure out period ratios
         if(A1[i].period<A2[j].period)
           {
           prat = rint(A2[j].period / A1[i].period);
           trat = rint((A1[i].tzero  - A2[j].tzero)/A1[i].period);

           pdiff1 = fabs(prat*A1[i].period - A2[j].period)/PERR(A1[i].period);
           tdiff1 = fabs(A1[i].tzero  - A2[j].tzero - trat*A1[i].period)/(sqrt(2)*TERR(A1[i].period));

           pdiff2=pdiff1;
           tdiff2=tdiff1;
           }
         else
           {
           prat = rint(A1[i].period / A2[j].period);
           trat = rint((A2[j].tzero  - A1[i].tzero)/A2[j].period);

           pdiff1 = fabs(prat*A2[j].period - A1[i].period)/PERR(A2[j].period);
           tdiff1 = fabs(A2[j].tzero  - A1[i].tzero - trat*A2[j].period)/(sqrt(2)*TERR(A2[j].period));

           pdiff2=pdiff1;
           tdiff2=tdiff1;
           }

        // Tried to do duration matching, but found it wasn't useful, so setting to 1, which basically disables duration matching
        ddiff1 = 1; //sqrt(2)*INVERFC(fabs((A1[i].duration-A2[j].duration)/A1[i].duration) - rint((A1[i].duration-A2[j].duration)/A1[i].duration));
        ddiff2 = 1; //sqrt(2)*INVERFC(fabs((A2[j].duration-A1[i].duration)/A2[j].duration) - rint((A2[j].duration-A1[i].duration)/A2[j].duration));

        // If it's a match, then print out results of match
        if(pdiff1<ptol && (tdiff1<ttol || A2[j].id=="RR-Lyr-pri"))  // P and T0 must match within tolerance. Except for RR-Lyr-pri, don't match epoch, it gets distorted by TPS too much.
          {
          outfile << fixed << setw(16)                                   << A1[i].id                               << " "
                          << setw(9)                     << setfill('0') << A1[i].kic                              << " "
                          << setw(11) << setprecision(6) << setfill(' ') << A1[i].period                           << " "
                          << setw(11)                                    << A1[i].tzero                            << " "
                          << setw(9)  << setprecision(4) << setfill(' ') << A1[i].duration                         << " "
                          << setw(8)  << setprecision(6) << setfill(' ') << A1[i].depth                            << " "
                          << setw(4)                                     << A1[i].disp                             << " "
                          << setw(6)  << setprecision(3)                 << A1[i].kepmag                           << " "
                          << setw(7)  << setprecision(3)                 << A1[i].period/A2[j].period              << " "
                          << setw(7)  << setprecision(1) << setfill(' ') << (A1[i].tzero-A2[j].tzero)/A1[i].period << " "
                          << setw(10) << setprecision(3) << setfill(' ') << pdiff1                                 << " "
                          << setw(10) << setprecision(3) << setfill(' ') << tdiff1                                 << " "
                          << setw(10) << setprecision(3) << setfill(' ') << ddiff1                                 << " "
                          << setw(7)  << setprecision(1) << setfill(' ') << posdiff                                << " ";
          for(l=0;l<4;l++)
            {
            outfile << setw(5) << setfill(' ') << A1[i].chan[l] << " "
                    << setw(4) << setfill(' ') << A1[i].mod[l]  << " "
                    << setw(4) << setfill(' ') << A1[i].out[l]  << " "
                    << setw(4) << setfill(' ') << A1[i].row[l]  << " "
                    << setw(4) << setfill(' ') << A1[i].col[l]  << " ";
            }
          outfile << endl;

          outfile << fixed << setw(16)                                   << A2[j].id                               << " "
                          << setw(9)                     << setfill('0') << A2[j].kic                              << " "
                          << setw(11) << setprecision(6) << setfill(' ') << A2[j].period                           << " "
                          << setw(11)                                    << A2[j].tzero                            << " "
                          << setw(9)  << setprecision(4) << setfill(' ') << A2[j].duration                         << " "
                          << setw(8)  << setprecision(6) << setfill(' ') << A2[j].depth                            << " "
                          << setw(4)                                     << A2[j].disp                             << " "
                          << setw(6)  << setprecision(3)                 << A2[j].kepmag                           << " "
                          << setw(7)  << setprecision(3)                 << A2[j].period/A1[i].period              << " "
                          << setw(7)  << setprecision(1) << setfill(' ') << (A2[j].tzero-A1[i].tzero)/A2[j].period << " "
                          << setw(10) << setprecision(3) << setfill(' ') << pdiff2                                 << " "
                          << setw(10) << setprecision(3) << setfill(' ') << tdiff2                                 << " "
                          << setw(10) << setprecision(3) << setfill(' ') << ddiff2                                 << " "
                          << setw(7)  << setprecision(1) << setfill(' ') << posdiff                                << " ";
          for(l=0;l<4;l++)
            {
            outfile << setw(5) << setfill(' ') << A2[j].chan[l] << " "
                    << setw(4) << setfill(' ') << A2[j].mod[l]  << " "
                    << setw(4) << setfill(' ') << A2[j].out[l]  << " "
                    << setw(4) << setfill(' ') << A2[j].row[l]  << " "
                    << setw(4) << setfill(' ') << A2[j].col[l]  << " ";
            }
          outfile << endl << endl;


          if(A1[i].period/A2[j].period > 1.0)
            k = rint(A1[i].period/A2[j].period);
          else
            k = rint(A2[j].period/A1[i].period);

          posfile << k << " " << setprecision(10) << A1[i].ra << " " << A1[i].dec << endl << setprecision(0) << k << " " << setprecision(10) << A2[j].ra << " " << A2[j].dec << endl;
          plotfile << "set arrow from " << A1[i].ra << "," << A1[i].dec << " to " << A2[j].ra << "," << A2[j].dec << " nohead lc ";
          if(k==1)
            plotfile << "-1" << endl;
          if(k==2)
            plotfile << "1" << endl;
          if(k>2)
            plotfile << "3" << endl;
          }
        }
      }
outfile.close();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void READTCEFILE() {

int z;

infile.open("TCEInput.txt");   // TCEInput.txt should have 6 columns: TCE Number, KIC, Period, T0, Width, Depth

z=0;
infile >> tmpdata.id;
while(!infile.eof())  {
  tcedat.push_back(tmpdata);
  infile >> tcedat[z].kic;
  infile >> tcedat[z].period;
  infile >> tcedat[z].tzero;
  infile >> tcedat[z].duration;  // Is in hours
  infile >> tcedat[z].depth;
  tcedat[z].disp="TCE";
  tcedat[z].tzero-=67.0;
  tcedat[z].depth/=1.0E6;  // Convert TCE depth from ppm to fractional
  if(tcedat[z].depth<=1E-6) tcedat[z].depth=1E-6;  // Some TCE depths are zero, so make the minimum 1E-6 so no dividing by zero occurs
  z++;
  infile >> tmpdata.id;
  }
infile.close();
ntce=z;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void READKOIFILE() {

int z;

infile.open("KOIInput.txt");   // KOIInput.txt should have 6 columns: KOI Number, KIC, Period, T0, Width, Depth

z=0;
infile >> tmpdata.id;
while(!infile.eof())  {
  koidat.push_back(tmpdata);
  infile >> koidat[z].kic;
  infile >> koidat[z].period;
  infile >> koidat[z].tzero;
  infile >> koidat[z].duration;  // Should be in hours
  infile >> koidat[z].depth;  // Should be fractioanl
  koidat[z].disp="KOI";
  z++;
  infile >> tmpdata.id;
  }
infile.close();
nkoi=z;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void READKEBFILE() {

int z;

infile.open("KEBInput.txt");  // KEBInput.txt should have 6 columns: EBID KICNumber  Period  Epoch  Width Depth

z=0;
infile >> tmpdata.id;
while(!infile.eof())
  {
  kebdat.push_back(tmpdata);
  infile >> kebdat[z].kic;
  infile >> kebdat[z].period;
  infile >> kebdat[z].tzero;
  infile >> kebdat[z].duration;
  infile >> kebdat[z].depth;
  kebdat[z].tzero-=54900.0;  // Offset from EB Catalog to Jason's time system
  kebdat[z].duration=kebdat[z].duration*kebdat[z].period*24.0;  // Converts from width in phase to width in hours
  kebdat[z].disp = "KEB";
  z++;
  infile >> tmpdata.id;
  }
infile.close();
nkeb = z;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void READGEBFILE() {

int z;

infile.open("GEBInput.txt");  // GEBinput.txt should have 7 columns: EBName Mag RA  Dec Period  Epoch  Depth

z=0;
infile >> tmpdata.id;
while(!infile.eof())
  {
  gebdat.push_back(tmpdata);
  infile >> gebdat[z].kepmag;
  infile >> gebdat[z].ra;
  infile >> gebdat[z].dec;
  infile >> gebdat[z].period;
  infile >> gebdat[z].tzero;
  infile >> gebdat[z].depth;
  gebdat[z].tzero-=2454900.0;  // Offset from Ground-based EBs to Jason's time system
  while(gebdat[z].tzero<0 && gebdat[z].ra>0)  // Make epoch positive -  Need the ra part else it will get in an infinite loop at end of file sometimes
    gebdat[z].tzero+=gebdat[z].period;
  gebdat[z].duration=0.0;  // Don't have duration info
  gebdat[z].kic = 0;  // No implicit KIC, but will do matching in other function
  gebdat[z].disp = "GEB";
  z++;
  infile >> tmpdata.id;
  }
infile.close();
ngeb = z;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void READCCPARAM() {

int i,j,k;

infile.open("KIC-CCD-Params.txt");  // 25 pipe-delimited fields: ID|KIC|RA|Dec|kepmag|Channel_0|Module_0|Output_0|Row_0|Column_0|Channel_1|Module_1|Output_1|Row_1|Column_1|Channel_2|Module_2|Output_2|Row_2|Column_2|Channel_3|Module_3|Output_3|Row_3|Column_3
while(!infile.eof())
  {
  getline(infile,tmpstr);
  iss << tmpstr;
  for(j=0;j<25;j++)
    getline(iss,tmpstr3[j],'|');

  for(i=0;i<ntce;i++)
    if(tcedat[i].kic == atoi(tmpstr3[1].c_str()) || tcedat[i].id == tmpstr3[0].c_str())
      {
      tcedat[i].ra  = atof(tmpstr3[2].c_str());
      tcedat[i].dec = atof(tmpstr3[3].c_str());
      tcedat[i].kepmag = atof(tmpstr3[4].c_str());
      for(j=0;j<4;j++)
        {
        tcedat[i].chan[j] = atof(tmpstr3[5*j+5].c_str());
        tcedat[i].mod[j]  = atof(tmpstr3[5*j+6].c_str());
        tcedat[i].out[j]  = atof(tmpstr3[5*j+7].c_str());
        tcedat[i].row[j]  = atof(tmpstr3[5*j+8].c_str());
        tcedat[i].col[j]  = atof(tmpstr3[5*j+9].c_str());
        }
      }

  for(i=0;i<nkoi;i++)
    if(koidat[i].kic == atoi(tmpstr3[1].c_str()) || koidat[i].id == tmpstr3[0].c_str())
      {
      koidat[i].ra  = atof(tmpstr3[2].c_str());
      koidat[i].dec = atof(tmpstr3[3].c_str());
      koidat[i].kepmag = atof(tmpstr3[4].c_str());
      for(j=0;j<4;j++)
        {
        koidat[i].chan[j] = atof(tmpstr3[5*j+5].c_str());
        koidat[i].mod[j]  = atof(tmpstr3[5*j+6].c_str());
        koidat[i].out[j]  = atof(tmpstr3[5*j+7].c_str());
        koidat[i].row[j]  = atof(tmpstr3[5*j+8].c_str());
        koidat[i].col[j]  = atof(tmpstr3[5*j+9].c_str());
        }
      }

  for(i=0;i<nkeb;i++)
    if(kebdat[i].kic == atoi(tmpstr3[1].c_str()) || kebdat[i].id == tmpstr3[0].c_str())
      {
      kebdat[i].ra  = atof(tmpstr3[2].c_str());
      kebdat[i].dec = atof(tmpstr3[3].c_str());
      kebdat[i].kepmag = atof(tmpstr3[4].c_str());
      for(j=0;j<4;j++)
        {
        kebdat[i].chan[j] = atof(tmpstr3[5*j+5].c_str());
        kebdat[i].mod[j]  = atof(tmpstr3[5*j+6].c_str());
        kebdat[i].out[j]  = atof(tmpstr3[5*j+7].c_str());
        kebdat[i].row[j]  = atof(tmpstr3[5*j+8].c_str());
        kebdat[i].col[j]  = atof(tmpstr3[5*j+9].c_str());
        }
      }

  for(i=0;i<ngeb;i++)
    if(gebdat[i].id == tmpstr3[0].c_str())
      {
      gebdat[i].kic=atoi(tmpstr3[1].c_str());
      gebdat[i].ra  = atof(tmpstr3[2].c_str());
      gebdat[i].dec = atof(tmpstr3[3].c_str());
      gebdat[i].kepmag = atof(tmpstr3[4].c_str());
      for(j=0;j<4;j++)
        {
        gebdat[i].chan[j] = atof(tmpstr3[5*j+5].c_str());
        gebdat[i].mod[j]  = atof(tmpstr3[5*j+6].c_str());
        gebdat[i].out[j]  = atof(tmpstr3[5*j+7].c_str());
        gebdat[i].row[j]  = atof(tmpstr3[5*j+8].c_str());
        gebdat[i].col[j]  = atof(tmpstr3[5*j+9].c_str());
        }
      }

  iss.clear();
  }
infile.close();
}

///////////////////////////////////////////////////////////////////////////////////////////

double INVERFC(double p)
  {
  double x, err, t, pp;

  if (p >= 2.) return -100.;
  if (p <= 0.0) return 100.;

  pp=(p < 1.0)? p:2.-p;
  t=sqrt(-2.*log(pp/2));

  x= -0.70711*((2.30753+t*0.27061)/(1+t*(0.99229+t*0.04481))-t);

  for (int j=0;j<2;j++)
    {
    err=erfc(x)-pp;
    x += err/(1.12837916709551257*exp(-x*x)-x*err);
    }
  return (p<1.0 ? x:-x);
  }
