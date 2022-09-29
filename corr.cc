#include "itensor/all.h"
#include "Corr.h"
#include "ReadInput.h"
using namespace itensor;
using namespace std;

int main (int argc, char* argv[])
{
    string infile = argv[1];
    if(argc != 2) { printfln("Usage: %s inputfile_exthubbard",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto mps_files = read_vector<string> (infile, "mps_file");

    MPS psi;
    for(auto mps_file : mps_files)
    {
        readFromFile (mps_file, psi);

        auto lamb_up = MeasureCorr (psi, "Cdagup", "Cup");
        auto lamb_dn = MeasureCorr (psi, "Cdagdn", "Cdn");
        auto corr_dag = MeasureCorr (psi, "Cdagup", "Cdagdn");
        auto corr_CC = MeasureCorr (psi, "Cup", "Cdn");

        cout << mps_file << endl;
        cout << "Cdagup Cup =" << endl;
        cout << lamb_up << endl << endl;
        cout << "Cdagdn Cdn =" << endl;
        cout << lamb_dn << endl << endl;
        cout << "Cdagup Cdagdn =" << endl;
        cout << corr_dag << endl << endl;
        cout << "Cup Cdn =" << endl;
        cout << corr_CC << endl << endl;
    }
}
