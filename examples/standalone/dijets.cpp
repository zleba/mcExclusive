//#include "dijetAnal.h"

#include "Pythia8/Pythia.h"
#include "ExclusiveHooks.h"

#include <iomanip>

using namespace std;

//ExclusiveHooks *exclHooks;

void Fill(Pythia8::Pythia *pythia);

int main()
{

    Pythia8::Pythia pythia;
    pythia.readString("Random:setSeed = on");
    //pythia.readString(TString::Format("Random:seed = %d", seed ).Data()  );

    pythia.readString("PDF:pSet =  LHAPDF6:MMHT2014lo68cl");
    pythia.readString("Beams:FrameType = 2");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");

    pythia.readString( "Beams:eA = 6500");
    pythia.readString( "Beams:eB = 6500");

    pythia.readString( "SigmaProcess:alphaSvalue= 0.127");
    pythia.readString( "SpaceShower:alphaSvalue = 0.127"); //0.135
    pythia.readString( "TimeShower:alphaSvalue  = 0.127");


    //QCD
    //pythia.readString("HardQCD:all = on");
    pythia.readString("HardQCD:gg2gg = on");


    pythia.readString("PartonLevel:MPI = off");

    pythia.readString("PartonLevel:FSR = on");
    pythia.readString("PartonLevel:ISR = off");

    pythia.readString("HadronLevel:all = on");



    pythia.readString("SpaceShower:pTmin = 1.5");


    pythia.readString("PhaseSpace:pTHatMin = 30.");
    //pythia.readString("PhaseSpace:mHatMin = 60");



    ExclusiveHooks *exclHooks = new ExclusiveHooks(&pythia);
    pythia.setUserHooksPtr(exclHooks);


    const double radius = 0.7;
    const double pTjetMin = 10;
    const double etaMax = 4;
    const int nSel = 2;
    Pythia8::SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);


    pythia.init();

    const int N = 40;

    double sumWgt = 0;
    double sum2Wgt = 0;

    for (int iEvent = 0; iEvent < N; ++iEvent) {
        if (!pythia.next()) continue;

        if(!exclHooks->IsExclusive()) continue;

        if(exclHooks->GetnMPI() == 0) {
            Double weight = exclHooks->ExclusiveWeight() * pythia.info.weight();
            sumWgt += weight;
            sum2Wgt += weight*weight;
            cout << "Weight is " <<setprecision(5)<< weight << endl;

            if(iEvent < 3) pythia.event.list();
            slowJet.analyze(pythia.event);
            if (iEvent < 5) slowJet.list();
        }

    }

    //Inclusive luminosity in nb
    Double lumi =  1e-6*pythia.info.weightSum() /  pythia.info.sigmaGen() ;

    double err = sqrt(sum2Wgt) / sumWgt;
    cout << "Cross section is     " << sumWgt / lumi <<" nb"<< endl;
    cout << "Cross section err is " <<err* sumWgt / lumi<<" nb"<< endl;


    //SaveAs( TString::Format("temp/higgsHistos/dijets_%d.root",seed).Data(), lumi);


    pythia.stat();

    exclHooks->PrintBestParams();

}


