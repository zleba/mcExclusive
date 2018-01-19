#include "Pythia8/Pythia.h"
#include "ExclusiveHooks.h"

#include <iomanip>

using namespace std;

//Standalone example code for the mcExclusive Add-On to Pythia, see [arXiv:1608.03765] for explanation of the method

////The currently implemented processes are (see http://home.thep.lu.se/~torbjorn/pythia82html/ProcessSelection.html for process IDs)
////Hard QCD
//ID=111 : HardQCD:gg2gg
//ID=112 : HardQCD:gg2qqbar
//ID=113 : HardQCD:qg2qg
//ID=114 : HardQCD:qq2qq
//ID=115 : HardQCD:qqbar2gg
//ID=116 : HardQCD:qqbar2qqbarNew
//
////charm production
//ID=121 : HardQCD:gg2ccbar
//ID=122 : HardQCD:qqbar2ccbar
//
////bottom production
//ID=123 : HardQCD:gg2bbbar
//ID=124 : HardQCD:qqbar2bbbar
//
////top production
//ID=601 : Top:gg2ttbar
//ID=602 : Top:qqbar2ttbar
//
////electroweak
//ID=204 : PromptPhoton:ffbar2gammagamma
//ID=205 : PromptPhoton:gg2gammagamma
//ID=224 : WeakSingleBoson:ffbar2ffbar(s:gmZ) (only to mu)
//
////Higgs to bb
//ID=901 : HiggsSM:ffbar2H (ony bb decay)
//ID=902 : HiggsSM:gg2H    (ony bb decay)



int main()
{

    Pythia8::Pythia pythia;
    pythia.readString("Next:numberShowEvent    = 0");
    pythia.readString("Next:numberShowProcess  = 0");
    pythia.readString("Random:setSeed = on");

    //Use LO sets which have regular low-scale behaviour, alternatively e.g  "LHAPDF6:MMHT2014lo68cl"
    pythia.readString("PDF:pSet =  13"); 
    pythia.readString("Beams:FrameType = 2");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");

    pythia.readString( "Beams:eA = 6500");
    pythia.readString( "Beams:eB = 6500");

    pythia.readString( "SigmaProcess:alphaSvalue= 0.127");
    pythia.readString( "SpaceShower:alphaSvalue = 0.127"); //0.135
    pythia.readString( "TimeShower:alphaSvalue  = 0.127");


    //Process selection
    //pythia.readString("HardQCD:all = on");
    pythia.readString("HardQCD:gg2gg = on");
    //pythia.readString("Top:gg2ttbar = on");
    //pythia.readString("HardQCD:gg2bbbar = on");

    //if MPI switched to on, the Soft rescattering probability would be simulated using MPI, 
    //othewise one needs to aplly own S2, i.e. from Superchic2 generator
    pythia.readString("PartonLevel:MPI = off"); 
                                                

    pythia.readString("PartonLevel:FSR = off"); 

    //if ISR switched to on, emissions in the Space-like shower are considered up to scale SpaceShower:pTmin, see arXiv:1608.03765, 
    //as it open new production chanel it fills the lower values in M12/MX distribution
    pythia.readString("PartonLevel:ISR = on"); 
                                                

    pythia.readString("HadronLevel:all = off");


    //The minimal treashold for the ISR, if it's switched on
    //is a free parameter of the model to be tuned on M12/MX distribution
    pythia.readString("SpaceShower:pTmin = 1.5"); 


    pythia.readString("PhaseSpace:pTHatMin = 30.");
    //pythia.readString("PhaseSpace:mHatMin = 60");



    ExclusiveHooks *exclHooks = new ExclusiveHooks(&pythia); //needs to be called before pythia.init()
    pythia.setUserHooksPtr(exclHooks);


    const double radius = 0.7;
    const double pTjetMin = 10;
    const double etaMax = 4;
    const int nSel = 2;
    Pythia8::SlowJet slowJet(-1, radius, pTjetMin, etaMax, nSel, 1); //AKT0.7 jets as an example


    pythia.init();

    const int N = 100;

    double sumWgt = 0;
    double sum2Wgt = 0;

    for (int iEvent = 0; iEvent < N; ++iEvent) {
        if (!pythia.next()) continue;

        if(!exclHooks->IsExclusive()) continue; //if not exclusive

        if(exclHooks->GetnMPI() == 0) {
            Double weight = exclHooks->ExclusiveWeight() * pythia.info.weight(); //the Pythia's weight multiplied by the probability of being exclusive
            sumWgt += weight;
            sum2Wgt += weight*weight;
            cout << "Weight is " <<setprecision(5)<< weight << endl;

            //if(iEvent < 3) pythia.process.list();
            if(iEvent < 20) pythia.event.list();
            slowJet.analyze(pythia.event);
            if (iEvent < 5) slowJet.list();
        }

    }

    //Inclusive luminosity in nb
    Double lumi =  1e-6*pythia.info.weightSum() /  pythia.info.sigmaGen() ;

    double err = sqrt(sum2Wgt) / sumWgt;
    cout << "Exclusive cross section is     " << sumWgt / lumi <<" nb"<< endl;
    cout << "Exclusive cross section err is " <<err* sumWgt / lumi<<" nb"<< endl;


    pythia.stat();

    //Important in case of ISR=on, then one can tune the parameters of the importance sampling to improve event generation efficiency
    exclHooks->PrintBestParams();
    //by method
    //void ExclusiveHooks::ModifyParameters(int procId, int pRest, int _p0, int _p1, int _p2)

}


