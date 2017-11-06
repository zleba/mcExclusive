#ifndef TopoSelector_H_
#define TopoSelector_H_

#include "Basics.h"
//#include "TH1D.h"

class TopoSelector {

	Double p0, p1, p2, pRest;
	Double prob0, prob1, prob2, probRest;
  int type;

	Double sum;

	Double weights[4];
	Double weights2[4];
	Double weights4[4];
	Double entries[4];
	Double aSum, aSum2;

	void CalcSum(){
		Double ExpA = exp(-a);

		sum = pRest + (p0-pRest)*ExpA + (p1-pRest) *(1-ExpA)/a + (p2-pRest) * (a-1 + ExpA)/ (a*a/2);


		prob0 = (p0-pRest)     *ExpA / sum;
		prob1 = (p1-pRest)/a *(1-ExpA)   / sum;
		prob2 = (p2-pRest)  * (a-1 + ExpA)/(a*a/2)  / sum;

		probRest = 1. - prob0 - prob1 - prob2;

    nEm = 0;
		
    //cout <<"probabilities "<< prob0<<" "<<prob1<<" "<<prob2<<" "<<probRest<<endl;
		//cout << "prob0, prob1 " << prob0 <<" "<< prob1 <<" "<<prob2<< endl;

	}
public:
	TopoSelector() { for(int i=0; i<4;++i) weights[i]=weights2[i]=entries[i]=0; }
	
  int nEm;
	Double a;


	void SetParam(Double _a, Double _pRest, Double _p0=0, Double _p1=0, Double _p2=0)
		{  
      a=_a; 
			pRest=_pRest;
			p0 = max(_p0, _pRest);
			p1 = max(_p1, _pRest);
			p2 = max(_p2, _pRest);

      CalcSum();
      type = CalcType();
    }
	void SetPs(Double _pRest, Double _p0=0, Double _p1=0, Double _p2=0) {
		pRest=_pRest;
		p0 = max(_p0, _pRest);
		p1 = max(_p1, _pRest);
		p2 = max(_p2, _pRest);
	}
	void SetA(Double _a) {
		a = _a;
		CalcSum();
		type = CalcType();
	}

	int CalcType() const {
		Double r = rand() / (RAND_MAX+0.);
    int tt;
		if( r <= probRest )
			tt= -1;
		else if( r <= probRest+prob0 )
			tt= 0;
		else if( r <= probRest+prob0+prob1 )
			tt= 1;
		else if( r <= probRest+prob0+prob1+prob2 )
			tt= 2;
    else {
      cout << "Strange type of r: "<< r << endl;
		cout << "ProbRest,prob0,prob1,prob2: "<<probRest<<" "<<prob0<<" "<<prob1<<" "<<prob2 << endl;
      exit(1);
    }
    return tt;
	}
  int Type() const { return type; }

	Double Weight() const {
    Double w;

		if     (type == 0 && nEm == 0)
			w = 1./prob0 * (1 - pRest/p0);
		else if(type == 1 && nEm == 1)
			w = 1./prob1 * (1 - pRest/p1);
		else if(type == 2 && nEm == 2)
			w = 1./prob2 * (1 - pRest/p2);

		else if(type == -1 && nEm == 0)
			w = 1./probRest * ( pRest/p0);
		else if(type == -1 && nEm == 1)
			w = 1./probRest * ( pRest/p1);
		else if(type == -1 && nEm == 2)
			w = 1./probRest * ( pRest/p2);

		else if(type == -1 && nEm > 2)
			w = 1./probRest;
		
		else
			w = 0;

    return w;

	}
	void FillEvent(Double weight) {
		int indx;
		if(nEm > 2) 
			indx = 0;
		else
			indx = nEm+1;

		weights[indx]  += weight;
		weights2[indx] += weight*weight;
		weights4[indx] += weight*weight*weight*weight;
		++entries[indx];
		aSum  += a;
		aSum2 += a*a;

	}
	void RemoveEvent(Double weight) {
		int indx;
		if(nEm > 2) 
			indx = 0;
		else
			indx = nEm+1;

		weights[indx]  -= weight;
		weights2[indx] -= weight*weight;
		weights4[indx] -= weight*weight*weight*weight;
		--entries[indx];
		aSum  -= a;
		aSum2 -= a*a;


	}

	int GetEntries() { return entries[0]+entries[1]+entries[2]+entries[3]; }

	void PrintResult() {
		Double ppRest, pp0, pp1, pp2;
		Double qqRest, qq0, qq1, qq2;

		int Ntot    = entries[0]+entries[1]+entries[2]+entries[3];
		a = aSum/Ntot;

		Double ExpA = exp(-a);
		Double scale = (1-ExpA - a*ExpA - a*a/2.*ExpA )/weights[0];

		ppRest = 1;
		pp0 = weights[1]*scale / (   ExpA );
		pp1 = weights[2]*scale / ( a*ExpA );
		pp2 = weights[3]*scale / ( a*a/2.*ExpA );

		cout <<"Used parameters : " << pRest<<" "<<p0<<" "<< p1<<" "<< p2<<  endl;

		cout << "Statistics " << Ntot << endl;

		cout <<"Params from events : "<< ppRest <<" "<<pp0<<" "<<pp1<<" "<<pp2 << endl;

		scale = (1-ExpA - a*ExpA - a*a/2.*ExpA )/weights2[0];
		pp0 = weights2[1]*scale / (   ExpA );
		pp1 = weights2[2]*scale / ( a*ExpA );
		pp2 = weights2[3]*scale / ( a*a/2.*ExpA );

		cout <<"Params from weights : " << ppRest <<" "<<pp0<<" "<<pp1<<" "<<pp2 << endl;
		
		Double wTot = weights[0]+weights[1]+weights[2]+weights[3];
		for(int i = 0; i < 4; ++i)
			weights[i] = (weights2[i] > 0) ? weights[i] : 1;

		qqRest = ( sqrt(weights2[0])/weights[0] )  /  (  1./sqrt(weights[0]/wTot *Ntot) );
		qq0    = ( sqrt(weights2[1])/weights[1] )  /  (  1./sqrt(weights[1]/wTot *Ntot) );
		qq1    = ( sqrt(weights2[2])/weights[2] )  /  (  1./sqrt(weights[2]/wTot *Ntot) );
		qq2    = ( sqrt(weights2[3])/weights[3] )  /  (  1./sqrt(weights[3]/wTot *Ntot) );

		cout <<"Qualities : "  << qqRest <<" "<<qq0<<" "<<qq1<<" "<<qq2 << endl;
		cout <<"The worst quality : "<< max( max(qqRest,qq0), max(qq1,qq2) ) << endl;

		//for(int i =0; i<4;++i) {
			//qqEr2[i] = sqrt(weights4[i]*weights4[i]/
		//}


	}

};


#if 0
class TH1my : public TObject {
public:
	TH1D *h;     //->
	TH1D *hN;    //->
	TH1my(const char *name, const char *title, int Nbins, Double low, Double high) {
		h  = new TH1D(name, title, Nbins, low, high);
		hN = new TH1D(TString::Format("%sN",name), title, Nbins, low, high);
	}

	TH1my(TH1D *_h, TH1D *_hN) {
		h  = _h;
		hN = _hN;
	}

	void Fill(Double val, Double w=1) { 
		h ->Fill(val, w);
		hN->Fill(val);
	}
	void Write() {
		h ->Write();
		hN->Write();
	}

	void DrawQuality(const char *parm="") {
		TH1D *hTemp = (TH1D *) h->Clone();
		int nTot = hN->Integral();
		Double intTot = h->Integral();
		for(int i = 1; i <= h->GetNbinsX(); ++i) {
			if(hN->GetBinContent(i) > 0) {
				Double nIdeal = hN->Integral() * h->GetBinContent(i)/ h->Integral();
				Double res = ( h->GetBinError(i) / h->GetBinContent(i) ) / ( 1. / sqrt(nIdeal) );
				//res = 1/res;
				hTemp->SetBinContent(i, res);
				hTemp->SetBinError(i, 0);
			}
			else {
				hTemp->SetBinContent(i, 1);
				hTemp->SetBinError(i, 0);
			}
		}
		hTemp->GetYaxis()->SetTitle("#frac{#Delta^{i}}{#Delta^{i}_{ideal}}");
		hTemp->SetMinimum(0);
		hTemp->SetMaximum(8.1);
		hTemp->Draw(parm);
	}
	static Double pow2(Double x) { return x*x;}

	void DrawErrSources(const char *parm="")
	{
		TH1D *hTemp = (TH1D *) h->Clone();
		int nTot = hN->Integral();
		Double intTot = h->Integral();

		for(int i = 1; i <= h->GetNbinsX(); ++i) {
			if(hN->GetBinContent(i) > 0) {
				Double nIdeal = nTot * h->GetBinContent(i)/ intTot;
				//#define pow2(x) ((x)*(x))
				Double res =   ( pow2(h->GetBinError(i)) - pow2(h->GetBinContent(i)/sqrt(nIdeal)) )/ (  pow2(intTot/sqrt(hN->Integral())) );

				hTemp->SetBinContent(i, res);
				hTemp->SetBinError(i, 0);
			}
			else {
				hTemp->SetBinContent(i, 0);
				hTemp->SetBinError(i, 0);
			}
		}

		hTemp->SetMinimum(-0.1);
		hTemp->SetMaximum(2.5);
		hTemp->GetYaxis()->SetTitle("#frac{#Delta^{i2} - #Delta^{i2}_{ideal}}{#Delta^{tot2}_{ideal}}");
		hTemp->Draw(parm);

	}


};
#endif 


#endif
