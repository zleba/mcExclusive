#include "ColorData.h"
#include "ColorTensor.h"
#include <vector>
using std::vector;

bool GeneralColorData::isGGtype(const ColorTensor *colL, const ColorTensor *colR)
{
	if(colL->GetNow() == "G" && colR->GetNow() == "G")
		return true;
	else if( (colL->GetNow() == "Q"  && colR->GetNow() == "Qx")
		|| (colL->GetNow() == "Qx" && colR->GetNow() == "Q" ) )
		return false;
	else {
		cout << "Something wrong with incoming flavour type" << endl;
		exit(1);
	}

}


void GeneralColorData::DataQQ( bool isColorSinglet, const ColorTensor *colL, const ColorTensor *colR, 
                        Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1)
{
	typedef void (*FunPtr)(const ColorTensor *, const ColorTensor *, Double &, Double &, Double &);

	FunPtr funQQ[2][2]={ {&SingletGG::DataQQ,   &SingletQQ::DataQQ  },
	                     {&InclusiveGG::DataQQ, &InclusiveQQ::DataQQ} };

	if( !( (colL->GetOrg() == "Q" || colL->GetOrg() == "Qx") &&
	       (colR->GetOrg() == "Q" || colR->GetOrg() == "Qx") ) ) {
		cout << "Incompatible flavour type of hard process" << endl;
		cout << "Q(x)Q(x) expected, in reality " << 
		        colL->GetOrg() <<" "<<colR->GetOrg() << endl;
		exit(1);
	}

	bool isGG =  isGGtype(colL, colR);

	funQQ[!isColorSinglet][!isGG](colL, colR, Dl1l1c, Dl1r1c, Dl1r1);

}

void GeneralColorData::DataGQ( bool isColorSinglet, const ColorTensor *colL, const ColorTensor *colR,
                   Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c)
{

	typedef void (*FunPtr)(const ColorTensor *, const ColorTensor *, Double &, Double &, Double &);

	FunPtr funGQ[2][2]={ {&SingletGG::DataGQ,   &SingletQQ::DataGQ  },
                       {&InclusiveGG::DataGQ, &InclusiveQQ::DataGQ} };

	if( !( ( (colL->GetOrg() == "Q" || colL->GetOrg() == "Qx") && colR->GetOrg() == "G") ||
	       ( (colR->GetOrg() == "Q" || colR->GetOrg() == "Qx") && colL->GetOrg() == "G") ) ) {
		cout << "Incompatible flavour type of hard process" << endl;
		cout << "QG or QxG expected, in reality " << 
		        colL->GetOrg() <<" "<<colR->GetOrg() << endl;
		exit(1);
	}

	bool isGG =  isGGtype(colL, colR);

	funGQ[!isColorSinglet][!isGG](colL, colR, Delta, Tl1cl1C, Tl1l1c);

}


void GeneralColorData::DataGG( bool isColorSinglet, const ColorTensor *colL, const ColorTensor *colR,
                Double &A12, Double &A13, Double &A14,
                Double &B234, Double &B243, Double &B342, Double &B432)
{

	typedef void (*FunPtr)(const ColorTensor *, const ColorTensor *,
	        Double &, Double &, Double &, Double &, Double &, Double &, Double &);

	FunPtr funGG[2][2]={ {&SingletGG::DataGG,   &SingletQQ::DataGG  },
	                     {&InclusiveGG::DataGG, &InclusiveQQ::DataGG} };

	if( !( colL->GetOrg() == "G"  && colR->GetOrg() == "G") ) {
		cout << "Incompatible flavour type of hard process" << endl;
		cout << "GG expected, in reality " << 
		        colL->GetOrg() <<" "<<colR->GetOrg() << endl;
		exit(1);
	}

	bool isGG =  isGGtype(colL, colR);

	funGG[!isColorSinglet][!isGG](colL, colR, A12, A13, A14,
	                                    B234, B243, B342, B432);
}
