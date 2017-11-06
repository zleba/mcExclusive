#ifndef _ColorData_
#define _ColorData_

//#include "ColorTensor.h"
#include "Basics.h"

class ColorTensor;


class GeneralColorData {
public:
static void DataQQ( bool isColorSinglet, const ColorTensor *colL, const ColorTensor *colR, 
                    Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1);

static void DataGQ( bool isColorSinglet, const ColorTensor *colL, const ColorTensor *colR,
                   Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c);

static void DataGG( bool isColorSinglet, const ColorTensor *colL, const ColorTensor *colR,
                Double &A12, Double &A13, Double &A14,
                Double &B234, Double &B243, Double &B342, Double &B432);

private:
static bool isGGtype(const ColorTensor *colL, const ColorTensor *colR);

};


//input partons from protons are GG
class SingletGG{
private:
static void DataQQ( const ColorTensor *colL, const ColorTensor *colR, 
                    Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1);

static void DataGQ( const ColorTensor *colL, const ColorTensor *colR,
                   Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c);

static void DataGG( const ColorTensor *colL, const ColorTensor *colR,
                Double &A12, Double &A13, Double &A14,
                Double &B234, Double &B243, Double &B342, Double &B432);

friend class GeneralColorData;
};

class InclusiveGG{
private:
static void DataQQ( const ColorTensor *colL, const ColorTensor *colR, 
                    Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1);

static void DataGQ( const ColorTensor *colL, const ColorTensor *colR,
                   Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c);

static void DataGG( const ColorTensor *colL, const ColorTensor *colR,
                Double &A12, Double &A13, Double &A14,
                Double &B234, Double &B243, Double &B342, Double &B432);
friend class GeneralColorData;
};


//input partons from protons are QQ
class SingletQQ{
private:
static void DataQQ( const ColorTensor *colL, const ColorTensor *colR, 
                    Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1);

static void DataGQ( const ColorTensor *colL, const ColorTensor *colR,
                   Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c);

static void DataGG( const ColorTensor *colL, const ColorTensor *colR,
                Double &A12, Double &A13, Double &A14,
                Double &B234, Double &B243, Double &B342, Double &B432);
friend class GeneralColorData;
};

class InclusiveQQ{
private:
static void DataQQ( const ColorTensor *colL, const ColorTensor *colR, 
                    Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1);

static void DataGQ( const ColorTensor *colL, const ColorTensor *colR,
                   Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c);

static void DataGG( const ColorTensor *colL, const ColorTensor *colR,
                Double &A12, Double &A13, Double &A14,
                Double &B234, Double &B243, Double &B342, Double &B432);
friend class GeneralColorData;
};



#endif
