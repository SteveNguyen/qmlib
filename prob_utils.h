//*****************************************************************************
//
// File Name	: 'prob_utils.h'
// Author	: Steve NGUYEN
// Contact      : steve.nguyen@college-de-france.fr
// Created	: jeudi, mai  2 2013
// Revised	:
// Version	:
// Target MCU	:
//
// This code is distributed under the GNU Public License
//		which can be found at http://www.gnu.org/licenses/gpl.txt
//
//
// Notes:	notes
//
//*****************************************************************************

#if !defined(PROB_UTILS_H)
#define PROB_UTILS_H

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <tbb/tbb.h>
#include <tbb/mutex.h>

#define INF 10e10
#define _P_EPSILON 0.0001 //Value considered to be small enough not be be stored
#define _PROBS_THRES 0.001 //Minimum proba


using namespace boost::numeric::ublas;
using namespace std;
using namespace tbb;

/**
 *
 * Definition of a Variable class containing tools for a probabilistic variable
 *
 *
 */

class Variable{
  public:
    std::string label;
    int _NumDimensions;
    std::vector<int> _Shape;
    int _NumElements;

    std::vector<double> _Min;
    std::vector<double> _Max;
    std::vector<bool> _Circular;

    Variable(){
        label="";
        _NumDimensions=0;
        _NumElements=0;
    }


  Variable(const std::string l, const std::vector<int> Shape)
      :label(l), _Shape(Shape)
    {
        _NumDimensions=_Shape.size();

        int total=1;
        for(std::vector<int>::const_iterator it = _Shape.begin() ; it != _Shape.end(); ++it)
            total *= *it;

        _NumElements=total;

    }


  Variable(const std::string l, const std::vector<int> Shape, const std::vector<double> minimums, const std::vector<double> maximums, const std::vector<bool> circular)
      :label(l), _Shape(Shape), _Min(minimums), _Max(maximums), _Circular(circular)
    {
        _NumDimensions=_Shape.size();

        int total=1;
        for(std::vector<int>::const_iterator it = _Shape.begin() ; it != _Shape.end(); ++it)
            total *= *it;

        _NumElements=total;


    }


    Variable(const Variable & a)
    {
        label=a.label;
        _NumDimensions=a._NumDimensions;
        _Shape=a._Shape;
        _NumElements=a._NumElements;
        _Min=a._Min;
        _Max=a._Max;
        _Circular=a._Circular;

    }


    ~Variable(){}

    std::string GetLabel(){return label;}
    int GetNumElements() const {return _NumElements;}
    int GetNumDimensions(){return _NumDimensions;}
    std::vector<int> GetShape(){return _Shape;}
    bool isCircular(int dim){return _Circular.at(dim);}
    double GetMax(int dim){return _Max.at(dim);}
    double GetMin(int dim){return _Min.at(dim);}

    Variable GetOnlyDims(std::vector<int> dims)
    {
        std::vector<int> shape;
        std::vector<double> mins;
        std::vector<double> maxs;
        std::vector<bool> circular;


        for(int i=0;i<dims.size();i++)
        {
            shape.push_back(_Shape.at(dims[i]));
            mins.push_back(_Min.at(dims[i]));
            maxs.push_back(_Max.at(dims[i]));
            circular.push_back(_Circular.at(dims[i]));

        }

        return Variable("",shape,mins,maxs,circular);

    }

    std::vector<int> Discretize(std::vector<double> v);

        //Discretize the value
    int Discretize(std::vector<double> v, int dim);

        //Continuize the value
    std::vector<double> Continuize(const std::vector<int> v);

    double Continuize(const std::vector<int> v, int dim);

    int GetFlatIdx(std::vector<int> i);


    std::vector<int> GetIdx(int i);


        //Serialization stuff
    friend class boost::serialization::access;

    template <class Archive>
        void save(Archive &ar, const unsigned int version) const
    {

        ar<<label;
            //ar << boost::serialization::base_object<Distribution>(*this);

        ar << _NumDimensions;
        ar << _NumElements;


        for(std::vector<int>::const_iterator it = _Shape.begin() ; it != _Shape.end(); ++it)
            ar<<*it;


    }

    template <class Archive>
        void load(Archive &ar, const unsigned int version)
    {

        ar>>label;

            //ar >> boost::serialization::base_object<Distribution>(*this);
        ar >> _NumDimensions;
        ar >> _NumElements;

        int card=0;
        for(int i=0; i<_NumDimensions;i++)
        {
            ar >> card;
            _Shape.push_back(card);
        }


    }


    BOOST_SERIALIZATION_SPLIT_MEMBER();

};


/**
 *
 * Definition of a class with tools for handeling a distribution of probability
 *
 *
 */

class Distribution{
  public:


    Variable Var;

        //std::vector<double> _FlatArray;
    mapped_vector<double> _FlatArray;

    Distribution(){}

  Distribution(const std::string l, const std::vector<int> Shape)
      :Var(l,Shape)
    {
            //_FlatArray= std::vector<double>(_NumElements,0.0);
        _FlatArray= mapped_vector<double>(Var._NumElements,0.0);


    }


  Distribution(const std::string l, const std::vector<int> Shape, const std::vector<double> minimums, const std::vector<double> maximums, const std::vector<bool> circular)
      :Var(l,Shape,minimums,maximums,circular)
    {
            //_FlatArray= std::vector<double>(_NumElements,0.0);
        _FlatArray= mapped_vector<double>(Var._NumElements,0.0);


    }


    Distribution(const Distribution & a)
    {
        Var=a.Var;

        _FlatArray=a._FlatArray;


    }


    ~Distribution(){}
    std::string GetLabel(){return Var.label;}
    int GetNumElements() const {return Var._NumElements;}

    std::vector<int> Discretize(std::vector<double> v)
    {
        return Var.Discretize(v);
    }

    std::vector<double> Continuize(std::vector<int> v)
    {
        return Var.Continuize(v);
    }


    int Discretize(std::vector<double> v, int dim)
    {
        return Var.Discretize(v,dim);
    }

    double Continuize(std::vector<int> v, int dim)
    {
        return Var.Continuize(v,dim);
    }


    int GetFlatIdx(std::vector<int> i)
    {
        return Var.GetFlatIdx(i);
    }


    std::vector<int> GetIdx(int i)
    {
        return Var.GetIdx(i);
    }

    double GetProb(std::vector<int> i)
    {
        return _FlatArray(Var.GetFlatIdx(i));
    }

    double GetProbIdx(int i)
    {
        return _FlatArray(i);
    }

    void SetProb(std::vector<int> i, double p){
        if(p>_P_EPSILON)
            _FlatArray(Var.GetFlatIdx(i))=p;
    }
    void SetProbIdx(int i, double p){
        if(p>_P_EPSILON)
            _FlatArray(i)=p;
    }


    int GetNumDimensions(){return Var._NumDimensions;}
    std::vector<int> GetShape(){return Var._Shape;}
    bool isCircular(int dim){return Var._Circular.at(dim);}
    double GetMax(int dim){return Var._Max.at(dim);}
    double GetMin(int dim){return Var._Min.at(dim);}


    friend class boost::serialization::access;

    template <class Archive>
        void save(Archive &ar, const unsigned int version) const
    {


        ar<<Var;


            //Standard vector
            /*
              for(std::vector<double>::const_iterator it = _FlatArray.begin() ; it != _FlatArray.end(); ++it)
              ar<<*it;
            */
            //Sparse vector

        ar<<_FlatArray;

    }

    template <class Archive>
        void load(Archive &ar, const unsigned int version)
    {

        ar >> Var;

            //Standard vector
            /*
              double prob=0.0;
              for(int i=0; i<_NumElements;i++)
              {
              ar >> prob;
              _FlatArray.push_back(prob);
              }
            */

        ar >> _FlatArray;

    }


    BOOST_SERIALIZATION_SPLIT_MEMBER();



};



/**
 *
 * Child class from Distribution for conditional distribution with a left and a right member
 *
 */

class CondDistribution: public Distribution{
  public:


    std::vector<int> _LShape;
    std::vector<int> _RShape;


        //to store the dimension index for states and actions on the right part of the distrib
    std::vector<int> _RStateDim;
    std::vector<int> _RActionDim;


    Variable Cond;
    Variable VarCond;

    int TotalElements;

    tbb::mutex Mutex;


  CondDistribution(const std::string ll,const std::string rl, const std::vector<int> LShape, const std::vector<int> RShape)
      :VarCond(ll,LShape), Cond(rl,RShape)
    {
            //_FlatArray= std::vector<double>(_NumElements,0.0);
        _FlatArray= mapped_vector<double>(VarCond._NumElements*Cond._NumElements,0.0);
        TotalElements=Var._NumElements*Cond._NumElements;

            //In fact it is much more convenient to invert the storage order
            //std::vector<int> S=LShape;
            //S.insert(S.end(),RShape.begin(),RShape.end());

        std::vector<int> S=RShape;
        S.insert(S.end(),LShape.begin(),LShape.end());

        std::string l=ll+"|"+rl;
        Var=Variable(l,S);


    }



  CondDistribution(const std::string ll,const std::string rl, const std::vector<int> LShape, const std::vector<int> RShape, const std::vector<int> RStateDim, const std::vector<int> RActionDim)
      :VarCond(ll,LShape), Cond(rl,RShape)
    {
            //_FlatArray= std::vector<double>(_NumElements,0.0);
        _FlatArray= mapped_vector<double>(VarCond._NumElements*Cond._NumElements,0.0);
        TotalElements=Var._NumElements*Cond._NumElements;

            //In fact it is much more convenient to invert the storage order
            //std::vector<int> S=LShape;
            //S.insert(S.end(),RShape.begin(),RShape.end());

        _RStateDim=RStateDim;
        _RActionDim=RActionDim;


        std::vector<int> S=RShape;
        S.insert(S.end(),LShape.begin(),LShape.end());

        std::string l=ll+"|"+rl;
        Var=Variable(l,S);


    }




  CondDistribution(const std::string ll, const std::string rl, const std::vector<int> LShape, const std::vector<int> RShape, const std::vector<double> Lminimums, const std::vector<double> Rminimums, const std::vector<double> Lmaximums, const std::vector<double> Rmaximums, const std::vector<bool> Lcircular, const std::vector<bool> Rcircular)
      :VarCond(ll,LShape,Lminimums,Lmaximums,Lcircular), Cond(rl,RShape,Rminimums,Rmaximums,Rcircular)
    {
            //_FlatArray= std::vector<double>(_NumElements,0.0);
        _FlatArray= mapped_vector<double>(VarCond._NumElements*Cond._NumElements,0.0);

        TotalElements=Var._NumElements*Cond._NumElements;


            //In fact it is much more convenient to invert the storage order
            //std::vector<int> S=LShape;
            //S.insert(S.end(),RShape.begin(),RShape.end());

        std::vector<int> S=RShape;
        S.insert(S.end(),LShape.begin(),LShape.end());

        std::string l=ll+"|"+rl;


            /* std::vector<double> mi=Lminimums; */
            /* mi.insert(mi.end(),Rminimums.begin(),Rminimums.end()); */

            /* std::vector<double> ma=Lmaximums; */
            /* ma.insert(ma.end(),Rmaximums.begin(),Rmaximums.end()); */

            /* std::vector<bool> b=Lcircular; */
            /* b.insert(b.end(),Rcircular.begin(),Rcircular.end()); */


        std::vector<double> mi=Rminimums;
        mi.insert(mi.end(),Lminimums.begin(),Lminimums.end());

        std::vector<double> ma=Rmaximums;
        ma.insert(ma.end(),Lmaximums.begin(),Lmaximums.end());

        std::vector<bool> b=Rcircular;
        b.insert(b.end(),Lcircular.begin(),Lcircular.end());


        Var=Variable(l,S,mi,ma,b);

    }



  CondDistribution(const std::string ll, const std::string rl, const std::vector<int> LShape, const std::vector<int> RShape, const std::vector<double> Lminimums, const std::vector<double> Rminimums, const std::vector<double> Lmaximums, const std::vector<double> Rmaximums, const std::vector<bool> Lcircular, const std::vector<bool> Rcircular, const std::vector<int> RStateDim, const std::vector<int> RActionDim)
      :VarCond(ll,LShape,Lminimums,Lmaximums,Lcircular), Cond(rl,RShape,Rminimums,Rmaximums,Rcircular)
    {
            //_FlatArray= std::vector<double>(_NumElements,0.0);
        _FlatArray= mapped_vector<double>(VarCond._NumElements*Cond._NumElements,0.0);

        TotalElements=Var._NumElements*Cond._NumElements;


            //In fact it is much more convenient to invert the storage order
            //std::vector<int> S=LShape;
            //S.insert(S.end(),RShape.begin(),RShape.end());


        _RStateDim=RStateDim;
        _RActionDim=RActionDim;


        std::vector<int> S=RShape;
        S.insert(S.end(),LShape.begin(),LShape.end());

        std::string l=ll+"|"+rl;


            /* std::vector<double> mi=Lminimums; */
            /* mi.insert(mi.end(),Rminimums.begin(),Rminimums.end()); */

            /* std::vector<double> ma=Lmaximums; */
            /* ma.insert(ma.end(),Rmaximums.begin(),Rmaximums.end()); */

            /* std::vector<bool> b=Lcircular; */
            /* b.insert(b.end(),Rcircular.begin(),Rcircular.end()); */


        std::vector<double> mi=Rminimums;
        mi.insert(mi.end(),Lminimums.begin(),Lminimums.end());

        std::vector<double> ma=Rmaximums;
        ma.insert(ma.end(),Lmaximums.begin(),Lmaximums.end());

        std::vector<bool> b=Rcircular;
        b.insert(b.end(),Lcircular.begin(),Lcircular.end());


        Var=Variable(l,S,mi,ma,b);

    }








        //Todo also from Distribution?
    CondDistribution(const CondDistribution & a)
    {
        Var=a.Var;
        VarCond=a.VarCond;
        Cond=a.Cond;

        TotalElements=a.TotalElements;

        _RStateDim=a._RStateDim;
        _RActionDim=a._RActionDim;


        _FlatArray=a._FlatArray;


    }


    ~CondDistribution(){}
    std::string GetLeftLabel(){return VarCond.label;}
    std::string GetRightLabel(){return Cond.label;}

    int GetNumLeftElements() const {return VarCond._NumElements;}
    int GetNumRightElements() const {return Cond._NumElements;}

    int GetNumElements() const {return TotalElements;}

    std::vector<int> GetRStateDim() const {return _RStateDim;}
    std::vector<int> GetRActionDim() const {return _RActionDim;}


    std::vector<int> DiscretizeLeft(std::vector<double> v)
    {
        return VarCond.Discretize(v);
    }

    std::vector<double> ContinuizeLeft(std::vector<int> v)
    {
        return VarCond.Continuize(v);
    }

    std::vector<int> DiscretizeRight(std::vector<double> v)
    {
        return Cond.Discretize(v);
    }

    std::vector<double> ContinuizeRight(std::vector<int> v)
    {
        return Cond.Continuize(v);
    }



        //Only one dim

    int DiscretizeLeft(std::vector<double> v, int dim)
    {
        return VarCond.Discretize(v,dim);
    }

    double ContinuizeLeft(std::vector<int> v, int dim)
    {
        return VarCond.Continuize(v,dim);
    }

    int DiscretizeRight(std::vector<double> v, int dim)
    {
        return Cond.Discretize(v,dim);
    }

    double ContinuizeRight(std::vector<int> v, int dim)
    {
        return Cond.Continuize(v,dim);
    }





    int GetLeftFlatIdx(std::vector<int> i)
    {
        return VarCond.GetFlatIdx(i);
    }

    std::vector<int> GetLeftIdx(int i)
    {
        return VarCond.GetIdx(i);
    }

    int GetRightFlatIdx(std::vector<int> i)
    {
        return Cond.GetFlatIdx(i);
    }

    std::vector<int> GetRightIdx(int i)
    {
        return Cond.GetIdx(i);
    }




        ///////
        //TODO: check



        //Get global index knowing left and right
    int GetFlatIdx(std::vector<int> left, std::vector<int> right)
    {
            /* std::vector<int> i=left; */
            /* i.insert(i.end(),right.begin(),right.end());  */

        std::vector<int> i=right;
        i.insert(i.end(),left.begin(),left.end());


        return Var.GetFlatIdx(i);
    }


    std::vector<int> GetIdx(int left, int right)
    {
            /* std::vector<int> i=GetLeftIdx(left); */
            /* std::vector<int> tmp=GetRightIdx(right); */
            /* i.insert(i.end(),tmp.begin(),tmp.end());  */

        std::vector<int> i=GetRightIdx(right);
        std::vector<int> tmp=GetLeftIdx(left);
        i.insert(i.end(),tmp.begin(),tmp.end());


        return i;
    }



        //Knowing the right part
    double GetProb(std::vector<int> left, std::vector<int> right)
    {
            /* std::vector<int> i=left; */
            /* i.insert(i.end(),right.begin(),right.end());  */

        std::vector<int> i=right;
        i.insert(i.end(),left.begin(),left.end());

        return _FlatArray(Var.GetFlatIdx(i));
    }

    double GetProbIdx(int left, int right)
    {
            /* std::vector<int> i=GetLeftIdx(left); */
            /* std::vector<int> tmp=GetRightIdx(right); */
            /* i.insert(i.end(),tmp.begin(),tmp.end());  */

        std::vector<int> i=GetRightIdx(right);
        std::vector<int> tmp=GetLeftIdx(left);
        i.insert(i.end(),tmp.begin(),tmp.end());

        return _FlatArray(Var.GetFlatIdx(i));
    }

    void SetProb(std::vector<int> left, std::vector<int> right, double p)
    {

        if(p>_P_EPSILON){
                /* std::vector<int> i=left; */
                /* i.insert(i.end(),right.begin(),right.end()); */
            std::vector<int> i=right;
            i.insert(i.end(),left.begin(),left.end());

            {

                tbb::mutex::scoped_lock lock(Mutex);
                _FlatArray(Var.GetFlatIdx(i))=p;
            }

        }
    }

    void SetProbIdx(int left, int right, double p)
    {
        if(p>_P_EPSILON){
                /* std::vector<int> i=GetLeftIdx(left); */
                /* std::vector<int> tmp=GetRightIdx(right); */
                /* i.insert(i.end(),tmp.begin(),tmp.end());  */

            std::vector<int> i=GetRightIdx(right);
            std::vector<int> tmp=GetLeftIdx(left);
            i.insert(i.end(),tmp.begin(),tmp.end());

            {
                    //These mutex seems to be important...
                tbb::mutex::scoped_lock lock(Mutex);
                _FlatArray(Var.GetFlatIdx(i))=p;
            }

        }

    }




    int GetLeftNumDimensions(){return VarCond._NumDimensions;}
    std::vector<int> GetLeftShape(){return VarCond._Shape;}
    bool isCircularLeft(int dim){return VarCond._Circular.at(dim);}
    double GetLeftMax(int dim){return VarCond._Max.at(dim);}
    double GetLeftMin(int dim){return VarCond._Min.at(dim);}


    int GetRightNumDimensions(){return Cond._NumDimensions;}
    std::vector<int> GetRightShape(){return Cond._Shape;}
    bool isCircularRight(int dim){return Cond._Circular.at(dim);}
    double GetRightMax(int dim){return Cond._Max.at(dim);}
    double GetRightMin(int dim){return Cond._Min.at(dim);}



    friend class boost::serialization::access;

    template <class Archive>
        void save(Archive &ar, const unsigned int version) const
    {


        ar<<Var;
        ar<<VarCond;
        ar<<Cond;
        ar<<_RStateDim;
        ar<<_RActionDim;

            //Standard vector
            /*
              for(std::vector<double>::const_iterator it = _FlatArray.begin() ; it != _FlatArray.end(); ++it)
              ar<<*it;
            */
            //Sparse vector

        ar<<_FlatArray;

    }

    template <class Archive>
        void load(Archive &ar, const unsigned int version)
    {

        ar >> Var;
        ar >> VarCond;
        ar >> Cond;
        ar >>_RStateDim;
        ar >>_RActionDim;

            //Standard vector
            /*
              double prob=0.0;
              for(int i=0; i<_NumElements;i++)
              {
              ar >> prob;
              _FlatArray.push_back(prob);
              }
            */

        ar >> _FlatArray;

    }


    BOOST_SERIALIZATION_SPLIT_MEMBER();



};



/**
 *
 * Definition of a class for handeling transition probabilities with states and actions
 *
 */

class TransitionProbs{
  public:

    std::vector<Distribution *> Distrib;
    unsigned int _NumDistrib;

    Variable _X;
    Variable _U;
    CondDistribution *Transition;

    TransitionProbs(){_NumDistrib=0;}

    ~TransitionProbs(){
        delete Transition;
    }

  TransitionProbs(Variable X, Variable U):_X(X),_U(U){_NumDistrib=0;}

    void AddDistrib(Distribution *d){
        Distrib.push_back(d);
        _NumDistrib++;

    }



        //as this vector will usually be short...
    Distribution* Get(std::string l){
        for(std::vector<Distribution*>::iterator it = Distrib.begin() ; it != Distrib.end(); ++it)
            if((*it)->GetLabel()==l)
                return (*it);
        return NULL;
    }


    virtual void MakeDistribution();



    friend class boost::serialization::access;

    template <class Archive>
        void save(Archive &ar,  const unsigned int version) const
    {
        ar<<_X;
        ar<<_U;
        ar<<_NumDistrib;


        for(std::vector<Distribution*>::const_iterator it = Distrib.begin() ; it != Distrib.end(); ++it)
            ar << *it;


        ar<<Transition;


    }

    template <class Archive>
        void load(Archive &ar,  const unsigned int version)
    {

        ar>>_X;
        ar>>_U;

        ar>>_NumDistrib;

        Distribution *d;

        for(int i=0;i<_NumDistrib;i++)
        {

            ar >> d;
            Distrib.push_back(d);


        }

        ar>>Transition;


    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();



};












#endif
