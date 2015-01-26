//*****************************************************************************
//
// File Name	: 'prob_utils.cc'
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


#include "prob_utils.h"



std::vector<int> Variable::Discretize(std::vector<double> v)
{
    std::vector<int> ret(this->_NumDimensions);
    for(int i = 0 ; i<this->_NumDimensions;i++)
    {
        double val=v.at(i);
        if(_Circular.at(i)){
            // while(val>_Max.at(i))
            //     val-=_Max.at(i);
            // while(val<_Min.at(i))
            //     val+=_Max.at(i);

            if(val>_Max.at(i))
                val=_Min.at(i)+(val-_Max.at(i));
            else if(val<_Min.at(i))
                val=_Max.at(i)-(_Min.at(i)-val);


        }


        int disc=0;
        bool done=false;

        for(double it=_Min.at(i);it<=_Max.at(i);it+=(_Max.at(i)-_Min.at(i))/double(_Shape.at(i)))
        {

            if(val <= it){
                if(disc==0)
                {
                    ret.at(i)=disc;
                    done=true;
                    break;

                }

                else if(!done){
                    ret.at(i)=disc-1;
                    done=true;
                    break;

                }

            }
            disc++;

        }
        if(!done)
            ret.at(i)=disc-2;

    }
    return ret;


}




    //Discretize the value
int Variable::Discretize(std::vector<double> v, int dim)
{
        //std::vector<int> ret(this->_NumDimensions);
    int ret=0;

    {
        double val=v.at(dim);
        if(_Circular.at(dim)){
            // while(val>_Max.at(dim))
            //     val-=_Max.at(dim);
            // while(val<_Min.at(dim))
            //     val+=_Max.at(dim);
            if(val>_Max.at(dim))
                val=_Min.at(dim)+(val-_Max.at(dim));
            else if(val<_Min.at(dim))
                val=_Max.at(dim)-(_Min.at(dim)-val);

        }


        int disc=0;
        bool done=false;

        for(double it=_Min.at(dim);it<=_Max.at(dim);it+=(_Max.at(dim)-_Min.at(dim))/double(_Shape.at(dim)))
        {

            if(val <= it){
                if(disc==0)
                {
                    ret=disc;
                    done=true;
                    break;

                }

                else if(!done){
                    ret=disc-1;
                    done=true;
                    break;

                }

            }
            disc++;

        }
        if(!done)
            ret=disc-2;

    }
    return ret;


}



    //Continuize the value
std::vector<double> Variable::Continuize(const std::vector<int> v)
{
    std::vector<int> vec=v;

    std::vector<double> ret(this->_NumDimensions);
    for(int i = 0 ; i<this->_NumDimensions;i++)
    {
            /*
              double val=v.at(i);
              if(_Circular.at(i)){
              while(val>_Max.at(i))
              val-=_Max.at(i);
              while(val<_Min.at(i))
              val+=_Max.at(i);

              }
            */

        int disc=0;
        bool done=false;
        double val=0.0;

        if(isCircular(i))
            vec.at(i)%=_Shape.at(i);


        for(double it=_Min.at(i);it<=_Max.at(i);it+=(_Max.at(i)-_Min.at(i))/double(_Shape.at(i)-1.0))
        {
                //std::cerr<<"DEBUG "<<i<<" "<<it<<" "<<disc<<" "<<v.at(i)<<std::endl;

            if(vec.at(i) == disc){
                ret.at(i)=it;
                done=true;
                break;

            }
            disc++;
            val=it;

        }
        if(!done)
            ret.at(i)=val;

    }


    return ret;

}




double Variable::Continuize(const std::vector<int> v, int dim)
{



    int vec=v.at(dim);

        //std::vector<double> ret(this->_NumDimensions);

    double ret=0.0;


        //for(int i = 0 ; i<this->_NumDimensions;i++)
    {
            /*
              double val=v.at(i);
              if(_Circular.at(i)){
              while(val>_Max.at(i))
              val-=_Max.at(i);
              while(val<_Min.at(i))
              val+=_Max.at(i);

              }
            */

        int disc=0;
        bool done=false;
        double val=0.0;

        if(isCircular(dim))
            vec%=_Shape.at(dim);



        for(double it=_Min.at(dim);it<=_Max.at(dim);it+=(_Max.at(dim)-_Min.at(dim))/double(_Shape.at(dim)-1.0))
        {


            if(vec == disc){
                ret=it;
                done=true;
                break;

            }
            disc++;
            val=it;

        }
        if(!done)
            ret=val;

    }


    return ret;

}





int Variable::GetFlatIdx(std::vector<int> i)
{
    int idx=0;

    for(unsigned int dim=1;dim<this->_NumDimensions;dim++)
    {

        int idxout=i[dim-1];

        for(unsigned int dimout=dim;dimout<this->_NumDimensions;dimout++){
            idxout*=_Shape[dimout];

        }

        idx+=idxout;


    }
    idx+=i[this->_NumDimensions-1];
    return idx;

}
    //TODO: blindage overshoot
std::vector<int> Variable::GetIdx(int i)
{
    std::vector<int> idx(this->_NumDimensions,0);
    int mod=1;
    int div=1;



    for(int dimout=this->_NumDimensions-1;dimout>=0;dimout--)
    {


        mod=_Shape[dimout];
        idx[dimout]=int(i/div)%mod;
        div*=_Shape[dimout];
    }


    return idx;

}
