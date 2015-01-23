//*****************************************************************************
//
// File Name	: 'test_pendule.cc'
// Author	: Steve NGUYEN
// Contact      : steve.nguyen@college-de-france.fr
// Created	: mercredi, mai  1 2013
// Revised	:
// Version	:
// Target MCU	:
//
// This code is distributed under the GNU Public License
//		which can be found at http://www.gnu.org/licenses/gpl.txt
//
//
// Notes:	clang++ -c test_pendule.cc -g -O3 -I/usr/lib/openmpi/include/
//          clang++ -o test_pendule test_pendule.o -lboost_serialization -ltbb -lboost_system -lboost_mpi -lmpi -lmpi_cxx -lboost_graph_parallel -g -O3
//
//*****************************************************************************


#include "qmlib.h"
#include "prob_utils.h"

//Parameters of the problem

#define SIGMA_XPOS 0.2 //Gaussian sigma
#define _XPOS_THRES 0.001 //Thresholds are optional

#define SIGMA_XVEL 0.2
#define _XVEL_THRES 0.001

//Dimensions cardinality
#define _XPOS_CARD 51
#define _XVEL_CARD 51
#define _U_CARD 21

//Dimensions ranges
#define _XPOS_MIN -M_PI
#define _XPOS_MAX M_PI

#define _XVEL_MIN -2.5
#define _XVEL_MAX 2.5

#define _U_MIN -0.25
#define _U_MAX 0.25








/**
 *
 * Model specific: Should be redefined
 *
 */


//We define a custom Build class
class BuildPend: public Build
{
public:

    CondDistribution* _XposU;
    CondDistribution* _XvelU;
    CondDistribution* _XposvelU;


    void Simulate(std::vector<double>& Xres, std::vector<double> U, double dt) const{};

    BuildPend(CondDistribution* XposU, CondDistribution* XvelU, CondDistribution* XposvelU)
        :Build(XposvelU),_XposU(XposU),_XvelU(XvelU),_XposvelU(XposvelU){}


    void operator() (const blocked_range<size_t> &r) const {

        // std::cerr<<r.begin()<<" - "<<r.end()<<std::endl;


        for(size_t i=r.begin(); i!=r.end(); i++)
        {



            std::vector<int> Rindex=_XposU->GetRightIdx(i);

            std::vector<double> Xres;

                //Continueized
            Xres.push_back(_XposU->ContinuizeRight(Rindex,0));
            Xres.push_back(_XvelU->ContinuizeRight(Rindex,1));
            std::vector<double> Uvec;
            Uvec.push_back(_XposU->ContinuizeRight(Rindex,2));


                //FIXME Can optimize this call
            Simulate(Xres,Uvec,0.25);


                //Now let's put gaussians
            double pxpos=0.0;
            double pxvel=0.0;



                //Iterate through the left member
            for(int ite=0;ite< _XposU->GetNumLeftElements();ite++)
            {
                    //Same dimension for Xpos and Xvel
                std::vector<int> Lindex=_XposU->GetLeftIdx(ite);

                    //FIXME Circularity
                pxpos=Gaussian(_XposU->ContinuizeLeft(Lindex,0),Xres[0],SIGMA_XPOS);
                pxvel=Gaussian(_XvelU->ContinuizeLeft(Lindex,0),Xres[1],SIGMA_XVEL);



                {

                        //tbb::mutex::scoped_lock lock(_XposU->Mutex);

                        //Just to keep it sparse
                    if(pxpos>_XPOS_THRES)
                        _XposU->SetProbIdx(ite,i,pxpos);

                    if(pxvel>_XVEL_THRES)
                        _XvelU->SetProbIdx(ite,i,pxvel);

                }



            }



                //Set the full proba
                //P(Xpos,Xvel|prevXpos,prevXvel,U)
                //    =P(Xpos|prevXpos,prevXvel,U)
                //        x P(Xvel|prevXpos,prevXvel,U)


                //Project the partial prob to a simple vector


            std::vector<int> Lindex=_XposU->GetLeftIdx(0);
            std::vector<int> LindexEnd=_XposU->GetLeftIdx(_XposU->GetNumLeftElements()-1);

                //Why using slice provoke an error?
                //slice spos(_XposU->GetFlatIdx(Lindex,Rindex),1,_XposU->GetFlatIdx(LindexEnd,Rindex));

                //Test with subrange

            boost::numeric::ublas::range rpos(_XposU->GetFlatIdx(Lindex,Rindex),_XposU->GetFlatIdx(LindexEnd,Rindex));

                //Same dimensions for both
                //mapped_vector<double> xposv=project(_XposU->_FlatArray,spos);
            mapped_vector<double> xposv=project(_XposU->_FlatArray,rpos);

            mapped_vector<double> xvelv=project(_XvelU->_FlatArray,rpos);




            //iterate through non zero elements
            for(mapped_vector<double>::const_iterator itp=xposv.begin();itp!=xposv.end();++itp)
            {
                for(mapped_vector<double>::const_iterator itv=xvelv.begin();itv!=xvelv.end();++itv)
                {
                    std::vector<int> Lindexfull;

                    Lindexfull.push_back(itp.index());
                    Lindexfull.push_back(itv.index());

                        //std::cerr<<*it<<" "<<it.index()<<std::endl;
                        //std::cerr<<*itp<<" "<<*itv<<" "<<itp.index()<<" "<<itv.index()<<std::endl;



                    _XposvelU->SetProb(Lindexfull,Rindex,(*itp)*(*itv));

                }


            }

        }
    }



};



//Test with the tbbexecutor style. Not very convincing...

class BuildPendEx: public Build
{
public:

    CondDistribution* _XposU;
    CondDistribution* _XvelU;
    CondDistribution* _XposvelU;


    void Simulate(std::vector<double>& Xres, std::vector<double> U, double dt) const{};

    BuildPendEx(CondDistribution* XposU, CondDistribution* XvelU, CondDistribution* XposvelU)
        :Build(XposvelU),_XposU(XposU),_XvelU(XvelU),_XposvelU(XposvelU){


        tbb::task_scheduler_init init;
        TbbExecutor tbbExec(this);
            //std::cerr<<"DEBUG pf "<<(P->Transition->GetNumLeftElements())<<std::endl;

            //std::cerr<<"DEBUG br "<<tbb::blocked_range<size_t>(0,(P->Transition->GetNumLeftElements()))<<std::endl;


        tbb::parallel_for(tbb::blocked_range<size_t>(0,(_XposU->GetNumRightElements())),tbbExec);


    }



    void init_probs_step(size_t i)
    {

            std::vector<int> Rindex=_XposU->GetRightIdx(i);

            std::vector<double> Xres;

                //Continueized
            Xres.push_back(_XposU->ContinuizeRight(Rindex,0));
            Xres.push_back(_XvelU->ContinuizeRight(Rindex,1));
            std::vector<double> Uvec;
            Uvec.push_back(_XposU->ContinuizeRight(Rindex,2));


                //FIXME Can optimize this call
            Simulate(Xres,Uvec,0.25);


                //Now let's put gaussians
            double pxpos=0.0;
            double pxvel=0.0;



                //Iterate through the left member
            for(int ite=0;ite< _XposU->GetNumLeftElements();ite++)
            {
                    //Same dimension for Xpos and Xvel
                std::vector<int> Lindex=_XposU->GetLeftIdx(ite);

                pxpos=Gaussian(_XposU->ContinuizeLeft(Lindex,0),Xres[0],SIGMA_XPOS);
                pxvel=Gaussian(_XvelU->ContinuizeLeft(Lindex,0),Xres[1],SIGMA_XVEL);



                {

                        //tbb::mutex::scoped_lock lock(_XposU->Mutex);

                        //Just to keep it sparse
                    if(pxpos>_XPOS_THRES)
                        _XposU->SetProbIdx(ite,i,pxpos);

                    if(pxvel>_XVEL_THRES)
                        _XvelU->SetProbIdx(ite,i,pxvel);

                }



            }



                //Set the full proba
                //P(Xpos,Xvel|prevXpos,prevXvel,U)
                //    =P(Xpos|prevXpos,prevXpos,prevXvel,U)
                //        xP(Xvel|prevXpos,prevXpos,prevXvel,U)


                //Project the partial prob to a simple vector


            std::vector<int> Lindex=_XposU->GetLeftIdx(0);
            std::vector<int> LindexEnd=_XposU->GetLeftIdx(_XposU->GetNumLeftElements()-1);

                //Why using slice provoke an error?
                //slice spos(_XposU->GetFlatIdx(Lindex,Rindex),1,_XposU->GetFlatIdx(LindexEnd,Rindex));

                //Test with subrange

            boost::numeric::ublas::range rpos(_XposU->GetFlatIdx(Lindex,Rindex),_XposU->GetFlatIdx(LindexEnd,Rindex));

                //Same dimensions for both
                //mapped_vector<double> xposv=project(_XposU->_FlatArray,spos);
            mapped_vector<double> xposv=project(_XposU->_FlatArray,rpos);

            mapped_vector<double> xvelv=project(_XvelU->_FlatArray,rpos);




            //iterate through non zero elements
            for(mapped_vector<double>::const_iterator itp=xposv.begin();itp!=xposv.end();++itp)
            {
                for(mapped_vector<double>::const_iterator itv=xvelv.begin();itv!=xvelv.end();++itv)
                {
                    std::vector<int> Lindexfull;

                    Lindexfull.push_back(itp.index());
                    Lindexfull.push_back(itv.index());

                        //std::cerr<<*it<<" "<<it.index()<<std::endl;
                        //std::cerr<<*itp<<" "<<*itv<<" "<<itp.index()<<" "<<itv.index()<<std::endl;



                    _XposvelU->SetProb(Lindexfull,Rindex,(*itp)*(*itv));

                }


            }

    }






        //TBB Trick
    struct TbbExecutor
    {
    public:

        BuildPendEx *_build;


        TbbExecutor(BuildPendEx* build) : _build(build) {}

        void operator() (const tbb::blocked_range<size_t> r) const {
            for (size_t i = r.begin(); i!=r.end(); ++i) {

                    //std::cerr<<"DEBUG "<<i<<std::endl;
                _build->init_probs_step(i);

            }
        }

    };



};







/**
 *
 * Model specific: Should be redefined
 *
 */


void Simulate(std::vector<double>& X, std::vector<double> U, double dt)
{

    double dtheta=X[1]+(U[0]+sin(X[0]))*dt;
    double theta=X[0] +  X[1]*dt+(dt*dt)/2.0*(U[0]+sin(X[0]));



    if(theta>M_PI)
        theta=theta-2.0*M_PI;
    else if(theta<-M_PI)
        theta=2.0*M_PI+theta;

    X[0]=theta;
    X[1]=dtheta;


}




/**
 *
 * Model specific: Should be redefined
 *
 */

void ParallelBuild(CondDistribution* Xpos, CondDistribution* Xvel, CondDistribution* Xposvel){



        //Iterate over the right member of the probs
        // parallel_for( blocked_range<size_t>(0,Xpos->GetNumRightElements()),BuildPend(Xpos,Xvel,Xposvel),tbb::auto_partitioner());


    //Specifying the grainsize? Hard to find the optimal value...
    parallel_for( blocked_range<size_t>(0,Xpos->GetNumRightElements(),1000000),BuildPend(Xpos,Xvel,Xposvel));





}



/**
 *
 * Model specific: Should be redefined
 *
 */


void TransitionProbs::MakeDistribution()
{

        //With Conditional distributions

                //Xpos
    std::vector<int> Xcard;
    Xcard.push_back(_XPOS_CARD);

    std::vector<double> Xma;
    Xma.push_back(_XPOS_MAX);

    std::vector<double> Xmi;
    Xmi.push_back(_XPOS_MIN);

    std::vector<bool> Xb;
    Xb.push_back(true);


            //Conditional Xpos
    std::vector<int> CondXUcard;
    CondXUcard.push_back(_XPOS_CARD); //Knowing Xpos t-1
    CondXUcard.push_back(_XVEL_CARD); //Xvel t-1
    CondXUcard.push_back(_U_CARD); //U t-1

    std::vector<double> CondXUma;
    CondXUma.push_back(_XPOS_MAX);
    CondXUma.push_back(_XVEL_MAX);
    CondXUma.push_back(_U_MAX);

    std::vector<double> CondXUmi;
    CondXUmi.push_back(_XPOS_MIN);
    CondXUmi.push_back(_XVEL_MIN);
    CondXUmi.push_back(_U_MIN);

    std::vector<bool> CondXUb;
    CondXUb.push_back(true);
    CondXUb.push_back(false);
    CondXUb.push_back(false);

    CondDistribution *Xpos= new CondDistribution("Xpos_t+1", "Xpos_t,Xvel_t,U_t",Xcard,CondXUcard,Xmi,CondXUmi,Xma,CondXUma,Xb,CondXUb);



        //Xvel
    std::vector<int> Xvcard;
    Xvcard.push_back(_XVEL_CARD);

    std::vector<double> Xvma;
    Xvma.push_back(_XVEL_MAX);

    std::vector<double> Xvmi;
    Xvmi.push_back(_XVEL_MIN);

    std::vector<bool> Xvb;
    Xvb.push_back(false);


            //Conditional Xvel same as Xpos

    CondDistribution *Xvel= new CondDistribution("Xvel_t+1", "Xpos_t,Xvel_t,U_t",Xvcard,CondXUcard,Xvmi,CondXUmi,Xvma,CondXUma,Xvb,CondXUb);


                    //Xposvel
    std::vector<int> Xpvcard;
    Xpvcard.push_back(_XPOS_CARD);
    Xpvcard.push_back(_XVEL_CARD);

    std::vector<double> Xpvma;
    Xpvma.push_back(_XPOS_MAX);
    Xpvma.push_back(_XVEL_MAX);

    std::vector<double> Xpvmi;
    Xpvmi.push_back(_XPOS_MIN);
    Xpvmi.push_back(_XVEL_MIN);

    std::vector<bool> Xpvb;
    Xpvb.push_back(true);
    Xpvb.push_back(false);


    Transition= new CondDistribution("Xpos_t+1,Xvel_t+1", "Xpos_t,Xvel_t,U_t",Xpvcard,CondXUcard,Xpvmi,CondXUmi,Xpvma,CondXUma,Xpvb,CondXUb);




    std::cerr<<"Building distributions..."<<std::endl;
    ParallelBuild(Xpos,Xvel,Transition);

        //BuildPendEx Builder(Xpos,Xvel,Transition);

    std::cerr<<"\t Done."<<std::endl;


    AddDistrib(Xpos);
    AddDistrib(Xvel);

}


int main(int argc, char* argv[])
{



        //First define the variables (cardinality and range of the spaces)

    std::vector<int> Xcard;
    Xcard.push_back(_XPOS_CARD);
    Xcard.push_back(_XVEL_CARD);

    std::vector<double> Xma;
    Xma.push_back(_XPOS_MAX);
    Xma.push_back(_XVEL_MAX);

    std::vector<double> Xmi;
    Xmi.push_back(_XPOS_MIN);
    Xmi.push_back(_XVEL_MIN);

    std::vector<bool> Xb;
    Xb.push_back(true);
    Xb.push_back(false);


        //X=state, contains position and velocity
    Variable X("X",Xcard,Xmi,Xma,Xb);




    std::vector<int> Ucard;
    Ucard.push_back(_U_CARD);

    std::vector<double> Uma;
    Uma.push_back(_U_MAX);

    std::vector<double> Umi;
    Umi.push_back(_U_MIN);

    std::vector<bool> Ub;
    Ub.push_back(false);


        //U=actions
    Variable U("U",Ucard,Umi,Uma,Ub);



        //A transition probability depending on X and U
    TransitionProbs tr(X,U);


        //Init for the parallel stuff
    task_scheduler_init init;


        //Build the distrib
    tr.MakeDistribution();



        //Save the distrib

    std::ofstream file("pend_test.dat");
    boost::archive::text_oarchive oa(file);

    oa << tr;

    file.close();

        ///GRAPH

        //Init for mpi
    boost::mpi::environment env(argc,argv);

    mpi_process_group pg;


        //Build the class for the quasimetric
    QuasiMetric QM(&tr);
        //Init the graph
    QM.Init();


        //Define the objective state
    std::vector<int> obj;
    obj.push_back(_XPOS_CARD/2);
    obj.push_back(_XVEL_CARD/2);

        //Compute the quasimetric
    QM.ComputeQM(obj);



    return 0;
}
