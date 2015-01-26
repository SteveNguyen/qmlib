//*****************************************************************************
//
// File Name	: 'qmlib.h'
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
//
//*****************************************************************************

#if !defined(QMLIB_H)
#define QMLIB_H

#include <iostream>
#include <vector>
#include <string>
#include <cfloat>
#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/split_free.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>



/* #include <boost/graph/graph_traits.hpp> */
/* #include <boost/graph/adjacency_list.hpp> */

/* #include <boost/graph/dijkstra_shortest_paths.hpp> */
/* #include <boost/graph/adjacency_matrix.hpp> */
/* #include <boost/graph/graph_utility.hpp> */
/* #include <boost/graph/copy.hpp> */

/* #include <boost/mpi.hpp> */



/////

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/copy.hpp>

// Enable PBGL interfaces to BGL algorithms
#include <boost/graph/use_mpi.hpp>

// Communication via MPI
#include <boost/graph/distributed/mpi_process_group.hpp>

// Dijkstra's single-source shortest paths algorithm
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/graph/distributed/delta_stepping_shortest_paths.hpp>

// Distributed adjacency list
#include <boost/graph/distributed/adjacency_list.hpp>


#include <boost/property_map/parallel/distributed_property_map.hpp>



#include <tbb/tbb.h>
#include <tbb/mutex.h>


#include "prob_utils.h"





using namespace boost::numeric::ublas;
using namespace tbb;
using namespace std;



/**
 *
 * Definition of the graph stuff
 *
 */

using namespace boost;
using boost::graph::distributed::mpi_process_group;



/**
 *
 * A vertex containing the quasidistance
 *
 */


class Vertex
{
  public:
        //tuple<int ,int> pos;

    std::vector<int> pos;

    double distance;



    friend class boost::serialization::access;

    template <class Archive>
        void save(Archive &ar,  const unsigned int version) const
    {

        int s=pos.size(); //pfff for some f***ing reason I need to do that using a temporary int...
        ar << s;

        for(int i=0;i<pos.size();i++)
            ar << pos[i];

    }

    template <class Archive>
        void load(Archive &ar,  const unsigned int version)
    {

        int s=0;

        ar >> s;
        pos=std::vector<int>(s);

        for(int i=0;i<s;i++)
            ar >> pos[i];



    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();

};



/**
 *
 * Definition of the edge class for the graph
 *
 */


class Edge
{
  public:

    double dist;
    double weight;


    friend class boost::serialization::access;

    template <class Archive>
        void save(Archive &ar,  const unsigned int version) const
    {

        ar << dist;
    }

    template <class Archive>
        void load(Archive &ar,  const unsigned int version)
    {


        ar >> dist;

    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();


};



/**
 * Spécific for parallel graph
 */

typedef boost::adjacency_list<boost::listS, distributedS<mpi_process_group, vecS>, boost::bidirectionalS, Vertex, Edge > graph_t;


typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
typedef boost::graph_traits<graph_t>::edge_descriptor edge_t;



// Todo:
// definir v_array: vecteur de vectex_t. Class avec getflatidx.
// definie QuasiDist: Class pour faire joli avec la gestion du graph
//



/**
 *
 * Class for the Quasimetric.
 *
 */


class QuasiMetric
{
  public:

    TransitionProbs *P;
    graph_t G;

    std::vector<vertex_t> Vertices;
    Variable VerticesVar;

    std::vector<double> QD;

    QuasiMetric(){}

  QuasiMetric(TransitionProbs *tr):P(tr){
            //TODO: init VerticesVar

    }

    ~QuasiMetric(){}


    virtual double cost(std::vector<int> right){return 1.0;}



    void Init()
    {
            //Number of vertices
        Vertices = std::vector<vertex_t> (P->Transition->GetNumLeftElements());

        BuildGraph();

        std::cerr<<"\tDone."<<std::endl;

    }

    void BuildGraph();

    double quasi_dist_init(Vertex orig, Vertex dest);

    void init_graph_step(size_t i);

    void ComputeQM(std::vector<int> obj);


        //TBB Trick
    struct TbbExecutor
    {
    public:

        QuasiMetric *_qm;


        TbbExecutor(QuasiMetric* qm) : _qm(qm) {}

        void operator() (const tbb::blocked_range<size_t> r) const {
            for (size_t i = r.begin(); i!=r.end(); ++i) {

                    //std::cerr<<"DEBUG "<<i<<std::endl;
                _qm->init_graph_step(i);

            }
        }

    };






};







/**
 *
 * Tool class used to build the distribution.
 * Usually needs to be derived for model specific usage.
 *
 *
 *
 */


class Build
{
public:

    CondDistribution* _XU;

    std::vector<double> _sigmas;


    Build(CondDistribution* XU)
        :_XU(XU){}

  Build(CondDistribution* XU, std::vector<double> sigmas)
      :_XU(XU), _sigmas(sigmas){}



    double Gaussian(double x, double mu, double sigma, double nbsigma=3.0) const
    {
        double dist_x_mu=x-mu;

        if(abs(dist_x_mu)>sigma*nbsigma)
            return 0.0;
        else
            return (1.0/(sigma*sqrt(2.0*M_PI)))*exp(-((dist_x_mu)*(dist_x_mu))/(2.0*(sigma*sigma)));

    }


    virtual void Simulate(std::vector<double>& Xres, std::vector<double> U, double dt)const{};



        /**
         *
         * We receive r as a range for the iteration among proba.
         * It relates to the right part of the conditional distribution.
         *
         */

    void operator() (const blocked_range<size_t> &r) const {


        for(size_t i=r.begin(); i!=r.end(); i++)
        {

                //Transform the iterator into a usable index for the right part of the proba
            std::vector<int> Rindex=_XU->GetRightIdx(i);

                //to be used to store the continuous values for the state
            std::vector<double> Xres;

                //Continuize the value corresponding to the discrete Rindex.
                //State part.

            for(std::vector<int>::const_iterator it = _XU->GetRStateDim().begin() ; it != _XU->GetRStateDim().end(); ++it)
                Xres.push_back(_XU->ContinuizeRight(Rindex,*it));

                //to be used to store the continuous values for the action
            std::vector<double> Uvec;

                //Continuize the value corresponding to the discrete Rindex.
                //Action part.

            for(std::vector<int>::const_iterator it = _XU->GetRActionDim().begin() ; it != _XU->GetRActionDim().end(); ++it)
                Uvec.push_back(_XU->ContinuizeRight(Rindex,*it));



                //FIXME Can optimize this call
            Simulate(Xres,Uvec,0.25);


                //Now let's put gaussians around the simulated values

            double gaussres=0.0;

                //Iterate through the left member
            for(int ite=0;ite< _XU->GetNumLeftElements();ite++)
            {

                std::vector<int> Lindex=_XU->GetLeftIdx(ite);

                    //For each dimension of the left part
                for(uint ii=0; ii<_XU->GetLeftNumDimensions(); ii++)
                {
                    gaussres=Gaussian(_XU->ContinuizeLeft(Lindex,ii),Xres[ii],_sigmas[ii]);





                        //tbb::mutex::scoped_lock lock(_XposU->Mutex);

                        //Just to keep it sparse
                    if(gaussres>_PROBS_THRES)
                        _XU->SetProbIdx(ite,i,gaussres);

                }



            }



                //Project the partial prob to a simple vector


            std::vector<int> Lindex=_XU->GetLeftIdx(0);
            std::vector<int> LindexEnd=_XU->GetLeftIdx(_XU->GetNumLeftElements()-1);

                //Why using slice provoke an error?
                //slice spos(_XposU->GetFlatIdx(Lindex,Rindex),1,_XposU->GetFlatIdx(LindexEnd,Rindex));

                //Test with subrange

            boost::numeric::ublas::range rpos(_XU->GetFlatIdx(Lindex,Rindex),_XU->GetFlatIdx(LindexEnd,Rindex));


            mapped_vector<double> xall=project(_XU->_FlatArray,rpos);

                //iterate through non zero elements
            for(mapped_vector<double>::const_iterator itp=xall.begin();itp!=xall.end();++itp)
            {
                std::vector<int> Lindexfull;

                Lindexfull.push_back(itp.index());


                        //std::cerr<<*it<<" "<<it.index()<<std::endl;
                        //std::cerr<<*itp<<" "<<*itv<<" "<<itp.index()<<" "<<itv.index()<<std::endl;



                _XU->SetProb(Lindexfull,Rindex,(*itp));

            }




        }
    }



};





#endif
