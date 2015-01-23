//*****************************************************************************
//
// File Name	: 'qmlib.cc'
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


#include "qmlib.h"


void QuasiMetric::BuildGraph()
{
        //Iterate over the left elements of probs and create a vertex for each.

    std::cerr<<"Building the graph..."<<std::endl;

    Vertex orig;

        //MPI stuff. I have no idea what I'm doing
        //Ok it's because the graph structure is distributed...
    if (process_id(G.process_group()) == 0)
    {
        for(int i=0;i<P->Transition->GetNumLeftElements();i++)
        {

            orig.pos = P->Transition->GetLeftIdx(i);
            vertex_t start = boost::add_vertex(G);


            G[start].pos=orig.pos;
            Vertices[i]=start;

        }
    }

    synchronize(G.process_group());//???


        //WARNING! DISTRIBUTED
    std::cerr<<"\tnb vertices: "<<num_vertices(G)<<endl;

        //Then iterate through all pair to initialize distance.

    tbb::task_scheduler_init init;
    TbbExecutor tbbExec(this);
        //std::cerr<<"DEBUG pf "<<(P->Transition->GetNumLeftElements())<<std::endl;

        //std::cerr<<"DEBUG br "<<tbb::blocked_range<size_t>(0,(P->Transition->GetNumLeftElements()))<<std::endl;


    tbb::parallel_for(tbb::blocked_range<size_t>(0,(P->Transition->GetNumLeftElements())),tbbExec);



}


double QuasiMetric::quasi_dist_init(Vertex orig, Vertex dest)
{
    if(orig.pos==dest.pos)
        return 0.0;

        //Iterate for U


        //Consider inversing origin-destination!
        /* std::vector<int> left=dest.pos; */
        /* std::vector<int> right=orig.pos; */

    std::vector<int> left=orig.pos;
    std::vector<int> right=dest.pos;


    std::vector<int> tmpright;

    std::vector<int> index;

    double p=0.0;
    double mind=INF;
    double tmp;


        //std::cerr<<"DEBUG u "<<P->_U.GetNumElements()<<std::endl;


        //FIXME iterate through non null only!!
    for(int u=0;u<P->_U.GetNumElements();u++)
    {
        tmpright=right;

            //Just add the U idx on the vector.
            //as we usually have P(X|prevX,U)

        std::vector<int> tmpu=P->_U.GetIdx(u);

        tmpright.insert(tmpright.end(),tmpu.begin(),tmpu.end());

            //std::cerr<<"TEST\n";

        p=P->Transition->GetProb(left,tmpright);


        if(p>_P_EPSILON)
        {
            tmp=(cost(tmpright)/p);

            if(tmp<mind)
                mind=tmp;

        }


    }

        //return *min_element(d.begin(),d.end());
    if(mind<INF)
        return mind;
    else
        return -1.0;
}


void QuasiMetric::init_graph_step(size_t i)
{
    double dist;
    Vertex tmp;
    Vertex orig;

        //std::cerr<<"DEBUG ";

        //std::cerr<<i<<" "<<P->Transition->GetLeftIdx(i)[0]<<" "<<P->Transition->GetLeftIdx(i)[1]<<" "<<P->Transition->GetNumLeftElements()<<std::endl;

    orig.pos = P->Transition->GetLeftIdx(i);

    for(int j=0;j<P->Transition->GetNumLeftElements();j++)
    {
            //std::cerr<<"DEBUG STEP "<<j<<std::endl;

        tmp.pos = P->Transition->GetLeftIdx(j);

        if(tmp.pos==orig.pos)
            dist=0.0;

        else
            dist=quasi_dist_init(orig,tmp);




        if(dist!=-1.0)
        {

                //boost::tie(e,b) = boost::add_edge(Vertices[i],Vertices[j],G);
                //Once again some obscure boost f*ck around that one have to summon Satan to discover
                //Find this: libs/graph_parallel/test/adjlist_build_test.cpp
            graph_t::lazy_add_edge lazy=boost::add_edge(Vertices[i],Vertices[j],G); //lazy = Developper didn't want to document


            std::pair<graph_traits<graph_t>::edge_descriptor, bool> result(lazy);
            if(result.second)
            {

                    //FIXME
                G[Vertices[j]]=tmp;

                G[result.first].weight =dist;

                    /* G[result.first].second =dist; */

            }


        }


    }

}



void QuasiMetric::ComputeQM(std::vector<int> obj)
{

    std::cerr<<"Computing QuasiMetric..."<<std::endl;

    int obj_state=P->Transition->GetLeftFlatIdx(obj);


        //WARNING! DISTRIBUTED
    QD=std::vector<double>(num_vertices(G));
    std::vector<vertex_t> pre(num_vertices(G));
    std::vector<double> w(num_vertices(G));

        //dijkstra_shortest_paths(G, Vertices[obj_state], predecessor_map(&pre[0]).weight_map(get(&Edge::dist, G)).distance_map(&QD[0]));





    vertex_t s = Vertices[obj_state];



        /* dijkstra_shortest_paths(G, s, */
        /*                         distance_map(get(vertex_distance, G))); */



        /* dijkstra_shortest_paths */
        /*     (G, s, */
        /*      predecessor_map( */
        /*          make_iterator_property_map(pre.begin(), get(vertex_index, G))). */
        /*      distance_map( */
        /*          make_iterator_property_map(QD.begin(), get(vertex_index, G))) */
        /*      ); */




        /* dijkstra_shortest_paths */
        /*     (G, s, */
        /*      predecessor_map(&pre[0]). */
        /*      weight_map(get(&Edge::dist, G)). */
        /*      distance_map(&QD[0]) */
        /*      ); */



    delta_stepping_shortest_paths(G, s,
                                  dummy_property_map(),
                                  get(&Vertex::distance, G),
                                  get(&Edge::weight, G));

    property_map<graph_t, double Vertex::*>::type
        distance_map = get(&Vertex::distance, G);



        //TEST, output
    for(int i=0; i<Vertices.size();i++)
    {
            //std::cout<<"VERTEX "<<G[Vertices[i]].pos[0]<<" "<<G[Vertices[i]].pos[1]<<" ";
        std::cout<<get(distance_map,Vertices[i])<<" ";


    }

    std::cout<<std::endl;
    std::cerr<<"\tDone."<<std::endl;



}
